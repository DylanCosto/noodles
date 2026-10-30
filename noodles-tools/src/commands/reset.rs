//! `ntools reset`: turn alignments back into unaligned reads (samtools reset).
//!
//! This follows samtools 1.24's `reset.c`: secondary and supplementary records are dropped;
//! the rest become unmapped, in their original orientation (reverse-strand reads are reverse
//! complemented back), with mate, pair, and (by default) duplicate information cleared and
//! alignment tags removed. The header keeps `@RG` lines and `@PG` lines (up to `--reject-PG`).

use std::{collections::HashSet, path::PathBuf};

use anyhow::bail;
use noodles_sam as sam;

use crate::{
    alignment::{AlignmentReader, RawEncoder, ReadOptions},
    output::{self, AlignmentWriter, OutputFormat, WriteOptions},
    program,
};

/// The 4-bit code of each base's complement (`seqi_rc`).
const COMPLEMENT: [u8; 16] = [0, 8, 4, 12, 2, 10, 6, 14, 1, 9, 5, 13, 3, 11, 7, 15];

/// Tags removed by default: alignment-specific ones.
const ALIGNMENT_TAGS: [&[u8; 2]; 15] = [
    b"AS", b"CC", b"CG", b"CP", b"H1", b"H2", b"HI", b"H0", b"IH", b"MC", b"MD", b"MQ", b"NM",
    b"SA", b"TS",
];

const FLAG_PAIRED: u16 = 0x1;
const FLAG_PROPER_PAIR: u16 = 0x2;
const FLAG_UNMAPPED: u16 = 0x4;
const FLAG_MATE_UNMAPPED: u16 = 0x8;
const FLAG_REVERSE: u16 = 0x10;
const FLAG_MATE_REVERSE: u16 = 0x20;
const FLAG_SECONDARY: u16 = 0x100;
const FLAG_DUPLICATE: u16 = 0x400;
const FLAG_SUPPLEMENTARY: u16 = 0x800;

#[derive(Debug, clap::Args)]
pub struct Args {
    /// Input file [default: stdin].
    #[arg(default_value = "-")]
    input: String,

    /// Output file.
    #[arg(short = 'o', default_value = "-")]
    output: PathBuf,

    /// Tags to remove (`XX,YY`), or with `^` (`^XX,YY`) the only tags to keep; repeatable.
    #[arg(short = 'x', long = "remove-tag")]
    remove_tags: Vec<String>,

    /// Tags to keep (as `-x ^STR`); repeatable.
    #[arg(long = "keep-tag")]
    keep_tags: Vec<String>,

    /// Remove the @PG line with this ID and those after it.
    #[arg(long = "reject-PG")]
    reject_pg: Option<String>,

    /// Remove @RG lines and RG tags.
    #[arg(long = "no-RG")]
    no_rg: bool,

    /// Do not add an @PG line to the header.
    #[arg(long = "no-PG")]
    no_pg: bool,

    /// Keep duplicate flags.
    #[arg(long = "dupflag")]
    dupflag: bool,

    /// Output format [default: from the output extension, else SAM].
    #[arg(short = 'O', long = "output-fmt", ignore_case = true)]
    output_format: Option<OutputFormat>,

    /// Reference FASTA (indexed) for CRAM input or output.
    #[arg(short = 'T', long = "reference")]
    reference: Option<PathBuf>,

    /// Additional threads for compression and decompression.
    #[arg(short = '@', long = "threads", default_value_t = 0)]
    threads: usize,
}

/// Which tags records keep (`update_aux_conf`): a keep list takes priority over a remove list.
enum TagFilter {
    Keep(HashSet<[u8; 2]>),
    Remove(HashSet<[u8; 2]>),
}

impl TagFilter {
    fn keeps(&self, tag: [u8; 2]) -> bool {
        match self {
            Self::Keep(tags) => tags.contains(&tag),
            Self::Remove(tags) => !tags.contains(&tag),
        }
    }
}

pub fn run(args: Args) -> anyhow::Result<()> {
    let mut keep: Option<HashSet<[u8; 2]>> = None;
    let mut remove: Option<HashSet<[u8; 2]>> = None;

    for list in &args.remove_tags {
        match list.strip_prefix('^') {
            Some(list) => parse_tag_list(keep.get_or_insert_default(), list)?,
            None => parse_tag_list(remove.get_or_insert_default(), list)?,
        }
    }

    for list in &args.keep_tags {
        parse_tag_list(keep.get_or_insert_default(), list)?;
    }

    let filter = match keep {
        Some(mut keep) => {
            if args.no_rg {
                keep.remove(b"RG");
            }

            TagFilter::Keep(keep)
        }
        None => {
            let mut remove = remove.unwrap_or_default();

            if args.no_rg {
                remove.insert(*b"RG");
            }

            remove.extend(ALIGNMENT_TAGS.iter().map(|&&tag| tag));
            TagFilter::Remove(remove)
        }
    };

    let read_options = ReadOptions {
        threads: args.threads,
        reference: args.reference.clone(),
    };

    let (mut reader, _) = AlignmentReader::open(&args.input, &read_options)?;
    let input_header = reader.read_header()?;
    let mut header = output_header(&input_header, &args)?;

    if !args.no_pg {
        program::add(&mut header)?;
    }

    let format = args
        .output_format
        .or_else(|| OutputFormat::from_path(&args.output))
        .unwrap_or(OutputFormat::Sam);

    let write_options = WriteOptions {
        format,
        level: output::parse_level(6)?,
        threads: args.threads,
        reference: args.reference.clone(),
    };

    let mut writer = AlignmentWriter::create(Some(&args.output), &write_options, &header)?;
    writer.write_header(&header)?;

    let mut encoder = RawEncoder::new();
    let mut record = Vec::new();

    let duplicates = if args.dupflag {
        Duplicates::Keep
    } else {
        Duplicates::Clear
    };

    reader.for_each(&input_header, |rec| {
        let raw: &[u8] = rec.encode_raw(&input_header, &mut encoder)?;

        if reset(raw, &filter, duplicates, &mut record)? {
            writer.write_raw(&header, &record)?;
        }

        Ok(())
    })?;

    writer.finish(&header)?;

    Ok(())
}

/// The output header: `@HD VN:1.6`, the input's `@RG` lines (unless `--no-RG`), and its `@PG`
/// lines up to `--reject-PG`.
fn output_header(input: &sam::Header, args: &Args) -> anyhow::Result<sam::Header> {
    let mut text = String::from("@HD\tVN:1.6\n");
    let mut formatted = sam::io::Writer::new(Vec::new());
    formatted.write_header(input)?;
    let formatted = String::from_utf8(formatted.into_inner())?;

    if !args.no_rg {
        for line in formatted.lines().filter(|l| l.starts_with("@RG\t")) {
            text.push_str(line);
            text.push('\n');
        }
    }

    for line in formatted.lines().filter(|l| l.starts_with("@PG\t")) {
        let id = line
            .split('\t')
            .find_map(|field| field.strip_prefix("ID:"))
            .unwrap_or_default();

        if args
            .reject_pg
            .as_deref()
            .is_some_and(|reject| !reject.is_empty() && reject == id)
        {
            break;
        }

        text.push_str(line);
        text.push('\n');
    }

    Ok(text.parse()?)
}

/// Resets one record into `dst`; returns `false` for records that are dropped (secondary and
/// supplementary).
/// What happens to the duplicate flag (`--dupflag` keeps it).
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
enum Duplicates {
    Clear,
    Keep,
}

fn reset(
    src: &[u8],
    filter: &TagFilter,
    duplicates: Duplicates,
    dst: &mut Vec<u8>,
) -> anyhow::Result<bool> {
    let flags = u16::from_le_bytes([src[14], src[15]]);

    if flags & (FLAG_SECONDARY | FLAG_SUPPLEMENTARY) != 0 {
        return Ok(false);
    }

    let mut new_flags = (flags & !FLAG_PROPER_PAIR) | FLAG_UNMAPPED;

    if flags & FLAG_PAIRED != 0 {
        new_flags |= FLAG_MATE_UNMAPPED;
    }

    new_flags &= !FLAG_MATE_REVERSE;

    if duplicates == Duplicates::Clear {
        new_flags &= !FLAG_DUPLICATE;
    }

    let l_name = usize::from(src[8]);
    let n_cigar = usize::from(u16::from_le_bytes([src[12], src[13]]));
    let len = u32::from_le_bytes([src[16], src[17], src[18], src[19]]) as usize;
    let seq_start = 32 + l_name + 4 * n_cigar;
    let qual_start = seq_start + len.div_ceil(2);
    let aux_start = qual_start + len;

    let code = |i: usize| (src[seq_start + i / 2] >> (4 * (1 - i % 2))) & 0xf;
    let reverse = flags & FLAG_REVERSE != 0;

    if reverse {
        new_flags &= !FLAG_REVERSE;
    }

    dst.clear();
    dst.extend_from_slice(&(-1i32).to_le_bytes());
    dst.extend_from_slice(&(-1i32).to_le_bytes());
    dst.push(src[8]);
    dst.push(0);
    dst.extend_from_slice(&4680u16.to_le_bytes());
    dst.extend_from_slice(&0u16.to_le_bytes());
    dst.extend_from_slice(&new_flags.to_le_bytes());
    dst.extend_from_slice(&(len as u32).to_le_bytes());
    dst.extend_from_slice(&(-1i32).to_le_bytes());
    dst.extend_from_slice(&(-1i32).to_le_bytes());
    dst.extend_from_slice(&0i32.to_le_bytes());
    dst.extend_from_slice(&src[32..32 + l_name]);

    let bases: Vec<u8> = if reverse {
        (0..len)
            .rev()
            .map(|i| COMPLEMENT[usize::from(code(i))])
            .collect()
    } else {
        (0..len).map(code).collect()
    };

    for pair in bases.chunks(2) {
        dst.push((pair[0] << 4) | pair.get(1).copied().unwrap_or(0));
    }

    let qualities = &src[qual_start..aux_start];

    if reverse {
        dst.extend(qualities.iter().rev());
    } else {
        dst.extend_from_slice(qualities);
    }

    let mut aux = Vec::with_capacity(src.len() - aux_start);
    crate::bam_raw::strip_fields(src, &mut aux, |tag| !filter.keeps(tag))?;
    dst.extend_from_slice(&aux[aux_start.min(aux.len())..]);

    Ok(true)
}

/// `parse_aux_list`: two-character tags separated by commas.
fn parse_tag_list(tags: &mut HashSet<[u8; 2]>, list: &str) -> anyhow::Result<()> {
    let mut b = list.as_bytes();

    while b.len() >= 2 {
        tags.insert([b[0], b[1]]);
        b = &b[2..];

        match b.first() {
            Some(b',') => b = &b[1..],
            Some(_) => break,
            None => {}
        }
    }

    if !b.is_empty() {
        bail!(
            "[reset] Error parsing option, auxiliary tags should be exactly two characters long."
        );
    }

    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_parse_tag_list() -> anyhow::Result<()> {
        let mut tags = HashSet::new();
        parse_tag_list(&mut tags, "XA,YB")?;
        assert!(tags.contains(b"XA") && tags.contains(b"YB"));
        assert!(parse_tag_list(&mut tags, "XAYB").is_err());
        assert!(parse_tag_list(&mut tags, "X").is_err());
        Ok(())
    }
}
