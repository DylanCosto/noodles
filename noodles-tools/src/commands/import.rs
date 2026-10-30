//! `ntools import`: convert FASTQ (or FASTA) to unaligned SAM/BAM/CRAM (samtools import).
//!
//! This follows `import_fastq` in samtools 1.24 and HTSlib's FASTQ reader (`fastq_parse1`):
//! read names lose `/1` and `/2` suffixes (setting the read flags), `-i` parses CASAVA
//! comments, `-U` moves Illumina UMIs from names to a tag, and `-T` parses SAM-style tags from
//! comments with HTSlib's lenient rules. Files given with `-0`/`-1`/`-2`/`--i1`/`--i2`/`-s` are
//! read in step, index reads becoming `BC`/`QT` tags of the reads after them. Unlike HTSlib,
//! tags after a CASAVA barcode are parsed.

use std::{
    collections::HashSet,
    io::{self, BufRead, BufReader},
    path::PathBuf,
};

use anyhow::{Context, bail};
use flate2::bufread::MultiGzDecoder;
use noodles_pileup::baq::SEQ_NT16_TABLE;
use noodles_sam as sam;
use regex::bytes::Regex;

use crate::{
    bam_raw::{self, BASES},
    output::{self, AlignmentWriter, OutputFormat, WriteOptions},
    source::{Compression, Input},
};

const FLAG_PAIRED: u16 = 0x1;
const FLAG_UNMAPPED: u16 = 0x4;
const FLAG_MATE_UNMAPPED: u16 = 0x8;
const FLAG_READ1: u16 = 0x40;
const FLAG_READ2: u16 = 0x80;
const FLAG_QC_FAIL: u16 = 0x200;

/// HTSlib's default Illumina read name, whose eighth field is the UMI.
const ILLUMINA_UMI: &str = "^[^:]+:[^:]+:[^:]+:[^:]+:[^:]+:[^:]+:[^:]+:([^:#/]+)";

#[derive(Debug, clap::Args)]
pub struct Args {
    /// FASTQ files: one is read as -s, two as -1 and -2.
    inputs: Vec<String>,

    /// Read paired-end data from a single file.
    #[arg(short = 's')]
    single: Option<String>,

    /// Read single-end data from this file.
    #[arg(short = '0')]
    unpaired: Option<String>,

    /// Read 1 from this file.
    #[arg(short = '1', long = "r1")]
    read1: Option<String>,

    /// Read 2 from this file.
    #[arg(short = '2', long = "r2")]
    read2: Option<String>,

    /// Index 1 from this file.
    #[arg(long = "i1")]
    index1: Option<String>,

    /// Index 2 from this file.
    #[arg(long = "i2")]
    index2: Option<String>,

    /// Add index sequences to read 2 too, not just read 1.
    #[arg(short = 'b')]
    index_both: bool,

    /// Parse CASAVA identifiers (read number, QC failure, and barcode) from comments.
    #[arg(short = 'i')]
    casava: bool,

    /// Parse the UMI from read names.
    #[arg(short = 'U', long = "umi")]
    umi: bool,

    /// Tag for UMI sequences.
    #[arg(long = "UMI-tag", default_value = "RX")]
    umi_tag: String,

    /// Tag for barcode sequences.
    #[arg(long = "barcode-tag", default_value = "BC")]
    barcode_tag: String,

    /// Tag for barcode qualities.
    #[arg(long = "quality-tag", default_value = "QT")]
    quality_tag: String,

    /// Use the second field of names as the read name (SRA format).
    #[arg(short = 'N', long = "name2")]
    name2: bool,

    /// Build up an @RG line (repeatable; fields are joined with tabs).
    #[arg(short = 'r', long = "rg-line")]
    rg_line: Vec<String>,

    /// Add an @RG line with this ID.
    #[arg(short = 'R', long = "rg")]
    rg: Option<String>,

    /// Parse SAM-style tags from comments: a comma-separated list, or `*` for all.
    #[arg(short = 'T')]
    tags: Option<String>,

    /// Output file.
    #[arg(short = 'o', default_value = "-")]
    output: PathBuf,

    /// Output format [default: from the output extension, else SAM].
    #[arg(short = 'O', long = "output-fmt", ignore_case = true)]
    output_format: Option<OutputFormat>,

    /// Uncompressed output.
    #[arg(short = 'u')]
    uncompressed: bool,

    /// Store each record's number in this tag (`TAG:N` for a zero-padded string of width N).
    #[arg(long = "order")]
    order: Option<String>,

    /// Reference FASTA (indexed) for CRAM output.
    #[arg(long = "reference")]
    reference: Option<PathBuf>,

    /// Additional threads for compression.
    #[arg(short = '@', long = "threads", default_value_t = 0)]
    threads: usize,

    /// Do not add an @PG line to the header.
    #[arg(long = "no-PG")]
    no_pg: bool,
}

/// The input files, in the order samtools reads them (index reads before the reads they
/// annotate).
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
enum Slot {
    Index1,
    Index2,
    Unpaired,
    Read1,
    Read2,
    Single,
}

/// How the FASTQ reader parses names and comments (HTSlib's `fastq_state` options).
struct ParseOptions {
    casava: bool,
    name2: bool,
    umi_tag: Option<[u8; 2]>,
    umi_regex: Regex,
    barcode_tag: [u8; 2],
    /// `-T`: `Some(None)` parses every tag, `Some(Some(set))` the listed ones.
    aux: Option<Option<HashSet<[u8; 2]>>>,
}

pub fn run(args: Args) -> anyhow::Result<()> {
    let mut files: Vec<(Slot, String)> = Vec::new();
    let mut single = args.single.clone();
    let (mut read1, mut read2) = (args.read1.clone(), args.read2.clone());

    match args.inputs.as_slice() {
        [] => {}
        [only] => single = Some(only.clone()),
        [first, second, ..] => {
            read1 = Some(first.clone());
            read2 = Some(second.clone());
        }
    }

    for (slot, path) in [
        (Slot::Index1, &args.index1),
        (Slot::Index2, &args.index2),
        (Slot::Unpaired, &args.unpaired),
        (Slot::Read1, &read1),
        (Slot::Read2, &read2),
        (Slot::Single, &single),
    ] {
        if let Some(path) = path {
            files.push((slot, path.clone()));
        }
    }

    if files.is_empty() {
        eprintln!("Usage: ntools import [options] [file.fastq ...]");
        return Ok(());
    }

    let options = ParseOptions {
        casava: args.casava,
        name2: args.name2,
        umi_tag: args.umi.then(|| parse_umi_tags(&args.umi_tag)).flatten(),
        umi_regex: Regex::new(ILLUMINA_UMI).expect("valid regex"),
        barcode_tag: tag(&args.barcode_tag),
        aux: args.tags.as_deref().map(parse_tag_list),
    };

    let mut readers = files
        .iter()
        .map(|(slot, path)| FastqReader::open(path).map(|reader| (*slot, reader)))
        .collect::<anyhow::Result<Vec<_>>>()?;

    let (text, rg) = header_text(&args, &files)?;
    let header: sam::Header = String::from_utf8_lossy(&text)
        .parse()
        .context("invalid header")?;

    let format = args
        .output_format
        .or_else(|| OutputFormat::from_path(&args.output))
        .unwrap_or(OutputFormat::Sam);

    let write_options = WriteOptions {
        format,
        level: output::parse_level(if args.uncompressed { 0 } else { 6 })?,
        threads: args.threads,
        reference: args.reference.clone(),
    };

    let mut writer = AlignmentWriter::create(Some(&args.output), &write_options, &header)?;
    writer.write_header_text(&header, &text)?;

    let order = args.order.as_deref().map(|order| {
        let width = (order.len() > 3 && order.as_bytes()[2] == b':')
            .then(|| atoi(&order[3..]))
            .filter(|&n| n != 0);
        (tag(order), width)
    });

    let index_tags = (tag(&args.barcode_tag), tag(&args.quality_tag));

    let n = readers.len();
    let mut record = Vec::new();
    let mut index_bases = Vec::new();
    let mut index_qualities = Vec::new();
    let mut read_number: u64 = 0;
    let mut eof = 0;

    // samtools reads one record from each file in turn until the last file ends; files that
    // end earlier are counted at each turn.
    loop {
        index_bases.clear();
        index_qualities.clear();
        let mut last = Read::Record;

        for i in 0..n {
            let slot = readers[i].0;

            last = readers[i].1.read(&options, &mut record)?;

            match last {
                Read::Eof => {
                    eof += 1;
                    continue;
                }
                Read::Truncated => break,
                Read::Record => {}
            }

            if matches!(slot, Slot::Index1 | Slot::Index2) {
                append_index(&mut index_bases, &mut index_qualities, &record);
                continue;
            }

            if !index_bases.is_empty()
                && (args.index_both || matches!(slot, Slot::Single | Slot::Unpaired | Slot::Read1))
            {
                append_string_tag(&mut record, index_tags.0, &index_bases);
                append_string_tag(&mut record, index_tags.1, &index_qualities);
            }

            let mut flags = record_flags(&record);

            match slot {
                Slot::Read1 => {
                    if flags & (FLAG_READ1 | FLAG_READ2) == 0 {
                        flags |= FLAG_READ1;
                    }

                    flags |= FLAG_PAIRED;

                    if i + 1 < n && readers[i + 1].0 == Slot::Read2 {
                        flags |= FLAG_MATE_UNMAPPED;
                    }
                }
                Slot::Read2 => {
                    flags |= FLAG_PAIRED | FLAG_READ2;

                    if i > 0 && readers[i - 1].0 == Slot::Read1 {
                        flags |= FLAG_MATE_UNMAPPED;
                    }
                }
                _ => {}
            }

            record[14..16].copy_from_slice(&flags.to_le_bytes());

            if let Some(rg) = &rg {
                append_string_tag(&mut record, *b"RG", rg);
            }

            if let Some((order_tag, width)) = order {
                match width {
                    Some(width) => {
                        let value = format!("{read_number:0width$}");
                        bam_raw::update_string_field(&mut record, order_tag, value.as_bytes())?;
                    }
                    None => update_integer_tag(&mut record, order_tag, read_number)?,
                }

                read_number += 1;
            }

            writer.write_raw(&header, &record)?;
        }

        if last != Read::Record {
            if last == Read::Truncated {
                writer.finish(&header)?;
                bail!("[import] truncated file. Aborting");
            }

            break;
        }
    }

    writer.finish(&header)?;

    if eof != n {
        bail!("[import] input files with differing number of records");
    }

    Ok(())
}

/// The header samtools builds: `@HD`, a `@CO` line saying how to convert back, `@PG`, and
/// `@RG`, in that order. Also returns the read group ID to tag records with.
fn header_text(
    args: &Args,
    files: &[(Slot, String)],
) -> anyhow::Result<(Vec<u8>, Option<Vec<u8>>)> {
    let mut reverse = String::new();
    let mut index_format = String::new();

    for (slot, _) in files {
        reverse.push_str(match slot {
            Slot::Index1 => " --i1 I1.fastq",
            Slot::Index2 => " --i2 I2.fastq",
            Slot::Unpaired => " -0 unpaired.fastq",
            Slot::Read1 => " -1 R1.fastq",
            Slot::Read2 => " -2 R2.fastq",
            Slot::Single => " -n -o paired.fastq",
        });

        if matches!(slot, Slot::Index1 | Slot::Index2) {
            index_format.push_str("i*");
        }
    }

    if args.casava {
        reverse.push_str(" -i");

        if index_format.is_empty() {
            reverse.push_str(" --index-format 'i*i*'");
        }
    }

    if args.umi {
        reverse.push_str(&format!(" -U --UMI-tag {}", args.umi_tag));
    }

    let mut text = b"@HD\tVN:1.6\tSO:unsorted\tGO:query\n".to_vec();

    // samtools formats this into a 2100-byte buffer.
    let mut comment = if index_format.is_empty() {
        format!("@CO\tReverse with: samtools fastq{reverse}\n")
    } else {
        format!("@CO\tReverse with: samtools fastq{reverse} --index-format=\"{index_format}\"\n")
    };

    if comment.len() > 2099 {
        comment.truncate(2099);
    }

    text.extend_from_slice(comment.as_bytes());

    if !args.no_pg {
        let command_line = std::env::args().collect::<Vec<_>>().join(" ");
        text.extend_from_slice(
            format!(
                "@PG\tID:ntools\tPN:ntools\tVN:{}\tCL:{command_line}\n",
                env!("CARGO_PKG_VERSION")
            )
            .as_bytes(),
        );
    }

    let mut rg_line = String::new();

    if !args.rg_line.is_empty() {
        for (i, part) in args.rg_line.iter().enumerate() {
            if i == 0 && !part.starts_with('@') {
                rg_line.push_str("@RG");
            }

            if !rg_line.is_empty() {
                rg_line.push('\t');
            }

            rg_line.push_str(part);
        }
    } else if let Some(rg) = &args.rg {
        rg_line = format!("@RG\tID:{rg}");
    }

    let mut rg = None;

    if !rg_line.is_empty() {
        let Some(start) = rg_line.find("\tID:") else {
            bail!("\"-r RG-LINE\" option contained no ID field");
        };

        let id = &rg_line[start + 4..];
        let id = &id[..id.find('\t').unwrap_or(id.len())];
        rg = Some(id.as_bytes().to_vec());

        text.extend_from_slice(rg_line.as_bytes());
        text.push(b'\n');
    }

    Ok((text, rg))
}

/// The result of reading one record.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
enum Read {
    Record,
    Eof,
    /// A malformed or truncated record (HTSlib's -2).
    Truncated,
}

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
enum SequenceFormat {
    Fastq,
    Fasta,
}

/// HTSlib's FASTQ/FASTA record reader.
struct FastqReader {
    reader: Box<dyn BufRead>,
    format: SequenceFormat,
    /// For FASTA, the name line of the next record, read while reading the sequence.
    pending_name: Option<Vec<u8>>,
    line: Vec<u8>,
    name: Vec<u8>,
    sequence: Vec<u8>,
    qualities: Vec<u8>,
}

impl FastqReader {
    fn open(src: &str) -> anyhow::Result<Self> {
        let input = Input::open(src)?;

        let mut reader: Box<dyn BufRead> = match input.compression {
            Compression::Bgzf => Box::new(BufReader::new(MultiGzDecoder::new(input.reader))),
            Compression::None => Box::new(input.reader),
        };

        let format = match reader.fill_buf()?.first() {
            Some(b'>') => SequenceFormat::Fasta,
            Some(b'@') | None => SequenceFormat::Fastq,
            Some(_) => bail!("{src}: not a FASTQ or FASTA file"),
        };

        Ok(Self {
            reader,
            format,
            pending_name: None,
            line: Vec::new(),
            name: Vec::new(),
            sequence: Vec::new(),
            qualities: Vec::new(),
        })
    }

    /// `hts_getline`: the next line without its line ending; `false` at the end of input.
    fn next_line(reader: &mut dyn BufRead, dst: &mut Vec<u8>) -> io::Result<bool> {
        dst.clear();

        if reader.read_until(b'\n', dst)? == 0 {
            return Ok(false);
        }

        if dst.last() == Some(&b'\n') {
            dst.pop();

            if dst.last() == Some(&b'\r') {
                dst.pop();
            }
        }

        Ok(true)
    }

    /// `fastq_parse1`: reads the next record as an encoded BAM record into `dst`.
    fn read(&mut self, options: &ParseOptions, dst: &mut Vec<u8>) -> io::Result<Read> {
        let fasta = self.format == SequenceFormat::Fasta;

        if let (true, Some(name)) = (fasta, self.pending_name.take()) {
            self.name = name;
        } else if !Self::next_line(&mut self.reader, &mut self.name)? {
            return Ok(Read::Eof);
        }

        let prefix = if fasta { b'>' } else { b'@' };

        if self.name.first() != Some(&prefix) {
            return Ok(Read::Truncated);
        }

        // The name runs from `name_start` to the first whitespace; the comment follows the
        // whitespace after it.
        let mut name_start = 1;
        let mut i = 0;

        if options.name2
            && let Some(cp) = self.name.iter().position(|&b| b == b' ' || b == b'\t')
        {
            let mut cp = cp;

            while matches!(self.name.get(cp), Some(b' ' | b'\t')) {
                cp += 1;
            }

            // HTSlib writes a fake '@' before the second field.
            cp -= 1;
            self.name[cp] = b'@';
            i = cp;
            name_start = cp + 1;
        }

        let l = self.name.len();

        while i < l && !is_space(self.name[i]) {
            i += 1;
        }

        let mut name_end = i;

        if i < l {
            i += 1;
        }

        while i < l && is_space(self.name[i]) {
            i += 1;
        }

        // The comment, NUL terminated, for HTSlib's string scanning.
        let mut comment = self.name[i..].to_vec();
        let comment_len = comment.len();
        comment.push(0);

        // Sequence lines, up to the '+' line (or the next '>' line or end, for FASTA).
        self.sequence.clear();

        loop {
            let more = Self::next_line(&mut self.reader, &mut self.line)?;

            if !more {
                if !fasta {
                    return Ok(Read::Truncated);
                }

                break;
            }

            if self.line.first() == Some(if fasta { &b'>' } else { &b'+' }) {
                if fasta {
                    self.pending_name = Some(std::mem::take(&mut self.line));
                }

                break;
            }

            self.sequence.extend_from_slice(&self.line);
        }

        self.qualities.clear();

        if !fasta {
            let mut remainder = self.sequence.len();

            loop {
                if !Self::next_line(&mut self.reader, &mut self.line)? {
                    return Ok(Read::Truncated);
                }

                if self.line.len() > remainder {
                    return Ok(Read::Truncated);
                }

                self.qualities
                    .extend(self.line.iter().map(|q| q.wrapping_sub(b'!')));
                remainder -= self.line.len();

                if remainder == 0 {
                    break;
                }
            }
        }

        let pair_flags = FLAG_MATE_UNMAPPED | FLAG_PAIRED;
        let mut flags = FLAG_UNMAPPED;

        if name_end > 2
            && self.name[name_end - 2] == b'/'
            && self.name[name_end - 1].is_ascii_digit()
        {
            flags |= match self.name[name_end - 1] {
                b'1' => FLAG_READ1 | pair_flags,
                b'2' => FLAG_READ2 | pair_flags,
                _ => FLAG_READ1 | FLAG_READ2 | pair_flags,
            };

            name_end -= 2;
        }

        // Strip an Illumina UMI off the name.
        let mut umi = None;
        let mut name: Vec<u8> = self.name[name_start..name_end].to_vec();

        if let Some(umi_tag) = options.umi_tag
            && let Some(m) = options
                .umi_regex
                .captures(&self.name[..name_end])
                .and_then(|c| c.get(1))
        {
            if m.len() > 255 {
                return Ok(Read::Truncated);
            }

            if !m.is_empty() {
                let bases: Vec<u8> = m
                    .as_bytes()
                    .iter()
                    .map(|&b| if b.is_ascii_alphabetic() { b } else { b'-' })
                    .collect();

                let mut end = m.start();

                if end > 0 && self.name[end - 1] == b':' {
                    end -= 1;
                }

                let mut stripped = self.name[..end].to_vec();
                stripped.extend_from_slice(&self.name[m.end()..name_end]);
                name = stripped[name_start.min(stripped.len())..].to_vec();
                umi = Some((umi_tag, bases));
            }
        }

        if name.len() > 254 {
            return Ok(Read::Truncated);
        }

        let qualities = (!fasta).then_some(self.qualities.as_slice());
        encode_unmapped(dst, &name, flags, &self.sequence, qualities);

        if let Some((umi_tag, bases)) = umi {
            append_string_tag(dst, umi_tag, &bases);
        }

        // Tags are parsed from after a CASAVA barcode. (HTSlib starts at the barcode's length
        // into the comment, inside the CASAVA field, so tags after a barcode are usually lost.)
        let mut aux_start = 0;

        if options.casava
            && let Some((barcode, casava_flags)) =
                parse_casava(&mut comment, comment_len, pair_flags)
        {
            let flags = record_flags(dst) | casava_flags;
            dst[14..16].copy_from_slice(&flags.to_le_bytes());

            if let Some(barcode) = barcode {
                aux_start = (barcode.end + 1).min(comment_len);
                append_string_tag(dst, options.barcode_tag, &comment[barcode]);
            }
        }

        if let Some(whitelist) = &options.aux
            && parse_aux(&comment, aux_start, comment_len, dst, whitelist.as_ref()).is_err()
        {
            return Ok(Read::Truncated);
        }

        Ok(Read::Record)
    }
}

/// HTSlib's `isspace_c`.
fn is_space(b: u8) -> bool {
    matches!(b, b' ' | b'\t' | b'\n' | 0x0b | 0x0c | b'\r')
}

/// The byte at `i` of a NUL-terminated buffer, or NUL past its end.
fn at(buf: &[u8], i: usize) -> u8 {
    buf.get(i).copied().unwrap_or(0)
}

/// Parses a CASAVA comment (`\d:[YN]:\d+:[ACGTN]+`), returning the barcode's range, if it
/// has one (a numeric one is skipped), and the flags to add. The barcode is NUL terminated in
/// place.
fn parse_casava(
    comment: &mut [u8],
    len: usize,
    pair_flags: u16,
) -> Option<(Option<std::ops::Range<usize>>, u16)> {
    // HTSlib tests `(s[1] | s[3]) == ':'`.
    if len <= 6 || (comment[1] | comment[3]) != b':' || !comment[0].is_ascii_digit() {
        return None;
    }

    let (n, end) = strtol(comment, 4);

    if n < 0 || end == 4 || at(comment, end) != b':' {
        return None;
    }

    let mut flags = match comment[0] {
        b'1' => FLAG_READ1 | pair_flags,
        b'2' => FLAG_READ2 | pair_flags,
        _ => FLAG_READ1 | FLAG_READ2 | pair_flags,
    };

    if comment[2] == b'Y' {
        flags |= FLAG_QC_FAIL;
    }

    let mut barcode = None;

    if !at(comment, end + 1).is_ascii_digit() {
        let start = end + 1;
        let mut i = start;

        while i < len && !is_space(comment[i]) {
            i += 1;
        }

        comment[i] = 0;
        barcode = Some(start..i);
    }

    Some((barcode, flags))
}

/// C's `strtol(s + start, &end, 10)`: the value and the end, which is `start` if no digits.
fn strtol(s: &[u8], start: usize) -> (i64, usize) {
    let mut i = start;

    while is_space(at(s, i)) {
        i += 1;
    }

    let negative = at(s, i) == b'-';

    if matches!(at(s, i), b'-' | b'+') {
        i += 1;
    }

    let digits = i;
    let mut n: i64 = 0;

    while at(s, i).is_ascii_digit() {
        n = n
            .saturating_mul(10)
            .saturating_add(i64::from(at(s, i) - b'0'));
        i += 1;
    }

    if i == digits {
        return (0, start);
    }

    (if negative { -n } else { n }, i)
}

/// `aux_parse` in lenient mode: SAM-style tags in a comment, from `start` to `end`. Malformed
/// tags are skipped; malformed `B` arrays are errors.
fn parse_aux(
    buf: &[u8],
    start: usize,
    end: usize,
    dst: &mut Vec<u8>,
    whitelist: Option<&HashSet<[u8; 2]>>,
) -> Result<(), ()> {
    let p = end as isize;
    let mut q = start as isize;

    let skip_token = |q: &mut isize| {
        while *q < p && !is_space(at(buf, *q as usize)) {
            *q += 1;
        }

        while *q < p && is_space(at(buf, *q as usize)) {
            *q += 1;
        }
    };

    while q < p {
        let checkpoint = dst.len();

        if p - q < 5 {
            break;
        }

        let c = |i: isize| at(buf, i as usize);

        if c(q) < b'!' || c(q + 1) < b'!' {
            skip_token(&mut q);
            continue;
        }

        if (c(q + 2) | c(q + 4)) != b':' {
            skip_token(&mut q);
            continue;
        }

        if let Some(whitelist) = whitelist
            && !whitelist.contains(&[c(q), c(q + 1)])
        {
            while q < p && c(q) != b'\t' {
                q += 1;
            }

            continue;
        }

        dst.extend_from_slice(&[c(q), c(q + 1)]);
        let ty = c(q + 3);
        q += 5;

        if ty != b'Z' && ty != b'H' && c(q) <= b'\t' {
            skip_token(&mut q);
            dst.truncate(checkpoint);
            continue;
        }

        let mut i = q as usize;

        match ty {
            b'A' | b'a' | b'c' | b'C' => {
                dst.push(b'A');
                dst.push(at(buf, i));
                i += 1;
            }
            b'i' | b'I' => {
                if at(buf, i) == b'-' {
                    let (n, end, _) = str2int(buf, i, 32);
                    i = end;
                    let n = n as i32;

                    if let Ok(n) = i8::try_from(n) {
                        dst.push(b'c');
                        dst.push(n as u8);
                    } else if let Ok(n) = i16::try_from(n) {
                        dst.push(b's');
                        dst.extend_from_slice(&n.to_le_bytes());
                    } else {
                        dst.push(b'i');
                        dst.extend_from_slice(&n.to_le_bytes());
                    }
                } else {
                    let (n, end, _) = str2uint(buf, i, 32);
                    i = end;
                    let n = n as u32;

                    if let Ok(n) = u8::try_from(n) {
                        dst.push(b'C');
                        dst.push(n);
                    } else if let Ok(n) = u16::try_from(n) {
                        dst.push(b'S');
                        dst.extend_from_slice(&n.to_le_bytes());
                    } else {
                        dst.push(b'I');
                        dst.extend_from_slice(&n.to_le_bytes());
                    }
                }
            }
            b'f' => {
                let (x, end) = strtod(buf, i);
                i = end;
                dst.push(b'f');
                dst.extend_from_slice(&(x as f32).to_le_bytes());
            }
            b'd' => {
                let (x, end) = strtod(buf, i);
                i = end;
                dst.push(b'd');
                dst.extend_from_slice(&x.to_le_bytes());
            }
            b'Z' | b'H' => {
                let value_end = buf[i.min(buf.len())..]
                    .iter()
                    .position(|&b| b == b'\t' || b == 0)
                    .map_or(buf.len(), |n| i + n);

                if ty == b'H' && !(value_end - i).is_multiple_of(2) {
                    skip_token(&mut q);
                    dst.truncate(checkpoint);
                    continue;
                }

                dst.push(ty);
                dst.extend_from_slice(&buf[i..value_end]);
                dst.push(0);
                i = value_end;
            }
            b'B' => {
                let subtype = at(buf, i);
                i += 1;

                if !matches!(at(buf, i), 0 | b',' | b'\t') {
                    skip_token(&mut q);
                    dst.truncate(checkpoint);
                    continue;
                }

                let mut retries = 0;
                i = parse_b_values(subtype, buf, i, dst, &mut retries)?;
            }
            _ => {
                skip_token(&mut q);
                dst.truncate(checkpoint);
                continue;
            }
        }

        // Skip any junk to the next tab.
        while at(buf, i) > b'\t' {
            i += 1;
        }

        q = i as isize + 1;
    }

    Ok(())
}

/// `sam_parse_B_vals_r`: a `B` array's values from the comma before the first. Values out of
/// the given type's range retype the array to the smallest type that holds them.
fn parse_b_values(
    subtype: u8,
    buf: &[u8],
    start: usize,
    dst: &mut Vec<u8>,
    retries: &mut u32,
) -> Result<usize, ()> {
    *retries += 1;

    if *retries > 2 {
        return Err(());
    }

    let (bits, signed) = match subtype {
        b'A' => (0, false),
        b'c' => (8, true),
        b'C' => (8, false),
        b's' => (16, true),
        b'S' => (16, false),
        b'i' => (32, true),
        b'I' => (32, false),
        b'f' => (32, true),
        _ => return Err(()),
    };

    let orig = dst.len();
    dst.extend_from_slice(&[b'B', subtype, 0, 0, 0, 0]);

    let mut q = start;
    let mut n: u32 = 0;
    let mut overflow = false;

    let skip_to_comma = |q: &mut usize| {
        while at(buf, *q) > b'\t' && at(buf, *q) != b',' {
            *q += 1;
        }
    };

    if subtype != b'A' {
        while at(buf, q) == b',' {
            n += 1;

            if subtype == b'f' {
                let (x, end) = strtod(buf, q + 1);
                q = end;
                dst.extend_from_slice(&(x as f32).to_le_bytes());
            } else if signed {
                let (x, end, failed) = str2int(buf, q + 1, bits);
                q = end;
                overflow |= failed;
                dst.extend_from_slice(&x.to_le_bytes()[..bits as usize / 8]);
            } else if at(buf, q + 1) != b'-' {
                let (x, end, failed) = str2uint(buf, q + 1, bits);
                q = end;
                overflow |= failed;
                dst.extend_from_slice(&x.to_le_bytes()[..bits as usize / 8]);
            } else {
                overflow = true;
                q += 1;
                skip_to_comma(&mut q);
            }
        }
    }

    if !matches!(at(buf, q), b'\t' | 0) {
        return Err(());
    }

    dst[orig + 2..orig + 6].copy_from_slice(&n.to_le_bytes());

    if !overflow {
        return Ok(q);
    }

    // Find the range of values, and retry with a type that holds it.
    let end = q;
    let mut q = start;
    let (mut min, mut max) = (0i64, 0i64);
    let mut failed = false;
    dst.truncate(orig);

    while q < end {
        let (x, next, f) = str2int(buf, q + 1, 64);
        failed |= f;
        q = next;
        max = max.max(x);
        min = min.min(x);
        skip_to_comma(&mut q);
    }

    if !failed {
        let retype = if min < 0 {
            if min >= i64::from(i8::MIN) && max <= i64::from(i8::MAX) {
                Some(b'c')
            } else if min >= i64::from(i16::MIN) && max <= i64::from(i16::MAX) {
                Some(b's')
            } else if min >= i64::from(i32::MIN) && max <= i64::from(i32::MAX) {
                Some(b'i')
            } else {
                None
            }
        } else if max < i64::from(u8::MAX) {
            // HTSlib's test is strict here.
            Some(b'C')
        } else if max <= i64::from(u16::MAX) {
            Some(b'S')
        } else if max <= i64::from(u32::MAX) {
            Some(b'I')
        } else {
            None
        };

        if let Some(subtype) = retype {
            return parse_b_values(subtype, buf, start, dst, retries);
        }
    }

    Err(())
}

/// `hts_str2int(s + start, &end, bits, &failed)`: an optionally signed decimal, saturating at
/// the type's limits. Returns the value, the end, and whether it overflowed.
fn str2int(s: &[u8], start: usize, bits: u32) -> (i64, usize, bool) {
    let mut i = start;
    let mut limit: u64 = (1 << (bits - 1)) - 1;
    let negative = at(s, i) == b'-';

    if negative {
        limit += 1;
        i += 1;
    } else if at(s, i) == b'+' {
        i += 1;
    }

    let (n, end, failed) = digits(s, i, limit);
    let n = n as i64;

    (if negative { n.wrapping_neg() } else { n }, end, failed)
}

/// `hts_str2uint(s + start, &end, bits, &failed)`
fn str2uint(s: &[u8], start: usize, bits: u32) -> (u64, usize, bool) {
    let mut i = start;
    let limit: u64 = if bits < 64 { (1 << bits) - 1 } else { u64::MAX };

    if at(s, i) == b'+' {
        i += 1;
    }

    digits(s, i, limit)
}

/// Decimal digits from `start`, saturating at `limit`.
fn digits(s: &[u8], start: usize, limit: u64) -> (u64, usize, bool) {
    let mut i = start;
    let mut n: u64 = 0;

    while at(s, i).is_ascii_digit() {
        let d = u64::from(at(s, i) - b'0');

        match n.checked_mul(10).and_then(|n| n.checked_add(d)) {
            Some(m) if m <= limit => {
                n = m;
                i += 1;
            }
            _ => {
                while at(s, i).is_ascii_digit() {
                    i += 1;
                }

                return (limit, i, true);
            }
        }
    }

    (n, i, false)
}

/// C's `strtod(s + start, &end)`: the longest prefix that is a decimal float (or `inf`/`nan`).
fn strtod(s: &[u8], start: usize) -> (f64, usize) {
    let mut i = start;

    while is_space(at(s, i)) {
        i += 1;
    }

    let number_start = i;

    if matches!(at(s, i), b'+' | b'-') {
        i += 1;
    }

    let rest = &s[i.min(s.len())..];

    for word in ["infinity", "inf", "nan"] {
        if rest.len() >= word.len() && rest[..word.len()].eq_ignore_ascii_case(word.as_bytes()) {
            let end = i + word.len();
            let text = String::from_utf8_lossy(&s[number_start..end]);
            return (text.parse().unwrap_or(0.0), end);
        }
    }

    let mut mantissa_digits = 0;

    while at(s, i).is_ascii_digit() {
        i += 1;
        mantissa_digits += 1;
    }

    if at(s, i) == b'.' {
        i += 1;

        while at(s, i).is_ascii_digit() {
            i += 1;
            mantissa_digits += 1;
        }
    }

    if mantissa_digits == 0 {
        return (0.0, start);
    }

    if matches!(at(s, i), b'e' | b'E') {
        let mut j = i + 1;

        if matches!(at(s, j), b'+' | b'-') {
            j += 1;
        }

        if at(s, j).is_ascii_digit() {
            while at(s, j).is_ascii_digit() {
                j += 1;
            }

            i = j;
        }
    }

    let text = String::from_utf8_lossy(&s[number_start..i]);
    (text.parse().unwrap_or(0.0), i)
}

/// `bam_set1` for an unmapped record without a position.
fn encode_unmapped(
    dst: &mut Vec<u8>,
    name: &[u8],
    flags: u16,
    sequence: &[u8],
    qualities: Option<&[u8]>,
) {
    const UNPLACED_BIN: u16 = 4680;

    dst.clear();
    dst.extend_from_slice(&(-1i32).to_le_bytes());
    dst.extend_from_slice(&(-1i32).to_le_bytes());
    dst.push(u8::try_from(name.len() + 1).expect("name length checked"));
    dst.push(0);
    dst.extend_from_slice(&UNPLACED_BIN.to_le_bytes());
    dst.extend_from_slice(&0u16.to_le_bytes());
    dst.extend_from_slice(&flags.to_le_bytes());
    dst.extend_from_slice(&(sequence.len() as u32).to_le_bytes());
    dst.extend_from_slice(&(-1i32).to_le_bytes());
    dst.extend_from_slice(&(-1i32).to_le_bytes());
    dst.extend_from_slice(&0i32.to_le_bytes());
    dst.extend_from_slice(name);
    dst.push(0);

    for pair in sequence.chunks(2) {
        let hi = SEQ_NT16_TABLE[usize::from(pair[0])];
        let lo = pair.get(1).map_or(0, |&b| SEQ_NT16_TABLE[usize::from(b)]);
        dst.push((hi << 4) | lo);
    }

    match qualities {
        Some(qualities) => dst.extend_from_slice(qualities),
        None => dst.resize(dst.len() + sequence.len(), 0xff),
    }
}

fn record_flags(record: &[u8]) -> u16 {
    u16::from_le_bytes([record[14], record[15]])
}

/// `append_index`: adds an index read's bases and qualities to the `BC`/`QT` values, separated
/// by `-` and ` `.
fn append_index(bases: &mut Vec<u8>, qualities: &mut Vec<u8>, record: &[u8]) {
    let l_name = usize::from(record[8]);
    let n_cigar = usize::from(u16::from_le_bytes([record[12], record[13]]));
    let len = u32::from_le_bytes([record[16], record[17], record[18], record[19]]) as usize;
    let seq_start = 32 + l_name + 4 * n_cigar;
    let qual_start = seq_start + len.div_ceil(2);

    if !bases.is_empty() {
        bases.push(b'-');
    }

    if !qualities.is_empty() {
        qualities.push(b' ');
    }

    for i in 0..len {
        let code = (record[seq_start + i / 2] >> (4 * (1 - i % 2))) & 0xf;
        bases.push(BASES[usize::from(code)]);
        qualities.push(record[qual_start + i].wrapping_add(b'!'));
    }
}

/// `bam_aux_append` of a `Z` value.
fn append_string_tag(record: &mut Vec<u8>, tag: [u8; 2], value: &[u8]) {
    record.extend_from_slice(&tag);
    record.push(b'Z');
    record.extend_from_slice(value);
    record.push(0);
}

/// `bam_aux_update_int`: sets an integer tag in the smallest unsigned type, in place if it
/// exists.
fn update_integer_tag(record: &mut Vec<u8>, tag: [u8; 2], n: u64) -> io::Result<()> {
    let mut field = tag.to_vec();

    if let Ok(n) = u8::try_from(n) {
        field.extend([b'C', n]);
    } else if let Ok(n) = u16::try_from(n) {
        field.push(b'S');
        field.extend(n.to_le_bytes());
    } else {
        field.push(b'I');
        field.extend((n as u32).to_le_bytes());
    }

    match bam_raw::find_field(record, tag)? {
        Some(range) => {
            record.splice(range, field);
        }
        None => record.extend_from_slice(&field),
    }

    Ok(())
}

/// A tag from the first two bytes of an option (NUL padded).
fn tag(s: &str) -> [u8; 2] {
    let b = s.as_bytes();
    [
        b.first().copied().unwrap_or(0),
        b.get(1).copied().unwrap_or(0),
    ]
}

/// `FASTQ_OPT_UMI`: the first of up to five comma-separated tags; `None` if the list is empty.
fn parse_umi_tags(s: &str) -> Option<[u8; 2]> {
    let b = s.as_bytes();

    if b.len() >= 2 && b[0].is_ascii_alphabetic() && b[1].is_ascii_alphanumeric() {
        Some([b[0], b[1]])
    } else {
        if !b.is_empty() {
            eprintln!("[W::fastq_state_set] Bad UMI tag list '{s}'");
        }

        None
    }
}

/// `FASTQ_OPT_AUX`: `*` (or empty) for every tag, else a list of tags (`XX,YY`), up to the
/// first malformed one.
fn parse_tag_list(s: &str) -> Option<HashSet<[u8; 2]>> {
    if s.is_empty() || s.starts_with('*') || s == "1" {
        return None;
    }

    let b = s.as_bytes();
    let mut tags = HashSet::new();
    let mut i = 0;

    while i + 3 <= b.len() + 1 {
        let next = b.get(i + 2).copied().unwrap_or(0);

        if b[i] == b',' || b.get(i + 1) == Some(&b',') || !(next == b',' || next == 0) {
            eprintln!(
                "[W::fastq_state_set] Bad tag format '{}'; skipping option",
                String::from_utf8_lossy(&b[i..(i + 3).min(b.len())])
            );
            break;
        }

        tags.insert([b[i], b.get(i + 1).copied().unwrap_or(0)]);
        i += 3;
    }

    Some(tags)
}

/// C's `atoi`.
fn atoi(s: &str) -> usize {
    let s = s.trim_start();
    let digits = s.find(|c: char| !c.is_ascii_digit()).unwrap_or(s.len());
    s[..digits].parse().unwrap_or(0)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn aux(comment: &str) -> Vec<u8> {
        let mut buf = comment.as_bytes().to_vec();
        let len = buf.len();
        buf.push(0);
        let mut dst = Vec::new();
        parse_aux(&buf, 0, len, &mut dst, None).expect("parsed");
        dst
    }

    #[test]
    fn test_parse_aux() {
        assert_eq!(aux("XX:Z:foo\tXI:i:-300"), b"XXZfoo\0XIs\xd4\xfe");
        assert_eq!(aux("1:N:0:ACGT\tXA:A:x"), b"XAAx");
        // Space-separated tags after a value are skipped up to the next tab.
        assert_eq!(aux("XA:i:1 XB:i:2"), b"XAC\x01");
        // An out-of-range B:C array is retyped.
        assert_eq!(aux("ML:B:C,1,300"), b"MLBS\x02\0\0\0\x01\0\x2c\x01");
    }

    #[test]
    fn test_parse_casava() {
        let mut comment = b"1:Y:0:ACGT rest\0".to_vec();
        let (barcode, flags) = parse_casava(&mut comment, 15, FLAG_PAIRED).expect("CASAVA");
        assert_eq!(barcode, Some(6..10));
        assert_eq!(flags, FLAG_READ1 | FLAG_PAIRED | FLAG_QC_FAIL);
        assert_eq!(&comment[6..11], b"ACGT\0");
    }

    #[test]
    fn test_number_parsing() {
        assert_eq!(str2int(b"-129", 0, 8), (-128, 4, true));
        assert_eq!(str2uint(b"256,", 0, 8), (255, 3, true));
        assert_eq!(strtod(b"1.5e2x", 0), (150.0, 5));
        assert_eq!(strtod(b"x", 0), (0.0, 0));
    }
}
