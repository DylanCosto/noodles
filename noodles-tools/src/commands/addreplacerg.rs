//! `ntools addreplacerg`: add or replace read group tags and header lines
//! (samtools addreplacerg).

use std::path::PathBuf;

use anyhow::{Context, bail};
use bstr::BString;
use noodles_bam as bam;
use noodles_sam::{self as sam, alignment::io::Write as _};

use crate::{
    alignment::{AlignmentReader, ReadOptions, Rec},
    bam_raw,
    output::{self, AlignmentWriter, OutputFormat, WriteOptions},
    program,
};

#[derive(Clone, Copy, Debug, Eq, PartialEq, clap::ValueEnum)]
#[clap(rename_all = "snake_case")]
enum Mode {
    /// Set the read group of every record.
    OverwriteAll,
    /// Set the read group only of records without one.
    OrphanOnly,
}

#[derive(Debug, clap::Args)]
pub struct Args {
    /// Input SAM, BAM, or CRAM file (`-` for stdin).
    input: String,

    /// @RG line text (`\t` escapes allowed); may be repeated to build one line.
    #[arg(short = 'r', conflicts_with = "read_group_id")]
    read_group_line: Vec<String>,

    /// ID of an @RG line already in the header.
    #[arg(short = 'R')]
    read_group_id: Option<String>,

    /// Mode of operation.
    #[arg(short = 'm', value_enum, default_value_t = Mode::OverwriteAll)]
    mode: Mode,

    /// Write output to this file instead of stdout.
    #[arg(short = 'o')]
    output: Option<PathBuf>,

    /// Replace an existing @RG line with the same ID.
    #[arg(short = 'w')]
    overwrite_header_line: bool,

    /// Output uncompressed BAM.
    #[arg(short = 'u')]
    uncompressed: bool,

    /// Output format [default: from the -o extension, else SAM].
    #[arg(short = 'O', long = "output-fmt", ignore_case = true)]
    output_format: Option<OutputFormat>,

    /// Reference FASTA (indexed) for CRAM input or output.
    #[arg(long = "reference")]
    reference: Option<PathBuf>,

    /// Additional threads for compression and decompression.
    #[arg(short = '@', long = "threads", default_value_t = 0)]
    threads: usize,

    /// Do not add an @PG line to the header.
    #[arg(long = "no-PG")]
    no_pg: bool,
}

pub fn run(args: Args) -> anyhow::Result<()> {
    let read_options = ReadOptions {
        threads: args.threads,
        reference: args.reference.clone(),
    };

    let (mut reader, _) = AlignmentReader::open(&args.input, &read_options)?;
    let input_header = reader.read_header()?;
    let mut header = input_header.clone();

    let read_group_id = set_read_group_line(&mut header, &args)?;

    if !args.no_pg {
        program::add(&mut header)?;
    }

    let format = args
        .output_format
        .or_else(|| args.output.as_deref().and_then(OutputFormat::from_path))
        .unwrap_or(OutputFormat::Sam);

    let write_options = WriteOptions {
        format,
        level: output::parse_level(if args.uncompressed { 0 } else { 6 })?,
        threads: args.threads,
        reference: args.reference,
    };

    let mut writer = AlignmentWriter::create(args.output.as_deref(), &write_options, &header)?;
    writer.write_header(&header)?;

    let mut rg_field = b"RGZ".to_vec();
    rg_field.extend_from_slice(read_group_id.as_ref());
    rg_field.push(0);

    let mut encoder = bam::io::Writer::from(Vec::new());
    let mut edited = Vec::new();

    reader.for_each(&input_header, |record| {
        let raw = match &record {
            Rec::Bam(r) => r.as_ref(),
            _ => {
                encoder.get_mut().clear();
                encoder.write_alignment_record(&input_header, record.as_dyn())?;
                &encoder.get_ref()[size_of::<u32>()..]
            }
        };

        bam_raw::strip_fields(raw, &mut edited, |tag| tag == *b"RG")?;
        let had_read_group = edited.len() != raw.len();

        match args.mode {
            Mode::OverwriteAll => edited.extend_from_slice(&rg_field),
            Mode::OrphanOnly if had_read_group => {
                edited.clear();
                edited.extend_from_slice(raw);
            }
            Mode::OrphanOnly => edited.extend_from_slice(&rg_field),
        }

        writer.write_raw(&header, &edited)?;
        Ok(())
    })?;

    writer.finish(&header)?;

    Ok(())
}

/// Applies -r/-R/-w to the header and returns the read group ID to assign.
fn set_read_group_line(header: &mut sam::Header, args: &Args) -> anyhow::Result<BString> {
    if args.read_group_line.is_empty() {
        let id = match &args.read_group_id {
            Some(id) => {
                if !header.read_groups().contains_key(id.as_bytes()) {
                    bail!(
                        "RG ID supplied does not exist in header. Supply full @RG line with -r instead?"
                    );
                }

                BString::from(id.as_str())
            }
            None => header
                .read_groups()
                .keys()
                .next()
                .cloned()
                .context("No RG specified on command line or in existing header.")?,
        };

        return Ok(id);
    }

    let line = read_group_line(&args.read_group_line)?;
    let id = line_id(&line).context("The supplied RG line lacks an ID tag.")?;

    if header.read_groups().contains_key(id.as_bytes()) {
        if !args.overwrite_header_line {
            bail!("RG line with ID:{id} already present in the header. Use -w to overwrite.");
        }

        header.read_groups_mut().shift_remove(id.as_bytes());
    }

    let parsed: sam::Header = format!("{line}\n").parse().context("invalid @RG line")?;

    for (rg_id, map) in parsed.read_groups() {
        header.read_groups_mut().insert(rg_id.clone(), map.clone());
    }

    if args.mode == Mode::OverwriteAll {
        header
            .read_groups_mut()
            .retain(|rg_id, _| rg_id.as_slice() == id.as_bytes());
    }

    Ok(BString::from(id))
}

/// Joins -r arguments into one @RG line and unescapes `\t` and `\\`, as samtools does.
fn read_group_line(parts: &[String]) -> anyhow::Result<String> {
    let mut line = String::new();

    for (i, part) in parts.iter().enumerate() {
        if i == 0 {
            // samtools prefixes "@RG\t" unless the text starts with @RG.
            if !part.starts_with("@RG") {
                line.push_str("@RG\t");
            }
        } else {
            line.push('\t');
        }

        line.push_str(part);
    }

    let mut unescaped = String::with_capacity(line.len());
    let mut chars = line.chars();

    while let Some(c) = chars.next() {
        if c != '\\' {
            unescaped.push(c);
            continue;
        }

        match chars.next() {
            Some('\\') => unescaped.push('\\'),
            Some('t') => unescaped.push('\t'),
            Some('n') => bail!("\\n in escape sequence is not supported."),
            Some(_) => bail!("Unsupported escape sequence."),
            None => bail!("Unterminated escape sequence."),
        }
    }

    Ok(unescaped)
}

/// The value of the `ID` field of an @RG line.
fn line_id(line: &str) -> Option<String> {
    let start = line.find("\tID:")? + 4;
    let rest = &line[start..];
    Some(rest.split('\t').next().unwrap_or(rest).to_owned())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_read_group_line() -> anyhow::Result<()> {
        assert_eq!(read_group_line(&["ID:a\\tSM:b".into()])?, "@RG\tID:a\tSM:b");
        assert_eq!(
            read_group_line(&["@RG\\tID:a".into(), "SM:b".into()])?,
            "@RG\tID:a\tSM:b"
        );
        assert!(read_group_line(&["ID:a\\n".into()]).is_err());
        assert_eq!(line_id("@RG\tID:a\tSM:b").as_deref(), Some("a"));
        assert_eq!(line_id("@RG\tSM:b"), None);
        Ok(())
    }
}
