//! `ntools view`: convert, filter, count, and query alignments (samtools view).

use std::{
    io::{self, Write},
    path::PathBuf,
};

use anyhow::bail;
use noodles_sam::alignment::record::Flags;

use crate::{
    alignment::{AlignmentReader, ReadOptions, Rec},
    flags,
    output::{self, AlignmentWriter, OutputFormat, WriteOptions},
    program,
};

/// Mapping quality samtools uses for "unavailable"; it passes every `-q` threshold.
const MISSING_MAPPING_QUALITY: u8 = 255;

#[derive(Debug, clap::Args)]
#[command(disable_help_flag = true)]
pub struct Args {
    /// Input SAM, BAM, or CRAM file (`-` for stdin).
    input: String,

    /// Regions to query (`chr1`, `chr1:100-200`, or `*` for unplaced unmapped reads); requires an
    /// index.
    regions: Vec<String>,

    /// Write output to this file instead of stdout.
    #[arg(short = 'o', long = "output")]
    output: Option<PathBuf>,

    /// Output format.
    #[arg(short = 'O', long = "output-fmt", ignore_case = true)]
    output_format: Option<OutputFormat>,

    /// Output BAM.
    #[arg(short = 'b', long = "bam")]
    bam: bool,

    /// Output CRAM (requires -T).
    #[arg(short = 'C', long = "cram")]
    cram: bool,

    /// Output uncompressed BAM.
    #[arg(short = 'u', long = "uncompressed")]
    uncompressed: bool,

    /// BGZF compression level for BAM output.
    #[arg(short = 'l', default_value_t = 6)]
    level: u8,

    /// Include the header in SAM output.
    #[arg(short = 'h', long = "with-header")]
    with_header: bool,

    /// Print the header only.
    #[arg(short = 'H', long = "header-only")]
    header_only: bool,

    /// Print only the count of matching records.
    #[arg(short = 'c', long = "count")]
    count: bool,

    /// Only include reads with all of these flags set.
    #[arg(short = 'f', long = "require-flags", value_parser = flags::parse)]
    require_flags: Option<Flags>,

    /// Exclude reads with any of these flags set.
    #[arg(short = 'F', long = "excl-flags", value_parser = flags::parse)]
    exclude_flags: Option<Flags>,

    /// Exclude reads with all of these flags set.
    #[arg(short = 'G', value_parser = flags::parse)]
    exclude_all_flags: Option<Flags>,

    /// Exclude reads with mapping quality below this value.
    #[arg(short = 'q', long = "min-MQ", default_value_t = 0)]
    min_mapping_quality: u8,

    /// Reference FASTA (indexed) for CRAM input or output.
    #[arg(short = 'T', long = "reference")]
    reference: Option<PathBuf>,

    /// Additional threads for BGZF compression and decompression.
    #[arg(short = '@', long = "threads", default_value_t = 0)]
    threads: usize,

    /// Do not add an @PG line to the header.
    #[arg(long = "no-PG")]
    no_pg: bool,

    /// Print help.
    #[arg(long, action = clap::ArgAction::Help)]
    help: Option<bool>,
}

/// A samtools-style record filter.
#[derive(Clone, Copy, Debug)]
pub struct Filter {
    pub require: Flags,
    pub exclude: Flags,
    pub exclude_all: Option<Flags>,
    pub min_mapping_quality: u8,
}

impl Filter {
    pub fn matches(&self, record: &Rec<'_>) -> io::Result<bool> {
        let flags = record.flags()?;

        if !flags.contains(self.require) || flags.intersects(self.exclude) {
            return Ok(false);
        }

        if let Some(all) = self.exclude_all
            && flags.contains(all)
        {
            return Ok(false);
        }

        if self.min_mapping_quality > 0 {
            let mapq = record.mapping_quality()?.unwrap_or(MISSING_MAPPING_QUALITY);

            if mapq < self.min_mapping_quality {
                return Ok(false);
            }
        }

        Ok(true)
    }
}

pub fn run(args: Args) -> anyhow::Result<()> {
    let format = output_format(&args)?;

    let read_options = ReadOptions {
        threads: args.threads,
        reference: args.reference.clone(),
    };

    let (mut reader, path) = AlignmentReader::open(&args.input, &read_options)?;
    let mut header = reader.read_header()?;

    let filter = Filter {
        require: args.require_flags.unwrap_or_else(Flags::empty),
        exclude: args.exclude_flags.unwrap_or_else(Flags::empty),
        exclude_all: args.exclude_all_flags,
        min_mapping_quality: args.min_mapping_quality,
    };

    if args.count {
        let mut n = 0u64;
        visit(
            &mut reader,
            path.as_deref(),
            &header,
            &args.regions,
            |rec| {
                if filter.matches(&rec)? {
                    n += 1;
                }
                Ok(())
            },
        )?;

        let mut stdout = io::stdout().lock();
        writeln!(stdout, "{n}")?;
        return Ok(());
    }

    if !args.no_pg {
        program::add(&mut header)?;
    }

    let write_options = WriteOptions {
        format,
        level: output::parse_level(if args.uncompressed { 0 } else { args.level })?,
        threads: args.threads,
        reference: args.reference.clone(),
    };

    let mut writer = AlignmentWriter::create(args.output.as_deref(), &write_options)?;

    if format != OutputFormat::Sam || args.with_header || args.header_only {
        writer.write_header(&header)?;
    }

    if !args.header_only {
        visit(
            &mut reader,
            path.as_deref(),
            &header,
            &args.regions,
            |rec| {
                if filter.matches(&rec)? {
                    writer.write(&header, &rec)?;
                }
                Ok(())
            },
        )?;
    }

    writer.finish(&header)?;

    Ok(())
}

fn visit<F>(
    reader: &mut AlignmentReader,
    path: Option<&std::path::Path>,
    header: &noodles_sam::Header,
    regions: &[String],
    f: F,
) -> anyhow::Result<()>
where
    F: FnMut(Rec<'_>) -> anyhow::Result<()>,
{
    if regions.is_empty() {
        reader.for_each(header, f)
    } else {
        let Some(path) = path else {
            bail!("region queries need an indexed file, not stdin");
        };

        reader.for_each_in_regions(path, header, regions, f)
    }
}

fn output_format(args: &Args) -> anyhow::Result<OutputFormat> {
    let explicit = [
        args.output_format,
        (args.bam || args.uncompressed).then_some(OutputFormat::Bam),
        args.cram.then_some(OutputFormat::Cram),
    ];

    let mut chosen = explicit.into_iter().flatten();

    match (chosen.next(), chosen.next()) {
        (Some(a), Some(b)) if a != b => bail!("conflicting output formats: {a:?} and {b:?}"),
        (Some(format), _) => Ok(format),
        (None, _) => Ok(args
            .output
            .as_deref()
            .and_then(OutputFormat::from_path)
            .unwrap_or(OutputFormat::Sam)),
    }
}
