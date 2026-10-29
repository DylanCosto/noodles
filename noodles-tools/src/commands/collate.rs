//! `ntools collate`: group records by read name (samtools collate).
//!
//! Records are grouped by sorting on a hash of the name, so the records of a template are
//! adjacent, in input order. The order of templates differs from samtools', which depends on its
//! internal hashing; tools that consume collated input (fastq, fixmate) only need the grouping.

use std::path::{Path, PathBuf};

use anyhow::bail;
use noodles_sam::{
    self as sam,
    alignment::record::Flags,
    header::record::value::{Map, map::header::tag::GROUP_ORDER},
};

use crate::{
    alignment::{AlignmentReader, ReadOptions},
    commands::sort::{SortedRecords, create_temp_dir, parse_memory, set_sort_order},
    output::{self, AlignmentWriter, OutputFormat, WriteOptions},
    program,
    sort_order::Order,
};

/// Memory for buffering records before spilling to temporary files.
const DEFAULT_MEMORY: &str = "768M";

#[derive(Debug, clap::Args)]
pub struct Args {
    /// Input SAM, BAM, or CRAM file (`-` for stdin).
    input: String,

    /// Write `PREFIX.bam` (when neither -o nor -O is given).
    prefix: Option<PathBuf>,

    /// Write output to stdout.
    #[arg(short = 'O')]
    stdout: bool,

    /// Write output to this file.
    #[arg(short = 'o')]
    output: Option<PathBuf>,

    /// Output uncompressed BAM.
    #[arg(short = 'u')]
    uncompressed: bool,

    /// Only keep primary alignments (samtools' fast mode).
    #[arg(short = 'f')]
    fast: bool,

    /// BGZF compression level.
    #[arg(short = 'l', default_value_t = 1)]
    level: u8,

    /// Memory for buffering records before spilling to temporary files (K/M/G suffixes).
    #[arg(short = 'm', default_value = DEFAULT_MEMORY, value_parser = parse_memory)]
    memory: usize,

    /// Prefix for temporary files.
    #[arg(short = 'T')]
    temp_prefix: Option<PathBuf>,

    /// Output format.
    #[arg(long = "output-fmt", ignore_case = true, default_value = "bam")]
    output_format: OutputFormat,

    /// Reference FASTA (indexed) for CRAM input or output.
    #[arg(long = "reference")]
    reference: Option<PathBuf>,

    /// Additional threads for sorting, compression, and decompression.
    #[arg(short = '@', long = "threads", default_value_t = 0)]
    threads: usize,

    /// Do not add an @PG line to the header.
    #[arg(long = "no-PG")]
    no_pg: bool,
}

pub fn run(args: Args) -> anyhow::Result<()> {
    let output: Option<PathBuf> = match (args.stdout, &args.output, &args.prefix) {
        (true, _, _) => None,
        (false, Some(path), _) => Some(path.clone()),
        (false, None, Some(prefix)) => {
            let mut path = prefix.as_os_str().to_owned();
            path.push(".");
            path.push(match args.output_format {
                OutputFormat::Sam => "sam",
                OutputFormat::Bam => "bam",
                OutputFormat::Cram => "cram",
            });
            Some(PathBuf::from(path))
        }
        (false, None, None) => bail!("give an output prefix, -o FILE, or -O for stdout"),
    };

    let options = ReadOptions {
        threads: args.threads,
        reference: args.reference.clone(),
    };

    let (mut reader, _) = AlignmentReader::open(&args.input, &options)?;
    let mut header = reader.read_header()?;

    let temp_dir = create_temp_dir(
        args.temp_prefix.as_deref(),
        output.as_deref().filter(|p| *p != Path::new("-")),
    )?;

    let not_primary = Flags::SECONDARY | Flags::SUPPLEMENTARY;
    let fast = args.fast;

    let sorted = SortedRecords::collect(
        &mut reader,
        &header,
        Order::Collate,
        args.memory.saturating_mul(args.threads + 1),
        args.threads,
        temp_dir.path(),
        |record| Ok(!fast || !record.flags()?.intersects(not_primary)),
    )?;

    set_group_order(&mut header)?;

    if !args.no_pg {
        program::add(&mut header)?;
    }

    let write_options = WriteOptions {
        format: args.output_format,
        level: output::parse_level(if args.uncompressed { 0 } else { args.level })?,
        threads: args.threads,
        reference: args.reference,
    };

    let mut writer = AlignmentWriter::create(output.as_deref(), &write_options, &header)?;
    writer.write_header(&header)?;
    sorted.write_to(&mut writer, &header)?;
    writer.finish(&header)?;

    Ok(())
}

/// Sets `@HD SO:unsorted GO:query`, as samtools collate does.
fn set_group_order(header: &mut sam::Header) -> anyhow::Result<()> {
    set_sort_order(header, Order::Collate)?;

    let hd = header.header_mut().get_or_insert_with(Map::default);
    hd.other_fields_mut().insert(GROUP_ORDER, "query".into());

    Ok(())
}
