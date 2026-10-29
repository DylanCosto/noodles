//! `ntools bedcov`: summed read depth over BED regions (samtools bedcov).

use std::{
    fs,
    io::{self, BufWriter, Write},
    path::PathBuf,
};

use anyhow::{Context, bail};
use noodles_sam::alignment::record::Flags;

use crate::{
    alignment::{AlignmentIndex, AlignmentReader, Queries, ReadOptions},
    depth::{CountOptions, DEFAULT_EXCLUDE_FLAGS, DepthWindow},
    filter::Filter,
    flags, output,
    regions::record_interval,
};

#[derive(Debug, clap::Args)]
#[command(disable_help_flag = true)]
pub struct Args {
    /// BED file of regions.
    bed: PathBuf,

    /// Indexed BAM or CRAM files.
    #[arg(required = true)]
    inputs: Vec<String>,

    /// Exclude reads with mapping quality below this value.
    #[arg(short = 'Q', long = "min-MQ", default_value_t = 0)]
    min_mapping_quality: u8,

    /// Remove these flags from the default exclusion list (UNMAP,SECONDARY,QCFAIL,DUP).
    #[arg(short = 'g', value_parser = flags::parse)]
    unexclude_flags: Option<Flags>,

    /// Add these flags to the exclusion list.
    #[arg(short = 'G', value_parser = flags::parse)]
    exclude_flags: Option<Flags>,

    /// Do not count deletions and reference skips.
    #[arg(short = 'j')]
    no_deletions: bool,

    /// Also report the number of reference bases with at least this depth.
    #[arg(short = 'd')]
    depth_threshold: Option<u32>,

    /// Also report read counts.
    #[arg(short = 'c')]
    read_counts: bool,

    /// Print a header line.
    #[arg(short = 'H')]
    print_header: bool,

    /// Reference FASTA (indexed) for CRAM input.
    #[arg(long = "reference")]
    reference: Option<PathBuf>,

    /// Additional decompression threads.
    #[arg(short = '@', long = "threads", default_value_t = 0)]
    threads: usize,

    /// Print help.
    #[arg(long, action = clap::ArgAction::Help)]
    help: Option<bool>,
}

/// Per-input results for one region.
#[derive(Clone, Copy, Debug, Default)]
struct RegionCounts {
    depth_sum: u64,
    deep_bases: u64,
    reads: u64,
}

pub fn run(args: Args) -> anyhow::Result<()> {
    let options = ReadOptions {
        threads: args.threads,
        reference: args.reference.clone(),
    };

    let mut inputs = Vec::with_capacity(args.inputs.len());

    for src in &args.inputs {
        let (mut reader, path) = AlignmentReader::open(src, &options)?;
        let header = reader.read_header()?;
        let path = path.with_context(|| format!("{src}: bedcov needs indexed files, not stdin"))?;
        let index = AlignmentIndex::read(reader.format(), &path, None)?;
        inputs.push((reader, index, header));
    }

    let exclude = (DEFAULT_EXCLUDE_FLAGS - args.unexclude_flags.unwrap_or_else(Flags::empty))
        | args.exclude_flags.unwrap_or_else(Flags::empty);

    let filter = Filter {
        exclude,
        min_mapping_quality: args.min_mapping_quality,
        ..Filter::default()
    };

    let count_options = CountOptions {
        include_deletions: !args.no_deletions,
        include_skips: !args.no_deletions,
        ..CountOptions::default()
    };

    let text =
        fs::read_to_string(&args.bed).with_context(|| format!("reading {}", args.bed.display()))?;

    let mut out = BufWriter::new(output::raw_stdout()?);

    if args.print_header {
        write_header(&mut out, &args)?;
    }

    let mut counts = vec![RegionCounts::default(); inputs.len()];

    for (i, line) in text.lines().enumerate() {
        let line = line.trim_end_matches('\r');

        if line.is_empty() || line.starts_with('#') || line.starts_with("track") {
            continue;
        }

        let mut fields = line.split('\t');
        let (Some(name), Some(start), Some(end)) = (fields.next(), fields.next(), fields.next())
        else {
            bail!(
                "{}:{}: expected at least 3 fields",
                args.bed.display(),
                i + 1
            );
        };

        let start: usize = start
            .parse()
            .with_context(|| format!("{}:{}: invalid start", args.bed.display(), i + 1))?;
        let end: usize = end
            .parse()
            .with_context(|| format!("{}:{}: invalid end", args.bed.display(), i + 1))?;

        for ((reader, index, header), counts) in inputs.iter_mut().zip(counts.iter_mut()) {
            *counts = count_region(
                reader,
                index,
                header,
                &filter,
                count_options,
                args.depth_threshold,
                (name, start, end),
            )?;
        }

        out.write_all(line.as_bytes())?;

        for c in &counts {
            write!(out, "\t{}", c.depth_sum)?;
        }

        if args.depth_threshold.is_some() {
            for c in &counts {
                write!(out, "\t{}", c.deep_bases)?;
            }
        }

        if args.read_counts {
            for c in &counts {
                write!(out, "\t{}", c.reads)?;
            }
        }

        writeln!(out)?;
    }

    out.flush()?;

    Ok(())
}

/// Sums depth over the 0-based half-open region `[start, end)` of one input.
fn count_region(
    reader: &mut AlignmentReader,
    index: &AlignmentIndex,
    header: &noodles_sam::Header,
    filter: &Filter,
    count_options: CountOptions,
    depth_threshold: Option<u32>,
    (name, start, end): (&str, usize, usize),
) -> anyhow::Result<RegionCounts> {
    let mut counts = RegionCounts::default();

    let Some(id) = header.reference_sequences().get_index_of(name.as_bytes()) else {
        return Ok(counts);
    };

    if end <= start {
        return Ok(counts);
    }

    let (first, last) = (start + 1, end);
    let region = format!("{name}:{first}-{last}");
    let regions = [region];

    let mut window = DepthWindow::new(count_options);

    let mut emit = |rid: usize, position: usize, cell: crate::depth::Cell| -> io::Result<()> {
        if rid == id && (first..=last).contains(&position) {
            counts.depth_sum += u64::from(cell.depth);

            if depth_threshold.is_some_and(|d| cell.depth >= d) {
                counts.deep_bases += 1;
            }
        }

        Ok(())
    };

    let mut reads = 0;

    reader.for_each_in_indexed_queries(index, header, &Queries::Regions(&regions), |record| {
        if !filter.matches(&record, header)? {
            return Ok(());
        }

        if let Some((rid, s, e)) = record_interval(record.as_dyn(), header)?
            && rid == id
            && s <= last
            && e >= first
        {
            reads += 1;
        }

        window.add(&record, header, &mut emit)?;
        Ok(())
    })?;

    window.finish(&mut emit)?;
    counts.reads = reads;

    Ok(counts)
}

fn write_header<W: Write>(out: &mut W, args: &Args) -> io::Result<()> {
    write!(out, "#chrom\tchromStart\tchromEnd")?;

    for src in &args.inputs {
        write!(out, "\t{src}_cov")?;
    }

    if args.depth_threshold.is_some() {
        for src in &args.inputs {
            write!(out, "\t{src}_depth")?;
        }
    }

    if args.read_counts {
        for src in &args.inputs {
            write!(out, "\t{src}_count")?;
        }
    }

    writeln!(out)
}
