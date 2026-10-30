//! `ntools bedcov`: summed read depth over BED regions (samtools bedcov).

use std::{
    fs,
    io::{self, BufWriter, Write},
    path::PathBuf,
    process::ExitCode,
};

use anyhow::Context;
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

pub fn run(args: Args) -> anyhow::Result<ExitCode> {
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

    let mut status = ExitCode::SUCCESS;

    for line in text.lines() {
        let line = line.trim_end_matches('\r');

        if line.is_empty()
            || line.starts_with('#')
            || line.starts_with("track ")
            || line.starts_with("browser ")
        {
            continue;
        }

        let Some((name, start, end)) = parse_bed_line(line, &inputs[0].2) else {
            eprintln!("Errors in BED line '{line}'");
            status = ExitCode::from(2);
            continue;
        };

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

    Ok(status)
}

/// A BED line's reference sequence name and 0-based start and end, read as samtools does: the
/// name up to the first whitespace, then two integers. `None` if they are missing, the name is
/// not a reference sequence of `header`, or the end comes before the start.
fn parse_bed_line<'l>(
    line: &'l str,
    header: &noodles_sam::Header,
) -> Option<(&'l str, usize, usize)> {
    let (name, rest) = line.split_once(|c: char| c.is_ascii_whitespace())?;

    header.reference_sequences().get_index_of(name.as_bytes())?;

    let (start, rest) = scan_integer(rest)?;
    let (end, _) = scan_integer(rest)?;

    if end < start {
        return None;
    }

    // Like HTSlib's iterators, negative coordinates count from 0.
    let clamp = |n: i64| usize::try_from(n.max(0)).unwrap_or(usize::MAX);
    Some((name, clamp(start), clamp(end)))
}

/// `sscanf("%ld")`: skips leading whitespace and reads an optionally signed integer prefix;
/// returns it and the rest of the string.
fn scan_integer(s: &str) -> Option<(i64, &str)> {
    let s = s.trim_start();
    let sign_len = usize::from(s.starts_with(['+', '-']));
    let digits = s[sign_len..].bytes().take_while(u8::is_ascii_digit).count();

    if digits == 0 {
        return None;
    }

    let (number, rest) = s.split_at(sign_len + digits);
    let n = number.parse().ok()?;
    Some((n, rest))
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

            // As in samtools, the threshold applies to aligned bases, without deletions and
            // reference skips.
            let aligned = if count_options.include_deletions {
                cell.depth - cell.gaps
            } else {
                cell.depth
            };

            if depth_threshold.is_some_and(|d| aligned >= d) {
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
