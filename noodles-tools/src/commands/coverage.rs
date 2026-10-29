//! `ntools coverage`: per-reference coverage summary (samtools coverage).

use std::{
    fs::File,
    io::{BufWriter, Write},
    path::{Path, PathBuf},
};

use anyhow::Context;
use noodles_sam::alignment::record::Flags;

use crate::{
    alignment::{AlignmentReader, Queries, ReadOptions},
    depth::{CountOptions, DEFAULT_EXCLUDE_FLAGS, DepthWindow, format_g},
    filter::Filter,
    flags, output,
    regions::{RegionSet, record_interval},
};

#[derive(Debug, clap::Args)]
#[command(disable_help_flag = true)]
pub struct Args {
    /// Coordinate-sorted SAM, BAM, or CRAM file.
    input: String,

    /// Only summarize this region.
    #[arg(short = 'r', long = "region")]
    region: Option<String>,

    /// Exclude bases with base quality below this value.
    #[arg(short = 'q', long = "min-BQ", default_value_t = 0)]
    min_base_quality: u8,

    /// Exclude reads with mapping quality below this value.
    #[arg(short = 'Q', long = "min-MQ", default_value_t = 0)]
    min_mapping_quality: u8,

    /// Exclude reads shorter than this.
    #[arg(short = 'l', long = "min-read-len", default_value_t = 0)]
    min_read_length: usize,

    /// Exclude reads with any of these flags [UNMAP,SECONDARY,QCFAIL,DUP].
    #[arg(long = "ff", value_parser = flags::parse)]
    exclude_flags: Option<Flags>,

    /// Only include reads with all of these flags.
    #[arg(long = "rf", value_parser = flags::parse)]
    require_flags: Option<Flags>,

    /// Do not print the header line.
    #[arg(short = 'H', long = "no-header")]
    no_header: bool,

    /// Write output to this file instead of stdout.
    #[arg(short = 'o', long = "output")]
    output: Option<PathBuf>,

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

#[derive(Clone, Copy, Debug, Default)]
struct Stats {
    reads: u64,
    covered_bases: u64,
    depth_sum: u64,
    base_quality_sum: u64,
    mapping_quality_sum: u64,
}

pub fn run(args: Args) -> anyhow::Result<()> {
    let options = ReadOptions {
        threads: args.threads,
        reference: args.reference.clone(),
    };

    let (mut reader, path) = AlignmentReader::open(&args.input, &options)?;
    let header = reader.read_header()?;

    let filter = Filter {
        require: args.require_flags.unwrap_or_else(Flags::empty),
        exclude: args.exclude_flags.unwrap_or(DEFAULT_EXCLUDE_FLAGS),
        min_mapping_quality: args.min_mapping_quality,
        ..Filter::default()
    };

    let mut region = None;

    if let Some(raw) = &args.region {
        let mut set = RegionSet::new(&header);
        set.add_region(&header, raw)?;
        region = Some(set);
    }

    let lengths: Vec<usize> = header
        .reference_sequences()
        .values()
        .map(|rs| rs.length().get())
        .collect();

    let mut stats = vec![Stats::default(); lengths.len()];

    let mut window = DepthWindow::new(CountOptions {
        min_base_quality: args.min_base_quality,
        include_deletions: false,
        track_quality: true,
    });

    let in_region = |id: usize, position: usize| match &region {
        Some(set) => set
            .intervals(id)
            .iter()
            .any(|&(s, e)| (s..=e).contains(&position)),
        None => true,
    };

    let mut emit = |id: usize, position: usize, cell: crate::depth::Cell| {
        if cell.depth > 0 && in_region(id, position) {
            let s = &mut stats[id];
            s.covered_bases += 1;
            s.depth_sum += u64::from(cell.depth);
            s.base_quality_sum += cell.quality_sum;
        }
        Ok(())
    };

    let mut read_stats = vec![(0u64, 0u64); lengths.len()];

    let mut visit = |record: crate::alignment::Rec<'_>| -> anyhow::Result<()> {
        if !filter.matches(&record, &header)? {
            return Ok(());
        }

        if args.min_read_length > 0 && record.as_dyn().sequence().len() < args.min_read_length {
            return Ok(());
        }

        if let Some((id, start, end)) = record_interval(record.as_dyn(), &header)? {
            let counted = match &region {
                Some(set) => set.first_overlap(id, start, end).is_some(),
                None => true,
            };

            if counted {
                let mapq = record.mapping_quality()?.unwrap_or(u8::MAX);
                read_stats[id].0 += 1;
                read_stats[id].1 += u64::from(mapq);
            }
        }

        window.add(&record, &header, &mut emit)?;
        Ok(())
    };

    match (&region, path.as_deref()) {
        (Some(set), Some(path)) => {
            reader.for_each_in_queries(path, None, &header, &Queries::Set(set), &mut visit)?;
        }
        _ => reader.for_each(&header, &mut visit)?,
    }

    window.finish(&mut emit)?;

    for (s, &(reads, mapq_sum)) in stats.iter_mut().zip(&read_stats) {
        s.reads = reads;
        s.mapping_quality_sum = mapq_sum;
    }

    let sink: output::Sink = match args.output.as_deref().filter(|p| *p != Path::new("-")) {
        Some(path) => {
            Box::new(File::create(path).with_context(|| format!("creating {}", path.display()))?)
        }
        None => output::raw_stdout()?,
    };

    let mut out = BufWriter::new(sink);

    if !args.no_header {
        writeln!(
            out,
            "#rname\tstartpos\tendpos\tnumreads\tcovbases\tcoverage\tmeandepth\tmeanbaseq\tmeanmapq"
        )?;
    }

    for (id, name) in header.reference_sequences().keys().enumerate() {
        let (start, end) = match &region {
            Some(set) => match set.intervals(id).first() {
                Some(&(s, e)) => (s, e.min(lengths[id])),
                None => continue,
            },
            None => (1, lengths[id]),
        };

        let s = stats[id];
        let len = (end + 1 - start) as f64;

        writeln!(
            out,
            "{name}\t{start}\t{end}\t{}\t{}\t{}\t{}\t{}\t{}",
            s.reads,
            s.covered_bases,
            format_g(s.covered_bases as f64 / len * 100.0, 6),
            format_g(s.depth_sum as f64 / len, 6),
            format_g(ratio(s.base_quality_sum, s.depth_sum), 3),
            format_g(ratio(s.mapping_quality_sum, s.reads), 3),
        )?;
    }

    out.flush()?;

    Ok(())
}

fn ratio(numerator: u64, denominator: u64) -> f64 {
    if denominator == 0 {
        0.0
    } else {
        numerator as f64 / denominator as f64
    }
}
