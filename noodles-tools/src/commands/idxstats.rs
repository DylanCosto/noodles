//! `ntools idxstats`: per-reference mapped and unmapped counts (samtools idxstats).

use std::{
    io::{self, BufWriter, Write},
    path::PathBuf,
};

use anyhow::Context;
use noodles_csi::BinningIndex;
use noodles_sam as sam;

use crate::{
    alignment::{AlignmentIndex, AlignmentReader, ReadOptions},
    source::Format,
};

#[derive(Debug, clap::Args)]
pub struct Args {
    /// Input BAM (indexed), CRAM, or SAM file.
    input: String,

    /// The BAM file's index (with `-X`).
    #[arg(requires = "customized_index")]
    index: Option<PathBuf>,

    /// Use the index file given after the input.
    #[arg(short = 'X', requires = "index")]
    customized_index: bool,

    /// Reference FASTA (indexed) for CRAM input.
    #[arg(short = 'T', long = "reference")]
    reference: Option<PathBuf>,

    /// Additional threads for BGZF decompression when scanning.
    #[arg(short = '@', long = "threads", default_value_t = 0)]
    threads: usize,
}

/// Per-reference counts, plus records with no reference.
struct Stats {
    by_reference: Vec<(u64, u64)>,
    unplaced: u64,
}

pub fn run(args: Args) -> anyhow::Result<()> {
    let options = ReadOptions {
        threads: args.threads,
        reference: args.reference,
    };

    let (mut reader, path) = AlignmentReader::open(&args.input, &options)?;
    let header = reader.read_header()?;

    // samtools reads counts from the index for BAM only.
    let index_stats = match (reader.format(), &path) {
        (Format::Bam, Some(path)) => {
            match AlignmentIndex::read(Format::Bam, path, args.index.as_deref()) {
                Ok(AlignmentIndex::Bam(index)) => from_index(&index, &header),
                _ => {
                    eprintln!(
                        "samtools idxstats: fail to load index for \"{}\", reverting to slow \
                         method",
                        args.input
                    );
                    None
                }
            }
        }
        _ => None,
    };

    let stats = match index_stats {
        Some(stats) => stats,
        None => scan(&mut reader, &header).context("counting records")?,
    };

    let mut stdout = BufWriter::new(io::stdout().lock());

    for ((name, reference_sequence), (mapped, unmapped)) in
        header.reference_sequences().iter().zip(&stats.by_reference)
    {
        writeln!(
            stdout,
            "{name}\t{}\t{mapped}\t{unmapped}",
            reference_sequence.length()
        )?;
    }

    writeln!(stdout, "*\t0\t0\t{}", stats.unplaced)?;
    stdout.flush()?;

    Ok(())
}

/// Reads counts from index metadata, as samtools does for BAM. Returns `None` when the index has
/// no metadata, so the caller can fall back to scanning.
fn from_index<I: BinningIndex>(index: &I, header: &sam::Header) -> Option<Stats> {
    let mut by_reference = vec![(0, 0); header.reference_sequences().len()];

    for (counts, reference_sequence) in by_reference.iter_mut().zip(index.reference_sequences()) {
        if let Some(metadata) = reference_sequence.metadata() {
            *counts = (
                metadata.mapped_record_count(),
                metadata.unmapped_record_count(),
            );
        }
    }

    Some(Stats {
        by_reference,
        unplaced: index.unplaced_unmapped_record_count()?,
    })
}

fn scan(reader: &mut AlignmentReader, header: &sam::Header) -> anyhow::Result<Stats> {
    let mut stats = Stats {
        by_reference: vec![(0, 0); header.reference_sequences().len()],
        unplaced: 0,
    };

    reader.for_each(header, |record| {
        let unmapped = record.flags()?.is_unmapped();

        match record.reference_sequence_id(header)? {
            Some(id) => {
                let counts = stats
                    .by_reference
                    .get_mut(id)
                    .with_context(|| format!("invalid reference sequence ID: {id}"))?;

                if unmapped {
                    counts.1 += 1;
                } else {
                    counts.0 += 1;
                }
            }
            None => stats.unplaced += 1,
        }

        Ok(())
    })?;

    Ok(stats)
}
