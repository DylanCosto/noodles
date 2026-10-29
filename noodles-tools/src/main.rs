//! `ntools`: a samtools-compatible command-line toolkit built on noodles.

mod alignment;
mod bam_blocks;
mod bam_raw;
mod bam_sam;
mod bgzf_stream;
mod commands;
mod depth;
mod filter;
mod flags;
mod output;
mod parallel_cram;
mod parallel_sam;
mod program;
mod regions;
mod runs;
mod sequence;
mod sort_order;
mod source;

use std::{io, process::ExitCode};

use clap::{Parser, Subcommand};

#[derive(Debug, Parser)]
#[command(name = "ntools", version, about)]
struct Cli {
    #[command(subcommand)]
    command: Command,
}

#[derive(Debug, Subcommand)]
enum Command {
    /// Convert, filter, count, and query SAM/BAM/CRAM.
    View(commands::view::Args),
    /// Count records by flag category.
    Flagstat(commands::flagstat::Args),
    /// Report per-reference mapped and unmapped counts.
    Idxstats(commands::idxstats::Args),
    /// Build a BAI, CSI, or CRAI index.
    Index(commands::index::Args),
    /// Index a FASTA file or extract regions from it.
    Faidx(commands::faidx::Args),
    /// Sort alignments by coordinate or read name.
    Sort(commands::sort::Args),
    /// Merge sorted BAM files.
    Merge(commands::merge::Args),
    /// Convert alignments to FASTQ.
    Fastq(commands::fastq::Args),
    /// Convert alignments to FASTA.
    Fasta(commands::fastq::Args),
    /// Check that files look intact.
    Quickcheck(commands::quickcheck::Args),
    /// Replace a BAM header without recompressing records.
    Reheader(commands::reheader::Args),
    /// Concatenate BAM files without recompressing records.
    Cat(commands::cat::Args),
    /// Report per-position read depth.
    Depth(commands::depth::Args),
    /// Summarize coverage per reference sequence.
    Coverage(commands::coverage::Args),
}

fn main() -> ExitCode {
    let cli = Cli::parse();

    let result = match cli.command {
        Command::View(args) => commands::view::run(args).map(|()| ExitCode::SUCCESS),
        Command::Flagstat(args) => commands::flagstat::run(args).map(|()| ExitCode::SUCCESS),
        Command::Idxstats(args) => commands::idxstats::run(args).map(|()| ExitCode::SUCCESS),
        Command::Index(args) => commands::index::run(args).map(|()| ExitCode::SUCCESS),
        Command::Faidx(args) => commands::faidx::run(args).map(|()| ExitCode::SUCCESS),
        Command::Sort(args) => commands::sort::run(args).map(|()| ExitCode::SUCCESS),
        Command::Merge(args) => commands::merge::run(args).map(|()| ExitCode::SUCCESS),
        Command::Fastq(args) => commands::fastq::run(args, commands::fastq::ReadFormat::Fastq)
            .map(|()| ExitCode::SUCCESS),
        Command::Fasta(args) => commands::fastq::run(args, commands::fastq::ReadFormat::Fasta)
            .map(|()| ExitCode::SUCCESS),
        Command::Quickcheck(args) => commands::quickcheck::run(args),
        Command::Reheader(args) => commands::reheader::run(args).map(|()| ExitCode::SUCCESS),
        Command::Cat(args) => commands::cat::run(args).map(|()| ExitCode::SUCCESS),
        Command::Depth(args) => commands::depth::run(args).map(|()| ExitCode::SUCCESS),
        Command::Coverage(args) => commands::coverage::run(args).map(|()| ExitCode::SUCCESS),
    };

    match result {
        Ok(code) => code,
        // A closed downstream pipe (`| head`) is not an error worth reporting.
        Err(e) if is_broken_pipe(&e) => ExitCode::SUCCESS,
        Err(e) => {
            eprintln!("ntools: {e:#}");
            ExitCode::FAILURE
        }
    }
}

fn is_broken_pipe(e: &anyhow::Error) -> bool {
    e.chain()
        .filter_map(|cause| cause.downcast_ref::<io::Error>())
        .any(|e| e.kind() == io::ErrorKind::BrokenPipe)
}
