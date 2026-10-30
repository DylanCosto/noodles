//! `ntools`: a samtools-compatible command-line toolkit built on noodles.

mod alignment;
mod bam_blocks;
mod bam_raw;
mod bam_sam;
mod bgzf_stream;
mod commands;
mod cram_file;
mod cram_raw;
mod depth;
mod expr;
mod file_list;
mod filter;
mod flags;
mod hts_region;
mod khash;
mod output;
mod parallel_cram;
mod parallel_sam;
mod program;
mod reference;
mod regions;
mod runs;
mod sanitize;
mod sequence;
mod sort_order;
mod source;
mod template_coordinate;

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
    /// Index a FASTQ file or extract regions from it.
    Fqidx(commands::faidx::Args),
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
    /// Write one file per read group or tag value.
    Split(commands::split::Args),
    /// Group records by read name.
    Collate(commands::collate::Args),
    /// Fill in mate information from name-grouped alignments.
    Fixmate(commands::fixmate::Args),
    /// Mark duplicate alignments.
    Markdup(commands::markdup::Args),
    /// Recalculate MD and NM tags.
    Calmd(commands::calmd::Args),
    /// Text pileup of alignments.
    Mpileup(commands::mpileup::Args),
    /// Alignment statistics for plot-bamstats.
    Stats(commands::stats::Args),
    /// Write a SAM sequence dictionary for a FASTA file.
    Dict(commands::dict::Args),
    /// Sum read depth over BED regions.
    Bedcov(commands::bedcov::Args),
    /// Add or replace read group tags and header lines.
    Addreplacerg(commands::addreplacerg::Args),
    /// Clip primer regions from reads.
    Ampliconclip(commands::ampliconclip::Args),
    /// Print the header and first records as SAM.
    Head(commands::head::Args),
    /// Convert FASTQ or FASTA to unaligned SAM/BAM/CRAM.
    Import(commands::import::Args),
    /// Turn alignments back into unaligned reads.
    Reset(commands::reset::Args),
    /// Order-independent checksums of reads, qualities, and tags.
    Checksum(commands::checksum::Args),
    /// Convert between numeric and textual flags.
    Flags(commands::flags::Args),
    /// Report CRAM block sizes by content ID and data series.
    #[command(alias = "cram_size")]
    CramSize(commands::cram_size::Args),
    /// Rebuild a reference from MD tags or embedded CRAM references.
    Reference(commands::reference::Args),
    /// Remove PCR duplicates (superseded by markdup).
    Rmdup(commands::rmdup::Args),
    /// List the samples (read group SM values) of alignment files.
    Samples(commands::samples::Args),
    /// Generate a consensus sequence from aligned reads.
    #[command(disable_help_flag = true)]
    Consensus(commands::consensus::Args),
}

fn main() -> ExitCode {
    let cli = Cli::parse();

    let result = match cli.command {
        Command::View(args) => commands::view::run(args).map(|()| ExitCode::SUCCESS),
        Command::Flagstat(args) => commands::flagstat::run(args).map(|()| ExitCode::SUCCESS),
        Command::Idxstats(args) => commands::idxstats::run(args).map(|()| ExitCode::SUCCESS),
        Command::Index(args) => commands::index::run(args).map(|()| ExitCode::SUCCESS),
        Command::Faidx(args) => commands::faidx::run(args),
        Command::Fqidx(args) => commands::faidx::run_fqidx(args),
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
        Command::Split(args) => commands::split::run(args).map(|()| ExitCode::SUCCESS),
        Command::Collate(args) => commands::collate::run(args).map(|()| ExitCode::SUCCESS),
        Command::Fixmate(args) => commands::fixmate::run(args).map(|()| ExitCode::SUCCESS),
        Command::Markdup(args) => commands::markdup::run(args).map(|()| ExitCode::SUCCESS),
        Command::Calmd(args) => commands::calmd::run(args).map(|()| ExitCode::SUCCESS),
        Command::Mpileup(args) => commands::mpileup::run(args).map(|()| ExitCode::SUCCESS),
        Command::Stats(args) => commands::stats::run(args).map(|()| ExitCode::SUCCESS),
        Command::Dict(args) => commands::dict::run(args).map(|()| ExitCode::SUCCESS),
        Command::Bedcov(args) => commands::bedcov::run(args),
        Command::Addreplacerg(args) => {
            commands::addreplacerg::run(args).map(|()| ExitCode::SUCCESS)
        }
        Command::Ampliconclip(args) => {
            commands::ampliconclip::run(args).map(|()| ExitCode::SUCCESS)
        }
        Command::Reset(args) => commands::reset::run(args).map(|()| ExitCode::SUCCESS),
        Command::Checksum(args) => commands::checksum::run(args).map(|()| ExitCode::SUCCESS),
        Command::Import(args) => commands::import::run(args).map(|()| ExitCode::SUCCESS),
        Command::Head(args) => commands::head::run(args).map(|()| ExitCode::SUCCESS),
        Command::Flags(args) => commands::flags::run(args),
        Command::CramSize(args) => commands::cram_size::run(args).map(|()| ExitCode::SUCCESS),
        Command::Reference(args) => commands::reference::run(args).map(|()| ExitCode::SUCCESS),
        Command::Rmdup(args) => commands::rmdup::run(args).map(|()| ExitCode::SUCCESS),
        Command::Samples(args) => commands::samples::run(args),
        Command::Consensus(args) => commands::consensus::run(args),
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
