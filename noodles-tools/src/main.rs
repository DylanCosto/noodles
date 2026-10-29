//! `ntools`: a samtools-compatible command-line toolkit built on noodles.

mod alignment;
mod bgzf_stream;
mod commands;
mod flags;
mod output;
mod program;
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
}

fn main() -> ExitCode {
    let cli = Cli::parse();

    let result = match cli.command {
        Command::View(args) => commands::view::run(args),
        Command::Flagstat(args) => commands::flagstat::run(args),
        Command::Idxstats(args) => commands::idxstats::run(args),
    };

    match result {
        Ok(()) => ExitCode::SUCCESS,
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
