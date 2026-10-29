//! `ntools quickcheck`: check that alignment files look intact (samtools quickcheck).

use std::{fs::File, path::PathBuf, process::ExitCode};

use anyhow::Context;

use crate::{
    alignment::{AlignmentReader, ReadOptions},
    bam_blocks,
    source::{Format, Input},
};

/// The CRAM 3 EOF container.
const CRAM_3_EOF: [u8; 38] = [
    0x0f, 0x00, 0x00, 0x00, 0xff, 0xff, 0xff, 0xff, 0x0f, 0xe0, 0x45, 0x4f, 0x46, 0x00, 0x00, 0x00,
    0x00, 0x01, 0x00, 0x05, 0xbd, 0xd9, 0x4f, 0x00, 0x01, 0x00, 0x06, 0x06, 0x01, 0x00, 0x01, 0x00,
    0x01, 0x00, 0xee, 0x63, 0x01, 0x4b,
];

#[derive(Debug, clap::Args)]
pub struct Args {
    /// Files to check.
    #[arg(required = true)]
    inputs: Vec<PathBuf>,

    /// Print the names of files that fail to stdout; repeat for reasons on stderr.
    #[arg(short = 'v', action = clap::ArgAction::Count)]
    verbose: u8,

    /// Do not print messages about failures to stderr.
    #[arg(short = 'q')]
    quiet: bool,

    /// Do not require @SQ lines (for unmapped files).
    #[arg(short = 'u')]
    unmapped: bool,
}

pub fn run(args: Args) -> anyhow::Result<ExitCode> {
    let mut failed = false;

    for path in &args.inputs {
        if let Err(e) = check(path, args.unmapped) {
            failed = true;

            if args.verbose > 0 {
                println!("{}", path.display());
            }

            if args.verbose > 1 && !args.quiet {
                eprintln!("{}: {e:#}", path.display());
            }
        } else if args.verbose > 1 && !args.quiet {
            eprintln!("{} is ok", path.display());
        }
    }

    Ok(if failed {
        ExitCode::FAILURE
    } else {
        ExitCode::SUCCESS
    })
}

fn check(path: &std::path::Path, unmapped: bool) -> anyhow::Result<()> {
    let src = path.to_str().context("path is not valid UTF-8")?;
    let input = Input::open(src)?;
    let format = input.format;

    let mut reader = AlignmentReader::from_input(input, &ReadOptions::default())?;
    let header = reader.read_header()?;

    if !unmapped && header.reference_sequences().is_empty() {
        anyhow::bail!("no @SQ lines in the header");
    }

    let mut file = File::open(path)?;
    let len = file.metadata()?.len();

    match format {
        Format::Bam => {
            if !bam_blocks::ends_with_eof(&mut file, len)? {
                anyhow::bail!("missing BGZF EOF marker");
            }
        }
        Format::Cram => {
            use std::io::{Read, Seek, SeekFrom};

            let eof_len = CRAM_3_EOF.len() as u64;
            let mut buf = [0; CRAM_3_EOF.len()];

            if len < eof_len {
                anyhow::bail!("file too short for a CRAM EOF container");
            }

            file.seek(SeekFrom::Start(len - eof_len))?;
            file.read_exact(&mut buf)?;

            if buf != CRAM_3_EOF {
                anyhow::bail!("missing CRAM EOF container");
            }
        }
        Format::Sam => {}
    }

    Ok(())
}
