//! `ntools reheader`: replace a BAM file's header without recompressing records
//! (samtools reheader).

use std::{
    fs::File,
    io::{BufReader, BufWriter, Write},
    path::PathBuf,
};

use anyhow::{Context, bail};
use noodles_sam as sam;

use crate::{
    bam_blocks, output, program,
    source::{Format, Input},
};

#[derive(Debug, clap::Args)]
pub struct Args {
    /// SAM file whose header lines become the new header.
    header: PathBuf,

    /// Input BAM file.
    input: PathBuf,

    /// Do not add an @PG line to the header.
    #[arg(short = 'P', long = "no-PG")]
    no_pg: bool,
}

pub fn run(args: Args) -> anyhow::Result<()> {
    let mut header = read_sam_header(&args.header)?;

    let src = args
        .input
        .to_str()
        .context("input path is not valid UTF-8")?;

    if Input::open(src)?.format != Format::Bam {
        bail!("{src}: only BAM files can be reheadered");
    }

    if !args.no_pg {
        program::add(&mut header)?;
    }

    let (_, body) = bam_blocks::split(&args.input)?;

    let mut dst = BufWriter::new(output::raw_stdout()?);
    bam_blocks::write_header(&mut dst, &header)?;
    body.copy_to(&mut dst)?;
    bam_blocks::write_eof(&mut dst)?;
    dst.flush()?;

    Ok(())
}

/// Reads the header lines of a SAM file.
pub fn read_sam_header(path: &std::path::Path) -> anyhow::Result<sam::Header> {
    let file = File::open(path).with_context(|| format!("opening {}", path.display()))?;
    sam::io::Reader::new(BufReader::new(file))
        .read_header()
        .with_context(|| format!("reading header from {}", path.display()))
}
