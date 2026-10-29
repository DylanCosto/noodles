//! `ntools cat`: concatenate BAM files without recompressing records (samtools cat).

use std::{
    fs::File,
    io::{BufWriter, Write},
    path::{Path, PathBuf},
};

use anyhow::{Context, bail};

use crate::{
    bam_blocks,
    commands::reheader::read_sam_header,
    output, program,
    source::{Format, Input},
};

#[derive(Debug, clap::Args)]
#[command(disable_help_flag = true)]
pub struct Args {
    /// Input BAM files, concatenated in order.
    #[arg(required = true)]
    inputs: Vec<PathBuf>,

    /// Use the header lines of this SAM file instead of the first input's header.
    #[arg(short = 'h', long = "header")]
    header: Option<PathBuf>,

    /// Write output to this file instead of stdout.
    #[arg(short = 'o', long = "output")]
    output: Option<PathBuf>,

    /// Do not add an @PG line to the header.
    #[arg(long = "no-PG")]
    no_pg: bool,

    /// Print help.
    #[arg(long, action = clap::ArgAction::Help)]
    help: Option<bool>,
}

pub fn run(args: Args) -> anyhow::Result<()> {
    for path in &args.inputs {
        let src = path.to_str().context("input path is not valid UTF-8")?;

        if Input::open(src)?.format != Format::Bam {
            bail!("{src}: only BAM files can be concatenated");
        }
    }

    let sink: output::Sink = match args.output.as_deref().filter(|p| *p != Path::new("-")) {
        Some(path) => {
            Box::new(File::create(path).with_context(|| format!("creating {}", path.display()))?)
        }
        None => output::raw_stdout()?,
    };

    let mut dst = BufWriter::new(sink);

    for (i, path) in args.inputs.iter().enumerate() {
        let (header, body) = bam_blocks::split(path)?;

        if i == 0 {
            let mut header = match &args.header {
                Some(path) => read_sam_header(path)?,
                None => header,
            };

            if !args.no_pg {
                program::add(&mut header)?;
            }

            bam_blocks::write_header(&mut dst, &header)?;
        }

        body.copy_to(&mut dst)?;
    }

    bam_blocks::write_eof(&mut dst)?;
    dst.flush()?;

    Ok(())
}
