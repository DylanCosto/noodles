//! `ntools cat`: concatenate BAM or CRAM files without recompressing records (samtools cat).

mod cram;

use std::{
    fs::File,
    io::{BufWriter, Write},
    path::{Path, PathBuf},
};

use anyhow::{Context, bail};
use noodles_sam as sam;

use crate::{
    bam_blocks,
    commands::{merge::read_file_list, reheader::read_sam_header},
    output, program,
    source::{Format, Input},
};

#[derive(Debug, clap::Args)]
#[command(disable_help_flag = true)]
pub struct Args {
    /// Input BAM or CRAM files, concatenated in order after those listed with -b.
    inputs: Vec<PathBuf>,

    /// A file listing input files, one per line (repeatable).
    #[arg(short = 'b', value_name = "FILE")]
    lists: Vec<PathBuf>,

    /// Copy the header from this file instead of the first input's.
    #[arg(short = 'h', long = "header", value_name = "FILE")]
    header: Option<PathBuf>,

    /// Write output to this file instead of stdout.
    #[arg(short = 'o', long = "output")]
    output: Option<PathBuf>,

    /// CRAM only: keep only REG, a region (`chr:start-end`, `*`) or a range of container
    /// numbers (`#:start-end`, from 0). Needs an index.
    #[arg(short = 'r', value_name = "REG")]
    region: Option<String>,

    /// CRAM only: keep part N of M (1 to M) of the file's containers, or of the region's.
    /// Needs an index.
    #[arg(short = 'p', value_name = "N/M")]
    part: Option<String>,

    /// CRAM only: copy containers that overlap the region whole, without filtering their
    /// records.
    #[arg(short = 'f')]
    fast: bool,

    /// CRAM only: print the number of indexed containers (in the region) instead of
    /// concatenating.
    #[arg(short = 'q')]
    query: bool,

    /// Reference FASTA (indexed) for re-encoding CRAM containers that are filtered or need new
    /// read group numbers [default: the files named by the header's @SQ UR fields].
    #[arg(long = "reference")]
    reference: Option<PathBuf>,

    /// Do not add an @PG line to the header.
    #[arg(long = "no-PG")]
    no_pg: bool,

    /// Print help.
    #[arg(long, action = clap::ArgAction::Help)]
    help: Option<bool>,
}

pub fn run(args: Args) -> anyhow::Result<()> {
    let mut inputs = Vec::new();

    for list in &args.lists {
        inputs.extend(read_file_list(list)?);
    }

    inputs.extend(args.inputs.iter().cloned());

    let Some(first) = inputs.first() else {
        bail!("Usage: ntools cat [options] <in1.bam> [... <inN.bam>]");
    };

    let format = input_format(first)?;

    if format != Format::Bam && format != Format::Cram {
        bail!("input is not BAM or CRAM");
    }

    for path in &inputs[1..] {
        if input_format(path)? != format {
            bail!("File {} is of different format!", path.display());
        }
    }

    if format == Format::Cram && args.query {
        return cram::print_container_counts(&inputs, args.region.as_deref());
    }

    let header = args.header.as_deref().map(read_sam_header).transpose()?;

    let sink: output::Sink = match args.output.as_deref().filter(|p| *p != Path::new("-")) {
        Some(path) => {
            Box::new(File::create(path).with_context(|| format!("creating {}", path.display()))?)
        }
        None => output::raw_stdout()?,
    };

    let mut dst = BufWriter::new(sink);

    match format {
        Format::Cram => {
            let options = cram::Options {
                region: args.region.as_deref(),
                part: args.part.as_deref(),
                filtering: if args.fast {
                    cram::Filtering::Fast
                } else {
                    cram::Filtering::Exact
                },
                reference: args.reference.as_deref(),
            };

            let no_pg = args.no_pg;
            cram::cat(&inputs, header, &options, &mut dst, |header| {
                if no_pg { Ok(()) } else { program::add(header) }
            })?;
        }
        _ => cat_bam(&inputs, header, args.no_pg, &mut dst)?,
    }

    dst.flush()?;

    Ok(())
}

fn input_format(path: &Path) -> anyhow::Result<Format> {
    let src = path.to_str().context("input path is not valid UTF-8")?;
    Ok(Input::open(src)?.format)
}

fn cat_bam<W: Write>(
    inputs: &[PathBuf],
    header: Option<sam::Header>,
    no_pg: bool,
    dst: &mut W,
) -> anyhow::Result<()> {
    let mut merged = header;
    let mut bodies = Vec::with_capacity(inputs.len());

    for path in inputs {
        let (header, body) = bam_blocks::split(path)?;
        merge_read_groups(&mut merged, header);
        bodies.push(body);
    }

    let mut header = merged.context("no input headers")?;

    if !no_pg {
        program::add(&mut header)?;
    }

    bam_blocks::write_header(dst, &header)?;

    for body in bodies {
        body.copy_to(dst)?;
    }

    bam_blocks::write_eof(dst)?;

    Ok(())
}

/// samtools' `cat_check_merge_hdr`: the first header (or the `-h` header) gains the `@RG` lines
/// of later headers whose IDs it lacks.
fn merge_read_groups(merged: &mut Option<sam::Header>, header: sam::Header) {
    match merged {
        None => *merged = Some(header),
        Some(merged) => {
            for (id, read_group) in header.read_groups() {
                if !merged.read_groups().contains_key(id) {
                    merged
                        .read_groups_mut()
                        .insert(id.clone(), read_group.clone());
                }
            }
        }
    }
}
