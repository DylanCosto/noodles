//! `ntools index`: build BAI, CSI, or CRAI indexes (samtools index).

use std::{
    fs::File,
    io::{BufReader, BufWriter, Write},
    num::NonZero,
    path::{Path, PathBuf},
};

use anyhow::{Context, bail};
use noodles_bam::{self as bam, bai, index::Indexer};
use noodles_bgzf as bgzf;
use noodles_core::Position;
use noodles_cram::{self as cram, crai};
use noodles_csi::{self as csi, binning_index::index::reference_sequence::bin::Chunk};
use noodles_sam::{
    self as sam,
    alignment::Record as _,
    header::record::value::map::header::{sort_order::COORDINATE, tag::SORT_ORDER},
};

use crate::{
    bgzf_stream::BgzfReader,
    source::{Format, Input, Stream},
};

#[derive(Debug, clap::Args)]
pub struct Args {
    /// Coordinate-sorted BAM or CRAM file to index.
    input: PathBuf,

    /// Write the index to this path instead of next to the input.
    #[arg(short = 'o', long = "output")]
    output: Option<PathBuf>,

    /// Build a BAI index (the default for BAM).
    #[arg(short = 'b', long = "bai", conflicts_with = "csi")]
    bai: bool,

    /// Build a CSI index.
    #[arg(short = 'c', long = "csi")]
    csi: bool,

    /// Additional threads for BGZF decompression.
    #[arg(short = '@', long = "threads", default_value_t = 0)]
    threads: usize,
}

pub fn run(args: Args) -> anyhow::Result<()> {
    let src = args
        .input
        .to_str()
        .context("input path is not valid UTF-8")?;
    let input = Input::open(src)?;

    match input.format {
        Format::Bam => {
            let format = if args.csi {
                bam::index::Format::Csi
            } else {
                bam::index::Format::Bai
            };

            let extension = if args.csi { "csi" } else { "bai" };
            let dst = args
                .output
                .unwrap_or_else(|| with_appended_extension(&args.input, extension));

            let index = index_bam(input.reader, format, args.threads)
                .with_context(|| format!("indexing {}", args.input.display()))?;

            write_bam_index(&dst, &index)
        }
        Format::Cram => {
            let dst = args
                .output
                .unwrap_or_else(|| with_appended_extension(&args.input, "crai"));

            if args.csi {
                bail!("CRAM files are indexed with CRAI, not CSI");
            }

            let index = cram::fs::index(&args.input)
                .with_context(|| format!("indexing {}", args.input.display()))?;

            let file = File::create(&dst).with_context(|| format!("creating {}", dst.display()))?;
            let mut writer = crai::io::Writer::new(file);
            writer.write_index(&index)?;
            writer.finish()?;

            Ok(())
        }
        Format::Sam => bail!("only BAM and CRAM files can be indexed"),
    }
}

fn index_bam(
    reader: BufReader<Stream>,
    format: bam::index::Format,
    threads: usize,
) -> anyhow::Result<bam::Index> {
    let mut reader = bam::io::Reader::from(BgzfReader::new(reader, threads));
    let header = reader.read_header()?;

    if !is_coordinate_sorted(&header) {
        bail!("the input is not coordinate-sorted (@HD SO:coordinate)");
    }

    let max_reference_sequence_length = header
        .reference_sequences()
        .values()
        .map(|rs| rs.length())
        .max()
        .unwrap_or(NonZero::<usize>::MIN);

    let max_position = Position::new(max_reference_sequence_length.get())
        .context("reference sequence length overflows a position")?;

    let mut indexer = Indexer::builder()
        .set_format(format)
        .set_max_position_hint(max_position)
        .build()?;

    let mut record = bam::Record::default();
    let mut start = bgzf::io::Read::virtual_position(reader.get_ref());

    while reader.read_record(&mut record)? != 0 {
        let end = bgzf::io::Read::virtual_position(reader.get_ref());

        let context = alignment_context(&record)?
            .map(|(id, start, end)| (id, start, end, !record.flags().is_unmapped()));

        indexer.add_record(context, Chunk::new(start, end))?;

        start = end;
    }

    Ok(indexer.build(header.reference_sequences().len()))
}

fn write_bam_index(dst: &Path, index: &bam::Index) -> anyhow::Result<()> {
    let file = File::create(dst).with_context(|| format!("creating {}", dst.display()))?;

    match index {
        bam::Index::Bai(index) => {
            let mut writer = bai::io::Writer::new(BufWriter::new(file));
            writer.write_index(index)?;
            writer.get_mut().flush()?;
        }
        bam::Index::Csi(index) => {
            let mut writer = csi::io::Writer::new(file);
            writer.write_index(index)?;
            writer.get_mut().try_finish()?;
        }
    }

    Ok(())
}

fn alignment_context(record: &bam::Record) -> anyhow::Result<Option<(usize, Position, Position)>> {
    match (
        record.reference_sequence_id().transpose()?,
        record.alignment_start().transpose()?,
        record.alignment_end().transpose()?,
    ) {
        (Some(id), Some(start), Some(end)) => Ok(Some((id, start, end))),
        _ => Ok(None),
    }
}

fn is_coordinate_sorted(header: &sam::Header) -> bool {
    header
        .header()
        .and_then(|hdr| hdr.other_fields().get(&SORT_ORDER))
        .is_some_and(|sort_order| sort_order == COORDINATE)
}

/// Appends an extension, as samtools names indexes: `in.bam` becomes `in.bam.bai`.
fn with_appended_extension(path: &Path, extension: &str) -> PathBuf {
    let mut s = path.as_os_str().to_owned();
    s.push(".");
    s.push(extension);
    PathBuf::from(s)
}
