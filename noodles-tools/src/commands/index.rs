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
use noodles_csi::{
    self as csi,
    binning_index::index::reference_sequence::{bin::Chunk, index::BinnedIndex},
};
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
    /// Coordinate-sorted BAM or CRAM files to index (several with `-M`); without `-M`, a second
    /// argument that does not exist or is an index names the index to write.
    #[arg(required = true)]
    inputs: Vec<PathBuf>,

    /// Write the index to this path instead of next to the input.
    #[arg(short = 'o', long = "output")]
    output: Option<PathBuf>,

    /// Build a BAI index (the default for BAM).
    #[arg(short = 'b', long = "bai", overrides_with_all = ["csi", "min_shift"])]
    bai: bool,

    /// Build a CSI index.
    #[arg(short = 'c', long = "csi", overrides_with = "bai")]
    csi: bool,

    /// Build a CSI index with intervals of 2^INT bases at its finest level [14].
    #[arg(
        short = 'm',
        long = "min-shift",
        overrides_with = "bai",
        allow_negative_numbers = true
    )]
    min_shift: Option<i32>,

    /// Index every file argument.
    #[arg(short = 'M')]
    multiple: bool,

    /// Additional threads for BGZF decompression.
    #[arg(short = '@', long = "threads", default_value_t = 0)]
    threads: usize,
}

/// HTSlib's default CSI minimum shift (`BAM_LIDX_SHIFT`).
const DEFAULT_MIN_SHIFT: i32 = 14;

pub fn run(args: Args) -> anyhow::Result<()> {
    let csi = args.csi || args.min_shift.is_some();
    let min_shift = args.min_shift.unwrap_or(DEFAULT_MIN_SHIFT);

    // As in HTSlib, a minimum shift of 0 or less means BAI.
    let kind = if csi && min_shift > 0 {
        IndexKind::CsiWithMinShift(u8::try_from(min_shift).context("invalid min_shift")?)
    } else {
        IndexKind::Default
    };

    let mut inputs = args.inputs;
    let mut output = args.output;

    // The legacy synopsis: `index in.bam out.index`.
    if inputs.len() == 2 && output.is_none() && is_missing_or_index(&inputs[1]) {
        output = inputs.pop();
    }

    if inputs.len() > 1 && !args.multiple {
        bail!("use -M to enable indexing more than one alignment file");
    }

    if output.is_some() && inputs.len() > 1 {
        bail!("can't use -o with multiple input alignment files");
    }

    for input in &inputs {
        index_file(input, output.clone(), kind, args.threads)?;
    }

    Ok(())
}

/// Whether a path does not exist or holds an index (BAI, CSI, TBI, or CRAI), as HTSlib's
/// `nonexistent_or_index` decides whether a second argument names the index to write.
fn is_missing_or_index(path: &Path) -> bool {
    use std::io::Read as _;

    let Ok(mut file) = File::open(path) else {
        return true;
    };

    let mut buf = Vec::new();

    if (&mut file).take(4).read_to_end(&mut buf).is_err() {
        return false;
    }

    if buf == b"BAI\x01" {
        return true;
    }

    if !buf.starts_with(&[0x1f, 0x8b]) {
        return false;
    }

    // Compressed: CSI and TBI have magic numbers; a CRAI is lines of six integers.
    let Ok(file) = File::open(path) else {
        return false;
    };

    let mut text = Vec::new();
    let mut decoder = flate2::read::MultiGzDecoder::new(file).take(256);

    if decoder.read_to_end(&mut text).is_err() && text.is_empty() {
        return false;
    }

    text.starts_with(b"CSI\x01") || text.starts_with(b"TBI\x01") || is_crai_text(&text)
}

/// Whether text starts like a CRAI index: tab-separated integers on each line.
fn is_crai_text(text: &[u8]) -> bool {
    let Some(line) = text.split(|&b| b == b'\n').next() else {
        return false;
    };

    let fields: Vec<&[u8]> = line.split(|&b| b == b'\t').collect();

    fields.len() == 6
        && fields.iter().all(|field| {
            let digits = field.strip_prefix(b"-").unwrap_or(field);
            !digits.is_empty() && digits.iter().all(u8::is_ascii_digit)
        })
}

/// Which index to build for BAM; CRAM always gets CRAI.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub enum IndexKind {
    /// BAI for BAM.
    Default,
    /// CSI with HTSlib's default minimum shift.
    Csi,
    /// CSI with intervals of 2^n bases at its finest level.
    CsiWithMinShift(u8),
}

/// Indexes a BAM or CRAM file, writing next to it (`in.bam.bai`) unless `output` is given.
pub fn index_file(
    input: &Path,
    output: Option<PathBuf>,
    kind: IndexKind,
    threads: usize,
) -> anyhow::Result<()> {
    let src = input.to_str().context("input path is not valid UTF-8")?;
    let opened = Input::open(src)?;

    match opened.format {
        Format::Bam => {
            let (min_shift, extension) = match kind {
                IndexKind::Default => (None, "bai"),
                IndexKind::Csi => (Some(DEFAULT_MIN_SHIFT as u8), "csi"),
                IndexKind::CsiWithMinShift(min_shift) => (Some(min_shift), "csi"),
            };

            let dst = output.unwrap_or_else(|| with_appended_extension(input, extension));

            let index = index_bam(opened.reader, min_shift, threads)
                .with_context(|| format!("indexing {}", input.display()))?;

            write_bam_index(&dst, &index)
        }
        // As in samtools, CSI options do not apply to CRAM.
        Format::Cram => {
            let dst = output.unwrap_or_else(|| with_appended_extension(input, "crai"));

            let index =
                cram::fs::index(input).with_context(|| format!("indexing {}", input.display()))?;

            let file = File::create(&dst).with_context(|| format!("creating {}", dst.display()))?;
            let mut writer = crai::io::Writer::new(file);
            writer.write_index(&index)?;
            writer.finish()?;

            Ok(())
        }
        Format::Sam => bail!("only BAM and CRAM files can be indexed"),
    }
}

/// HTSlib's CSI layout for a longest reference sequence (`hts_adjust_csi_settings`): enough
/// levels above `min_shift` to cover it, or, past 9 levels, a larger minimum shift. Returns the
/// minimum shift and depth.
fn csi_layout(max_length: usize, min_shift: u8) -> (u8, u8) {
    const MAX_DEPTH: u8 = 9;

    let max_length = u128::try_from(max_length).unwrap_or(u128::MAX) + 256;
    let max_position = |min_shift: u8, depth: u8| {
        1u128
            .checked_shl(u32::from(min_shift) + 3 * u32::from(depth))
            .unwrap_or(u128::MAX)
    };

    if max_length <= max_position(min_shift, MAX_DEPTH) {
        let mut depth = 0;

        while max_length > max_position(min_shift, depth) {
            depth += 1;
        }

        (min_shift, depth)
    } else {
        let mut shift = min_shift;

        while max_length > max_position(shift, MAX_DEPTH) {
            shift += 1;
        }

        eprintln!(
            "[W::hts_adjust_csi_settings] Adjusted min_shift from {min_shift} to {shift} due to \
             longest reference of {} bases.",
            max_length - 256
        );

        (shift, MAX_DEPTH)
    }
}

/// A BAI or CSI indexer.
enum BamIndexer {
    Bai(Indexer),
    Csi(csi::binning_index::Indexer<BinnedIndex>),
}

impl BamIndexer {
    fn add_record(
        &mut self,
        context: Option<(usize, Position, Position, bool)>,
        chunk: Chunk,
    ) -> anyhow::Result<()> {
        match self {
            Self::Bai(indexer) => indexer.add_record(context, chunk)?,
            Self::Csi(indexer) => indexer.add_record(context, chunk)?,
        }

        Ok(())
    }

    fn build(self, reference_sequence_count: usize) -> bam::Index {
        match self {
            Self::Bai(indexer) => indexer.build(reference_sequence_count),
            Self::Csi(indexer) => bam::Index::Csi(indexer.build(reference_sequence_count)),
        }
    }
}

/// Builds a BAI, or with a minimum shift, a CSI.
fn index_bam(
    reader: BufReader<Stream>,
    min_shift: Option<u8>,
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

    let mut indexer = match min_shift {
        Some(min_shift) => {
            let (min_shift, depth) = csi_layout(max_reference_sequence_length.get(), min_shift);
            let indexer = csi::binning_index::Indexer::new(min_shift, depth)
                .with_context(|| format!("invalid CSI layout: min_shift {min_shift}"))?;
            BamIndexer::Csi(indexer)
        }
        None => BamIndexer::Bai(
            Indexer::builder()
                .set_format(bam::index::Format::Bai)
                .set_max_position_hint(max_position)
                .build()?,
        ),
    };

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
