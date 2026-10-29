//! Alignment writers for SAM, BAM, and CRAM output to a file or stdout.

use std::{
    fs::File,
    io::{self, BufWriter, Write},
    path::{Path, PathBuf},
};

use anyhow::{Context, bail};
use noodles_bam as bam;
use noodles_bgzf::io::writer::CompressionLevel;
use noodles_cram as cram;
use noodles_sam::{self as sam, alignment::io::Write as _};

use crate::{
    alignment::{Rec, reference_repository},
    bgzf_stream::BgzfWriter,
    parallel_sam::ParallelSamWriter,
    runs::write_raw_record,
};

const SAM_BUFFER_CAPACITY: usize = 1 << 17;

/// An output alignment format.
#[derive(Clone, Copy, Debug, Eq, PartialEq, clap::ValueEnum)]
pub enum OutputFormat {
    Sam,
    Bam,
    Cram,
}

impl OutputFormat {
    /// Infers a format from an output path's extension, as samtools does.
    pub fn from_path(path: &Path) -> Option<Self> {
        match path.extension()?.to_str()? {
            "sam" => Some(Self::Sam),
            "bam" => Some(Self::Bam),
            "cram" => Some(Self::Cram),
            _ => None,
        }
    }
}

/// Options for creating an alignment writer.
#[derive(Clone, Debug)]
pub struct WriteOptions {
    pub format: OutputFormat,
    /// The BGZF compression level for BAM output.
    pub level: CompressionLevel,
    /// Additional threads for BGZF compression, or for formatting SAM.
    pub threads: usize,
    /// A FASTA reference for CRAM encoding.
    pub reference: Option<PathBuf>,
}

pub type Sink = Box<dyn Write + Send>;

/// An open alignment writer.
pub enum AlignmentWriter {
    Sam(sam::io::Writer<BufWriter<Sink>>),
    /// SAM formatted on worker threads. Records that are not raw BAM are encoded first.
    ParallelSam {
        writer: ParallelSamWriter,
        encoder: bam::io::Writer<Vec<u8>>,
    },
    Bam(bam::io::Writer<BgzfWriter<Sink>>),
    Cram(Box<cram::io::Writer<Sink>>),
}

impl AlignmentWriter {
    /// Creates a writer to `dst`, or to stdout when `dst` is `None` or `-`, for records with
    /// reference sequence IDs relative to `header`.
    pub fn create(
        dst: Option<&Path>,
        options: &WriteOptions,
        header: &sam::Header,
    ) -> anyhow::Result<Self> {
        let sink: Sink = match dst {
            Some(path) if path != Path::new("-") => Box::new(
                File::create(path).with_context(|| format!("creating {}", path.display()))?,
            ),
            _ => raw_stdout()?,
        };

        Ok(match options.format {
            OutputFormat::Sam if options.threads > 0 => Self::ParallelSam {
                writer: ParallelSamWriter::new(sink, header, options.threads),
                encoder: bam::io::Writer::from(Vec::new()),
            },
            OutputFormat::Sam => Self::Sam(sam::io::Writer::new(BufWriter::with_capacity(
                SAM_BUFFER_CAPACITY,
                sink,
            ))),
            OutputFormat::Bam => Self::Bam(bam::io::Writer::from(BgzfWriter::new(
                sink,
                options.level,
                options.threads,
            ))),
            OutputFormat::Cram => {
                if options.reference.is_none() {
                    bail!("CRAM output requires a reference (-T/--reference)");
                }

                let repository = reference_repository(options.reference.as_deref())?;
                Self::Cram(Box::new(
                    cram::io::writer::Builder::default()
                        .set_reference_sequence_repository(repository)
                        .build_from_writer(sink),
                ))
            }
        })
    }

    pub fn write_header(&mut self, header: &sam::Header) -> io::Result<()> {
        match self {
            Self::Sam(w) => w.write_header(header),
            Self::ParallelSam { writer, .. } => {
                let mut text = sam::io::Writer::new(Vec::new());
                text.write_header(header)?;
                writer.write_text(text.into_inner())
            }
            Self::Bam(w) => w.write_header(header),
            Self::Cram(w) => w.write_header(header),
        }
    }

    /// Writes a record. BAM records going to BAM output are forwarded without re-encoding, so
    /// `header` must have the same reference sequences as the header the record was read with.
    pub fn write(&mut self, header: &sam::Header, record: &Rec<'_>) -> io::Result<()> {
        match (self, record) {
            (Self::Bam(w), Rec::Bam(r)) => write_raw_record(w.get_mut(), r.as_ref()),
            (Self::ParallelSam { writer, .. }, Rec::Bam(r)) => writer.push(r.as_ref()),
            (Self::ParallelSam { writer, encoder }, _) => {
                encoder.get_mut().clear();
                encoder.write_alignment_record(header, record.as_dyn())?;
                writer.push(&encoder.get_ref()[size_of::<u32>()..])
            }
            (Self::Sam(w), _) => w.write_alignment_record(header, record.as_dyn()),
            (Self::Bam(w), _) => w.write_alignment_record(header, record.as_dyn()),
            (Self::Cram(w), _) => w.write_alignment_record(header, record.as_dyn()),
        }
    }

    /// Writes an encoded BAM record (without `block_size`) whose reference sequence IDs are
    /// relative to `header`. BAM output copies it; other formats decode it first.
    pub fn write_raw(&mut self, header: &sam::Header, record: &[u8]) -> io::Result<()> {
        match self {
            Self::Bam(w) => write_raw_record(w.get_mut(), record),
            Self::ParallelSam { writer, .. } => writer.push(record),
            _ => {
                let record = bam::Record::try_from(record.to_vec())?;
                self.write(header, &Rec::Bam(&record))
            }
        }
    }

    /// Returns whether pre-formatted SAM text ([`Self::write_sam_text`]) or pre-encoded BAM
    /// records ([`Self::write_bam_records`]) can be written directly, and which.
    pub fn preformatted_kind(&self) -> Option<OutputFormat> {
        match self {
            Self::Sam(_) | Self::ParallelSam { .. } => Some(OutputFormat::Sam),
            Self::Bam(_) => Some(OutputFormat::Bam),
            Self::Cram(_) => None,
        }
    }

    /// Writes SAM record lines formatted elsewhere. Only valid for SAM output.
    pub fn write_sam_text(&mut self, text: Vec<u8>) -> io::Result<()> {
        match self {
            Self::Sam(w) => w.get_mut().write_all(&text),
            Self::ParallelSam { writer, .. } => writer.write_text(text),
            _ => Err(io::Error::new(
                io::ErrorKind::InvalidInput,
                "SAM text written to non-SAM output",
            )),
        }
    }

    /// Writes `block_size`-prefixed BAM records encoded elsewhere. Only valid for BAM output.
    pub fn write_bam_records(&mut self, records: &[u8]) -> io::Result<()> {
        match self {
            Self::Bam(w) => w.get_mut().write_all(records),
            _ => Err(io::Error::new(
                io::ErrorKind::InvalidInput,
                "BAM records written to non-BAM output",
            )),
        }
    }

    /// Flushes buffered records and writes any end-of-file marker.
    pub fn finish(self, header: &sam::Header) -> io::Result<()> {
        match self {
            Self::Sam(w) => w.into_inner().flush(),
            Self::ParallelSam { writer, .. } => writer.finish(),
            Self::Bam(w) => w.into_inner().finish()?.flush(),
            Self::Cram(mut w) => {
                w.finish(header)?;
                w.get_mut().flush()
            }
        }
    }
}

/// Returns stdout without Rust's line buffering, which scans every write for newlines; callers
/// buffer output themselves.
#[cfg(unix)]
pub fn raw_stdout() -> anyhow::Result<Sink> {
    use std::os::fd::AsFd;

    let fd = io::stdout()
        .as_fd()
        .try_clone_to_owned()
        .context("duplicating stdout")?;

    Ok(Box::new(File::from(fd)))
}

#[cfg(not(unix))]
pub fn raw_stdout() -> anyhow::Result<Sink> {
    Ok(Box::new(io::stdout()))
}

/// Parses a samtools compression level (`0`–`9`).
pub fn parse_level(level: u8) -> anyhow::Result<CompressionLevel> {
    CompressionLevel::new(level).with_context(|| format!("invalid compression level: {level}"))
}
