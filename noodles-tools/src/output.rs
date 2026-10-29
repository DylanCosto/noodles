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

use crate::{alignment::Rec, alignment::reference_repository, bgzf_stream::BgzfWriter};

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
    /// Additional compression threads.
    pub threads: usize,
    /// A FASTA reference for CRAM encoding.
    pub reference: Option<PathBuf>,
}

type Sink = Box<dyn Write + Send>;

/// An open alignment writer.
pub enum AlignmentWriter {
    Sam(sam::io::Writer<BufWriter<Sink>>),
    Bam(bam::io::Writer<BgzfWriter<Sink>>),
    Cram(Box<cram::io::Writer<Sink>>),
}

impl AlignmentWriter {
    /// Creates a writer to `dst`, or to stdout when `dst` is `None` or `-`.
    pub fn create(dst: Option<&Path>, options: &WriteOptions) -> anyhow::Result<Self> {
        let sink: Sink = match dst {
            Some(path) if path != Path::new("-") => Box::new(
                File::create(path).with_context(|| format!("creating {}", path.display()))?,
            ),
            _ => raw_stdout()?,
        };

        Ok(match options.format {
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
            Self::Bam(w) => w.write_header(header),
            Self::Cram(w) => w.write_header(header),
        }
    }

    /// Writes a record. BAM records going to BAM output are forwarded without re-encoding, so
    /// `header` must have the same reference sequences as the header the record was read with.
    pub fn write(&mut self, header: &sam::Header, record: &Rec<'_>) -> io::Result<()> {
        match (self, record) {
            (Self::Bam(w), Rec::Bam(r)) => crate::runs::write_raw_record(w.get_mut(), r.as_ref()),
            (Self::Sam(w), _) => w.write_alignment_record(header, record.as_dyn()),
            (Self::Bam(w), _) => w.write_alignment_record(header, record.as_dyn()),
            (Self::Cram(w), _) => w.write_alignment_record(header, record.as_dyn()),
        }
    }

    /// Writes an encoded BAM record (without `block_size`) whose reference sequence IDs are
    /// relative to `header`. BAM output copies it; other formats decode it first.
    pub fn write_raw(&mut self, header: &sam::Header, record: &[u8]) -> io::Result<()> {
        match self {
            Self::Bam(w) => crate::runs::write_raw_record(w.get_mut(), record),
            _ => {
                let record = bam::Record::try_from(record.to_vec())?;
                self.write(header, &Rec::Bam(&record))
            }
        }
    }

    /// Flushes buffered records and writes any end-of-file marker.
    pub fn finish(self, header: &sam::Header) -> io::Result<()> {
        match self {
            Self::Sam(w) => w.into_inner().flush(),
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
fn raw_stdout() -> anyhow::Result<Sink> {
    use std::os::fd::AsFd;

    let fd = io::stdout()
        .as_fd()
        .try_clone_to_owned()
        .context("duplicating stdout")?;

    Ok(Box::new(File::from(fd)))
}

#[cfg(not(unix))]
fn raw_stdout() -> anyhow::Result<Sink> {
    Ok(Box::new(io::stdout()))
}

/// Parses a samtools compression level (`0`–`9`).
pub fn parse_level(level: u8) -> anyhow::Result<CompressionLevel> {
    CompressionLevel::new(level).with_context(|| format!("invalid compression level: {level}"))
}
