//! Format-specific alignment readers behind one streaming interface.
//!
//! Each format keeps its fastest native path: BAM records stay lazy and can be forwarded as raw
//! bytes, SAM records are parsed lazily, and CRAM records are decoded into record buffers.

use std::{
    io::{self, BufRead, BufReader},
    path::{Path, PathBuf},
};

use anyhow::{Context, bail};
use noodles_bam as bam;
use noodles_core::Region;
use noodles_cram as cram;
use noodles_fasta as fasta;
use noodles_sam::{
    self as sam,
    alignment::{RecordBuf, record::Flags},
};

use crate::{
    bgzf_stream::BgzfReader,
    source::{Compression, Format, Input, Stream},
};

/// The region name samtools uses for unplaced, unmapped records.
const UNMAPPED_REGION: &str = "*";

/// Options shared by every command that reads alignments.
#[derive(Clone, Debug, Default)]
pub struct ReadOptions {
    /// Additional decompression threads.
    pub threads: usize,
    /// A FASTA reference for CRAM decoding.
    pub reference: Option<PathBuf>,
}

/// A borrowed alignment record from any supported format.
pub enum Rec<'a> {
    Bam(&'a bam::Record),
    Sam(&'a sam::Record),
    Buf(&'a RecordBuf),
}

impl Rec<'_> {
    /// Returns the record as a format-independent trait object.
    pub fn as_dyn(&self) -> &dyn sam::alignment::Record {
        match self {
            Self::Bam(r) => *r,
            Self::Sam(r) => *r,
            Self::Buf(r) => *r,
        }
    }

    pub fn flags(&self) -> io::Result<Flags> {
        match self {
            Self::Bam(r) => Ok(r.flags()),
            Self::Sam(r) => r.flags(),
            Self::Buf(r) => Ok(r.flags()),
        }
    }

    /// Returns the mapping quality, with `255` (missing) as `None`.
    pub fn mapping_quality(&self) -> io::Result<Option<u8>> {
        match self {
            Self::Bam(r) => Ok(r.mapping_quality().map(u8::from)),
            Self::Sam(r) => r.mapping_quality().transpose().map(|q| q.map(u8::from)),
            Self::Buf(r) => Ok(r.mapping_quality().map(u8::from)),
        }
    }

    pub fn reference_sequence_id(&self, header: &sam::Header) -> io::Result<Option<usize>> {
        self.as_dyn().reference_sequence_id(header).transpose()
    }

    pub fn mate_reference_sequence_id(&self, header: &sam::Header) -> io::Result<Option<usize>> {
        self.as_dyn().mate_reference_sequence_id(header).transpose()
    }
}

type BamReader = bam::io::Reader<BgzfReader<BufReader<Stream>>>;

/// An open alignment reader for SAM, BAM, or CRAM.
pub enum AlignmentReader {
    Sam(sam::io::Reader<Box<dyn BufRead + Send>>),
    Bam(BamReader),
    Cram(cram::io::Reader<BufReader<Stream>>),
}

impl AlignmentReader {
    /// Opens `src` (`-` for stdin), detecting its format.
    pub fn open(src: &str, options: &ReadOptions) -> anyhow::Result<(Self, Option<PathBuf>)> {
        let input = Input::open(src)?;
        let path = input.path.clone();
        let reader = Self::from_input(input, options)?;
        Ok((reader, path))
    }

    pub fn from_input(input: Input, options: &ReadOptions) -> anyhow::Result<Self> {
        let Input {
            reader,
            format,
            compression,
            ..
        } = input;

        Ok(match (format, compression) {
            (Format::Bam, _) => Self::Bam(bam::io::Reader::from(BgzfReader::new(
                reader,
                options.threads,
            ))),
            (Format::Sam, Compression::Bgzf) => {
                let inner: Box<dyn BufRead + Send> =
                    Box::new(BgzfReader::new(reader, options.threads));
                Self::Sam(sam::io::Reader::new(inner))
            }
            (Format::Sam, Compression::None) => {
                let inner: Box<dyn BufRead + Send> = Box::new(reader);
                Self::Sam(sam::io::Reader::new(inner))
            }
            (Format::Cram, _) => {
                let repository = reference_repository(options.reference.as_deref())?;
                Self::Cram(
                    cram::io::reader::Builder::default()
                        .set_reference_sequence_repository(repository)
                        .build_from_reader(reader),
                )
            }
        })
    }

    pub fn format(&self) -> Format {
        match self {
            Self::Sam(_) => Format::Sam,
            Self::Bam(_) => Format::Bam,
            Self::Cram(_) => Format::Cram,
        }
    }

    pub fn read_header(&mut self) -> anyhow::Result<sam::Header> {
        match self {
            Self::Sam(r) => r.read_header(),
            Self::Bam(r) => r.read_header(),
            Self::Cram(r) => r.read_header(),
        }
        .context("reading header")
    }

    /// Calls `f` with every record in file order.
    pub fn for_each<F>(&mut self, header: &sam::Header, mut f: F) -> anyhow::Result<()>
    where
        F: FnMut(Rec<'_>) -> anyhow::Result<()>,
    {
        match self {
            Self::Bam(reader) => {
                let mut record = bam::Record::default();
                while reader.read_record(&mut record)? != 0 {
                    f(Rec::Bam(&record))?;
                }
            }
            Self::Sam(reader) => {
                let mut record = sam::Record::default();
                while reader.read_record(&mut record)? != 0 {
                    f(Rec::Sam(&record))?;
                }
            }
            Self::Cram(reader) => {
                for result in reader.records(header) {
                    let record = result?;
                    f(Rec::Buf(&record))?;
                }
            }
        }

        Ok(())
    }

    /// Calls `f` with the records overlapping each region, region by region, as samtools does
    /// without `-M`. A record overlapping several regions is visited once per region.
    pub fn for_each_in_regions<F>(
        &mut self,
        path: &Path,
        header: &sam::Header,
        regions: &[String],
        mut f: F,
    ) -> anyhow::Result<()>
    where
        F: FnMut(Rec<'_>) -> anyhow::Result<()>,
    {
        match self {
            Self::Bam(reader) => {
                let index = bam::fs::read_associated_index(path)
                    .with_context(|| format!("reading index for {}", path.display()))?;

                for raw_region in regions {
                    if raw_region == UNMAPPED_REGION {
                        for result in reader.query_unmapped(&index)? {
                            let record = result?;
                            f(Rec::Bam(&record))?;
                        }
                        continue;
                    }

                    let region = parse_region(raw_region)?;
                    let mut query = reader.query(header, &index, &region)?;
                    let mut record = bam::Record::default();

                    while query.read_record(&mut record)? != 0 {
                        f(Rec::Bam(&record))?;
                    }
                }
            }
            Self::Cram(reader) => {
                let index = cram::fs::read_associated_index(path)
                    .with_context(|| format!("reading index for {}", path.display()))?;

                for raw_region in regions {
                    if raw_region == UNMAPPED_REGION {
                        for result in reader.query_unmapped(header, &index)? {
                            let record = result?;
                            f(Rec::Buf(&record))?;
                        }
                        continue;
                    }

                    let region = parse_region(raw_region)?;
                    let mut query = reader.query(header, &index, &region)?;
                    let mut record = RecordBuf::default();

                    while query.read_record_buf(&mut record)? != 0 {
                        f(Rec::Buf(&record))?;
                    }
                }
            }
            Self::Sam(_) => bail!("region queries on SAM input are not supported; use BAM or CRAM"),
        }

        Ok(())
    }
}

fn parse_region(s: &str) -> anyhow::Result<Region> {
    s.parse().with_context(|| format!("invalid region: {s}"))
}

/// Builds a reference sequence repository from an indexed FASTA, or an empty one.
pub fn reference_repository(src: Option<&Path>) -> anyhow::Result<fasta::Repository> {
    let Some(src) = src else {
        return Ok(fasta::Repository::default());
    };

    let reader = fasta::io::indexed_reader::Builder::default()
        .build_from_path(src)
        .with_context(|| format!("opening reference {}", src.display()))?;

    Ok(fasta::Repository::new(
        fasta::repository::adapters::IndexedReader::new(reader),
    ))
}
