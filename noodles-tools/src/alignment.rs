//! Format-specific alignment readers behind one streaming interface.
//!
//! Each format keeps its fastest native path: BAM records stay lazy and can be forwarded as raw
//! bytes, SAM records are parsed lazily, and CRAM records are decoded into record buffers.

use std::{
    io::{self, BufRead, BufReader},
    mem,
    path::{Path, PathBuf},
};

use anyhow::{Context, bail};
use noodles_bam as bam;
use noodles_core::{Position, Region};
use noodles_cram as cram;
use noodles_fasta as fasta;
use noodles_sam::{
    self as sam,
    alignment::{RecordBuf, record::Flags},
};

use crate::{
    bgzf_stream::BgzfReader,
    regions::{RegionSet, record_interval, region_string},
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
    /// Any other record, e.g., a CRAM record decoded on a worker.
    Dyn(&'a dyn sam::alignment::Record),
}

impl Rec<'_> {
    /// Returns the record as a format-independent trait object.
    pub fn as_dyn(&self) -> &dyn sam::alignment::Record {
        match self {
            Self::Bam(r) => *r,
            Self::Sam(r) => *r,
            Self::Buf(r) => *r,
            Self::Dyn(r) => *r,
        }
    }

    pub fn flags(&self) -> io::Result<Flags> {
        match self {
            Self::Bam(r) => Ok(r.flags()),
            Self::Sam(r) => r.flags(),
            Self::Buf(r) => Ok(r.flags()),
            Self::Dyn(r) => r.flags(),
        }
    }

    /// Returns the mapping quality, with `255` (missing) as `None`.
    pub fn mapping_quality(&self) -> io::Result<Option<u8>> {
        match self {
            Self::Bam(r) => Ok(r.mapping_quality().map(u8::from)),
            Self::Sam(r) => r.mapping_quality().transpose().map(|q| q.map(u8::from)),
            Self::Buf(r) => Ok(r.mapping_quality().map(u8::from)),
            Self::Dyn(r) => r.mapping_quality().transpose().map(|q| q.map(u8::from)),
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
    /// A CRAM reader, the reference repository it decodes with, and where the reference comes
    /// from.
    Cram(
        cram::io::Reader<BufReader<Stream>>,
        fasta::Repository,
        ReferenceSource,
    ),
}

/// Where a CRAM reader's reference sequences come from.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub enum ReferenceSource {
    /// A FASTA file given with `--reference`.
    Given,
    /// The FASTA file named by the header's `@SQ` `UR` fields, as HTSlib finds it.
    Header,
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
                let source = if options.reference.is_some() {
                    ReferenceSource::Given
                } else {
                    ReferenceSource::Header
                };

                Self::Cram(
                    cram::io::reader::Builder::default()
                        .set_reference_sequence_repository(repository.clone())
                        .build_from_reader(reader),
                    repository,
                    source,
                )
            }
        })
    }

    pub fn format(&self) -> Format {
        match self {
            Self::Sam(_) => Format::Sam,
            Self::Bam(_) => Format::Bam,
            Self::Cram(..) => Format::Cram,
        }
    }

    pub fn read_header(&mut self) -> anyhow::Result<sam::Header> {
        let header = match self {
            Self::Sam(r) => r.read_header(),
            Self::Bam(r) => r.read_header(),
            Self::Cram(r, _, _) => r.read_header(),
        }
        .context("reading header")?;

        if let Self::Cram(r, repository, ReferenceSource::Header) = self
            && let Some(path) = header_reference(&header)
        {
            *repository = reference_repository(Some(&path))?;
            r.set_reference_sequence_repository(repository.clone());
        }

        Ok(header)
    }

    /// Calls `f` with every record in file order.
    pub fn for_each<F>(&mut self, header: &sam::Header, mut f: F) -> anyhow::Result<()>
    where
        F: FnMut(Rec<'_>) -> anyhow::Result<()>,
    {
        self.for_each_while(header, |rec| f(rec).map(|()| true))
    }

    /// Calls `f` with records in file order until it returns `false`.
    pub fn for_each_while<F>(&mut self, header: &sam::Header, mut f: F) -> anyhow::Result<()>
    where
        F: FnMut(Rec<'_>) -> anyhow::Result<bool>,
    {
        match self {
            Self::Bam(reader) => {
                let mut record = bam::Record::default();
                while reader.read_record(&mut record)? != 0 {
                    if !f(Rec::Bam(&record))? {
                        return Ok(());
                    }
                }
            }
            Self::Sam(reader) => {
                let mut record = sam::Record::default();
                while reader.read_record(&mut record)? != 0 {
                    if !f(Rec::Sam(&record))? {
                        return Ok(());
                    }
                }
            }
            Self::Cram(reader, repository, _) => {
                // Decoding containers directly avoids converting each record to a `RecordBuf`, and
                // encoding records as BAM lets callers use their raw BAM paths.
                let mut container = cram::io::reader::Container::default();
                let mut bam_record = bam::Record::default();

                while reader.read_container(&mut container)? != 0 {
                    let compression_header = container.compression_header()?;

                    for result in container.slices() {
                        let slice = result?;
                        let (core_data_src, external_data_srcs) = slice.decode_blocks()?;

                        let records = slice.records(
                            repository.clone(),
                            header,
                            &compression_header,
                            &core_data_src,
                            &external_data_srcs,
                        )?;

                        for record in &records {
                            encode_cram_record(header, record, &mut bam_record)?;
                            if !f(Rec::Bam(&bam_record))? {
                                return Ok(());
                            }
                        }
                    }
                }
            }
        }

        Ok(())
    }

    /// Calls `f` with the records of each indexed query, in query order.
    ///
    /// `index` names an index file explicitly (`-X`); otherwise the index next to `path` is used.
    pub fn for_each_in_queries<F>(
        &mut self,
        path: &Path,
        index: Option<&Path>,
        header: &sam::Header,
        queries: &Queries<'_>,
        f: F,
    ) -> anyhow::Result<()>
    where
        F: FnMut(Rec<'_>) -> anyhow::Result<()>,
    {
        let index = AlignmentIndex::read(self.format(), path, index)?;
        self.for_each_in_indexed_queries(&index, header, queries, f)
    }

    /// Like [`Self::for_each_in_queries`], with an index already read, for callers that query
    /// the same file many times.
    pub fn for_each_in_indexed_queries<F>(
        &mut self,
        index: &AlignmentIndex,
        header: &sam::Header,
        queries: &Queries<'_>,
        f: F,
    ) -> anyhow::Result<()>
    where
        F: FnMut(Rec<'_>) -> anyhow::Result<()>,
    {
        let queries = queries.to_list(header)?;
        self.run_queries(index, header, &queries, f)
    }

    /// Calls `f` with the records overlapping the 0-based, half-open interval `[start, end)` of
    /// one reference sequence, as `sam_itr_queryi` selects them. An `end` of `None` means the
    /// end of the reference sequence.
    pub fn for_each_in_interval<F>(
        &mut self,
        index: &AlignmentIndex,
        header: &sam::Header,
        reference_sequence_id: usize,
        start: usize,
        end: Option<usize>,
        f: F,
    ) -> anyhow::Result<()>
    where
        F: FnMut(Rec<'_>) -> anyhow::Result<()>,
    {
        let Some(name) = header
            .reference_sequences()
            .get_index(reference_sequence_id)
            .map(|(name, _)| name)
        else {
            return Ok(());
        };

        let name = std::str::from_utf8(name).context("reference sequence name is not UTF-8")?;
        let start = Position::try_from(start + 1)?;

        let region = match end {
            Some(end) if end < start.get() => return Ok(()),
            Some(end) => Region::new(name, start..=Position::try_from(end)?),
            None if start == Position::MIN => Region::new(name, ..),
            None => Region::new(name, start..),
        };

        self.run_queries(index, header, &[Query::Region(region)], f)
    }

    fn run_queries<F>(
        &mut self,
        index: &AlignmentIndex,
        header: &sam::Header,
        queries: &[Query<'_>],
        mut f: F,
    ) -> anyhow::Result<()>
    where
        F: FnMut(Rec<'_>) -> anyhow::Result<()>,
    {
        // In a region set, a record is visited only while querying the first merged interval it
        // overlaps, so each record is visited once and in file order.
        let mut visit = |record: Rec<'_>, query: &Query| -> anyhow::Result<()> {
            if let Query::Interval {
                set,
                reference_sequence_id,
                index,
                ..
            } = query
            {
                match record_interval(record.as_dyn(), header)? {
                    Some((id, start, end))
                        if id == *reference_sequence_id
                            && set.first_overlap(id, start, end) == Some(*index) => {}
                    _ => return Ok(()),
                }
            }

            f(record)
        };

        match (self, index) {
            (Self::Bam(reader), AlignmentIndex::Bam(index)) => {
                for query in queries {
                    match query {
                        Query::Unmapped => {
                            for result in reader.query_unmapped(index)? {
                                let record = result?;
                                visit(Rec::Bam(&record), query)?;
                            }
                        }
                        Query::Region(region) | Query::Interval { region, .. } => {
                            let mut query_reader = reader.query(header, index, region)?;
                            let mut record = bam::Record::default();

                            while query_reader.read_record(&mut record)? != 0 {
                                visit(Rec::Bam(&record), query)?;
                            }
                        }
                    }
                }
            }
            (Self::Cram(reader, _, _), AlignmentIndex::Cram(index)) => {
                for query in queries {
                    match query {
                        Query::Unmapped => {
                            for result in reader.query_unmapped(header, index)? {
                                let record = result?;
                                visit(Rec::Buf(&record), query)?;
                            }
                        }
                        Query::Region(region) | Query::Interval { region, .. } => {
                            let mut query_reader = reader.query(header, index, region)?;
                            let mut record = RecordBuf::default();

                            while query_reader.read_record_buf(&mut record)? != 0 {
                                visit(Rec::Buf(&record), query)?;
                            }
                        }
                    }
                }
            }
            (Self::Sam(_), _) => {
                bail!("region queries on SAM input are not supported; use BAM or CRAM")
            }
            _ => bail!("the index does not match the input format"),
        }

        Ok(())
    }
}

/// Encodes a CRAM record as a BAM record into `dst`, reusing its buffer.
pub fn encode_cram_record(
    header: &sam::Header,
    record: &cram::Record<'_>,
    dst: &mut bam::Record,
) -> io::Result<()> {
    let mut buf = Vec::from(mem::take(dst));
    record.encode_bam_record(header, &mut buf)?;
    *dst = bam::Record::try_from(buf)?;
    Ok(())
}

/// A BAM (BAI or CSI) or CRAM (CRAI) index.
pub enum AlignmentIndex {
    Bam(bam::Index),
    Cram(cram::crai::Index),
}

impl AlignmentIndex {
    /// Reads the index named by `explicit`, or the one next to `path`, for an input of `format`.
    pub fn read(format: Format, path: &Path, explicit: Option<&Path>) -> anyhow::Result<Self> {
        match format {
            Format::Bam => read_bam_index(path, explicit).map(Self::Bam),
            Format::Cram => match explicit {
                Some(index) => cram::crai::fs::read(index),
                None => cram::fs::read_associated_index(path),
            }
            .map(Self::Cram)
            .with_context(|| format!("reading index for {}", path.display())),
            Format::Sam => bail!("region queries on SAM input are not supported; use BAM or CRAM"),
        }
    }
}

/// What to query through an index.
pub enum Queries<'a> {
    /// samtools regions, queried in the order given. Without `-M`, a record overlapping several
    /// regions is visited once per region.
    Regions(&'a [String]),
    /// Merged intervals (`-M`, `-L`), visiting each overlapping record once, in file order.
    Set(&'a RegionSet),
}

enum Query<'a> {
    Unmapped,
    Region(Region),
    Interval {
        region: Region,
        set: &'a RegionSet,
        reference_sequence_id: usize,
        index: usize,
    },
}

impl<'a> Queries<'a> {
    fn to_list(&self, header: &sam::Header) -> anyhow::Result<Vec<Query<'a>>> {
        match *self {
            Self::Regions(raw_regions) => raw_regions
                .iter()
                .map(|raw| {
                    if raw == UNMAPPED_REGION {
                        Ok(Query::Unmapped)
                    } else {
                        parse_region(raw).map(Query::Region)
                    }
                })
                .collect(),
            Self::Set(set) => {
                let mut queries = Vec::new();

                for (id, name) in header.reference_sequences().keys().enumerate() {
                    for (i, &(start, end)) in set.intervals(id).iter().enumerate() {
                        queries.push(Query::Interval {
                            region: region_string(name, start, end)?,
                            set,
                            reference_sequence_id: id,
                            index: i,
                        });
                    }
                }

                Ok(queries)
            }
        }
    }
}

/// Reads an explicit BAI or CSI index (by extension), or the one associated with `path`.
fn read_bam_index(path: &Path, index: Option<&Path>) -> anyhow::Result<bam::Index> {
    let result = match index {
        Some(index) if index.extension().is_some_and(|ext| ext == "csi") => {
            noodles_csi::fs::read(index).map(bam::Index::Csi)
        }
        Some(index) => bam::bai::fs::read(index).map(bam::Index::Bai),
        None => bam::fs::read_associated_index(path),
    };

    result.with_context(|| format!("reading index for {}", path.display()))
}

fn parse_region(s: &str) -> anyhow::Result<Region> {
    s.parse().with_context(|| format!("invalid region: {s}"))
}

/// The local FASTA file a CRAM header's `@SQ` lines name in their `UR` fields (the first that
/// exists, with an index), as HTSlib uses them when no reference is given.
fn header_reference(header: &sam::Header) -> Option<PathBuf> {
    use sam::header::record::value::map::reference_sequence::tag::URI;

    header
        .reference_sequences()
        .values()
        .filter_map(|reference_sequence| reference_sequence.other_fields().get(&URI))
        .map(|uri| {
            let uri = uri.to_string();
            PathBuf::from(uri.strip_prefix("file://").unwrap_or(&uri))
        })
        .find(|path| {
            let mut index = path.clone().into_os_string();
            index.push(".fai");
            path.is_file() && Path::new(&index).is_file()
        })
}

/// Builds a reference sequence repository from an indexed FASTA, or an empty one.
pub fn reference_repository(src: Option<&Path>) -> anyhow::Result<fasta::Repository> {
    let Some(src) = src else {
        return Ok(fasta::Repository::default());
    };

    crate::reference::ensure_index(src)?;

    let reader = fasta::io::indexed_reader::Builder::default()
        .build_from_path(src)
        .with_context(|| format!("opening reference {}", src.display()))?;

    Ok(fasta::Repository::new(
        fasta::repository::adapters::IndexedReader::new(reader),
    ))
}
