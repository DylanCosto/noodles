//! `samtools cat` for CRAM (`cram_cat`): containers are copied as they are, except those that
//! must be filtered to a region or need new read group numbers, which are decoded and encoded
//! again.

use std::{
    fs::File,
    io::{self, Read, Seek, SeekFrom, Write},
    path::{Path, PathBuf},
    sync::Arc,
};

use anyhow::{Context, bail};
use bstr::BString;
use noodles_bam as bam;
use noodles_cram as cram;
use noodles_pileup::md::{self, FillMdOptions};
use noodles_sam::{self as sam, alignment::io::Write as _};

use crate::{
    alignment::{AlignmentReader, ReadOptions, encode_cram_record},
    bam_raw,
    cram_file::{ContainerHeader, FileDefinition, MULTI_REFERENCE, encode_header_container},
    hts_region,
    regions::record_interval,
};

/// The CRAM 3 EOF container.
const EOF_3: [u8; 38] = [
    0x0f, 0x00, 0x00, 0x00, 0xff, 0xff, 0xff, 0xff, 0x0f, 0xe0, 0x45, 0x4f, 0x46, 0x00, 0x00, 0x00,
    0x00, 0x01, 0x00, 0x05, 0xbd, 0xd9, 0x4f, 0x00, 0x01, 0x00, 0x06, 0x06, 0x01, 0x00, 0x01, 0x00,
    0x01, 0x00, 0xee, 0x63, 0x01, 0x4b,
];

/// The CRAM 2.1 EOF container.
const EOF_2: [u8; 30] = [
    0x0b, 0x00, 0x00, 0x00, 0xff, 0xff, 0xff, 0xff, 0xff, 0xe0, 0x45, 0x4f, 0x46, 0x00, 0x00, 0x00,
    0x00, 0x01, 0x00, 0x00, 0x01, 0x00, 0x06, 0x06, 0x01, 0x00, 0x01, 0x00, 0x01, 0x00,
];

/// HTSlib's `HTS_POS_MAX`-like "after every reference" ID for unplaced records (`*`).
const UNPLACED: i64 = i64::MAX;

/// How containers that partly overlap a region are handled.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub enum Filtering {
    /// Keep only the records in the region, as `samtools view` would.
    Exact,
    /// Copy the containers whole (`-f`).
    Fast,
}

pub struct Options<'a> {
    pub region: Option<&'a str>,
    pub part: Option<&'a str>,
    pub filtering: Filtering,
    pub reference: Option<&'a Path>,
}

/// An input file positioned after its header container.
struct CramInput {
    path: PathBuf,
    file: File,
    definition: FileDefinition,
    header: sam::Header,
}

impl CramInput {
    fn open(path: &Path) -> anyhow::Result<Self> {
        let mut file =
            File::open(path).with_context(|| format!("fail to open file '{}'", path.display()))?;

        let header = cram::io::Reader::new(&mut file)
            .read_header()
            .with_context(|| format!("fail to read the header of file '{}'", path.display()))?;

        file.seek(SeekFrom::Start(0))?;
        let definition = FileDefinition::read(&mut file)?;
        let container = ContainerHeader::read(&mut file, &definition)?
            .context("missing CRAM header container")?;
        file.seek(SeekFrom::Current(to_i64(container.data_len)))?;

        Ok(Self {
            path: path.to_path_buf(),
            file,
            definition,
            header,
        })
    }
}

/// Concatenates CRAM files (`cram_cat`). `add_pg` adds the `@PG` line to the output header.
pub fn cat<W, F>(
    inputs: &[PathBuf],
    header: Option<sam::Header>,
    options: &Options<'_>,
    dst: &mut W,
    add_pg: F,
) -> anyhow::Result<()>
where
    W: Write,
    F: FnOnce(&mut sam::Header) -> anyhow::Result<()>,
{
    let mut inputs = inputs
        .iter()
        .map(|path| CramInput::open(path))
        .collect::<anyhow::Result<Vec<_>>>()?;

    let definition = inputs[0].definition;

    if inputs
        .iter()
        .any(|input| input.definition.0[4..6] != definition.0[4..6])
    {
        bail!("input files have differing version numbers");
    }

    let mut merged = header;

    for input in &inputs {
        merge_read_groups(&mut merged, &input.header)?;
    }

    let mut header = merged.context("no input headers")?;

    // The read group number each input's single read group gets, if it changes.
    let new_read_groups = inputs
        .iter()
        .map(|input| new_read_group(&header, &input.header))
        .collect::<anyhow::Result<Vec<_>>>()?;

    add_pg(&mut header)?;

    dst.write_all(&definition.0)?;
    dst.write_all(&encode_header_container(&header, &definition)?)?;

    for (input, &new_read_group) in inputs.iter_mut().zip(&new_read_groups) {
        let done = cat_input(input, &header, new_read_group, options, dst)?;

        if done == Progress::Stop {
            break;
        }
    }

    if definition.major_version() >= 3 {
        dst.write_all(&EOF_3)?;
    } else {
        dst.write_all(&EOF_2)?;
    }

    Ok(())
}

/// Adds the read groups of `header` that `merged` lacks; for CRAM, inputs whose read groups
/// are as many as the merged ones must list them in the same order.
fn merge_read_groups(merged: &mut Option<sam::Header>, header: &sam::Header) -> anyhow::Result<()> {
    let Some(merged) = merged else {
        *merged = Some(header.clone());
        return Ok(());
    };

    for (id, read_group) in header.read_groups() {
        if !merged.read_groups().contains_key(id) {
            merged
                .read_groups_mut()
                .insert(id.clone(), read_group.clone());
        }
    }

    let count = header.read_groups().len();

    if count > 1
        && merged.read_groups().len() == count
        && !merged.read_groups().keys().eq(header.read_groups().keys())
    {
        bail!("Same size @RG lists but differing order / contents");
    }

    Ok(())
}

/// The index of an input's only read group in the merged header, when it is not the first.
fn new_read_group(merged: &sam::Header, header: &sam::Header) -> anyhow::Result<Option<usize>> {
    if header.read_groups().len() != 1 {
        return Ok(None);
    }

    let Some((id, _)) = header.read_groups().first() else {
        return Ok(None);
    };

    let index = merged
        .read_groups()
        .get_index_of(id)
        .with_context(|| format!("fail to find @RG line '{id}' in the new header"))?;

    Ok((index != 0).then_some(index))
}

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
enum Progress {
    Continue,
    /// `-p` selected nothing, which ends the whole concatenation in samtools.
    Stop,
}

/// What `-r` selects in one input.
enum Selection {
    /// Containers by number (`#:start-end`).
    Containers,
    /// A region: a reference sequence ID (or [`UNPLACED`]) and a 0-based, half-open interval.
    Region { id: i64, start: i64, end: i64 },
}

fn cat_input<W: Write>(
    input: &mut CramInput,
    header: &sam::Header,
    new_read_group: Option<usize>,
    options: &Options<'_>,
    dst: &mut W,
) -> anyhow::Result<Progress> {
    let Some(Scope { bounds, selection }) = scope(input, options)? else {
        return Ok(Progress::Stop);
    };

    if let Some((start, _)) = bounds {
        input.file.seek(SeekFrom::Start(start))?;
    }

    let mut reencoder: Option<Reencoder> = None;
    let mut last_reference_sequence_id: i64 = -99;
    let mut data = Vec::new();

    loop {
        let offset = input.file.stream_position()?;

        let Some(container) = ContainerHeader::read(&mut input.file, &input.definition)? else {
            break;
        };

        if container.is_eof() {
            input
                .file
                .seek(SeekFrom::Current(to_i64(container.data_len)))?;
            continue;
        }

        if let Some((_, end)) = bounds
            && offset > end
        {
            break;
        }

        if new_read_group.is_some() {
            if options.region.is_some() {
                bail!("Cannot specify a region while transcoding RG lines");
            }

            input
                .file
                .seek(SeekFrom::Current(to_i64(container.data_len)))?;
            let reencoder = Reencoder::get(&mut reencoder, input, options)?;
            reencoder.reencode(offset, &input.header, header, None, dst)?;
            continue;
        }

        let filter = container_filter(
            &container,
            selection.as_ref(),
            options.filtering,
            last_reference_sequence_id,
        );

        match (filter, selection.as_ref()) {
            (ContainerFilter::End, _) => break,
            (
                ContainerFilter::Skip | ContainerFilter::Decode,
                Some(&Selection::Region { id, start, end }),
            ) => {
                input
                    .file
                    .seek(SeekFrom::Current(to_i64(container.data_len)))?;

                // `cram_filter_container`: containers of one other reference, or ending before
                // the region, are skipped without decoding.
                // HTSlib's range for `*` is reference -1 from position 0.
                let (range_id, range_start) = if id == UNPLACED {
                    (-1, 0)
                } else {
                    (id, start + 1)
                };
                let container_id = i64::from(container.reference_sequence_id);

                if container.reference_sequence_id != MULTI_REFERENCE
                    && (range_id != container_id
                        || range_start > i64::from(container.start) + i64::from(container.span) - 1)
                {
                    last_reference_sequence_id = container_id;
                    continue;
                }

                let region = RecordRange {
                    id: range_id,
                    start: range_start,
                    end,
                };

                let reencoder = Reencoder::get(&mut reencoder, input, options)?;
                last_reference_sequence_id =
                    reencoder.reencode(offset, &input.header, header, Some(&region), dst)?;
            }
            _ => {
                data.resize(container.data_len, 0);
                input.file.read_exact(&mut data)?;
                dst.write_all(&container.raw)?;
                dst.write_all(&data)?;
            }
        }
    }

    Ok(Progress::Continue)
}

/// The part of an input `-r` and `-p` select: the byte range of its containers and what `-r`
/// selects in them.
struct Scope {
    bounds: Option<(u64, u64)>,
    selection: Option<Selection>,
}

/// Finds an input's [`Scope`], or `None` when `-p` selects nothing.
fn scope(input: &CramInput, options: &Options<'_>) -> anyhow::Result<Option<Scope>> {
    let index = if options.region.is_some() || options.part.is_some() {
        Some(ContainerIndex::read(&input.path)?)
    } else {
        None
    };

    let mut bounds: Option<(u64, u64)> = None;
    let mut selection = None;

    if let (Some(region), Some(index)) = (options.region, &index) {
        if let Some(range) = region.strip_prefix("#:") {
            bounds = Some(index.container_number_range(range)?);
            selection = Some(Selection::Containers);
        } else {
            let (id, start, end) = parse_region(region, &input.header)?;
            bounds = Some(index.extents(id, start, end)?);
            selection = Some(Selection::Region { id, start, end });
        }
    }

    if let (Some(part), Some(index)) = (options.part, &index) {
        match index.part(part, bounds)? {
            Some(range) => bounds = Some(range),
            None => return Ok(None),
        }
    }

    Ok(Some(Scope { bounds, selection }))
}

/// How a container is handled for a `-r` region: copied as is, skipped, or decoded and filtered
/// record by record.
fn container_filter(
    container: &ContainerHeader,
    selection: Option<&Selection>,
    filtering: Filtering,
    last_reference_sequence_id: i64,
) -> ContainerFilter {
    let Some(&Selection::Region { id, start, end }) = selection else {
        return ContainerFilter::Copy;
    };

    let container_id = match i64::from(container.reference_sequence_id) {
        -1 => UNPLACED,
        n => n,
    };

    let container_start = i64::from(container.start);
    let container_end = container_start + i64::from(container.span);

    let filter = if container_id > id || container_start > end {
        return ContainerFilter::End;
    } else if container.reference_sequence_id == MULTI_REFERENCE {
        ContainerFilter::Decode
    } else if container_id < id || container_end < start {
        ContainerFilter::Skip
    } else if container_start < start || container_end > end {
        match filtering {
            Filtering::Exact => ContainerFilter::Decode,
            Filtering::Fast => ContainerFilter::Copy,
        }
    } else {
        ContainerFilter::Copy
    };

    if last_reference_sequence_id == -1 && id == UNPLACED {
        ContainerFilter::Copy
    } else {
        filter
    }
}

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
enum ContainerFilter {
    Copy,
    Skip,
    Decode,
    /// The region ends before this container, so the input is done.
    End,
}

/// Parses a region against an input's header (`sam_itr_querys`): `*` is the unplaced records.
fn parse_region(region: &str, header: &sam::Header) -> anyhow::Result<(i64, i64, i64)> {
    if region == "*" {
        return Ok((UNPLACED, 0, 0));
    }

    let reference_sequences = header.reference_sequences();

    let (id, start, end) =
        hts_region::parse_region(region, |name| reference_sequences.get_index_of(name))
            .with_context(|| format!("Unable to parse region {region}"))?;

    Ok((to_i64(id), start, end))
}

/// The records `cram_filter_container` keeps: on reference `id` (-1 for unplaced), overlapping
/// the 1-based closed interval `[start, end]`.
struct RecordRange {
    id: i64,
    start: i64,
    end: i64,
}

/// Decodes containers and encodes the records again, with the output header.
struct Reencoder {
    reader: cram::io::Reader<std::io::BufReader<crate::source::Stream>>,
    repository: noodles_fasta::Repository,
    /// The last reference sequence used for `MD`/`NM`, by name.
    reference: Option<(BString, Arc<noodles_fasta::record::Sequence>)>,
}

impl Reencoder {
    fn get<'r>(
        slot: &'r mut Option<Self>,
        input: &CramInput,
        options: &Options<'_>,
    ) -> anyhow::Result<&'r mut Self> {
        if input.definition.major_version() < 3 {
            bail!(
                "ntools cannot re-encode CRAM {}.x containers",
                input.definition.major_version()
            );
        }

        if slot.is_none() {
            let src = input
                .path
                .to_str()
                .context("input path is not valid UTF-8")?;
            let read_options = ReadOptions {
                threads: 0,
                reference: options.reference.map(Path::to_path_buf),
            };

            let (mut reader, _) = AlignmentReader::open(src, &read_options)?;
            reader.read_header()?;

            let AlignmentReader::Cram(reader, repository, _) = reader else {
                bail!("{src}: not a CRAM file");
            };

            *slot = Some(Self {
                reader,
                repository,
                reference: None,
            });
        }

        Ok(slot.as_mut().expect("reencoder"))
    }

    /// Writes the container at `offset` again, keeping the records in `range` (all records
    /// without one). Returns the reference sequence ID all examined records share, -2 if they
    /// differ, or -3 if there were none (`cram_filter_container`'s `ref_id`).
    fn reencode<W: Write>(
        &mut self,
        offset: u64,
        input_header: &sam::Header,
        output_header: &sam::Header,
        range: Option<&RecordRange>,
        dst: &mut W,
    ) -> anyhow::Result<i64> {
        self.reader.get_mut().seek(SeekFrom::Start(offset))?;

        let mut container = cram::io::reader::Container::default();
        self.reader.read_container(&mut container)?;

        let mut writer = cram::io::writer::Builder::default()
            .set_reference_sequence_repository(self.repository.clone())
            .build_from_writer(Vec::new());

        let compression_header = container.compression_header()?;
        let mut bam_record = bam::Record::default();
        let mut fixed_reference_sequence_id: i64 = -3;

        'slices: for result in container.slices() {
            let slice = result?;
            let (core_data_src, external_data_srcs) = slice.decode_blocks()?;

            let records = slice.records(
                self.repository.clone(),
                input_header,
                &compression_header,
                &core_data_src,
                &external_data_srcs,
            )?;

            for record in &records {
                encode_cram_record(input_header, record, &mut bam_record)?;

                if let Some(range) = range {
                    let record_id = bam_record
                        .reference_sequence_id()
                        .transpose()?
                        .map_or(-1, to_i64);

                    fixed_reference_sequence_id = match fixed_reference_sequence_id {
                        -3 => record_id,
                        id if id == record_id => id,
                        _ => -2,
                    };

                    match keep(range, record_id, &bam_record, input_header)? {
                        Keep::Yes => {}
                        Keep::Skip => continue,
                        Keep::Stop => break 'slices,
                    }
                }

                self.drop_regenerated_md_nm(input_header, &mut bam_record)?;
                writer.write_alignment_record(output_header, &bam_record)?;
            }
        }

        writer.try_finish(output_header)?;
        copy_containers(&writer.into_inner(), dst)?;

        Ok(fixed_reference_sequence_id)
    }

    /// Removes `MD` and `NM` tags that equal the values computed from the reference, as
    /// HTSlib's CRAM writer does, so readers regenerate them where HTSlib's would.
    fn drop_regenerated_md_nm(
        &mut self,
        header: &sam::Header,
        record: &mut bam::Record,
    ) -> anyhow::Result<()> {
        const UNMAPPED: u16 = 0x04;

        let raw: &[u8] = record.as_ref();
        let [md, nm] = bam_raw::find_fields(raw, [*b"MD", *b"NM"])?;

        if (md.is_none() && nm.is_none()) || bam_raw::Record::new(raw)?.flags() & UNMAPPED != 0 {
            return Ok(());
        }

        let Some(id) = record.reference_sequence_id().transpose()? else {
            return Ok(());
        };

        let Some((name, _)) = header.reference_sequences().get_index(id) else {
            return Ok(());
        };

        if self.reference.as_ref().is_none_or(|(n, _)| n != name) {
            let Some(sequence) = self.repository.get(name).transpose()? else {
                return Ok(());
            };

            self.reference = Some((name.clone(), sequence));
        }

        let Some((_, sequence)) = &self.reference else {
            return Ok(());
        };

        let mut read = noodles_pileup::Read::new(raw.to_vec())?;

        let Some(report) = md::fill_md(
            &mut read,
            sequence.as_ref().as_ref(),
            &FillMdOptions::default(),
        ) else {
            return Ok(());
        };

        let drop_md = md.is_some() && report.md.is_none();
        let drop_nm = nm.is_some() && report.nm.is_none();

        if drop_md || drop_nm {
            let mut stripped = Vec::with_capacity(raw.len());
            bam_raw::strip_fields(raw, &mut stripped, |tag| {
                (drop_md && tag == *b"MD") || (drop_nm && tag == *b"NM")
            })?;
            *record = bam::Record::try_from(stripped)?;
        }

        Ok(())
    }
}

enum Keep {
    Yes,
    Skip,
    /// Past the range: no later record in the container is kept.
    Stop,
}

fn keep(
    range: &RecordRange,
    record_id: i64,
    record: &bam::Record,
    header: &sam::Header,
) -> io::Result<Keep> {
    if range.id != record_id {
        return Ok(if range.id > record_id || range.id == -1 {
            Keep::Skip
        } else {
            Keep::Stop
        });
    }

    let (start, end) = match record_interval(record, header)? {
        Some((_, start, end)) => (to_i64(start), to_i64(end)),
        None => (0, 0),
    };

    if end < range.start {
        Ok(Keep::Skip)
    } else if start > range.end {
        Ok(Keep::Stop)
    } else {
        Ok(Keep::Yes)
    }
}

/// Copies the containers of a CRAM stream without its EOF containers. The stream has no file
/// definition or header container.
fn copy_containers<W: Write>(mut src: &[u8], dst: &mut W) -> io::Result<()> {
    let definition = FileDefinition(*b"CRAM\x03\x00\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0");

    while let Some(container) = ContainerHeader::read(&mut src, &definition)? {
        let data = src
            .split_off(..container.data_len)
            .ok_or_else(|| io::Error::from(io::ErrorKind::UnexpectedEof))?;

        if !container.is_eof() {
            dst.write_all(&container.raw)?;
            dst.write_all(data)?;
        }
    }

    Ok(())
}

/// Prints the number of indexed containers of each input, or of those covering a region (`-q`).
pub fn print_container_counts(inputs: &[PathBuf], region: Option<&str>) -> anyhow::Result<()> {
    let mut stdout = io::stdout().lock();

    for path in inputs {
        let index = ContainerIndex::read(path).context("No index found")?;

        let (start, end) = match region {
            Some(region) => {
                let header = cram::io::Reader::new(File::open(path)?)
                    .read_header()
                    .context("Unable to read header")?;
                let (id, start, end) = parse_region(region, &header)?;
                index.extents(id, start, end)?
            }
            None => (0, 0),
        };

        let (first, last) = index.numbers_between(start, end);
        writeln!(
            stdout,
            "{}\t{}\t{first}\t{last}",
            path.display(),
            last - first + 1
        )?;
    }

    Ok(())
}

fn to_i64<N: TryInto<i64>>(n: N) -> i64 {
    n.try_into().unwrap_or(i64::MAX)
}

/// A slice's entry in a CRAM index.
#[derive(Clone, Copy, Debug)]
struct IndexEntry {
    /// 1-based start and end (inclusive).
    start: i64,
    end: i64,
    /// The offset of the slice's container, from the start of the file.
    offset: u64,
}

/// A CRAM index, as HTSlib numbers its containers.
struct ContainerIndex {
    /// Entries by reference sequence ID + 1 (0 for unplaced slices), in index order.
    by_reference: Vec<Vec<IndexEntry>>,
    /// Container offsets in HTSlib's order: by reference sequence, unplaced last.
    offsets: Vec<u64>,
}

impl ContainerIndex {
    fn read(path: &Path) -> anyhow::Result<Self> {
        let index = cram::fs::read_associated_index(path)
            .with_context(|| format!("failed to load index for {}", path.display()))?;

        let mut by_reference: Vec<Vec<IndexEntry>> = Vec::new();

        for record in &index {
            let slot = record.reference_sequence_id().map_or(0, |id| id + 1);

            if by_reference.len() <= slot {
                by_reference.resize_with(slot + 1, Vec::new);
            }

            let start = record
                .alignment_start()
                .map_or(0, |p| to_i64(usize::from(p)));

            by_reference[slot].push(IndexEntry {
                start,
                end: start + to_i64(record.alignment_span()) - 1,
                offset: record.offset(),
            });
        }

        // `cram_container_num2offset`: references in order, then unplaced; repeated offsets
        // (multi-reference containers) count once.
        let mut offsets: Vec<u64> = Vec::new();
        let order = (1..by_reference.len()).chain(by_reference.first().map(|_| 0));

        for slot in order {
            for entry in &by_reference[slot] {
                if offsets.last() != Some(&entry.offset) {
                    offsets.push(entry.offset);
                }
            }
        }

        Ok(Self {
            by_reference,
            offsets,
        })
    }

    fn entries(&self, id: i64) -> &[IndexEntry] {
        let slot = if id == UNPLACED { 0 } else { id + 1 };

        usize::try_from(slot)
            .ok()
            .and_then(|slot| self.by_reference.get(slot))
            .map_or(&[], Vec::as_slice)
    }

    /// `cram_index_query`: the first entry overlapping `pos`, or the one before it.
    fn query(&self, id: i64, pos: i64) -> Option<usize> {
        let entries = self.entries(id);
        let pos = if id == UNPLACED { 0 } else { pos };

        if entries.is_empty() {
            return None;
        }

        let (mut i, mut j) = (0, entries.len() - 1);
        let mut k = j / 2;

        while k != i {
            if entries[k].start >= pos {
                j = k;
            } else {
                i = k;
            }

            k = (j - i) / 2 + i;
        }

        if entries[j].start < pos {
            i = j;
        }

        while i > 0 && entries[i - 1].end >= pos {
            i -= 1;
        }

        while i + 1 < entries.len() && entries[i].end < pos {
            i += 1;
        }

        Some(i)
    }

    /// `cram_index_extents`: the offsets of the first and last containers covering a region.
    fn extents(&self, id: i64, start: i64, end: i64) -> anyhow::Result<(u64, u64)> {
        let entries = self.entries(id);

        let first = self.query(id, start).context("Failed to query index")?;
        let mut last = self.query(id, end).context("Failed to query index")?;

        // `cram_index_query_last`: later entries starting by the end, then any with the same
        // offset.
        while last + 1 < entries.len()
            && entries[last + 1].start <= if id == UNPLACED { 0 } else { end }
        {
            last += 1;
        }

        while last + 1 < entries.len() && entries[last + 1].offset == entries[last].offset {
            last += 1;
        }

        Ok((entries[first].offset, entries[last].offset))
    }

    /// `cram_num_containers_between`: the first and last numbers of containers starting in
    /// `[start, end]` (no end if 0), or -1 for none.
    fn numbers_between(&self, start: u64, end: u64) -> (i64, i64) {
        let mut numbers = self
            .offsets
            .iter()
            .enumerate()
            .filter(|&(_, &offset)| offset >= start && (end == 0 || offset <= end))
            .map(|(i, _)| to_i64(i));

        let first = numbers.next().unwrap_or(-1);
        let last = numbers.next_back().unwrap_or(first);

        (first, last)
    }

    /// `cram_container_offset2num`: the number of the first container at or after `offset`.
    fn number_at(&self, offset: u64) -> Option<usize> {
        self.offsets.iter().position(|&o| o >= offset)
    }

    /// Parses `start-end` or `start` container numbers (`#:` regions) into offsets.
    fn container_number_range(&self, range: &str) -> anyhow::Result<(u64, u64)> {
        let (start, end) = match range.split_once('-') {
            Some((start, end)) => (start, Some(end)),
            None => (range, None),
        };

        let start: usize =
            leading_int(start).with_context(|| format!("ERROR: Malformed region: #:{range}"))?;
        let end = match end.and_then(leading_int) {
            Some(end) => end,
            None => start,
        };

        if end >= self.offsets.len() {
            bail!(
                "Too many containers.  The end range should be < {}",
                self.offsets.len()
            );
        }

        let offset = |n: usize| {
            self.offsets
                .get(n)
                .copied()
                .context("container number out of range")
        };

        Ok((offset(start)?, offset(end)?))
    }

    /// `cram_subdivide_part`: the containers of part N of M (`N/M`) of `bounds` (or the whole
    /// file), or `None` if the part is empty.
    fn part(&self, part: &str, bounds: Option<(u64, u64)>) -> anyhow::Result<Option<(u64, u64)>> {
        let malformed = || format!("malformed region {part}. Should be e.g. '1/10'");
        let (a, b) = part.split_once('/').with_context(malformed)?;
        let a: i64 = leading_int(a).with_context(malformed)?;
        let mut b: i64 = leading_int(b).with_context(malformed)?;

        let (first, last) = match bounds {
            Some((start, end)) => (
                self.number_at(start).map_or(-1, to_i64),
                self.number_at(end).map_or(-1, to_i64),
            ),
            None => (0, to_i64(self.offsets.len()) - 1),
        };

        let n = last - first + 1;
        b = b.min(n);

        #[expect(
            clippy::cast_possible_truncation,
            clippy::cast_precision_loss,
            reason = "samtools' arithmetic"
        )]
        let (start, end) = (
            ((a - 1) as f64 * n as f64 / b as f64) as i64,
            (a as f64 * n as f64 / b as f64) as i64 - 1,
        );

        if start < 0 || end >= n {
            return Ok(None);
        }

        let offset = |k: i64| {
            usize::try_from(k + first)
                .ok()
                .and_then(|k| self.offsets.get(k).copied())
                .context("container number out of range")
        };

        Ok(Some((offset(start)?, offset(end)?)))
    }
}

/// Parses a leading decimal integer, as `sscanf("%d")` does.
fn leading_int<N: std::str::FromStr>(s: &str) -> Option<N> {
    let s = s.trim_start();
    let end = s
        .char_indices()
        .find(|&(i, c)| !(c.is_ascii_digit() || (i == 0 && (c == '-' || c == '+'))))
        .map_or(s.len(), |(i, _)| i);

    s[..end].parse().ok()
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::cram_file::FILE_DEFINITION_LEN;

    #[test]
    fn test_leading_int() {
        assert_eq!(leading_int::<i64>("12"), Some(12));
        assert_eq!(leading_int::<i64>("3x"), Some(3));
        assert_eq!(leading_int::<i64>("-2"), Some(-2));
        assert_eq!(leading_int::<i64>("x"), None);
    }

    #[test]
    fn test_eof_containers() -> io::Result<()> {
        for (eof, major) in [(&EOF_3[..], 3), (&EOF_2[..], 2)] {
            let mut definition = [0; FILE_DEFINITION_LEN];
            definition[..5].copy_from_slice(b"CRAM\0");
            definition[4] = major;

            let mut src = eof;
            let container =
                ContainerHeader::read(&mut src, &FileDefinition(definition))?.expect("container");

            assert!(container.is_eof());
            assert_eq!(src.len(), container.data_len);
        }

        Ok(())
    }
}
