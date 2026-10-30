//! `ntools merge`: merge sorted alignment files (samtools merge).

mod header;

use std::{
    collections::HashMap,
    fs::File,
    io::{self, BufRead, BufReader, Read},
    path::{Path, PathBuf},
    thread,
    time::{SystemTime, UNIX_EPOCH},
};

use anyhow::{Context, bail};
use bstr::{BStr, BString, ByteSlice};
use crossbeam_channel::{Receiver, Sender};
use flate2::read::MultiGzDecoder;
use noodles_sam as sam;

use crate::{
    alignment::{AlignmentIndex, AlignmentReader, Queries, RawEncoder, ReadOptions},
    bam_raw,
    bgzf_stream::BgzfReader,
    commands::sort::parse_tag,
    hts_region::{self, POS_MAX},
    output::{self, AlignmentWriter, OutputFormat, WriteOptions},
    program,
    regions::RegionSet,
    runs::{self, RecordSource},
    sort_order::{Order, TagTieBreak, remap_reference_sequence_ids},
    source::{Format, Stream},
};

use self::header::{Collisions, HeaderMerger, Translation};

/// Raw records are sent from reading threads in batches of about this many bytes.
const BATCH_SIZE: usize = 1 << 18;

#[derive(Debug, clap::Args)]
#[command(disable_help_flag = true)]
pub struct Args {
    /// Output file, then sorted input files (or give the output with -o). With -X, the inputs
    /// are followed by their index files.
    files: Vec<PathBuf>,

    /// Write output to this file; all positional arguments are then inputs.
    #[arg(short = 'o', long = "output")]
    output: Option<PathBuf>,

    /// Overwrite the output file if it exists.
    #[arg(short = 'f', long = "force")]
    force: bool,

    /// Output format [default: from the output extension, else BAM].
    #[arg(short = 'O', long = "output-fmt", ignore_case = true)]
    output_format: Option<OutputFormat>,

    /// BGZF compression level for BAM output.
    #[arg(short = 'l', default_value_t = 6)]
    level: u8,

    /// Compress at level 1.
    #[arg(short = '1')]
    fast: bool,

    /// Output uncompressed BAM.
    #[arg(short = 'u')]
    uncompressed: bool,

    /// Inputs are sorted by read name, compared naturally.
    #[arg(short = 'n', conflicts_with = "name_ascii")]
    name_natural: bool,

    /// Inputs are sorted by read name, bytewise.
    #[arg(short = 'N')]
    name_ascii: bool,

    /// Inputs are sorted by the value of TAG (then by position, or by name with -n/-N).
    #[arg(short = 't', value_name = "TAG", value_parser = parse_tag)]
    tag: Option<[u8; 2]>,

    /// Attach an RG tag to each record, named after its input file.
    #[arg(short = 'r')]
    attach_read_groups: bool,

    /// Merge only the records in this region (inputs must be indexed).
    #[arg(short = 'R', value_name = "STR", conflicts_with = "bed")]
    region: Option<String>,

    /// Merge only the records in the regions of this BED file (inputs must be indexed).
    #[arg(short = 'L', value_name = "FILE")]
    bed: Option<PathBuf>,

    /// Start the merged header with the header of this file.
    #[arg(short = 'h', value_name = "FILE")]
    header: Option<PathBuf>,

    /// Combine @RG lines with colliding IDs instead of renaming them.
    #[arg(short = 'c')]
    combine_read_groups: bool,

    /// Combine @PG lines with colliding IDs instead of renaming them.
    #[arg(short = 'p')]
    combine_programs: bool,

    /// Seed for the random suffixes of renamed IDs [default: the current time].
    #[arg(short = 's', value_name = "VALUE")]
    seed: Option<i64>,

    /// A file listing input files, one per line (repeatable); these come before the other
    /// inputs.
    #[arg(short = 'b', value_name = "FILE", conflicts_with = "index_files")]
    lists: Vec<PathBuf>,

    /// Each input file is followed by its index file, after all the inputs.
    #[arg(short = 'X')]
    index_files: bool,

    /// Reference FASTA (indexed) for CRAM input and output.
    #[arg(long = "reference")]
    reference: Option<PathBuf>,

    /// Additional threads for compression and decompression.
    #[arg(short = '@', long = "threads", default_value_t = 0)]
    threads: usize,

    /// Do not add an @PG line to the header.
    #[arg(long = "no-PG")]
    no_pg: bool,

    /// Print help.
    #[arg(long, action = clap::ArgAction::Help)]
    help: Option<bool>,
}

/// One input file and how to read it.
struct InputFile {
    path: PathBuf,
    index: Option<PathBuf>,
}

pub fn run(args: Args) -> anyhow::Result<()> {
    let order = order(&args);
    let (output, inputs) = files(&args)?;

    if !args.force && output != Path::new("-") && output.is_file() {
        bail!(
            "File '{}' exists. Please apply '-f' to overwrite. Abort.",
            output.display()
        );
    }

    let seed = args.seed.unwrap_or_else(|| {
        SystemTime::now()
            .duration_since(UNIX_EPOCH)
            .map_or(0, |d| i64::try_from(d.as_secs()).unwrap_or(0))
    });

    let collisions = |combine| {
        if combine {
            Collisions::Combine
        } else {
            Collisions::Rename
        }
    };

    let mut merger = HeaderMerger::new(
        collisions(args.combine_read_groups),
        collisions(args.combine_programs),
        seed,
    );

    let read_options = ReadOptions {
        threads: args.threads / inputs.len(),
        reference: args.reference.clone(),
    };

    if let Some(path) = &args.header {
        let src = path.to_str().context("header path is not valid UTF-8")?;
        let (mut reader, _) = AlignmentReader::open(src, &read_options)?;
        let header = reader
            .read_header()
            .with_context(|| format!("couldn't read headers from \"{src}\""))?;
        merger.add(&header, None);
    }

    let mut opened = Vec::with_capacity(inputs.len());
    let mut reset_order = false;

    for input in &inputs {
        let src = input
            .path
            .to_str()
            .context("input path is not valid UTF-8")?;
        let (mut reader, _) = AlignmentReader::open(src, &read_options)
            .with_context(|| format!("fail to open \"{src}\""))?;
        let header = reader
            .read_header()
            .with_context(|| format!("failed to read header from \"{src}\""))?;

        let read_group = args
            .attach_read_groups
            .then(|| read_group_name(&input.path));

        let translation = merger.add(&header, read_group.as_ref().map(|s| s.as_bstr()));

        if translation.lost_coordinate_order && order == Order::Coordinate {
            eprintln!(
                "[bam_merge_core] Order of targets in file {src} caused coordinate sort to be lost"
            );
            reset_order = true;
        }

        opened.push((reader, header, translation, read_group));
    }

    let mut header = merger.finish(reset_order);

    if header.header().is_none() {
        eprintln!("[W::bam_merge_core2] No @HD tag found.");
    }

    let selection = Selection::new(&args, &header)?;

    let format = args
        .output_format
        .or_else(|| OutputFormat::from_path(&output))
        .unwrap_or(OutputFormat::Bam);

    let level = if args.uncompressed {
        0
    } else if args.fast {
        1
    } else {
        args.level
    };

    let write_options = WriteOptions {
        format,
        level: output::parse_level(level)?,
        threads: args.threads,
        reference: args.reference.clone(),
    };

    if !args.no_pg {
        program::add(&mut header)?;
    }

    thread::scope(|scope| -> anyhow::Result<()> {
        let mut sources = Vec::with_capacity(opened.len());
        let mut readers = Vec::new();

        for ((reader, input_header, translation, read_group), input) in
            opened.into_iter().zip(&inputs)
        {
            let records = match (&selection, reader) {
                (Selection::All, AlignmentReader::Bam(reader)) => Records::Bam(reader.into_inner()),
                (selection, reader) => {
                    let query =
                        selection.query(input, reader.format(), &input_header, &translation)?;
                    let (tx, rx) = crossbeam_channel::bounded(2);
                    readers
                        .push(scope.spawn(move || read_records(reader, &input_header, query, &tx)));
                    Records::Channel {
                        receiver: rx,
                        batch: Vec::new(),
                        position: 0,
                    }
                }
            };

            sources.push(MergeInput {
                records,
                translation: RecordTranslation::new(translation, read_group),
            });
        }

        let mut writer = AlignmentWriter::create(Some(&output), &write_options, &header)?;
        writer.write_header(&header)?;

        runs::merge(&mut sources, order, |record| {
            writer.write_raw(&header, record)?;
            Ok(())
        })?;

        writer.finish(&header)?;

        for reader in readers {
            reader
                .join()
                .map_err(|_| anyhow::anyhow!("a reading thread panicked"))??;
        }

        Ok(())
    })
}

/// The merge order, as samtools picks it.
fn order(args: &Args) -> Order {
    let name_order = if args.name_natural {
        Some(TagTieBreak::NameNatural)
    } else if args.name_ascii {
        Some(TagTieBreak::NameAscii)
    } else {
        None
    };

    match (args.tag, name_order) {
        (Some(tag), Some(then)) => Order::Tag { tag, then },
        (Some(tag), None) => Order::Tag {
            tag,
            then: TagTieBreak::Coordinate,
        },
        (None, Some(TagTieBreak::NameAscii)) => Order::NameAscii,
        (None, Some(_)) => Order::NameNatural,
        (None, None) => Order::Coordinate,
    }
}

/// Splits the arguments into the output and the inputs (files from `-b` lists first).
fn files(args: &Args) -> anyhow::Result<(PathBuf, Vec<InputFile>)> {
    let mut positional = args.files.iter().cloned();

    let output = match &args.output {
        Some(output) => output.clone(),
        None => positional
            .next()
            .context("You must at least specify the output file")?,
    };

    let positional: Vec<PathBuf> = positional.collect();
    let mut inputs = Vec::new();

    for list in &args.lists {
        for path in read_file_list(list)? {
            inputs.push(InputFile { path, index: None });
        }
    }

    if args.index_files {
        if !positional.len().is_multiple_of(2) {
            bail!("Odd number of filenames detected! Each BAM file should have an index file");
        }

        let (paths, indexes) = positional.split_at(positional.len() / 2);

        inputs.extend(paths.iter().zip(indexes).map(|(path, index)| InputFile {
            path: path.clone(),
            index: Some(index.clone()),
        }));
    } else {
        inputs.extend(
            positional
                .into_iter()
                .map(|path| InputFile { path, index: None }),
        );
    }

    if inputs.is_empty() {
        bail!("You must specify at least one (and usually two or more) input files");
    }

    Ok((output, inputs))
}

/// Reads a list of file names, one per line, skipping empty lines (`hts_readlines`). A gzipped
/// list is decompressed; a name that is not a file but starts with `:` is a comma-separated
/// list.
pub fn read_file_list(src: &Path) -> anyhow::Result<Vec<PathBuf>> {
    let file = match File::open(src) {
        Ok(file) => file,
        Err(e) => {
            if let Some(list) = src.to_str().and_then(|s| s.strip_prefix(':')) {
                return Ok(list.split(',').map(PathBuf::from).collect());
            }

            return Err(e).with_context(|| format!("Invalid file list \"{}\"", src.display()));
        }
    };

    let mut reader = BufReader::new(file);
    let is_gzipped = reader.fill_buf()?.starts_with(&[0x1f, 0x8b]);

    let mut text = Vec::new();

    if is_gzipped {
        MultiGzDecoder::new(reader).read_to_end(&mut text)?;
    } else {
        reader.read_to_end(&mut text)?;
    }

    // Lines end at `\n`, with any `\r` before it removed.
    let paths = ByteSlice::lines(text.as_slice())
        .filter(|line| !line.is_empty())
        .map(|line| line.to_path().map(Path::to_path_buf))
        .collect::<Result<_, _>>()
        .with_context(|| format!("Invalid file list \"{}\"", src.display()))?;

    Ok(paths)
}

/// The read group name `-r` gives an input: its file name without the directory and without a
/// `.bam`, `.sam`, or `.cram` extension.
fn read_group_name(path: &Path) -> BString {
    let s = path.as_os_str().as_encoded_bytes();
    let mut len = s.len();

    if len > 4 && (s.ends_with(b".bam") || s.ends_with(b".sam")) {
        len -= 4;
    }

    // samtools compares the rest of the whole name, so this only applies if nothing was removed.
    if len > 5 && len == s.len() && s.ends_with(b".cram") {
        len -= 5;
    }

    let start = s[..len]
        .iter()
        .rposition(|&b| b == b'/')
        .map_or(0, |i| i + 1);

    BString::from(&s[start..len])
}

/// Which records of each input to merge.
enum Selection {
    All,
    /// `-R`: a reference sequence of the merged header and a 0-based, half-open interval.
    Region(usize, i64, i64),
    /// `-L`: the regions of a BED file.
    Bed(PathBuf),
}

/// What a reading thread reads from its input.
enum Query {
    All,
    Interval(AlignmentIndex, Option<(usize, usize, Option<usize>)>),
    Set(AlignmentIndex, RegionSet),
}

impl Selection {
    fn new(args: &Args, header: &sam::Header) -> anyhow::Result<Self> {
        if let Some(region) = &args.region {
            let reference_sequences = header.reference_sequences();

            let (id, start, end) =
                hts_region::parse_region(region, |name| reference_sequences.get_index_of(name))
                    .with_context(|| {
                        format!("Badly formatted region or unknown reference name: \"{region}\"")
                    })?;

            Ok(Self::Region(id, start, end))
        } else if let Some(bed) = &args.bed {
            if !has_bed_records(bed)? {
                bail!("Null or empty region list");
            }

            Ok(Self::Bed(bed.clone()))
        } else {
            Ok(Self::All)
        }
    }

    /// The query for one input, loading its index.
    fn query(
        &self,
        input: &InputFile,
        format: Format,
        input_header: &sam::Header,
        translation: &Translation,
    ) -> anyhow::Result<Query> {
        let index = || {
            AlignmentIndex::read(format, &input.path, input.index.as_deref()).with_context(|| {
                format!(
                    "failed to load index for {}. Random alignment retrieval only works for \
                     indexed BAM or CRAM files.",
                    input.path.display()
                )
            })
        };

        Ok(match self {
            Self::All => Query::All,
            Self::Region(id, start, end) => {
                // The input's reference sequence that maps to the merged one, if any.
                let input_id = translation
                    .reference_sequence_ids
                    .iter()
                    .position(|&merged_id| usize::try_from(merged_id).ok() == Some(*id));

                let start = usize::try_from(*start).unwrap_or(0);
                let end = (*end != POS_MAX).then(|| usize::try_from(*end).unwrap_or(usize::MAX));

                Query::Interval(index()?, input_id.map(|id| (id, start, end)))
            }
            Self::Bed(path) => {
                let mut set = RegionSet::new(input_header);
                set.add_bed(input_header, path)?;
                set.merge();
                Query::Set(index()?, set)
            }
        })
    }
}

/// Whether a BED file has any region lines.
fn has_bed_records(src: &Path) -> anyhow::Result<bool> {
    let text = std::fs::read_to_string(src)
        .with_context(|| format!("Could not read BED file: \"{}\"", src.display()))?;

    Ok(text.lines().any(|line| {
        let line = line.trim_end();
        !(line.is_empty()
            || line.starts_with('#')
            || line.starts_with("track")
            || line.starts_with("browser"))
    }))
}

/// Reads an input's selected records on its own thread, sending them as batches of
/// `block_size`-prefixed raw records.
fn read_records(
    mut reader: AlignmentReader,
    header: &sam::Header,
    query: Query,
    tx: &Sender<Vec<u8>>,
) -> anyhow::Result<()> {
    let mut encoder = RawEncoder::new();
    let mut batch = Vec::with_capacity(BATCH_SIZE);

    // A closed channel means the merge stopped early; the merge reports why.
    let mut send = |record: crate::alignment::Rec<'_>| -> anyhow::Result<()> {
        let raw = record.encode_raw(header, &mut encoder)?;
        runs::write_raw_record(&mut batch, raw)?;

        if batch.len() >= BATCH_SIZE {
            let full = std::mem::replace(&mut batch, Vec::with_capacity(BATCH_SIZE));
            tx.send(full)
                .map_err(|_| anyhow::anyhow!("merge stopped"))?;
        }

        Ok(())
    };

    match query {
        Query::All => reader.for_each(header, &mut send),
        Query::Interval(_, None) => Ok(()),
        Query::Interval(index, Some((id, start, end))) => {
            reader.for_each_in_interval(&index, header, id, start, end, &mut send)
        }
        Query::Set(index, set) => {
            reader.for_each_in_indexed_queries(&index, header, &Queries::Set(&set), &mut send)
        }
    }?;

    if !batch.is_empty() {
        tx.send(batch).ok();
    }

    Ok(())
}

/// The records of one input.
enum Records {
    /// A BAM file read directly.
    Bam(BgzfReader<BufReader<Stream>>),
    /// Batches of raw records from a reading thread.
    Channel {
        receiver: Receiver<Vec<u8>>,
        batch: Vec<u8>,
        position: usize,
    },
}

impl Records {
    fn next_record(&mut self, buf: &mut Vec<u8>) -> io::Result<bool> {
        match self {
            Self::Bam(reader) => runs::read_raw_record(reader, buf),
            Self::Channel {
                receiver,
                batch,
                position,
            } => loop {
                let mut rest = &batch[*position..];

                if runs::read_raw_record(&mut rest, buf)? {
                    *position = batch.len() - rest.len();
                    return Ok(true);
                }

                match receiver.recv() {
                    Ok(next) => {
                        *batch = next;
                        *position = 0;
                    }
                    Err(_) => return Ok(false),
                }
            },
        }
    }
}

/// A sorted input, with its records translated to the merged header.
struct MergeInput {
    records: Records,
    translation: RecordTranslation,
}

impl RecordSource for MergeInput {
    fn next_record(&mut self, buf: &mut Vec<u8>) -> io::Result<bool> {
        if !self.records.next_record(buf)? {
            return Ok(false);
        }

        self.translation.apply(buf)?;

        Ok(true)
    }
}

/// How an input's records change in the merged output (`bam_translate`).
struct RecordTranslation {
    reference_sequence_ids: Option<Vec<i32>>,
    /// Merged read group IDs by input ID; `None` for IDs not in the header (already reported).
    read_groups: HashMap<BString, Option<BString>>,
    programs: HashMap<BString, Option<BString>>,
    /// The `-r` read group of every record.
    read_group: Option<BString>,
    tail: Vec<u8>,
}

impl RecordTranslation {
    fn new(translation: Translation, read_group: Option<BString>) -> Self {
        let Translation {
            reference_sequence_ids,
            read_groups,
            programs,
            ..
        } = translation;

        let is_identity = reference_sequence_ids
            .iter()
            .enumerate()
            .all(|(i, &id)| usize::try_from(id).is_ok_and(|id| id == i));

        let to_options =
            |map: HashMap<BString, BString>| map.into_iter().map(|(k, v)| (k, Some(v))).collect();

        Self {
            reference_sequence_ids: (!is_identity).then_some(reference_sequence_ids),
            read_groups: to_options(read_groups),
            programs: to_options(programs),
            read_group,
            tail: Vec::new(),
        }
    }

    /// Translates reference sequence IDs and the `RG` and `PG` tags, which samtools moves to
    /// the end (`RG` before `PG`); with `-r`, `RG` is replaced and goes last.
    fn apply(&mut self, record: &mut Vec<u8>) -> io::Result<()> {
        if let Some(map) = &self.reference_sequence_ids {
            remap_reference_sequence_ids(record, map);
        }

        let [rg, pg] = bam_raw::find_fields(record, [*b"RG", *b"PG"])?;

        if rg.is_none() && pg.is_none() && self.read_group.is_none() {
            return Ok(());
        }

        self.tail.clear();

        if let Some(range) = &rg {
            let id = translate_tag(&mut self.read_groups, record, range, "RG");

            if self.read_group.is_none()
                && let Some(id) = id
            {
                push_string_field(&mut self.tail, *b"RG", &id);
            }
        }

        if let Some(range) = &pg {
            let id = translate_tag(&mut self.programs, record, range, "PG");

            if let Some(id) = id {
                push_string_field(&mut self.tail, *b"PG", &id);
            }
        }

        if let Some(id) = &self.read_group {
            push_string_field(&mut self.tail, *b"RG", id);
        }

        let mut ranges: Vec<_> = [rg, pg].into_iter().flatten().collect();
        ranges.sort_by_key(|range| std::cmp::Reverse(range.start));

        for range in ranges {
            record.drain(range);
        }

        record.extend_from_slice(&self.tail);

        Ok(())
    }
}

/// Looks up a record's `RG` or `PG` value, warning once per ID that the header lacks.
fn translate_tag(
    map: &mut HashMap<BString, Option<BString>>,
    record: &[u8],
    range: &std::ops::Range<usize>,
    tag: &str,
) -> Option<BString> {
    let field = &record[range.clone()];

    let value: &BStr = match field[2] {
        b'Z' | b'H' => field[3..field.len() - 1].as_bstr(),
        _ => b"".as_bstr(),
    };

    if let Some(id) = map.get(value) {
        return id.clone();
    }

    eprintln!(
        "[bam_translate] {tag} tag \"{value}\" on read \"{}\" encountered with no \
         corresponding entry in header, tag lost. Unknown tags are only reported once per \
         input file for each tag ID.",
        crate::sort_order::name(record).as_bstr()
    );

    map.insert(value.into(), None);

    None
}

fn push_string_field(dst: &mut Vec<u8>, tag: [u8; 2], value: &[u8]) {
    dst.extend_from_slice(&tag);
    dst.push(b'Z');
    dst.extend_from_slice(value);
    dst.push(0);
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_read_group_name() {
        assert_eq!(read_group_name(Path::new("dir/a.bam")), "a");
        assert_eq!(read_group_name(Path::new("/x/y/b.cram")), "b");
        assert_eq!(read_group_name(Path::new("c.sam.cram")), "c.sam");
        assert_eq!(read_group_name(Path::new("d.txt")), "d.txt");
    }
}
