//! `ntools markdup`: mark duplicate alignments (samtools markdup).
//!
//! This follows `bam_mark_duplicates` in samtools 1.24: reads are held in a coordinate window,
//! keyed by the unclipped 5' positions of the read (and its mate, from `MC`), and the read with
//! the highest base quality score (plus the `ms` mate score from `fixmate -m`) is kept. Hash
//! entries are removed only by the read that inserted them, as in samtools.
//!
//! Not yet supported: `-S` (supplementary marking), `-t` and `-d` (duplicate chains, optical
//! duplicates), and barcodes.

use std::{
    collections::{HashMap, VecDeque},
    fs::File,
    io::{self, Write},
    path::PathBuf,
};

use anyhow::{Context, bail};
use noodles_bam as bam;
use noodles_sam::{
    self as sam,
    alignment::{
        io::Write as _,
        record::{
            Flags,
            cigar::op::Kind,
            data::field::{Tag, Value},
        },
    },
    header::record::value::map::header::{sort_order::QUERY_NAME, tag::SORT_ORDER},
};

use crate::{
    alignment::{AlignmentReader, ReadOptions, Rec},
    bam_raw::{self, field, le_u16, le_u32},
    output::{self, AlignmentWriter, OutputFormat, WriteOptions},
    program,
};

/// Base qualities below this are not scored (`MD_MIN_QUALITY`).
const MIN_SCORED_QUALITY: u8 = 15;

// Orientations and read sides, the values samtools hashes.
const O_FF: i8 = 2;
const O_RR: i8 = 3;
const O_FR: i8 = 5;
const O_RF: i8 = 7;
const R_LE: i8 = 11;
const R_RI: i8 = 13;

const MATE_CIGAR: Tag = Tag::new(b'M', b'C');
const MATE_SCORE: Tag = Tag::new(b'm', b's');
const DUPLICATE_COUNT: Tag = Tag::new(b'd', b'c');

#[derive(Clone, Copy, Debug, Eq, PartialEq, clap::ValueEnum)]
enum Mode {
    /// Template start and end positions.
    T,
    /// Sequence start positions.
    S,
}

#[derive(Debug, clap::Args)]
pub struct Args {
    /// Coordinate-sorted input, after `fixmate -m`.
    input: String,

    /// Output file.
    output: PathBuf,

    /// Remove duplicate reads.
    #[arg(short = 'r')]
    remove: bool,

    /// Maximum read length.
    #[arg(short = 'l', default_value_t = 300)]
    max_length: i64,

    /// Report statistics on stderr.
    #[arg(short = 's')]
    stats: bool,

    /// Write statistics to this file (implies -s).
    #[arg(short = 'f')]
    stats_file: Option<PathBuf>,

    /// Clear existing duplicate flags and dt/do tags.
    #[arg(short = 'c')]
    clear: bool,

    /// Duplicate decision method for pairs.
    #[arg(short = 'm', long = "mode", value_enum, default_value_t = Mode::T)]
    mode: Mode,

    /// Output uncompressed BAM.
    #[arg(short = 'u')]
    uncompressed: bool,

    /// Include QC-failed reads.
    #[arg(long = "include-fails")]
    include_fails: bool,

    /// Use read groups in duplicate matching.
    #[arg(long = "use-read-groups")]
    use_read_groups: bool,

    /// Record each kept read's duplicate count (including itself) in a dc tag.
    #[arg(long = "duplicate-count")]
    duplicate_count: bool,

    /// Output format [default: from the output extension, else BAM].
    #[arg(short = 'O', long = "output-fmt", ignore_case = true)]
    output_format: Option<OutputFormat>,

    /// Reference FASTA (indexed) for CRAM input or output.
    #[arg(long = "reference")]
    reference: Option<PathBuf>,

    /// Additional threads for compression and decompression.
    #[arg(short = '@', long = "threads", default_value_t = 0)]
    threads: usize,

    /// Do not add an @PG line to the header.
    #[arg(long = "no-PG")]
    no_pg: bool,
}

/// `key_data_t`. Single-read keys leave the mate fields zero, so derived equality matches
/// samtools' `key_equal`.
#[derive(Clone, Copy, Debug, Default, Eq, Hash, PartialEq)]
struct Key {
    this_coord: i64,
    other_coord: i64,
    this_ref: i32,
    other_ref: i32,
    read_group: i32,
    single: bool,
    leftmost: i8,
    orientation: i8,
}

/// `read_queue_t`
struct Read {
    record: bam::Record,
    /// The window coordinate: the 5' unclipped position for examined reads, else `pos`.
    pos: i64,
    tid: i32,
    /// Keys this read inserted, which it removes when written.
    pair_key: Option<Key>,
    single_key: Option<Key>,
    duplicate_count: i64,
    read_group: usize,
}

#[derive(Debug, Default, Clone, Copy)]
struct Stats {
    reading: i64,
    writing: i64,
    excluded: i64,
    examined: i64,
    pair: i64,
    single: i64,
    duplicate: i64,
    single_duplicate: i64,
}

struct Markdup {
    options: Options,
    reads: VecDeque<Read>,
    /// ID of `reads[0]`; IDs increase with input order.
    front_id: u64,
    pair_hash: HashMap<Key, u64>,
    single_hash: HashMap<Key, u64>,
    stats: Vec<Stats>,
}

#[derive(Clone, Copy, Debug)]
struct Options {
    remove: bool,
    max_length: i64,
    mode: Mode,
    duplicate_count: bool,
}

pub fn run(args: Args) -> anyhow::Result<()> {
    let read_options = ReadOptions {
        threads: args.threads,
        reference: args.reference.clone(),
    };

    let (mut reader, _) = AlignmentReader::open(&args.input, &read_options)?;
    let mut header = reader.read_header()?;

    let is_queryname_sorted = header
        .header()
        .and_then(|hd| hd.other_fields().get(&SORT_ORDER))
        .is_some_and(|so| so == QUERY_NAME);

    if is_queryname_sorted {
        bail!("[markdup] error, queryname sorted, must be sorted by coordinate.");
    }

    let read_groups: HashMap<Vec<u8>, usize> = if args.use_read_groups {
        header
            .read_groups()
            .keys()
            .enumerate()
            .map(|(i, id)| (id.to_vec(), i + 1))
            .collect()
    } else {
        HashMap::new()
    };

    let group_count = read_groups.len();

    if !args.no_pg {
        program::add(&mut header)?;
    }

    let format = args
        .output_format
        .or_else(|| OutputFormat::from_path(&args.output))
        .unwrap_or(OutputFormat::Bam);

    let write_options = WriteOptions {
        format,
        level: output::parse_level(if args.uncompressed { 0 } else { 6 })?,
        threads: args.threads,
        reference: args.reference.clone(),
    };

    let mut writer = AlignmentWriter::create(Some(&args.output), &write_options, &header)?;
    writer.write_header(&header)?;

    let mut state = Markdup {
        options: Options {
            remove: args.remove,
            max_length: args.max_length,
            mode: args.mode,
            duplicate_count: args.duplicate_count,
        },
        reads: VecDeque::new(),
        front_id: 0,
        pair_hash: HashMap::new(),
        single_hash: HashMap::new(),
        stats: vec![Stats::default(); group_count + 1],
    };

    let exclude = if args.include_fails {
        Flags::SECONDARY | Flags::SUPPLEMENTARY | Flags::UNMAPPED
    } else {
        Flags::SECONDARY | Flags::SUPPLEMENTARY | Flags::UNMAPPED | Flags::QC_FAIL
    };

    let mut prev_tid: i32 = 0;
    let mut prev_coord: i64 = 0;
    let mut encoder = bam::io::Writer::from(Vec::new());
    let mut cleared = Vec::new();

    reader.for_each(&header.clone(), |rec| {
        let raw: Vec<u8> = match &rec {
            Rec::Bam(r) => r.as_ref().to_vec(),
            _ => {
                encoder.get_mut().clear();
                encoder.write_alignment_record(&header, rec.as_dyn())?;
                encoder.get_ref()[size_of::<u32>()..].to_vec()
            }
        };

        let mut raw = raw;

        if args.clear && flags_of(&raw).is_duplicate() {
            let flags = flags_of(&raw) - Flags::DUPLICATE;
            set_flags(&mut raw, flags);
            bam_raw::strip_fields(&raw, &mut cleared, |tag| tag == *b"dt" || tag == *b"do")?;
            std::mem::swap(&mut raw, &mut cleared);
        }

        let record = bam::Record::try_from(raw)?;

        let tid = i32_field(record.as_ref(), field::REFERENCE_SEQUENCE_ID);
        let pos = i64::from(i32_field(record.as_ref(), field::POSITION));

        if tid >= 0 && (tid < prev_tid || (tid == prev_tid && pos < prev_coord)) {
            bail!("[markdup] error, input is not coordinate sorted.");
        }

        prev_coord = pos;
        prev_tid = tid;

        let read_group = if args.use_read_groups {
            read_group_of(&record, &read_groups)?
        } else {
            0
        };

        state.stats[read_group].reading += 1;

        let id = state.front_id + state.reads.len() as u64;
        state.reads.push_back(Read {
            record,
            pos,
            tid,
            pair_key: None,
            single_key: None,
            duplicate_count: 1,
            read_group,
        });

        let flags = flags_of(state.reads.back().expect("just pushed").record.as_ref());

        if flags.intersects(exclude) {
            state.stats[read_group].excluded += 1;
        } else {
            state.stats[read_group].examined += 1;
            state.examine(id, read_group)?;
        }

        state.flush(prev_tid, prev_coord, &header, &mut writer)?;

        Ok(())
    })?;

    state.flush_all(&header, &mut writer)?;
    writer.finish(&header)?;

    if args.stats || args.stats_file.is_some() {
        let mut out: Box<dyn Write> = match &args.stats_file {
            Some(path) => Box::new(
                File::create(path).with_context(|| format!("creating {}", path.display()))?,
            ),
            None => Box::new(io::stderr()),
        };

        let command_line = std::env::args().collect::<Vec<_>>().join(" ");

        let total = state.stats.iter().fold(Stats::default(), |mut t, s| {
            t.reading += s.reading;
            t.writing += s.writing;
            t.excluded += s.excluded;
            t.examined += s.examined;
            t.pair += s.pair;
            t.single += s.single;
            t.duplicate += s.duplicate;
            t.single_duplicate += s.single_duplicate;
            t
        });

        write_stats(&mut out, "COMMAND: ", &command_line, &total)?;
        writeln!(out)?;

        if args.use_read_groups {
            if state.stats[0].reading > 0 {
                write_stats(&mut out, "READ GROUP: ", "ungrouped", &state.stats[0])?;
                writeln!(out)?;
            }

            for (id, stats) in header.read_groups().keys().zip(&state.stats[1..]) {
                write_stats(&mut out, "READ GROUP: ", &id.to_string(), stats)?;
                writeln!(out)?;
            }
        }
    }

    Ok(())
}

impl Markdup {
    fn read(&self, id: u64) -> Option<&Read> {
        id.checked_sub(self.front_id)
            .and_then(|i| self.reads.get(usize::try_from(i).ok()?))
    }

    fn read_mut(&mut self, id: u64) -> Option<&mut Read> {
        id.checked_sub(self.front_id)
            .and_then(|i| self.reads.get_mut(usize::try_from(i).ok()?))
    }

    /// Keys the read and resolves duplicates against earlier reads with the same key.
    fn examine(&mut self, id: u64, read_group: usize) -> anyhow::Result<()> {
        let read_group_key = i32::try_from(read_group).unwrap_or(i32::MAX);
        let record = &self.read(id).expect("examined read is buffered").record;

        if has_mate(record) {
            let pair_key = make_pair_key(record, self.options.mode, read_group_key)?;
            let single_key = make_single_key(record, read_group_key);

            self.stats[read_group].pair += 1;
            self.read_mut(id).expect("buffered").pos = single_key.this_coord;

            // A paired read replaces a single read at its position, which becomes a duplicate.
            match self.existing(&self.single_hash, &single_key) {
                None => {
                    self.single_hash.insert(single_key, id);
                    self.read_mut(id).expect("buffered").single_key = Some(single_key);
                }
                Some(other) => {
                    if !has_mate(&self.read(other).expect("buffered").record) {
                        self.single_hash.insert(single_key, id);
                        let other_count = self.read(other).expect("buffered").duplicate_count;
                        let _ = other_count;
                        self.read_mut(id).expect("buffered").duplicate_count += 1;
                        self.mark_duplicate(other);
                        self.stats[read_group].single_duplicate += 1;
                    }
                }
            }

            match self.existing(&self.pair_hash, &pair_key) {
                None => {
                    self.pair_hash.insert(pair_key, id);
                    self.read_mut(id).expect("buffered").pair_key = Some(pair_key);
                }
                Some(other) => {
                    let old = &self.read(other).expect("buffered").record;
                    let new = &self.read(id).expect("buffered").record;

                    let (old_score, new_score) = if flags_of(old.as_ref()).is_qc_fail()
                        != flags_of(new.as_ref()).is_qc_fail()
                    {
                        if flags_of(old.as_ref()).is_qc_fail() {
                            (0, 1)
                        } else {
                            (1, 0)
                        }
                    } else {
                        (
                            read_score(old.as_ref()) + mate_score(old)?,
                            read_score(new.as_ref()) + mate_score(new)?,
                        )
                    };

                    let tie_add = if new_score == old_score {
                        if name(new.as_ref()) < name(old.as_ref()) {
                            1
                        } else {
                            -1
                        }
                    } else {
                        0
                    };

                    if new_score + tie_add > old_score {
                        let other_count = self.read(other).expect("buffered").duplicate_count;
                        self.read_mut(id).expect("buffered").duplicate_count += other_count;
                        self.pair_hash.insert(pair_key, id);
                        self.mark_duplicate(other);
                    } else {
                        self.read_mut(other).expect("buffered").duplicate_count += 1;
                        self.mark_duplicate(id);
                    }

                    self.stats[read_group].duplicate += 1;
                }
            }
        } else {
            let single_key = make_single_key(record, read_group_key);

            self.stats[read_group].single += 1;
            self.read_mut(id).expect("buffered").pos = single_key.this_coord;

            match self.existing(&self.single_hash, &single_key) {
                None => {
                    self.single_hash.insert(single_key, id);
                    self.read_mut(id).expect("buffered").single_key = Some(single_key);
                }
                Some(other) => {
                    if has_mate(&self.read(other).expect("buffered").record) {
                        self.read_mut(other).expect("buffered").duplicate_count += 1;
                        self.mark_duplicate(id);
                    } else {
                        let old_score =
                            read_score(self.read(other).expect("buffered").record.as_ref());
                        let new_score =
                            read_score(self.read(id).expect("buffered").record.as_ref());

                        if new_score > old_score {
                            let other_count = self.read(other).expect("buffered").duplicate_count;
                            self.read_mut(id).expect("buffered").duplicate_count += other_count;
                            self.single_hash.insert(single_key, id);
                            self.mark_duplicate(other);
                        } else {
                            self.read_mut(other).expect("buffered").duplicate_count += 1;
                            self.mark_duplicate(id);
                        }
                    }

                    self.stats[read_group].single_duplicate += 1;
                }
            }
        }

        Ok(())
    }

    /// The read a key maps to, if that read is still buffered.
    fn existing(&self, hash: &HashMap<Key, u64>, key: &Key) -> Option<u64> {
        hash.get(key).copied().filter(|&id| self.read(id).is_some())
    }

    fn mark_duplicate(&mut self, id: u64) {
        let read = self.read_mut(id).expect("duplicate is buffered");
        let mut raw = std::mem::take(&mut read.record).as_ref().to_vec();
        let flags = flags_of(&raw) | Flags::DUPLICATE;
        set_flags(&mut raw, flags);
        read.record = bam::Record::try_from(raw).expect("flag edits keep records valid");
    }

    /// Writes reads that can no longer gain duplicates: those behind the window.
    fn flush(
        &mut self,
        prev_tid: i32,
        prev_coord: i64,
        header: &sam::Header,
        writer: &mut AlignmentWriter,
    ) -> anyhow::Result<()> {
        while let Some(front) = self.reads.front() {
            let in_window = front.pos + self.options.max_length > prev_coord
                && front.tid == prev_tid
                && (prev_tid != -1 || prev_coord != -1);

            if in_window {
                break;
            }

            self.write_front(header, writer)?;
        }

        Ok(())
    }

    fn flush_all(
        &mut self,
        header: &sam::Header,
        writer: &mut AlignmentWriter,
    ) -> anyhow::Result<()> {
        while !self.reads.is_empty() {
            self.write_front(header, writer)?;
        }

        Ok(())
    }

    fn write_front(
        &mut self,
        header: &sam::Header,
        writer: &mut AlignmentWriter,
    ) -> anyhow::Result<()> {
        let read = self.reads.pop_front().expect("front exists");
        self.front_id += 1;

        let is_duplicate = flags_of(read.record.as_ref()).is_duplicate();

        if !self.options.remove || !is_duplicate {
            if self.options.duplicate_count && !is_duplicate {
                let mut raw = Vec::new();
                bam_raw::strip_fields(read.record.as_ref(), &mut raw, |tag| tag == *b"dc")?;
                append_integer_tag(&mut raw, DUPLICATE_COUNT, read.duplicate_count);
                writer.write_raw(header, &raw)?;
            } else {
                writer.write_raw(header, read.record.as_ref())?;
            }

            self.stats[read.read_group].writing += 1;
        }

        if let Some(key) = read.pair_key {
            self.pair_hash.remove(&key);
        }

        if let Some(key) = read.single_key {
            self.single_hash.remove(&key);
        }

        Ok(())
    }
}

fn i32_field(record: &[u8], range: std::ops::Range<usize>) -> i32 {
    le_u32(&record[range]).cast_signed()
}

fn flags_of(record: &[u8]) -> Flags {
    Flags::from_bits_retain(le_u16(&record[field::FLAGS]))
}

fn set_flags(record: &mut [u8], flags: Flags) {
    record[field::FLAGS].copy_from_slice(&flags.bits().to_le_bytes());
}

fn name(record: &[u8]) -> &[u8] {
    crate::sort_order::name(record)
}

/// `calc_score`: the sum of base qualities of at least 15 (missing qualities count as 255).
fn read_score(record: &[u8]) -> i64 {
    let Ok(data) = bam_raw::data_range(record) else {
        return 0;
    };

    let base_count = usize::try_from(le_u32(&record[field::SEQUENCE_LENGTH])).unwrap_or(0);
    let quality_scores = &record[data.start - base_count..data.start];

    quality_scores
        .iter()
        .filter(|&&q| q >= MIN_SCORED_QUALITY)
        .map(|&q| i64::from(q))
        .sum()
}

/// `get_mate_score`, from the `ms` tag.
fn mate_score(record: &bam::Record) -> anyhow::Result<i64> {
    let data = record.data();

    match sam::alignment::record::Data::get(&data, &MATE_SCORE).transpose()? {
        Some(value) => value
            .as_int()
            .context("[markdup] error, ms tag is not an integer"),
        None => {
            bail!("[markdup] error, no ms score tag. Please run samtools fixmate on file first.")
        }
    }
}

fn read_group_of(record: &bam::Record, read_groups: &HashMap<Vec<u8>, usize>) -> io::Result<usize> {
    let data = record.data();

    Ok(
        match sam::alignment::record::Data::get(&data, &Tag::READ_GROUP).transpose()? {
            Some(Value::String(id)) => read_groups.get(id.as_ref() as &[u8]).copied().unwrap_or(0),
            _ => 0,
        },
    )
}

/// `has_mate`: the read's mate is present and mapped.
fn has_mate(record: &bam::Record) -> bool {
    let raw = record.as_ref();
    let flags = flags_of(raw);
    let mtid = i32_field(raw, field::MATE_REFERENCE_SEQUENCE_ID);
    let mpos = i32_field(raw, 24..28);

    flags.is_segmented() && !flags.is_mate_unmapped() && !(mtid == -1 && mpos == -1)
}

/// The CIGAR as `(kind, len)` pairs.
fn cigar_ops(record: &bam::Record) -> io::Result<Vec<(Kind, i64)>> {
    record
        .cigar()
        .iter()
        .map(|result| result.map(|op| (op.kind(), i64::try_from(op.len()).unwrap_or(i64::MAX))))
        .collect()
}

/// HTSlib's `bam_endpos`.
fn end_pos(record: &bam::Record, ops: &[(Kind, i64)]) -> i64 {
    let raw = record.as_ref();
    let pos = i64::from(i32_field(raw, field::POSITION));

    let span: i64 = if flags_of(raw).is_unmapped() {
        0
    } else {
        ops.iter()
            .filter(|(kind, _)| kind.consumes_reference())
            .map(|(_, len)| len)
            .sum()
    };

    pos + span.max(1)
}

fn is_clip(kind: Kind) -> bool {
    matches!(kind, Kind::SoftClip | Kind::HardClip)
}

/// `unclipped_start(b, 1)`
fn unclipped_start(record: &bam::Record, ops: &[(Kind, i64)]) -> i64 {
    let pos = i64::from(i32_field(record.as_ref(), field::POSITION));
    let clipped: i64 = ops
        .iter()
        .take_while(|(kind, _)| is_clip(*kind))
        .map(|(_, len)| len)
        .sum();

    pos - clipped + 1
}

/// `unclipped_end(b, 1)`
fn unclipped_end(record: &bam::Record, ops: &[(Kind, i64)]) -> i64 {
    let clipped: i64 = ops
        .iter()
        .rev()
        .take_while(|(kind, _)| is_clip(*kind))
        .map(|(_, len)| len)
        .sum();

    end_pos(record, ops) + clipped
}

/// Parses a CIGAR string as samtools' `unclipped_other_*` do: numbers via `strtol`, and an
/// operation without a number counts 1.
fn parse_cigar_string(s: &[u8]) -> Vec<(u8, i64)> {
    let mut ops = Vec::new();
    let mut i = 0;

    while i < s.len() && s[i] != b'*' {
        let start = i;

        while i < s.len() && s[i].is_ascii_digit() {
            i += 1;
        }

        let n = if i > start {
            std::str::from_utf8(&s[start..i])
                .ok()
                .and_then(|t| t.parse().ok())
                .unwrap_or(i64::MAX)
        } else {
            1
        };

        if i < s.len() {
            ops.push((s[i], n));
            i += 1;
        }
    }

    ops
}

/// `unclipped_other_start(mpos, MC, 1)`
fn unclipped_other_start(mpos: i64, cigar: &[(u8, i64)]) -> i64 {
    let mut clipped = 0;

    for &(op, n) in cigar {
        match op {
            b'S' | b'H' => clipped += n,
            _ => break,
        }
    }

    mpos - clipped + 1
}

/// `unclipped_other_end(mpos, MC, 1)`
fn unclipped_other_end(mpos: i64, cigar: &[(u8, i64)]) -> i64 {
    let mut reference_position = 0;
    let mut skip = true;

    for &(op, n) in cigar {
        match op {
            b'M' | b'D' | b'N' | b'=' | b'X' => {
                reference_position += n;
                skip = false;
            }
            b'S' | b'H' if !skip => reference_position += n,
            _ => {}
        }
    }

    mpos + reference_position
}

/// `make_single_key`
fn make_single_key(record: &bam::Record, read_group: i32) -> Key {
    let ops = cigar_ops(record).unwrap_or_default();
    let raw = record.as_ref();
    let this_ref = i32_field(raw, field::REFERENCE_SEQUENCE_ID) + 1;

    let (this_coord, orientation) = if flags_of(raw).is_reverse_complemented() {
        (unclipped_end(record, &ops), O_RR)
    } else {
        (unclipped_start(record, &ops), O_FF)
    };

    Key {
        this_coord,
        this_ref,
        read_group,
        single: true,
        orientation,
        ..Key::default()
    }
}

/// `make_pair_key`
fn make_pair_key(record: &bam::Record, mode: Mode, read_group: i32) -> anyhow::Result<Key> {
    let ops = cigar_ops(record)?;
    let raw = record.as_ref();
    let flags = flags_of(raw);

    let this_ref = i32_field(raw, field::REFERENCE_SEQUENCE_ID) + 1;
    let other_ref = i32_field(raw, field::MATE_REFERENCE_SEQUENCE_ID) + 1;
    let mpos = i64::from(i32_field(raw, 24..28));
    let pos = i64::from(i32_field(raw, field::POSITION));

    let mut this_coord = unclipped_start(record, &ops);
    let this_end = unclipped_end(record, &ops);

    let data = record.data();
    let mate_cigar = match sam::alignment::record::Data::get(&data, &MATE_CIGAR).transpose()? {
        Some(Value::String(s)) => parse_cigar_string(s.as_ref()),
        Some(_) => bail!(
            "[markdup] error, MC tag wrong type. Please use the MC tag provided by samtools fixmate."
        ),
        None => bail!("[markdup] error, no MC tag. Please run samtools fixmate on file first."),
    };

    let other_end = unclipped_other_end(mpos, &mate_cigar);
    let mut other_coord = unclipped_other_start(mpos, &mate_cigar);

    let is_rev = flags.is_reverse_complemented();
    let is_mrev = flags.is_mate_reverse_complemented();
    let is_read1 = flags.is_first_segment();

    let leftmost;
    let orientation;

    match mode {
        Mode::T => {
            leftmost = if this_ref != other_ref {
                this_ref < other_ref
            } else if is_rev == is_mrev {
                if is_rev {
                    this_end <= other_end
                } else {
                    this_coord <= other_coord
                }
            } else if is_rev {
                this_end <= other_coord
            } else {
                this_coord <= other_end
            };

            if leftmost {
                if is_rev == is_mrev {
                    other_coord = other_end;
                    orientation = match (is_rev, is_read1) {
                        (false, true) | (true, false) => O_FF,
                        (false, false) | (true, true) => O_RR,
                    };
                } else if !is_rev {
                    orientation = O_FR;
                    other_coord = other_end;
                } else {
                    orientation = O_RF;
                    this_coord = this_end;
                }
            } else if is_rev == is_mrev {
                this_coord = this_end;
                orientation = match (is_rev, is_read1) {
                    (false, true) | (true, false) => O_RR,
                    (false, false) | (true, true) => O_FF,
                };
            } else if !is_rev {
                orientation = O_RF;
                other_coord = other_end;
            } else {
                orientation = O_FR;
                this_coord = this_end;
            }
        }
        Mode::S => {
            let difference = if this_ref != other_ref {
                i64::from(this_ref - other_ref)
            } else if is_rev == is_mrev {
                if is_rev {
                    this_end - other_end
                } else {
                    this_coord - other_coord
                }
            } else if is_rev {
                this_end - other_coord
            } else {
                this_coord - other_end
            };

            leftmost = match difference.cmp(&0) {
                std::cmp::Ordering::Less => true,
                std::cmp::Ordering::Greater => false,
                std::cmp::Ordering::Equal => {
                    if pos == mpos {
                        is_read1
                    } else {
                        pos < mpos
                    }
                }
            };

            orientation = match (leftmost, is_rev == is_mrev, is_rev) {
                (true, true, false) => O_FF,
                (true, true, true) => O_RR,
                (true, false, false) => O_FR,
                (true, false, true) => O_RF,
                (false, true, false) => O_RR,
                (false, true, true) => O_FF,
                (false, false, false) => O_RF,
                (false, false, true) => O_FR,
            };

            this_coord = if is_rev {
                unclipped_end(record, &ops)
            } else {
                unclipped_start(record, &ops)
            };

            other_coord = if is_mrev {
                unclipped_other_end(mpos, &mate_cigar)
            } else {
                unclipped_other_start(mpos, &mate_cigar)
            };
        }
    }

    Ok(Key {
        this_coord,
        other_coord,
        this_ref,
        other_ref,
        read_group,
        single: false,
        leftmost: if leftmost { R_LE } else { R_RI },
        orientation,
    })
}

/// Appends an integer tag with the smallest type that holds it, as HTSlib does.
fn append_integer_tag(record: &mut Vec<u8>, tag: Tag, n: i64) {
    let tag: [u8; 2] = tag.into();
    record.extend_from_slice(&tag);

    if let Ok(n) = u8::try_from(n) {
        record.push(b'C');
        record.push(n);
    } else if let Ok(n) = i8::try_from(n) {
        record.push(b'c');
        record.extend_from_slice(&n.to_le_bytes());
    } else if let Ok(n) = u16::try_from(n) {
        record.push(b'S');
        record.extend_from_slice(&n.to_le_bytes());
    } else if let Ok(n) = i16::try_from(n) {
        record.push(b's');
        record.extend_from_slice(&n.to_le_bytes());
    } else if let Ok(n) = u32::try_from(n) {
        record.push(b'I');
        record.extend_from_slice(&n.to_le_bytes());
    } else {
        record.push(b'i');
        let n = i32::try_from(n).unwrap_or(i32::MAX);
        record.extend_from_slice(&n.to_le_bytes());
    }
}

/// The Picard library size estimate samtools reports.
fn estimate_library_size(paired_reads: i64, paired_duplicate_reads: i64, optical: i64) -> i64 {
    let coverage = |x: f64, c: f64, n: f64| c / x - 1.0 + (-n / x).exp();

    let non_optical_pairs = (paired_reads - optical) / 2;
    let unique_pairs = (paired_reads - paired_duplicate_reads) / 2;
    let duplicate_pairs = (paired_duplicate_reads - optical) / 2;

    if non_optical_pairs == 0
        || duplicate_pairs == 0
        || unique_pairs == 0
        || non_optical_pairs <= duplicate_pairs
    {
        eprintln!(
            "[markdup] warning, unable to calculate estimated library size. Read pairs \
             {non_optical_pairs} should be greater than duplicate pairs {duplicate_pairs}, which \
             should both be non zero."
        );
        return 0;
    }

    let (u, n) = (unique_pairs as f64, non_optical_pairs as f64);
    let (mut m, mut big_m) = (1.0, 100.0);

    if coverage(m * u, u, n) < 0.0 {
        eprintln!("[markdup] warning, unable to calculate estimated library size.");
        return 0;
    }

    while coverage(big_m * u, u, n) > 0.0 {
        big_m *= 10.0;
    }

    for _ in 0..40 {
        let r = (m + big_m) / 2.0;
        let value = coverage(r * u, u, n);

        if value > 0.0 {
            m = r;
        } else if value < 0.0 {
            big_m = r;
        } else {
            break;
        }
    }

    #[expect(
        clippy::cast_possible_truncation,
        reason = "samtools truncates the estimate to an unsigned long the same way"
    )]
    let size = (u * (m + big_m) / 2.0) as i64;
    size
}

/// `write_stats`. The library size is estimated (possibly warning) before the title is written.
fn write_stats<W: Write>(out: &mut W, title: &str, name: &str, stats: &Stats) -> io::Result<()> {
    let library_size = estimate_library_size(stats.pair, stats.duplicate, 0);

    writeln!(out, "{title}{name}")?;

    writeln!(out, "READ: {}", stats.reading)?;
    writeln!(out, "WRITTEN: {}", stats.writing)?;
    writeln!(out, "EXCLUDED: {}", stats.excluded)?;
    writeln!(out, "EXAMINED: {}", stats.examined)?;
    writeln!(out, "PAIRED: {}", stats.pair)?;
    writeln!(out, "SINGLE: {}", stats.single)?;
    writeln!(out, "DUPLICATE PAIR: {}", stats.duplicate)?;
    writeln!(out, "DUPLICATE SINGLE: {}", stats.single_duplicate)?;
    writeln!(out, "DUPLICATE PAIR OPTICAL: 0")?;
    writeln!(out, "DUPLICATE SINGLE OPTICAL: 0")?;
    writeln!(out, "DUPLICATE NON PRIMARY: 0")?;
    writeln!(out, "DUPLICATE NON PRIMARY OPTICAL: 0")?;
    writeln!(
        out,
        "DUPLICATE PRIMARY TOTAL: {}",
        stats.single_duplicate + stats.duplicate
    )?;
    writeln!(
        out,
        "DUPLICATE TOTAL: {}",
        stats.single_duplicate + stats.duplicate
    )?;
    writeln!(out, "ESTIMATED_LIBRARY_SIZE: {library_size}")
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_unclipped_other_positions() {
        let cigar = parse_cigar_string(b"5S10M2D3M4S2H");
        assert_eq!(unclipped_other_start(100, &cigar), 96);
        assert_eq!(unclipped_other_end(100, &cigar), 100 + 10 + 2 + 3 + 4 + 2);
        assert!(parse_cigar_string(b"*").is_empty());
    }
}
