//! `ntools rmdup`: remove PCR duplicates (samtools rmdup).
//!
//! This follows samtools 1.24's `bam_rmdup.c` and `bam_rmdupse.c`. Input must be sorted by
//! coordinate. Duplicates are compared within a library (the read group's `LB`), and the record
//! with the highest sum of base qualities is kept (the first one on ties).
//!
//! - Paired-end mode (default): pairs whose first-in-coordinate record ("head", positive TLEN)
//!   has the same position and TLEN are duplicates. The removed heads' mates are removed when
//!   they are reached, by name. Unpaired records, pairs with an unmapped read, and pairs across
//!   reference sequences are kept as they are.
//! - Single-end mode (`-s`, or `-S` to treat paired reads as single): forward reads with the
//!   same start, and reverse reads with the same end, are duplicates.
//!
//! samtools' single-end mode holds up to 2^20 records before writing any (and then may modify
//! records it has already freed); ntools writes each record as soon as no later read can
//! replace it, in the same order, so memory stays small.

use std::{
    collections::{HashMap, HashSet, VecDeque, hash_map::Entry},
    path::PathBuf,
};

use anyhow::bail;
use noodles_bam as bam;
use noodles_sam::{
    self as sam, alignment::io::Write as _, header::record::value::map::read_group::tag::LIBRARY,
};

use crate::{
    alignment::{AlignmentReader, ReadOptions, Rec},
    bam_raw::{self, field, le_u16, le_u32},
    output::{self, AlignmentWriter, OutputFormat, WriteOptions},
};

const FLAG_PAIRED: u16 = 0x1;
const FLAG_UNMAPPED: u16 = 0x4;
const FLAG_MATE_UNMAPPED: u16 = 0x8;
const FLAG_REVERSE: u16 = 0x10;

/// The library name samtools reports for records without one.
const NO_LIBRARY: &str = "\t";

#[derive(Debug, clap::Args)]
pub struct Args {
    /// Coordinate-sorted input.
    input: String,

    /// Output file (BAM unless the extension or -O says otherwise).
    output: PathBuf,

    /// Remove duplicates of single-end reads (paired reads are kept as they are).
    #[arg(short = 's')]
    single_end: bool,

    /// Treat paired reads as single-end reads (implies -s).
    #[arg(short = 'S')]
    force_single_end: bool,

    /// Output format.
    #[arg(short = 'O', long = "output-fmt", ignore_case = true)]
    output_format: Option<OutputFormat>,

    /// FASTA reference for CRAM input and output.
    #[arg(short = 'T', long = "reference")]
    reference: Option<PathBuf>,

    /// Additional threads for compression and decompression.
    #[arg(short = '@', long = "threads", default_value_t = 0)]
    threads: usize,
}

/// Which records `-s`/`-S` deduplicate.
#[derive(Clone, Copy, PartialEq, Eq)]
enum SingleEnd {
    /// Only unpaired reads.
    Unpaired,
    /// Every mapped read (`-S`).
    All,
}

pub fn run(args: Args) -> anyhow::Result<()> {
    let read_options = ReadOptions {
        threads: args.threads,
        reference: args.reference.clone(),
    };

    let (mut reader, _) = AlignmentReader::open(&args.input, &read_options)?;
    let header = reader.read_header()?;

    if header.reference_sequences().is_empty() {
        bail!("[bam_rmdup] input SAM does not have header. Abort!");
    }

    let write_options = WriteOptions {
        format: args
            .output_format
            .or_else(|| OutputFormat::from_path(&args.output))
            .unwrap_or(OutputFormat::Bam),
        level: output::parse_level(6)?,
        threads: args.threads,
        reference: args.reference.clone(),
    };

    let mut writer = AlignmentWriter::create(Some(&args.output), &write_options, &header)?;
    writer.write_header(&header)?;

    let libraries = Libraries::new(&header);
    let mut encoder = bam::io::Writer::from(Vec::new());

    let single_end = if args.force_single_end {
        Some(SingleEnd::All)
    } else {
        args.single_end.then_some(SingleEnd::Unpaired)
    };

    let mut write = |record: &[u8]| writer.write_raw(&header, record);

    let stats = match single_end {
        Some(mode) => {
            let mut dedup = SingleEndDedup::new(libraries, mode);
            reader.for_each(&header, |rec| {
                let record = raw_record(&header, &rec, &mut encoder)?;
                dedup.push(record, &mut write)
            })?;
            dedup.finish(&mut write)?
        }
        None => {
            let mut dedup = PairedDedup::new(libraries, &header);
            reader.for_each(&header, |rec| {
                let record = raw_record(&header, &rec, &mut encoder)?;
                dedup.push(record, &mut write)
            })?;
            dedup.finish(&mut write)?
        }
    };

    let tag = if single_end.is_some() {
        "bam_rmdupse_core"
    } else {
        "bam_rmdup_core"
    };

    for (name, s) in stats.names.iter().zip(&stats.counts) {
        eprintln!(
            "[{tag}] {} / {} = {:.4} in library '{name}'",
            s.removed,
            s.checked,
            s.removed as f64 / s.checked as f64
        );
    }

    writer.finish(&header)?;

    Ok(())
}

/// Returns a record's raw BAM encoding, encoding non-BAM records into `encoder`.
fn raw_record<'r>(
    header: &sam::Header,
    rec: &'r Rec<'_>,
    encoder: &'r mut bam::io::Writer<Vec<u8>>,
) -> anyhow::Result<&'r [u8]> {
    Ok(match rec {
        Rec::Bam(r) => r.as_ref(),
        _ => {
            encoder.get_mut().clear();
            encoder.write_alignment_record(header, rec.as_dyn())?;
            &encoder.get_ref()[size_of::<u32>()..]
        }
    })
}

/// Fixed fields of a raw BAM record.
struct Fields {
    reference_sequence_id: i32,
    position: i32,
    flags: u16,
    mate_reference_sequence_id: i32,
    template_length: i32,
}

impl Fields {
    fn new(record: &[u8]) -> Self {
        let i32_at =
            |i: usize| i32::from_le_bytes([record[i], record[i + 1], record[i + 2], record[i + 3]]);

        Self {
            reference_sequence_id: i32_at(field::REFERENCE_SEQUENCE_ID.start),
            position: i32_at(field::POSITION.start),
            flags: le_u16(&record[field::FLAGS]),
            mate_reference_sequence_id: i32_at(field::MATE_REFERENCE_SEQUENCE_ID.start),
            template_length: i32_at(28),
        }
    }
}

/// The read name, without its NUL terminator.
fn name(record: &[u8]) -> &[u8] {
    let len = usize::from(record[field::NAME_LENGTH]);
    &record[field::NAME..field::NAME + len - 1]
}

/// The sum of base qualities (samtools' `sum_qual`).
fn quality_sum(record: &[u8]) -> u64 {
    let name_len = usize::from(record[field::NAME_LENGTH]);
    let n_cigar = usize::from(le_u16(&record[field::CIGAR_OP_COUNT]));
    let l_seq = le_u32(&record[field::SEQUENCE_LENGTH]) as usize;
    let start = field::NAME + name_len + 4 * n_cigar + l_seq.div_ceil(2);
    record[start..start + l_seq]
        .iter()
        .map(|&q| u64::from(q))
        .sum()
}

/// The 0-based exclusive end position (HTSlib's `bam_endpos`): at least one past the start.
fn end_position(
    record: &[u8],
    flags: u16,
    position: i32,
    cigar: &mut Vec<u32>,
) -> anyhow::Result<i64> {
    let mut len = 0;

    if flags & FLAG_UNMAPPED == 0 {
        bam_raw::read_cigar(record, cigar)?;
        len = cigar
            .iter()
            .filter(|&&op| matches!(op & 0xf, 0 | 2 | 3 | 7 | 8))
            .map(|&op| i64::from(op >> 4))
            .sum();
    }

    Ok(i64::from(position) + len.max(1))
}

/// Libraries by read group, with duplicate counts.
struct Libraries {
    /// Library indices by read group ID.
    by_read_group: HashMap<Vec<u8>, usize>,
    /// Library indices by library name.
    by_name: HashMap<String, usize>,
    stats: Stats,
}

#[derive(Default)]
struct Stats {
    /// Library names, in the order they are first seen.
    names: Vec<String>,
    counts: Vec<Counts>,
}

#[derive(Clone, Copy, Default)]
struct Counts {
    checked: u64,
    removed: u64,
}

impl Libraries {
    fn new(header: &sam::Header) -> Self {
        let mut by_name = HashMap::new();
        let mut by_read_group = HashMap::new();
        let mut names = Vec::new();

        for (id, read_group) in header.read_groups() {
            if let Some(library) = read_group.other_fields().get(&LIBRARY) {
                let library = library.to_string();
                let next = by_name.len();
                let i = *by_name.entry(library.clone()).or_insert_with(|| {
                    names.push(library);
                    next
                });
                by_read_group.insert(id.to_vec(), i);
            }
        }

        let next = by_name.len();
        by_name.insert(NO_LIBRARY.into(), next);
        names.push(NO_LIBRARY.into());

        Self {
            by_read_group,
            by_name,
            stats: Stats {
                counts: vec![Counts::default(); names.len()],
                names,
            },
        }
    }

    /// Returns the library index of a record (HTSlib's `bam_get_library`), counting it as checked.
    fn check(&mut self, record: &[u8]) -> anyhow::Result<usize> {
        let library = bam_raw::find_field(record, *b"RG")?
            .filter(|range| record[range.start + 2] == b'Z')
            .and_then(|range| {
                self.by_read_group
                    .get(&record[range.start + 3..range.end - 1])
            })
            .copied()
            .unwrap_or(self.by_name[NO_LIBRARY]);

        self.stats.counts[library].checked += 1;

        Ok(library)
    }

    fn remove(&mut self, library: usize) {
        self.stats.counts[library].removed += 1;
    }

    /// The libraries with checked records, in the order they were defined.
    fn into_stats(self) -> Stats {
        let (names, counts) = self
            .stats
            .names
            .into_iter()
            .zip(self.stats.counts)
            .filter(|(_, c)| c.checked > 0)
            .unzip();

        Stats { names, counts }
    }
}

/// Paired-end deduplication (samtools' `bam_rmdup_core`).
struct PairedDedup<'h> {
    header: &'h sam::Header,
    libraries: Libraries,
    last_reference_sequence_id: i32,
    last_position: i32,
    /// Once unplaced unmapped records are reached, the rest are copied.
    copying: bool,
    /// The best head at the current position for each library and TLEN, in arrival order.
    best: Vec<Vec<u8>>,
    best_index: HashMap<(usize, i32), usize>,
    /// Names of removed heads whose mates are still to come.
    removed: HashSet<Vec<u8>>,
}

impl<'h> PairedDedup<'h> {
    fn new(libraries: Libraries, header: &'h sam::Header) -> Self {
        Self {
            header,
            libraries,
            last_reference_sequence_id: -1,
            last_position: -1,
            copying: false,
            best: Vec::new(),
            best_index: HashMap::new(),
            removed: HashSet::new(),
        }
    }

    fn push<W>(&mut self, record: &[u8], write: &mut W) -> anyhow::Result<()>
    where
        W: FnMut(&[u8]) -> std::io::Result<()>,
    {
        if self.copying {
            return Ok(write(record)?);
        }

        let fields = Fields::new(record);

        if fields.reference_sequence_id != self.last_reference_sequence_id
            || fields.position != self.last_position
        {
            self.flush(write)?;

            if fields.reference_sequence_id != self.last_reference_sequence_id {
                if !self.removed.is_empty() {
                    eprintln!("[bam_rmdup_core] {} unmatched pairs", self.removed.len());
                    self.removed.clear();
                }

                if fields.reference_sequence_id == -1 {
                    self.copying = true;
                    return Ok(write(record)?);
                }

                self.last_reference_sequence_id = fields.reference_sequence_id;

                let name = usize::try_from(fields.reference_sequence_id)
                    .ok()
                    .and_then(|i| self.header.reference_sequences().get_index(i))
                    .map(|(name, _)| name.to_string())
                    .unwrap_or_default();
                eprintln!("[bam_rmdup_core] processing reference {name}...");
            }
        }

        let flags = fields.flags;

        if flags & FLAG_PAIRED == 0
            || flags & (FLAG_UNMAPPED | FLAG_MATE_UNMAPPED) != 0
            || (fields.mate_reference_sequence_id >= 0
                && fields.reference_sequence_id != fields.mate_reference_sequence_id)
        {
            write(record)?;
        } else if fields.template_length > 0 {
            let library = self.libraries.check(record)?;

            match self.best_index.entry((library, fields.template_length)) {
                Entry::Occupied(entry) => {
                    self.libraries.remove(library);
                    let best = &mut self.best[*entry.get()];

                    let removed_name = if quality_sum(best) < quality_sum(record) {
                        let removed_name = name(best).to_vec();
                        best.clear();
                        best.extend_from_slice(record);
                        removed_name
                    } else {
                        name(record).to_vec()
                    };

                    if !self.removed.insert(removed_name) {
                        eprintln!(
                            "[bam_rmdup_core] inconsistent BAM file for pair '{}'. Continue anyway.",
                            String::from_utf8_lossy(name(record))
                        );
                    }
                }
                Entry::Vacant(entry) => {
                    entry.insert(self.best.len());
                    self.best.push(record.to_vec());
                }
            }
        } else if !self.removed.remove(name(record)) {
            write(record)?;
        }

        self.last_position = fields.position;

        Ok(())
    }

    /// Writes the best heads of the current position.
    fn flush<W>(&mut self, write: &mut W) -> anyhow::Result<()>
    where
        W: FnMut(&[u8]) -> std::io::Result<()>,
    {
        for record in self.best.drain(..) {
            write(&record)?;
        }

        self.best_index.clear();

        Ok(())
    }

    fn finish<W>(mut self, write: &mut W) -> anyhow::Result<Stats>
    where
        W: FnMut(&[u8]) -> std::io::Result<()>,
    {
        self.flush(write)?;
        Ok(self.libraries.into_stats())
    }
}

/// How a queued single-end record can still change.
#[derive(Clone, Copy)]
enum Pending {
    /// Not a duplicate candidate; written in turn.
    Fixed,
    /// A forward read, replaceable by a better read with the same start.
    Forward,
    /// A reverse read, discardable for a better read with the same end.
    Reverse,
    /// Discarded for a better reverse read.
    Discarded,
}

/// A duplicate key: library, reverse strand, and start (forward) or end (reverse).
type Key = (usize, bool, i64);

struct Queued {
    record: Vec<u8>,
    pending: Pending,
    /// The record's duplicate key, for candidates.
    key: Option<Key>,
    position: i32,
    end: i64,
    score: u64,
}

/// Single-end deduplication (samtools' `bam_rmdupse_core`).
struct SingleEndDedup {
    libraries: Libraries,
    mode: SingleEnd,
    last_reference_sequence_id: Option<i32>,
    /// Records in input order; replaced forward reads keep their place.
    queue: VecDeque<Queued>,
    /// The sequence number of the queue's first record.
    queue_start: u64,
    /// The best read for each library, strand, and start (forward) or end (reverse).
    best: HashMap<Key, u64>,
    cigar: Vec<u32>,
}

impl SingleEndDedup {
    fn new(libraries: Libraries, mode: SingleEnd) -> Self {
        Self {
            libraries,
            mode,
            last_reference_sequence_id: None,
            queue: VecDeque::new(),
            queue_start: 0,
            best: HashMap::new(),
            cigar: Vec::new(),
        }
    }

    fn push<W>(&mut self, record: &[u8], write: &mut W) -> anyhow::Result<()>
    where
        W: FnMut(&[u8]) -> std::io::Result<()>,
    {
        let fields = Fields::new(record);

        if self.last_reference_sequence_id != Some(fields.reference_sequence_id) {
            self.flush(None, write)?;
            self.best.clear();
            self.last_reference_sequence_id = Some(fields.reference_sequence_id);
        } else {
            self.flush(Some(fields.position), write)?;
        }

        let flags = fields.flags;
        let end = end_position(record, flags, fields.position, &mut self.cigar)?;
        let score = quality_sum(record);

        let candidate =
            flags & FLAG_UNMAPPED == 0 && (flags & FLAG_PAIRED == 0 || self.mode == SingleEnd::All);

        if !candidate {
            self.enqueue(record, Pending::Fixed, None, fields.position, end, score);
            return Ok(());
        }

        let library = self.libraries.check(record)?;
        let reverse = flags & FLAG_REVERSE != 0;
        let key = (
            library,
            reverse,
            if reverse {
                end
            } else {
                i64::from(fields.position)
            },
        );

        let Some(&seq) = self.best.get(&key) else {
            let pending = if reverse {
                Pending::Reverse
            } else {
                Pending::Forward
            };
            let seq = self.enqueue(record, pending, Some(key), fields.position, end, score);
            self.best.insert(key, seq);
            return Ok(());
        };

        self.libraries.remove(library);

        let i = usize::try_from(seq - self.queue_start)?;
        let best = &mut self.queue[i];

        if best.score < score {
            if reverse {
                best.pending = Pending::Discarded;
                let seq = self.enqueue(
                    record,
                    Pending::Reverse,
                    Some(key),
                    fields.position,
                    end,
                    score,
                );
                self.best.insert(key, seq);
            } else {
                best.record.clear();
                best.record.extend_from_slice(record);
                best.end = end;
                best.score = score;
            }
        }

        Ok(())
    }

    fn enqueue(
        &mut self,
        record: &[u8],
        pending: Pending,
        key: Option<Key>,
        position: i32,
        end: i64,
        score: u64,
    ) -> u64 {
        self.queue.push_back(Queued {
            record: record.to_vec(),
            pending,
            key,
            position,
            end,
            score,
        });

        self.queue_start + self.queue.len() as u64 - 1
    }

    /// Writes queued records from the front while no read at `position` or later can replace
    /// them (all of them for `None`).
    fn flush<W>(&mut self, position: Option<i32>, write: &mut W) -> anyhow::Result<()>
    where
        W: FnMut(&[u8]) -> std::io::Result<()>,
    {
        while let Some(front) = self.queue.front() {
            let settled = match (front.pending, position) {
                (_, None) | (Pending::Fixed | Pending::Discarded, _) => true,
                (Pending::Forward, Some(p)) => front.position < p,
                (Pending::Reverse, Some(p)) => front.end <= i64::from(p),
            };

            if !settled {
                break;
            }

            let front = self.queue.pop_front().expect("front exists");
            let seq = self.queue_start;
            self.queue_start += 1;

            // A settled record's entry can no longer match a later read.
            if let Some(key) = front.key
                && self.best.get(&key) == Some(&seq)
            {
                self.best.remove(&key);
            }

            if !matches!(front.pending, Pending::Discarded) {
                write(&front.record)?;
            }
        }

        Ok(())
    }

    fn finish<W>(mut self, write: &mut W) -> anyhow::Result<Stats>
    where
        W: FnMut(&[u8]) -> std::io::Result<()>,
    {
        self.flush(None, write)?;
        Ok(self.libraries.into_stats())
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    const HEADER: &str = "@SQ\tSN:chr\tLN:1000\n@RG\tID:a\tLB:x\n@RG\tID:b\tLB:y\n";

    fn encode(header: &sam::Header, line: &str) -> anyhow::Result<Vec<u8>> {
        let record = sam::Record::try_from(line.as_bytes())?;
        let mut writer = bam::io::Writer::from(Vec::new());
        writer.write_alignment_record(header, &record)?;
        Ok(writer.get_ref()[size_of::<u32>()..].to_vec())
    }

    /// Runs `lines` through single-end (`Some`) or paired-end deduplication and returns the
    /// names of the records written, in order.
    fn dedup(mode: Option<SingleEnd>, lines: &[&str]) -> anyhow::Result<Vec<String>> {
        let header: sam::Header = HEADER.parse()?;
        let records = lines
            .iter()
            .map(|line| encode(&header, line))
            .collect::<anyhow::Result<Vec<_>>>()?;

        let mut out = Vec::new();
        let mut write = |record: &[u8]| {
            out.push(String::from_utf8_lossy(name(record)).into_owned());
            Ok(())
        };

        match mode {
            Some(mode) => {
                let mut d = SingleEndDedup::new(Libraries::new(&header), mode);
                for record in &records {
                    d.push(record, &mut write)?;
                }
                d.finish(&mut write)?;
            }
            None => {
                let mut d = PairedDedup::new(Libraries::new(&header), &header);
                for record in &records {
                    d.push(record, &mut write)?;
                }
                d.finish(&mut write)?;
            }
        }

        Ok(out)
    }

    #[test]
    fn test_single_end() -> anyhow::Result<()> {
        let names = dedup(
            Some(SingleEnd::Unpaired),
            &[
                // Forward reads starting at 10: the better f2 takes f1's place.
                "f1\t0\tchr\t10\t60\t4M\t*\t0\t0\tACGT\t####",
                "f2\t0\tchr\t10\t60\t4M\t*\t0\t0\tACGT\tIIII",
                // Another library is compared separately.
                "f3\t0\tchr\t10\t60\t4M\t*\t0\t0\tACGT\tIIII\tRG:Z:b",
                // Reverse reads ending at 15: the better r2 is written in its own place.
                "r1\t16\tchr\t11\t60\t5M\t*\t0\t0\tACGTA\t#####",
                "p1\t1\tchr\t11\t60\t4M\t=\t11\t0\tACGT\t####",
                "r2\t16\tchr\t12\t60\t4M\t*\t0\t0\tACGT\tIIII",
                // An equal score keeps the first read.
                "r3\t16\tchr\t13\t60\t3M\t*\t0\t0\tCGT\tIII\t",
            ]
            .map(|line| line.trim_end_matches('\t')),
        )?;

        assert_eq!(names, ["f2", "f3", "p1", "r2"]);

        Ok(())
    }

    #[test]
    fn test_paired_end() -> anyhow::Result<()> {
        let names = dedup(
            None,
            &[
                "a\t99\tchr\t10\t60\t4M\t=\t20\t14\tACGT\t####",
                "b\t99\tchr\t10\t60\t4M\t=\t20\t14\tACGT\tIIII",
                "u\t0\tchr\t10\t60\t4M\t*\t0\t0\tACGT\tIIII",
                "c\t99\tchr\t10\t60\t4M\t=\t30\t24\tACGT\t####",
                "a\t147\tchr\t20\t60\t4M\t=\t10\t-14\tACGT\t####",
                "b\t147\tchr\t20\t60\t4M\t=\t10\t-14\tACGT\tIIII",
                "c\t147\tchr\t30\t60\t4M\t=\t10\t-24\tACGT\t####",
            ],
        )?;

        // Unpaired reads are written at once, ahead of the position's best heads.
        assert_eq!(names, ["u", "b", "c", "b", "c"]);

        Ok(())
    }
}
