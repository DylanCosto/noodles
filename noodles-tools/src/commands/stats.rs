//! `ntools stats`: alignment statistics for plot-bamstats (samtools stats).
//!
//! This is a port of `stats.c` from samtools 1.24, kept close to the original so that the
//! output matches it byte for byte: the same buffers and their growth rules, single-precision
//! (`float`) arithmetic where samtools uses it, and its quirks where they show in the output
//! (for example, the GC-depth table leaves out its last bin, and the coverage buffer keeps only
//! some of its contents when a long read enlarges it).

use std::{
    collections::HashMap,
    fs::{self, File},
    io::{self, BufWriter, Write},
    path::{Path, PathBuf},
};

use anyhow::{Context, bail};
use noodles_fasta::{self as fasta, fai};
use noodles_sam::{self as sam, alignment::record::Flags};

use crate::{
    alignment::{AlignmentReader, Queries, ReadOptions, Rec},
    bam_raw::{self, field, le_u16, le_u32},
    commands::faidx::read_bases,
    flags, output,
    regions::RegionSet,
};

/// The shortest read BWA trims (`BWA_MIN_RDLEN`).
const BWA_MIN_READ_LENGTH: usize = 35;
/// Read pairs held for overlap removal before stale ones are dropped (`DEFAULT_PAIR_MAX`).
const PAIR_MAX: u32 = 10_000;
/// The number of quality bins (`nquals`): every byte value, plus one for `qual + 1`.
const QUALITY_BINS: usize = 257;
/// The number of GC content bins (`ngc`).
const GC_BINS: usize = 200;
/// The initial per-cycle buffer length (`nbases`); also the longest indel counted.
const INITIAL_BASES: usize = 300;

/// The CIGAR operation codes.
const CIGAR_MATCH: u32 = 0;
const CIGAR_INSERTION: u32 = 1;
const CIGAR_DELETION: u32 = 2;
const CIGAR_SKIP: u32 = 3;
const CIGAR_SOFT_CLIP: u32 = 4;
const CIGAR_HARD_CLIP: u32 = 5;
const CIGAR_PAD: u32 = 6;
const CIGAR_SEQUENCE_MATCH: u32 = 7;
const CIGAR_SEQUENCE_MISMATCH: u32 = 8;

#[derive(Debug, clap::Args)]
#[command(disable_help_flag = true)]
pub struct Args {
    /// Input SAM, BAM, or CRAM file (`-` for stdin).
    input: String,

    /// Regions to restrict the statistics to (`chr`, `chr:from-to`); requires an index. With
    /// -X, the index file comes first.
    regions: Vec<String>,

    /// Coverage distribution min,max,step.
    #[arg(short = 'c', long = "coverage", default_value = "1,1000,1")]
    coverage: String,

    /// Exclude reads marked as duplicates.
    #[arg(short = 'd', long = "remove-dups")]
    remove_dups: bool,

    /// Use a customized index file, given after the input.
    #[arg(short = 'X', long = "customized-index")]
    customized_index: bool,

    /// Required flags.
    #[arg(short = 'f', long = "required-flag", value_parser = flags::parse)]
    required_flags: Option<Flags>,

    /// Filtering flags.
    #[arg(short = 'F', long = "filtering-flag", value_parser = flags::parse)]
    filtering_flags: Vec<Flags>,

    /// The size of GC-depth bins.
    #[arg(long = "GC-depth", default_value_t = 2e4)]
    gc_depth_bin_size: f64,

    /// Maximum insert size (0 for no limit).
    #[arg(short = 'i', long = "insert-size", default_value_t = 8000)]
    insert_size: i32,

    /// Include only this read group or sample name.
    #[arg(short = 'I', long = "id")]
    group_id: Option<String>,

    /// Include only reads of this length.
    #[arg(short = 'l', long = "read-length", default_value_t = -1)]
    read_length: i64,

    /// Report only the main part of inserts.
    #[arg(short = 'm', long = "most-inserts", default_value_t = 0.99)]
    most_inserts: f32,

    /// The BWA trimming parameter.
    #[arg(short = 'q', long = "trim-quality", default_value_t = 0)]
    trim_quality: i32,

    /// Indexed reference FASTA, for GC-depth and mismatches per cycle.
    #[arg(short = 'r', long = "ref-seq")]
    reference: Option<PathBuf>,

    /// Ignored (the input format is detected).
    #[arg(short = 's', long = "sam")]
    sam: bool,

    /// Also write statistics to separate files split by this tag.
    #[arg(short = 'S', long = "split")]
    split_tag: Option<String>,

    /// Path or prefix for the -S files [default: the input name].
    #[arg(short = 'P', long = "split-prefix")]
    split_prefix: Option<String>,

    /// Restrict the statistics to these regions (chr, from, to; 1-based, inclusive).
    #[arg(short = 't', long = "target-regions")]
    targets: Option<PathBuf>,

    /// Leave out insert size rows with no pairs.
    #[arg(short = 'x', long = "sparse")]
    sparse: bool,

    /// Remove overlaps of read pairs from coverage and base counts.
    #[arg(short = 'p', long = "remove-overlaps")]
    remove_overlaps: bool,

    /// Count only bases with coverage above this in the target percentage.
    #[arg(short = 'g', long = "cov-threshold", default_value_t = 0)]
    coverage_threshold: i32,

    /// Report statistics on the reference sequences.
    #[arg(long = "ref-stats")]
    ref_stats: bool,

    /// Reference retrieval chunk size for --ref-stats, in MiB.
    #[arg(long = "ref-stats-chunk", default_value_t = 1)]
    ref_stats_chunk: i64,

    /// Reference FASTA (indexed) for CRAM input.
    #[arg(long = "reference")]
    cram_reference: Option<PathBuf>,

    /// Additional decompression threads.
    #[arg(short = '@', long = "threads", default_value_t = 0)]
    threads: usize,

    /// Print help.
    #[arg(short = 'h', long, action = clap::ArgAction::Help)]
    help: Option<bool>,
}

/// `stats_info_t`: options shared by every set of statistics.
struct Info {
    flag_require: u16,
    flag_filter: u16,
    gcd_bin_size: i32,
    nisize: i32,
    trim_qual: i32,
    isize_main_bulk: f32,
    cov_min: i32,
    cov_max: i32,
    cov_step: i32,
    filter_readlen: i64,
    remove_overlaps: bool,
    cov_threshold: i32,
    ref_chunk_size: i64,
    reference: Option<Reference>,
    reference_names: Vec<Vec<u8>>,
    reference_lengths: Vec<u64>,
    /// Threads for per-cycle counts.
    workers: usize,
}

/// An indexed FASTA file, read through its index.
struct Reference {
    file: File,
    records: HashMap<Vec<u8>, fai::Record>,
}

impl Reference {
    fn open(path: &Path) -> anyhow::Result<Self> {
        let index = fasta::fs::read_associated_index(path)
            .with_context(|| format!("Could not load faidx: {}", path.display()))?;
        let records = index
            .as_ref()
            .iter()
            .map(|record| (record.name().to_vec(), record.clone()))
            .collect();
        let file = File::open(path).with_context(|| format!("opening {}", path.display()))?;
        Ok(Self { file, records })
    }

    /// `faidx_fetch_seq64` of the 0-based closed interval `[beg, end]`, clamped to the sequence
    /// as htslib does. Returns `None` if the sequence is not in the index.
    fn fetch(&self, name: &[u8], beg: i64, end: i64, dst: &mut Vec<u8>) -> anyhow::Result<bool> {
        let Some(record) = self.records.get(name) else {
            return Ok(false);
        };

        let length = i64::try_from(record.length())?;
        let mut beg = beg;
        let mut end = end;

        if end < beg {
            beg = end;
        }

        beg = beg.clamp(0, length);
        end = if end < 0 { 0 } else { end.min(length - 1) };

        dst.clear();

        if end >= beg {
            let mut raw = Vec::new();
            read_bases(
                &self.file,
                record,
                usize::try_from(beg + 1)?,
                usize::try_from(end + 1)?,
                &mut raw,
                dst,
            )?;
        }

        Ok(true)
    }
}

/// `gc_depth_t`
#[derive(Clone, Copy, Debug, Default)]
struct GcDepth {
    gc: f32,
    depth: u32,
}

/// `round_buffer_t`: per-position depth for the coverage distribution.
struct RoundBuffer {
    pos: i64,
    size: i32,
    start: i32,
    buffer: Vec<i32>,
}

/// Per-cycle counts of each 4-bit BAM base code (`=ACMGRSVTWYHKDBN`), indexed by the code so
/// that counting needs no branch. [`BaseCounts::from`] groups them as samtools does.
type CodeCounts = [u64; 16];

/// The complement of each base code for A, C, G, and T; other codes map to `=`, which is never
/// reported for read-oriented counts.
const COMPLEMENT: [u8; 16] = [0, 8, 4, 0, 2, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0];

/// The number of G and C bases in each byte of a 4-bit packed sequence.
const GC_IN_BYTE: [u8; 256] = {
    const fn is_gc(code: usize) -> u8 {
        (code == 2 || code == 4) as u8
    }

    let mut table = [0; 256];
    let mut b = 0;

    while b < 256 {
        table[b] = is_gc(b >> 4) + is_gc(b & 0xf);
        b += 1;
    }

    table
};

/// Bases of pending reads that trigger counting them (see [`PendingCycles`]).
const PENDING_BASES: usize = 1 << 24;

/// `acgtno_count_t`
#[derive(Clone, Copy, Debug, Default)]
struct BaseCounts {
    a: u64,
    c: u64,
    g: u64,
    t: u64,
    n: u64,
    other: u64,
}

impl From<&CodeCounts> for BaseCounts {
    fn from(counts: &CodeCounts) -> Self {
        let (a, c, g, t, n) = (counts[1], counts[2], counts[4], counts[8], counts[15]);
        let other = counts.iter().sum::<u64>() - a - c - g - t - n;
        Self {
            a,
            c,
            g,
            t,
            n,
            other,
        }
    }
}

/// A read whose per-cycle base and quality counts are pending, with its sequence and qualities in
/// [`PendingCycles::data`].
struct PendingRead {
    start: usize,
    len: usize,
    reverse: bool,
    first: bool,
}

/// Reads whose per-cycle counts are deferred so that worker threads can count them together,
/// each over its own range of cycles. The counts are sums, so this gives the same tables as
/// counting each read as it is read.
#[derive(Default)]
struct PendingCycles {
    data: Vec<u8>,
    reads: Vec<PendingRead>,
    max_len: usize,
}

/// The per-cycle tables for cycles `lo..lo + revcomp.len()`.
struct CycleTables<'a> {
    lo: usize,
    first: &'a mut [CodeCounts],
    last: &'a mut [CodeCounts],
    revcomp: &'a mut [CodeCounts],
    quals_first: &'a mut [u64],
    quals_last: &'a mut [u64],
    rows: usize,
}

/// Counts the bases and qualities of one first or last read at the cycles `tables` covers
/// (`collect_orig_read_stats`). Cycle `c` holds base `c` of the read as sequenced: from the end of
/// the stored sequence for a reverse-strand read.
fn count_cycles(seq: &[u8], qual: &[u8], reverse: bool, first: bool, tables: &mut CycleTables<'_>) {
    let len = qual.len();
    let lo = tables.lo;
    let hi = (lo + tables.revcomp.len()).min(len);

    if lo >= hi {
        return;
    }

    let (cycles, quals) = if first {
        (&mut *tables.first, &mut *tables.quals_first)
    } else {
        (&mut *tables.last, &mut *tables.quals_last)
    };

    let complement: [u8; 16] = if reverse {
        COMPLEMENT
    } else {
        std::array::from_fn(|i| i as u8)
    };

    let rows = tables.rows;

    for cycle in lo..hi {
        let i = if reverse { len - cycle - 1 } else { cycle };
        let code = usize::from((seq[i / 2] >> ((!i & 1) << 2)) & 0xf);
        let j = cycle - lo;

        cycles[j][code] += 1;
        tables.revcomp[j][usize::from(complement[code])] += 1;
        quals[j * rows + usize::from(qual[i])] += 1;
    }
}

/// Ends of up to `workers` cycle ranges covering `0..max_len` with about equal numbers of bases:
/// cycle `c` has a base from every read longer than `c`.
fn balanced_bounds(reads: &[PendingRead], max_len: usize, workers: usize) -> Vec<usize> {
    let mut ending = vec![0u64; max_len + 1];

    for read in reads {
        ending[read.len] += 1;
    }

    let total: u64 = reads.iter().map(|read| read.len as u64).sum();
    let target = total.div_ceil(workers as u64).max(1);

    let mut bounds = Vec::with_capacity(workers);
    let mut active = reads.len() as u64;
    let mut work = 0;

    for (cycle, &n) in ending.iter().enumerate().take(max_len) {
        active -= n;
        work += active;

        if work >= target && bounds.len() + 1 < workers {
            bounds.push(cycle + 1);
            work = 0;
        }
    }

    bounds.push(max_len);
    bounds.dedup();
    bounds
}

/// Target intervals of one reference sequence (`regions_t`), 1-based and closed.
#[derive(Clone, Debug, Default)]
struct Regions {
    intervals: Vec<(i64, i64)>,
    /// The first interval not yet passed.
    cpos: usize,
}

/// `isize_t`: insert size counts, inward, outward, and other.
enum InsertSizes {
    /// Sizes up to a bound (`-i`).
    Dense(Vec<[u64; 3]>),
    /// Unbounded (`-i 0`), and the largest size seen.
    Sparse(HashMap<i32, [u64; 3]>, i32),
}

impl InsertSizes {
    fn new(bound: i32) -> Self {
        if bound <= 0 {
            Self::Sparse(HashMap::new(), 0)
        } else {
            Self::Dense(vec![[0; 3]; usize::try_from(bound).unwrap_or(0)])
        }
    }

    fn len(&self) -> i32 {
        match self {
            Self::Dense(v) => i32::try_from(v.len()).unwrap_or(i32::MAX),
            Self::Sparse(_, max) => max + 1,
        }
    }

    fn get(&self, isize: i32) -> [u64; 3] {
        match self {
            Self::Dense(v) => v[usize::try_from(isize).expect("nonnegative")],
            Self::Sparse(map, _) => map.get(&isize).copied().unwrap_or_default(),
        }
    }

    fn increment(&mut self, isize: i32, orientation: usize) {
        match self {
            Self::Dense(v) => v[usize::try_from(isize).expect("nonnegative")][orientation] += 1,
            Self::Sparse(map, max) => {
                map.entry(isize).or_default()[orientation] += 1;
                *max = (*max).max(isize);
            }
        }
    }

    /// Halves every count (each pair was counted twice), truncating as samtools does.
    fn halve(&mut self) {
        let halve = |counts: &mut [u64; 3]| {
            for n in counts {
                #[expect(
                    clippy::cast_possible_truncation,
                    clippy::cast_sign_loss,
                    reason = "samtools stores `count * 0.5` back into a uint64_t"
                )]
                let half = (*n as f64 * 0.5) as u64;
                *n = half;
            }
        };

        match self {
            Self::Dense(v) => v.iter_mut().for_each(halve),
            Self::Sparse(map, _) => map.values_mut().for_each(halve),
        }
    }
}

const INWARD: usize = 0;
const OUTWARD: usize = 1;
const OTHER: usize = 2;

/// `pair_t`: the reference chunks a template's first-seen read covers, for `-p`.
struct Pair {
    first: u32,
    chunks: Vec<(i64, i64)>,
}

/// `refstats`
struct ReferenceStats {
    total_count: usize,
    count: i32,
    combined_length: i64,
    min_length: i64,
    max_length: i64,
    average_length: f32,
    average_gc: f32,
    rows: Vec<Option<(String, i64, f32, i64)>>,
}

/// `stats_t`
struct Stats {
    nbases: usize,
    ngc: usize,
    nindels: usize,

    /// Quality by cycle, `[cycle * quals_rows + quality]`, for `nbases` cycles.
    quals_1st: Vec<u64>,
    quals_2nd: Vec<u64>,
    quals_rows: usize,
    gc_1st: Vec<u64>,
    gc_2nd: Vec<u64>,
    acgtno_cycles_1st: Vec<CodeCounts>,
    acgtno_cycles_2nd: Vec<CodeCounts>,
    acgtno_revcomp: Vec<CodeCounts>,
    read_lengths: Vec<u64>,
    read_lengths_1st: Vec<u64>,
    read_lengths_2nd: Vec<u64>,
    insertions: Vec<u64>,
    deletions: Vec<u64>,
    ins_cycles_1st: Vec<u64>,
    ins_cycles_2nd: Vec<u64>,
    del_cycles_1st: Vec<u64>,
    del_cycles_2nd: Vec<u64>,
    isize: InsertSizes,
    mapping_qualities: [u64; 256],

    max_len: usize,
    max_len_1st: usize,
    max_len_2nd: usize,
    max_qual: usize,
    is_sorted: bool,

    total_len: u64,
    total_len_1st: u64,
    total_len_2nd: u64,
    total_len_dup: u64,
    nreads_1st: u64,
    nreads_2nd: u64,
    nreads_other: u64,
    nreads_filtered: u64,
    nreads_dup: u64,
    nreads_unmapped: u64,
    nreads_single_mapped: u64,
    nreads_paired_and_mapped: u64,
    nreads_properly_paired: u64,
    nreads_paired_tech: u64,
    nreads_anomalous: u64,
    nreads_mq0: u64,
    nbases_mapped: u64,
    nbases_mapped_cigar: u64,
    nbases_trimmed: u64,
    nmismatches: u64,
    nreads_qc_failed: u64,
    nreads_secondary: u64,
    nreads_supplementary: u64,
    checksum_names: u32,
    checksum_reads: u32,
    checksum_quals: u32,

    igcd: usize,
    gcd: Vec<GcDepth>,
    tid: i32,
    gcd_pos: i64,
    pos: i64,

    ncov: usize,
    cov: Vec<u64>,
    cov_rbuf: RoundBuffer,

    rseq_buf: Vec<u8>,
    mrseq_buf: i64,
    rseq_pos: i64,
    nrseq_buf: i64,
    mpc_buf: Option<Vec<u64>>,

    regions: Option<Vec<Regions>>,
    reg_from: i64,
    reg_to: i64,
    chunks: Vec<(i64, i64)>,
    nchunks: usize,
    target_count: u64,

    /// The sum of qualities; samtools sums in a double, which is exact for these integers.
    sum_qual: u64,
    pending: PendingCycles,
    /// Threads for per-cycle counts; 0 counts each read as it is read.
    workers: usize,
    read_groups: Option<Vec<Vec<u8>>>,
    split_name: Option<Vec<u8>>,

    pair_count: u32,
    last_pair_tid: i32,
    last_read_flush: u32,

    reference_stats: Option<ReferenceStats>,
}

/// The fixed fields and parts of an encoded BAM record.
struct Bam<'a> {
    tid: i32,
    pos: i64,
    mapq: u8,
    flag: u16,
    mtid: i32,
    mpos: i64,
    isize: i32,
    name: &'a [u8],
    cigar: Vec<u32>,
    seq: &'a [u8],
    qual: &'a [u8],
    l_seq: usize,
    raw: &'a [u8],
}

impl<'a> Bam<'a> {
    fn parse(raw: &'a [u8]) -> anyhow::Result<Self> {
        let i32_at = |offset: usize| le_u32(&raw[offset..offset + 4]).cast_signed();

        let name_len = usize::from(raw[field::NAME_LENGTH]);
        let n_cigar = usize::from(le_u16(&raw[field::CIGAR_OP_COUNT]));
        let l_seq = usize::try_from(le_u32(&raw[field::SEQUENCE_LENGTH]))?;

        let cigar_start = field::NAME + name_len;
        let seq_start = cigar_start + 4 * n_cigar;
        let qual_start = seq_start + l_seq.div_ceil(2);
        let data = bam_raw::data_range(raw)?;

        let mut cigar: Vec<u32> = raw[cigar_start..seq_start]
            .as_chunks::<4>()
            .0
            .iter()
            .map(|b| u32::from_le_bytes(*b))
            .collect();

        // As in htslib, a `kSmN` placeholder CIGAR is replaced by the CG tag's operations.
        if cigar.len() == 2
            && cigar[0] & 0xf == CIGAR_SOFT_CLIP
            && (cigar[0] >> 4) as usize == l_seq
            && cigar[1] & 0xf == CIGAR_SKIP
            && let Some(range) = bam_raw::find_field(raw, *b"CG")?
            && raw[range.start + 2] == b'B'
            && raw[range.start + 3] == b'I'
        {
            cigar = raw[range.start + 8..range.end]
                .as_chunks::<4>()
                .0
                .iter()
                .map(|b| u32::from_le_bytes(*b))
                .collect();
        }

        Ok(Self {
            tid: i32_at(0),
            pos: i64::from(i32_at(4)),
            mapq: raw[9],
            flag: le_u16(&raw[field::FLAGS]),
            mtid: i32_at(20),
            mpos: i64::from(i32_at(24)),
            isize: i32_at(28),
            name: &raw[field::NAME..cigar_start - 1],
            cigar,
            seq: &raw[seq_start..qual_start],
            qual: &raw[qual_start..data.start],
            l_seq,
            raw,
        })
    }

    fn has(&self, flag: Flags) -> bool {
        self.flag & flag.bits() != 0
    }

    fn base(&self, i: usize) -> u8 {
        (self.seq[i / 2] >> ((!i & 1) << 2)) & 0xf
    }

    fn is_paired_and_mapped(&self) -> bool {
        self.has(Flags::SEGMENTED) && !self.has(Flags::UNMAPPED) && !self.has(Flags::MATE_UNMAPPED)
    }

    fn is_original(&self) -> bool {
        !self.has(Flags::SECONDARY) && !self.has(Flags::SUPPLEMENTARY)
    }

    /// `READ_ORDER_*`: 1 for a first (or single) read, 2 for a last read, 0 or 3 otherwise.
    fn order(&self) -> u32 {
        if self.has(Flags::SEGMENTED) {
            u32::from(self.has(Flags::FIRST_SEGMENT)) + 2 * u32::from(self.has(Flags::LAST_SEGMENT))
        } else {
            READ_ORDER_FIRST
        }
    }

    /// htslib's `bam_endpos`.
    fn end_pos(&self) -> i64 {
        let span: i64 = if self.has(Flags::UNMAPPED) {
            0
        } else {
            self.cigar
                .iter()
                .filter(|&&op| consumes_reference(op & 0xf))
                .map(|&op| i64::from(op >> 4))
                .sum()
        };

        self.pos + span.max(1)
    }
}

const READ_ORDER_FIRST: u32 = 1;
const READ_ORDER_LAST: u32 = 2;

fn consumes_reference(op: u32) -> bool {
    matches!(
        op,
        CIGAR_MATCH | CIGAR_DELETION | CIGAR_SKIP | CIGAR_SEQUENCE_MATCH | CIGAR_SEQUENCE_MISMATCH
    )
}

/// htslib's `bam_aux2i`: the value of an integer field, or 0.
fn integer_field(raw: &[u8], tag: [u8; 2]) -> anyhow::Result<Option<i64>> {
    let Some(range) = bam_raw::find_field(raw, tag)? else {
        return Ok(None);
    };

    let value = &raw[range.start + 3..range.end];

    Ok(Some(match raw[range.start + 2] {
        b'c' => i64::from(value[0].cast_signed()),
        b'C' => i64::from(value[0]),
        b's' => i64::from(i16::from_le_bytes([value[0], value[1]])),
        b'S' => i64::from(le_u16(value)),
        b'i' => i64::from(le_u32(value).cast_signed()),
        b'I' => i64::from(le_u32(value)),
        _ => 0,
    }))
}

/// The value bytes of a field, as samtools passes `bam_aux_get(...) + 1` to string functions:
/// up to the first NUL.
fn string_value(raw: &[u8], tag: [u8; 2]) -> anyhow::Result<Option<&[u8]>> {
    let Some(range) = bam_raw::find_field(raw, tag)? else {
        return Ok(None);
    };

    let value = &raw[range.start + 3..range.end];
    let end = memchr::memchr(0, value).unwrap_or(value.len());
    Ok(Some(&value[..end]))
}

pub fn run(mut args: Args) -> anyhow::Result<()> {
    let index = if args.customized_index {
        if args.regions.is_empty() {
            bail!("No index file provided");
        }

        Some(PathBuf::from(args.regions.remove(0)))
    } else {
        None
    };

    let read_options = ReadOptions {
        threads: args.threads,
        reference: args.cram_reference.clone(),
    };

    let (mut reader, path) = AlignmentReader::open(&args.input, &read_options)?;
    let header = reader.read_header()?;

    let (cov_min, cov_max, cov_step) = parse_coverage(&args.coverage)?;

    if args.coverage_threshold < 0 || args.coverage_threshold == i32::MAX {
        bail!(
            "Unsupported value for coverage threshold {}",
            args.coverage_threshold
        );
    }

    let flag_filter = args
        .filtering_flags
        .iter()
        .fold(0, |f, flags| f | flags.bits())
        | if args.remove_dups {
            Flags::DUPLICATE.bits()
        } else {
            0
        };

    #[expect(
        clippy::cast_possible_truncation,
        reason = "samtools stores atof(--GC-depth) in an int"
    )]
    let gcd_bin_size = args.gc_depth_bin_size as i32;

    let reference = args.reference.as_deref().map(Reference::open).transpose()?;

    let mut info = Info {
        flag_require: args.required_flags.map_or(0, |f| f.bits()),
        flag_filter,
        gcd_bin_size,
        nisize: args.insert_size,
        trim_qual: args.trim_quality,
        isize_main_bulk: args.most_inserts,
        cov_min,
        cov_max,
        cov_step,
        filter_readlen: args.read_length,
        remove_overlaps: args.remove_overlaps,
        cov_threshold: args.coverage_threshold,
        ref_chunk_size: args.ref_stats_chunk.max(1) * 1024 * 1024,
        reference,
        reference_names: header
            .reference_sequences()
            .keys()
            .map(|name| name.to_vec())
            .collect(),
        reference_lengths: header
            .reference_sequences()
            .values()
            .map(|rs| rs.length().get() as u64)
            .collect(),
        workers: args.threads,
    };

    let mut all = Stats::new(
        &mut info,
        &header,
        args.group_id.as_deref(),
        args.targets.as_deref(),
    )?;
    let mut splits: Vec<(Vec<u8>, Stats)> = Vec::new();
    let mut split_index: HashMap<Vec<u8>, usize> = HashMap::new();
    let mut pairs: HashMap<Vec<u8>, Pair> = HashMap::new();

    let split_tag = args
        .split_tag
        .as_deref()
        .map(|tag| -> anyhow::Result<[u8; 2]> {
            tag.as_bytes()
                .try_into()
                .with_context(|| format!("invalid tag: {tag}"))
        })
        .transpose()?;

    let mut encoder = noodles_bam::io::Writer::from(Vec::new());

    let mut process = |rec: Rec<'_>,
                       all: &mut Stats,
                       info: &Info,
                       pairs: &mut HashMap<Vec<u8>, Pair>|
     -> anyhow::Result<()> {
        use noodles_sam::alignment::io::Write as _;

        let raw: &[u8] = match &rec {
            Rec::Bam(r) => r.as_ref(),
            _ => {
                encoder.get_mut().clear();
                encoder.write_alignment_record(&header, rec.as_dyn())?;
                &encoder.get_ref()[size_of::<u32>()..]
            }
        };

        let bam = Bam::parse(raw)?;

        if let Some(tag) = split_tag {
            let Some(value) = (match bam_raw::find_field(raw, tag)? {
                Some(range) if raw[range.start + 2] == b'Z' => {
                    Some(raw[range.start + 3..range.end - 1].to_vec())
                }
                _ => None,
            }) else {
                bail!(
                    "Tag '{}' not found in bam_line.",
                    String::from_utf8_lossy(&tag)
                );
            };

            let i = match split_index.get(&value) {
                Some(&i) => i,
                None => {
                    let mut stats = Stats::new_with_info(info, args.targets.as_deref())?;
                    stats.split_name = Some(value.clone());
                    splits.push((value.clone(), stats));
                    split_index.insert(value, splits.len() - 1);
                    splits.len() - 1
                }
            };

            splits[i].1.collect(&bam, info, pairs)?;
        }

        all.collect(&bam, info, pairs)
    };

    if !args.regions.is_empty() {
        let Some(path) = path else {
            bail!("Random alignment retrieval only works for indexed files");
        };

        let mut set = RegionSet::new(&header);

        for region in &args.regions {
            set.add_region(&header, region)?;
        }

        set.merge();

        if args.targets.is_none() {
            all.replicate_regions(&set, &info, args.regions.len());
        }

        if all.regions.as_ref().is_some_and(|r| !r.is_empty()) {
            reader.for_each_in_queries(
                &path,
                index.as_deref(),
                &header,
                &Queries::Set(&set),
                |rec| process(rec, &mut all, &info, &mut pairs),
            )?;
        }
    } else {
        if info.cov_threshold > 0 && args.targets.is_none() {
            bail!("Coverage percentage calculation requires a list of target regions");
        }

        reader.for_each(&header, |rec| process(rec, &mut all, &info, &mut pairs))?;
    }

    if args.ref_stats {
        all.collect_reference_stats(&info)?;
    }

    all.round_buffer_flush(-1, &info)?;

    let command_line = std::iter::once("stats".to_owned())
        .chain(std::env::args().skip(2))
        .collect::<Vec<_>>()
        .join(" ");

    let mut out = BufWriter::new(output::raw_stdout()?);
    all.output(&mut out, &info, &command_line, None, args.sparse)?;
    out.flush()?;

    for (name, mut stats) in splits {
        stats.round_buffer_flush(-1, &info)?;

        let prefix = args.split_prefix.as_deref().unwrap_or(&args.input);
        let path = format!("{prefix}_{}.bamstat", String::from_utf8_lossy(&name));
        let file =
            File::create(&path).with_context(|| format!("Could not open '{path}' for writing."))?;
        let mut out = BufWriter::new(file);
        stats.output(
            &mut out,
            &info,
            &command_line,
            args.split_tag.as_deref(),
            args.sparse,
        )?;
        out.flush()?;
    }

    Ok(())
}

/// Parses `-c min,max,step` as `sscanf("%d,%d,%d")`.
fn parse_coverage(s: &str) -> anyhow::Result<(i32, i32, i32)> {
    let parts: Vec<&str> = s.split(',').collect();

    let parse = |i: usize| -> Option<i32> {
        parts.get(i).and_then(|p| {
            let p = p.trim_start();
            let len = p
                .char_indices()
                .take_while(|&(j, c)| c.is_ascii_digit() || (j == 0 && (c == '-' || c == '+')))
                .count();
            p[..len].parse().ok()
        })
    };

    match (parse(0), parse(1), parse(2)) {
        (Some(min), Some(max), Some(step)) => Ok((min, max, step)),
        _ => bail!("Unable to parse -c {s}"),
    }
}

impl Stats {
    /// `stats_init` and `init_stat_structs`.
    fn new(
        info: &mut Info,
        header: &sam::Header,
        group_id: Option<&str>,
        targets: Option<&Path>,
    ) -> anyhow::Result<Self> {
        // Coverage bins
        if info.cov_step > info.cov_max - info.cov_min + 1 {
            info.cov_step = info.cov_max - info.cov_min;

            if info.cov_step <= 0 {
                info.cov_step = 1;
            }
        }

        let ncov = usize::try_from(3 + (info.cov_max - info.cov_min) / info.cov_step)
            .context("invalid coverage range")?;
        info.cov_max =
            info.cov_min + ((info.cov_max - info.cov_min) / info.cov_step + 1) * info.cov_step - 1;

        let mut stats = Self::new_with_info(info, targets)?;
        stats.ncov = ncov;
        stats.cov = vec![0; ncov];

        if let Some(id) = group_id {
            stats.read_groups = Some(read_groups_for(header, id));
        }

        Ok(stats)
    }

    /// `stats_init` and `init_stat_structs` with the coverage options already adjusted.
    fn new_with_info(info: &Info, targets: Option<&Path>) -> anyhow::Result<Self> {
        let nbases = INITIAL_BASES;
        let ncov = usize::try_from(3 + (info.cov_max - info.cov_min) / info.cov_step)
            .context("invalid coverage range")?;

        let mut stats = Self {
            nbases,
            ngc: GC_BINS,
            nindels: nbases,
            quals_1st: vec![0; 2 * nbases],
            quals_2nd: vec![0; 2 * nbases],
            quals_rows: 2,
            gc_1st: vec![0; GC_BINS],
            gc_2nd: vec![0; GC_BINS],
            acgtno_cycles_1st: vec![CodeCounts::default(); nbases],
            acgtno_cycles_2nd: vec![CodeCounts::default(); nbases],
            acgtno_revcomp: vec![CodeCounts::default(); nbases],
            read_lengths: vec![0; nbases],
            read_lengths_1st: vec![0; nbases],
            read_lengths_2nd: vec![0; nbases],
            insertions: vec![0; nbases],
            deletions: vec![0; nbases],
            ins_cycles_1st: vec![0; nbases + 1],
            ins_cycles_2nd: vec![0; nbases + 1],
            del_cycles_1st: vec![0; nbases + 1],
            del_cycles_2nd: vec![0; nbases + 1],
            isize: InsertSizes::new(if info.nisize != 0 { info.nisize + 1 } else { 0 }),
            mapping_qualities: [0; 256],
            max_len: 0,
            max_len_1st: 0,
            max_len_2nd: 0,
            max_qual: 0,
            is_sorted: true,
            total_len: 0,
            total_len_1st: 0,
            total_len_2nd: 0,
            total_len_dup: 0,
            nreads_1st: 0,
            nreads_2nd: 0,
            nreads_other: 0,
            nreads_filtered: 0,
            nreads_dup: 0,
            nreads_unmapped: 0,
            nreads_single_mapped: 0,
            nreads_paired_and_mapped: 0,
            nreads_properly_paired: 0,
            nreads_paired_tech: 0,
            nreads_anomalous: 0,
            nreads_mq0: 0,
            nbases_mapped: 0,
            nbases_mapped_cigar: 0,
            nbases_trimmed: 0,
            nmismatches: 0,
            nreads_qc_failed: 0,
            nreads_secondary: 0,
            nreads_supplementary: 0,
            checksum_names: 0,
            checksum_reads: 0,
            checksum_quals: 0,
            igcd: 0,
            gcd: Vec::new(),
            tid: -1,
            gcd_pos: -1,
            pos: 0,
            ncov,
            cov: vec![0; ncov],
            cov_rbuf: RoundBuffer {
                pos: 0,
                size: i32::try_from(nbases * 5).expect("small"),
                start: 0,
                buffer: vec![0; nbases * 5],
            },
            rseq_buf: Vec::new(),
            mrseq_buf: 0,
            rseq_pos: -1,
            nrseq_buf: 0,
            mpc_buf: info
                .reference
                .is_some()
                .then(|| vec![0; QUALITY_BINS * nbases]),
            regions: None,
            reg_from: 0,
            reg_to: 0,
            chunks: Vec::new(),
            nchunks: 0,
            target_count: 0,
            sum_qual: 0,
            pending: PendingCycles::default(),
            workers: info.workers,
            read_groups: None,
            split_name: None,
            pair_count: 0,
            last_pair_tid: -2,
            last_read_flush: 0,
            reference_stats: None,
        };

        stats.realloc_rseq_buffer(info);

        if let Some(path) = targets {
            stats.init_regions(path, info)?;
        }

        Ok(stats)
    }

    /// `init_regions`: target regions from a file of `chr from to` lines (1-based, closed).
    fn init_regions(&mut self, path: &Path, info: &Info) -> anyhow::Result<()> {
        let text = fs::read(path).with_context(|| format!("{}", path.display()))?;
        let mut regions: Vec<Regions> = Vec::new();
        let mut warned = false;
        let mut prev_tid: Option<usize> = None;
        let mut prev_pos = -1;
        let mut nchunks = 0;

        // Lines as `kfgetline` reads them: no empty line after a final newline, and no `\r`.
        let text = text.strip_suffix(b"\n").unwrap_or(&text);

        for line in text.split(|&b| b == b'\n') {
            let line = line.strip_suffix(b"\r").unwrap_or(line);

            if line.first() == Some(&b'#') {
                continue;
            }

            let name_end = line
                .iter()
                .position(|&b| b.is_ascii_whitespace())
                .with_context(|| {
                    format!(
                        "Could not parse the file: {} [{}]",
                        path.display(),
                        String::from_utf8_lossy(line)
                    )
                })?;

            let name = &line[..name_end];

            let Some(tid) = info.reference_names.iter().position(|n| n == name) else {
                if !warned {
                    eprintln!(
                        "Warning: Some sequences not present in the BAM, e.g. \"{}\". This message is printed only once.",
                        String::from_utf8_lossy(name)
                    );
                }

                warned = true;
                continue;
            };

            let rest = String::from_utf8_lossy(&line[name_end + 1..]);
            let mut fields = rest.split_ascii_whitespace();
            let (Some(Ok(beg)), Some(Ok(end))) = (
                fields.next().map(str::parse::<i64>),
                fields.next().map(str::parse::<i64>),
            ) else {
                bail!("Could not parse the region [{rest}]");
            };

            if regions.len() <= tid {
                regions.resize_with(tid + 100, Regions::default);
            }

            if prev_tid != Some(tid) {
                prev_tid = Some(tid);
                prev_pos = beg;
            }

            if prev_pos > beg {
                bail!(
                    "The positions are not in chromosomal order ({}:{beg} comes after {prev_pos})",
                    String::from_utf8_lossy(name)
                );
            }

            regions[tid].intervals.push((beg, end));
            nchunks = nchunks.max(regions[tid].intervals.len());
        }

        if regions.is_empty() {
            bail!("Unable to map the -t sequences to the BAM sequences.");
        }

        // Sort the intervals and merge overlapping ones (adjacent ones stay apart).
        for (tid, reg) in regions.iter_mut().enumerate() {
            if reg.intervals.len() > 1 {
                reg.intervals
                    .sort_by(|a, b| a.0.cmp(&b.0).then(a.1.cmp(&b.1)));

                let mut merged: Vec<(i64, i64)> = Vec::with_capacity(reg.intervals.len());

                for &(beg, end) in &reg.intervals {
                    match merged.last_mut() {
                        Some(last) if last.1 >= beg => last.1 = last.1.max(end),
                        _ => merged.push((beg, end)),
                    }
                }

                reg.intervals = merged;
            }

            for &(beg, end) in &reg.intervals {
                if end < i64::MAX {
                    self.target_count = self.target_count.wrapping_add_signed(end - beg + 1);
                } else if let Some(&len) = info.reference_lengths.get(tid)
                    && len > 0
                {
                    self.target_count = self
                        .target_count
                        .wrapping_add(len.wrapping_sub(beg as u64).wrapping_add(1));
                }
            }
        }

        self.regions = Some(regions);
        self.chunks = vec![(0, 0); nchunks];
        self.nchunks = 0;

        Ok(())
    }

    /// `replicate_regions`: target regions from the command-line regions.
    fn replicate_regions(&mut self, set: &RegionSet, info: &Info, region_count: usize) {
        let mut regions = Vec::new();

        for tid in 0..info.reference_names.len() {
            let intervals = set.intervals(tid);

            if intervals.is_empty() {
                continue;
            }

            if regions.len() <= tid {
                regions.resize_with(tid + 10, Regions::default);
            }

            // htslib ends a region without an end (or ending at the sequence end here) at
            // `HTS_POS_MAX`.
            let length = info.reference_lengths[tid];

            regions[tid].intervals = intervals
                .iter()
                .map(|&(beg, end)| {
                    let end = if end as u64 == length {
                        i64::MAX
                    } else {
                        end as i64
                    };
                    (beg as i64, end)
                })
                .collect();

            for &(beg, end) in &regions[tid].intervals {
                let end = if end == i64::MAX { length as i64 } else { end };
                self.target_count = self.target_count.wrapping_add_signed(end - beg + 1);
            }
        }

        self.regions = Some(regions);
        self.chunks = vec![(0, 0); region_count.max(1)];
    }

    /// `realloc_rseq_buffer`
    fn realloc_rseq_buffer(&mut self, info: &Info) {
        let n = i64::try_from(self.nbases * 10)
            .unwrap_or(i64::MAX)
            .max(i64::from(info.gcd_bin_size));

        if self.mrseq_buf < n {
            self.rseq_buf.resize(usize::try_from(n).unwrap_or(0), 0);
            self.mrseq_buf = n;
        }
    }

    /// The per-cycle tables for every cycle.
    fn cycle_tables(&mut self) -> CycleTables<'_> {
        CycleTables {
            lo: 0,
            first: &mut self.acgtno_cycles_1st,
            last: &mut self.acgtno_cycles_2nd,
            revcomp: &mut self.acgtno_revcomp,
            quals_first: &mut self.quals_1st,
            quals_last: &mut self.quals_2nd,
            rows: self.quals_rows,
        }
    }

    /// Counts the pending reads on `workers` threads, splitting the cycles among them.
    fn flush_cycles(&mut self) {
        let pending = std::mem::take(&mut self.pending);

        if pending.reads.is_empty() {
            self.pending = pending;
            return;
        }

        let workers = self.workers.max(1);
        let rows = self.quals_rows;
        let bounds = balanced_bounds(&pending.reads, pending.max_len, workers);

        let mut first = &mut self.acgtno_cycles_1st[..];
        let mut last = &mut self.acgtno_cycles_2nd[..];
        let mut revcomp = &mut self.acgtno_revcomp[..];
        let mut quals_first = &mut self.quals_1st[..];
        let mut quals_last = &mut self.quals_2nd[..];

        let mut parts = Vec::with_capacity(workers);
        let mut lo = 0;

        for &hi in &bounds {
            let n = hi - lo;

            let (a, rest) = first.split_at_mut(n);
            first = rest;
            let (b, rest) = last.split_at_mut(n);
            last = rest;
            let (c, rest) = revcomp.split_at_mut(n);
            revcomp = rest;
            let (d, rest) = quals_first.split_at_mut(n * rows);
            quals_first = rest;
            let (e, rest) = quals_last.split_at_mut(n * rows);
            quals_last = rest;

            parts.push(CycleTables {
                lo,
                first: a,
                last: b,
                revcomp: c,
                quals_first: d,
                quals_last: e,
                rows,
            });

            lo += n;
        }

        let reads = &pending.reads;
        let data = &pending.data;

        std::thread::scope(|scope| {
            for mut tables in parts {
                scope.spawn(move || {
                    for read in reads {
                        let seq_end = read.start + read.len.div_ceil(2);
                        let seq = &data[read.start..seq_end];
                        let qual = &data[seq_end..seq_end + read.len];
                        count_cycles(seq, qual, read.reverse, read.first, &mut tables);
                    }
                });
            }
        });

        let mut pending = pending;
        pending.data.clear();
        pending.reads.clear();
        pending.max_len = 0;
        self.pending = pending;
    }

    /// `realloc_qual_buffers`: grows the quality tables to `new_len` cycles and `new_max_qual + 2`
    /// qualities.
    fn realloc_qual_buffers(&mut self, old_len: usize, new_len: usize, new_max_qual: usize) {
        let old_rows = self.quals_rows;
        let rows = (new_max_qual + 2).max(old_rows);

        for quals in [&mut self.quals_1st, &mut self.quals_2nd] {
            if rows == old_rows {
                quals.resize(rows * new_len, 0);
                continue;
            }

            let mut grown = vec![0; rows * new_len];

            for (cycle, counts) in quals.chunks_exact(old_rows).take(old_len).enumerate() {
                grown[cycle * rows..cycle * rows + old_rows].copy_from_slice(counts);
            }

            *quals = grown;
        }

        self.quals_rows = rows;
    }

    /// `realloc_buffers`: grows the per-cycle buffers for a read of `seq_len` bases.
    fn realloc_buffers(&mut self, seq_len: usize, info: &Info) {
        let n = 2 * (1 + seq_len - self.nbases) + self.nbases;

        self.realloc_qual_buffers(self.nbases, n, self.max_qual);

        if let Some(mpc) = &mut self.mpc_buf {
            mpc.resize(n * QUALITY_BINS, 0);
        }

        for v in [
            &mut self.acgtno_cycles_1st,
            &mut self.acgtno_cycles_2nd,
            &mut self.acgtno_revcomp,
        ] {
            v.resize(n, CodeCounts::default());
        }

        for v in [
            &mut self.read_lengths,
            &mut self.read_lengths_1st,
            &mut self.read_lengths_2nd,
            &mut self.insertions,
            &mut self.deletions,
        ] {
            v.resize(n, 0);
        }

        for v in [
            &mut self.ins_cycles_1st,
            &mut self.ins_cycles_2nd,
            &mut self.del_cycles_1st,
            &mut self.del_cycles_2nd,
        ] {
            v.resize(n + 1, 0);
        }

        self.nbases = n;

        // As in samtools, the old contents are copied by byte counts that are meant to be
        // element counts, so only some of the pending depths survive.
        let new_size = seq_len * 5;
        let old = &self.cov_rbuf;
        let old_bytes: Vec<u8> = old.buffer.iter().flat_map(|d| d.to_ne_bytes()).collect();
        let mut new_bytes = vec![0u8; new_size * 4];
        let start = usize::try_from(old.start).unwrap_or(0);
        let tail = usize::try_from(old.size - old.start).unwrap_or(0);

        new_bytes[..tail].copy_from_slice(&old_bytes[start * 4..start * 4 + tail]);

        if old.start > 1 {
            new_bytes[tail * 4..tail * 4 + start].copy_from_slice(&old_bytes[..start]);
        }

        self.cov_rbuf.buffer = new_bytes
            .as_chunks::<4>()
            .0
            .iter()
            .map(|b| i32::from_ne_bytes(*b))
            .collect();
        self.cov_rbuf.start = 0;
        self.cov_rbuf.size = i32::try_from(new_size).unwrap_or(i32::MAX);

        self.realloc_rseq_buffer(info);
    }

    /// `coverage_idx`
    fn coverage_index(&self, info: &Info, depth: i32) -> usize {
        if depth < info.cov_min {
            0
        } else if depth > info.cov_max {
            self.ncov - 1
        } else {
            usize::try_from(1 + (depth - info.cov_min) / info.cov_step).expect("nonnegative")
        }
    }

    /// `round_buffer_flush`
    fn round_buffer_flush(&mut self, pos: i64, info: &Info) -> anyhow::Result<()> {
        if pos == self.cov_rbuf.pos {
            return Ok(());
        }

        let new_pos = pos;
        let size = i64::from(self.cov_rbuf.size);
        let mut pos = pos;

        if pos == -1 || pos - self.cov_rbuf.pos >= size {
            pos = self.cov_rbuf.pos + size - 1;
        }

        if pos < self.cov_rbuf.pos {
            bail!(
                "Expected coordinates in ascending order, got {pos} after {}",
                self.cov_rbuf.pos
            );
        }

        let mut ifrom = self.cov_rbuf.start;
        let ito = lidx2ridx(
            self.cov_rbuf.start,
            self.cov_rbuf.size,
            self.cov_rbuf.pos,
            pos - 1,
        );

        if ifrom > ito {
            for ibuf in ifrom..self.cov_rbuf.size {
                self.flush_depth(ibuf as usize, info);
            }

            ifrom = 0;
        }

        for ibuf in ifrom..=ito {
            self.flush_depth(ibuf as usize, info);
        }

        self.cov_rbuf.start = if new_pos == -1 {
            0
        } else {
            lidx2ridx(
                self.cov_rbuf.start,
                self.cov_rbuf.size,
                self.cov_rbuf.pos,
                pos,
            )
        };
        self.cov_rbuf.pos = new_pos;

        Ok(())
    }

    fn flush_depth(&mut self, ibuf: usize, info: &Info) {
        let depth = self.cov_rbuf.buffer[ibuf];

        if depth != 0 {
            let i = self.coverage_index(info, depth);
            self.cov[i] += 1;
            self.cov_rbuf.buffer[ibuf] = 0;
        }
    }

    /// `round_buffer_insert_read`: adds the 0-based half-open `[from, to)`.
    fn insert_read(&mut self, from: i64, to: i64) -> anyhow::Result<()> {
        let rbuf = &mut self.cov_rbuf;

        if to - from > i64::from(rbuf.size) {
            bail!(
                "The read length too big ({}), please increase the buffer length (currently {})",
                to - from,
                rbuf.size
            );
        }

        if from < rbuf.pos {
            bail!(
                "The reads are not sorted ({from} comes after {}).",
                rbuf.pos
            );
        }

        let mut ifrom = lidx2ridx(rbuf.start, rbuf.size, rbuf.pos, from);
        let ito = lidx2ridx(rbuf.start, rbuf.size, rbuf.pos, to);

        if ifrom > ito {
            for d in &mut rbuf.buffer[ifrom as usize..rbuf.size as usize] {
                *d += 1;
            }

            ifrom = 0;
        }

        for d in &mut rbuf.buffer[ifrom as usize..ito as usize] {
            *d += 1;
        }

        Ok(())
    }

    /// `read_ref_seq`
    fn read_ref_seq(&mut self, info: &Info, tid: i32, pos: i64, end: i64) -> anyhow::Result<()> {
        let reference = info.reference.as_ref().expect("reference is loaded");

        if end >= pos + self.mrseq_buf - 1 && self.mrseq_buf < end - pos {
            self.mrseq_buf = end - pos;
            self.rseq_buf.resize(usize::try_from(self.mrseq_buf)?, 0);
        }

        let name = &info.reference_names[usize::try_from(tid)?];
        let mut sequence = Vec::new();

        if !reference.fetch(name, pos, pos + self.mrseq_buf - 1, &mut sequence)? {
            bail!(
                "Failed to fetch the sequence \"{}\"",
                String::from_utf8_lossy(name)
            );
        }

        for (dst, &b) in self.rseq_buf.iter_mut().zip(&sequence) {
            *dst = match b {
                b'A' | b'a' => 1,
                b'C' | b'c' => 2,
                b'G' | b'g' => 4,
                b'T' | b't' => 8,
                _ => 0,
            };
        }

        self.rseq_buf[sequence.len()..].fill(0);
        self.nrseq_buf = i64::try_from(sequence.len())?;
        self.rseq_pos = pos;
        self.tid = tid;

        Ok(())
    }

    /// `fai_gc_content`
    fn gc_content(&self, pos: i64, len: i32) -> f32 {
        let i = usize::try_from(pos - self.rseq_pos).expect("in the buffer");
        let end = (pos - self.rseq_pos + i64::from(len)).min(self.nrseq_buf);
        let end = usize::try_from(end).unwrap_or(0);

        let (mut gc, mut count) = (0u32, 0u32);

        for &c in self.rseq_buf.get(i..end).unwrap_or(&[]) {
            match c {
                2 | 4 => {
                    gc += 1;
                    count += 1;
                }
                1 | 8 => count += 1,
                _ => {}
            }
        }

        if count > 0 {
            gc as f32 / count as f32
        } else {
            0.0
        }
    }

    /// `count_indels`
    fn count_indels(&mut self, bam: &Bam<'_>) -> anyhow::Result<()> {
        let is_fwd = !bam.has(Flags::REVERSE_COMPLEMENTED);
        let order = bam.order();
        let read_len = i64::try_from(bam.l_seq)?;
        let mut icycle: i64 = 0;

        for &op in &bam.cigar {
            let (cig, ncig) = (op & 0xf, i64::from(op >> 4));

            if ncig == 0 {
                continue;
            }

            if cig == CIGAR_INSERTION {
                let idx = if is_fwd {
                    icycle
                } else {
                    read_len - icycle - ncig
                };

                if idx < 0 {
                    bail!("FIXME: read_len={read_len} vs icycle={icycle}");
                }

                let idx = usize::try_from(idx)?;

                if idx >= self.nbases {
                    bail!("FIXME: {idx} vs {}", self.nbases);
                }

                if order == READ_ORDER_FIRST {
                    self.ins_cycles_1st[idx] += 1;
                }

                if order == READ_ORDER_LAST {
                    self.ins_cycles_2nd[idx] += 1;
                }

                icycle += ncig;

                if ncig as usize <= self.nindels {
                    self.insertions[ncig as usize - 1] += 1;
                }

                continue;
            }

            if cig == CIGAR_DELETION {
                let idx = if is_fwd {
                    icycle - 1
                } else {
                    read_len - icycle - 1
                };

                if idx < 0 {
                    continue;
                }

                let idx = usize::try_from(idx)?;

                if idx >= self.nbases {
                    bail!("FIXME: {idx} vs {}", self.nbases);
                }

                if order == READ_ORDER_FIRST {
                    self.del_cycles_1st[idx] += 1;
                }

                if order == READ_ORDER_LAST {
                    self.del_cycles_2nd[idx] += 1;
                }

                if ncig as usize <= self.nindels {
                    self.deletions[ncig as usize - 1] += 1;
                }

                continue;
            }

            if cig != CIGAR_SKIP && cig != CIGAR_HARD_CLIP && cig != CIGAR_PAD {
                icycle += ncig;
            }
        }

        Ok(())
    }

    /// `count_mismatches_per_cycle`
    fn count_mismatches_per_cycle(&mut self, bam: &Bam<'_>, read_len: usize) -> anyhow::Result<()> {
        let is_fwd = !bam.has(Flags::REVERSE_COMPLEMENTED);
        let mut iread = 0usize;
        let mut icycle = 0usize;
        let mut iref = bam.pos - self.rseq_pos;
        let max_len = self.max_len;
        let nbases = self.nbases;
        let mpc = self.mpc_buf.as_mut().expect("reference is loaded");

        for &op in &bam.cigar {
            let (cig, ncig) = (op & 0xf, (op >> 4) as usize);

            match cig {
                CIGAR_INSERTION => {
                    iread += ncig;
                    icycle += ncig;
                    continue;
                }
                CIGAR_DELETION => {
                    iref += ncig as i64;
                    continue;
                }
                CIGAR_SOFT_CLIP => {
                    icycle += ncig;
                    iread += ncig;
                    continue;
                }
                CIGAR_HARD_CLIP => {
                    icycle += ncig;
                    continue;
                }
                CIGAR_SKIP | CIGAR_PAD => continue,
                CIGAR_MATCH | CIGAR_SEQUENCE_MATCH | CIGAR_SEQUENCE_MISMATCH => {}
                _ => bail!("TODO: cigar {cig}"),
            }

            if ncig as i64 + iref > self.nrseq_buf {
                bail!(
                    "FIXME: {ncig}+{iref} > {}, {}",
                    self.nrseq_buf,
                    String::from_utf8_lossy(bam.name)
                );
            }

            for _ in 0..ncig {
                let cread = bam.base(iread);
                let cref = self.rseq_buf[usize::try_from(iref)?];

                let cycle_index = |icycle: usize| {
                    if is_fwd {
                        icycle
                    } else {
                        read_len - icycle - 1
                    }
                };

                if cread == 15 {
                    let idx = cycle_index(icycle);

                    if idx > max_len {
                        bail!("mpc: {idx}>{max_len}");
                    }

                    let idx = idx * QUALITY_BINS;

                    if idx >= QUALITY_BINS * nbases {
                        bail!("FIXME: mpc_buf overflow");
                    }

                    mpc[idx] += 1;
                } else if cref != 0 && cread != 0 && cref != cread {
                    let qual = usize::from(bam.qual[iread]) + 1;
                    let idx = cycle_index(icycle);

                    if idx > max_len {
                        bail!("mpc: {idx}>{max_len}");
                    }

                    let idx = idx * QUALITY_BINS + qual;

                    if idx >= QUALITY_BINS * nbases {
                        bail!("FIXME: mpc_buf overflow");
                    }

                    mpc[idx] += 1;
                }

                iref += 1;
                iread += 1;
                icycle += 1;
            }
        }

        Ok(())
    }

    /// `update_checksum`
    fn update_checksum(&mut self, bam: &Bam<'_>) {
        let crc = |bytes: &[u8]| {
            let mut crc = flate2::Crc::new();
            crc.update(bytes);
            crc.sum()
        };

        self.checksum_names = self.checksum_names.wrapping_add(crc(bam.name));

        if bam.l_seq == 0 {
            return;
        }

        self.checksum_reads = self.checksum_reads.wrapping_add(crc(bam.seq));
        self.checksum_quals = self.checksum_quals.wrapping_add(crc(bam.qual));
    }

    /// `collect_orig_read_stats`; returns the read's GC count.
    fn collect_orig_read_stats(&mut self, bam: &Bam<'_>, info: &Info) -> usize {
        let seq_len = bam.l_seq;
        self.total_len += seq_len as u64;

        if bam.has(Flags::QC_FAIL) {
            self.nreads_qc_failed += 1;
        }

        if bam.has(Flags::SEGMENTED) {
            self.nreads_paired_tech += 1;
        }

        let order = bam.order();
        let reverse = bam.has(Flags::REVERSE_COMPLEMENTED);
        let mut gc_count = 0;

        if order == READ_ORDER_FIRST || order == READ_ORDER_LAST {
            gc_count = bam
                .seq
                .iter()
                .map(|&b| usize::from(GC_IN_BYTE[usize::from(b)]))
                .sum();
        }

        let gc_idx_min = gc_count * (self.ngc - 1) / seq_len;
        let gc_idx_max = ((gc_count + 1) * (self.ngc - 1) / seq_len).min(self.ngc - 1);

        let quals = match order {
            READ_ORDER_FIRST => {
                self.nreads_1st += 1;
                self.total_len_1st += seq_len as u64;

                for n in &mut self.gc_1st[gc_idx_min..gc_idx_max.max(gc_idx_min)] {
                    *n += 1;
                }

                Some(READ_ORDER_FIRST)
            }
            READ_ORDER_LAST => {
                self.nreads_2nd += 1;
                self.total_len_2nd += seq_len as u64;

                for n in &mut self.gc_2nd[gc_idx_min..gc_idx_max.max(gc_idx_min)] {
                    *n += 1;
                }

                Some(READ_ORDER_LAST)
            }
            _ => {
                self.nreads_other += 1;
                None
            }
        };

        if info.trim_qual > 0 {
            self.nbases_trimmed += bwa_trim_read(info.trim_qual, bam.qual, seq_len, reverse) as u64;
        }

        if let Some(order) = quals {
            let max_qual = usize::from(bam.qual.iter().copied().max().unwrap_or(0));

            if max_qual > self.max_qual {
                self.flush_cycles();
                self.realloc_qual_buffers(self.nbases, self.nbases, max_qual);
                self.max_qual = max_qual;
            }

            self.sum_qual += bam.qual.iter().map(|&q| u64::from(q)).sum::<u64>();

            let first = order == READ_ORDER_FIRST;

            if self.workers == 0 {
                let mut tables = self.cycle_tables();
                count_cycles(bam.seq, bam.qual, reverse, first, &mut tables);
            } else {
                let pending = &mut self.pending;
                let start = pending.data.len();
                pending.data.extend_from_slice(bam.seq);
                pending.data.extend_from_slice(bam.qual);
                pending.reads.push(PendingRead {
                    start,
                    len: seq_len,
                    reverse,
                    first,
                });
                pending.max_len = pending.max_len.max(seq_len);

                if pending.data.len() >= PENDING_BASES {
                    self.flush_cycles();
                }
            }
        }

        if bam.has(Flags::UNMAPPED) {
            self.nreads_unmapped += 1;
        } else {
            self.nbases_mapped += seq_len as u64;

            if bam.mapq == 0 {
                self.nreads_mq0 += 1;
            }

            if !bam.is_paired_and_mapped() {
                self.nreads_single_mapped += 1;
            } else {
                self.nreads_paired_and_mapped += 1;

                if bam.has(Flags::PROPERLY_SEGMENTED) {
                    self.nreads_properly_paired += 1;
                }

                if bam.tid != bam.mtid {
                    self.nreads_anomalous += 1;
                }
            }
        }

        gc_count
    }

    /// `is_in_regions`
    fn is_in_regions(&mut self, bam: &Bam<'_>) -> anyhow::Result<bool> {
        let Some(regions) = &mut self.regions else {
            return Ok(true);
        };

        let Ok(tid) = usize::try_from(bam.tid) else {
            return Ok(false);
        };

        if tid >= regions.len() {
            return Ok(false);
        }

        if !self.is_sorted {
            bail!("The BAM must be sorted in order for -t to work.");
        }

        let reg = &mut regions[tid];

        if reg.cpos == reg.intervals.len() {
            return Ok(false);
        }

        let mut i = reg.cpos;

        while i < reg.intervals.len() && reg.intervals[i].1 <= bam.pos {
            i += 1;
        }

        if i >= reg.intervals.len() {
            reg.cpos = reg.intervals.len();
            return Ok(false);
        }

        let end_pos = bam.end_pos();

        if end_pos < reg.intervals[i].0 {
            return Ok(false);
        }

        reg.cpos = i;
        self.reg_from = reg.intervals[i].0;
        self.reg_to = reg.intervals[i].1;

        self.nchunks = 0;

        for &(beg, end) in &reg.intervals[i..] {
            if bam.pos < end && end_pos >= beg {
                if self.nchunks >= self.chunks.len() {
                    self.chunks.push((0, 0));
                }

                self.chunks[self.nchunks] = ((bam.pos + 1).max(beg), end_pos.min(end));
                self.nchunks += 1;
            }
        }

        Ok(true)
    }

    /// `remove_overlaps`: adds `[pmin, pmax)` to the coverage, less what the mate already covered.
    fn remove_overlaps(
        &mut self,
        bam: &Bam<'_>,
        pairs: &mut HashMap<Vec<u8>, Pair>,
        pmin: i64,
        pmax: i64,
    ) -> anyhow::Result<()> {
        let order =
            u32::from(bam.has(Flags::FIRST_SEGMENT)) + 2 * u32::from(bam.has(Flags::LAST_SEGMENT));

        if !bam.has(Flags::SEGMENTED)
            || bam.has(Flags::MATE_UNMAPPED)
            || i64::from(bam.isize).abs() >= 2 * bam.l_seq as i64
            || (order != READ_ORDER_FIRST && order != READ_ORDER_LAST)
        {
            if pmin >= 0 {
                self.insert_read(pmin, pmax)?;
            }

            return Ok(());
        }

        let mut pmin = pmin;

        match pairs.get_mut(bam.name) {
            None => {
                if pmin == -1 {
                    return Ok(());
                }

                pairs.insert(
                    bam.name.to_vec(),
                    Pair {
                        first: order,
                        chunks: vec![(pmin, pmax)],
                    },
                );
                self.pair_count += 1;
            }
            Some(pair) if order == pair.first => {
                if pmin == -1 {
                    return Ok(());
                }

                pair.chunks.push((pmin, pmax));
            }
            Some(pair) => {
                if pmin == -1 {
                    pairs.remove(bam.name);
                    self.pair_count = self.pair_count.wrapping_sub(1);
                    return Ok(());
                }

                for &(beg, end) in &pair.chunks {
                    if pmin >= end {
                        continue;
                    }

                    if pmax <= beg {
                        break;
                    }

                    if pmin < beg {
                        self.insert_read(pmin, beg)?;
                        pmin = beg;
                    }

                    if pmax <= end {
                        self.nbases_mapped_cigar =
                            self.nbases_mapped_cigar.wrapping_add_signed(-(pmax - pmin));
                        return Ok(());
                    }

                    self.nbases_mapped_cigar =
                        self.nbases_mapped_cigar.wrapping_add_signed(-(end - pmin));
                    pmin = end;
                }
            }
        }

        self.insert_read(pmin, pmax)
    }

    /// `collect_stats`
    fn collect(
        &mut self,
        bam: &Bam<'_>,
        info: &Info,
        pairs: &mut HashMap<Vec<u8>, Pair>,
    ) -> anyhow::Result<()> {
        if !self.is_in_regions(bam)? {
            return Ok(());
        }

        if let Some(read_groups) = &self.read_groups {
            let Some(rg) = string_value(bam.raw, *b"RG")? else {
                return Ok(());
            };

            if !read_groups.iter().any(|id| id == rg) {
                return Ok(());
            }
        }

        if info.flag_require != 0 && bam.flag & info.flag_require != info.flag_require {
            self.nreads_filtered += 1;
            return Ok(());
        }

        if info.flag_filter != 0 && bam.flag & info.flag_filter != 0 {
            self.nreads_filtered += 1;
            return Ok(());
        }

        if info.filter_readlen != -1 && bam.l_seq as i64 != info.filter_readlen {
            return Ok(());
        }

        self.update_checksum(bam);

        if bam.has(Flags::SECONDARY) {
            self.nreads_secondary += 1;
            return Ok(());
        }

        if bam.has(Flags::SUPPLEMENTARY) {
            self.nreads_supplementary += 1;
        }

        let seq_len = bam.l_seq;

        if seq_len == 0 {
            return Ok(());
        }

        if bam.has(Flags::DUPLICATE) {
            self.total_len_dup += seq_len as u64;
            self.nreads_dup += 1;
        }

        let order = bam.order();

        let read_len = seq_len
            + bam
                .cigar
                .iter()
                .filter(|&&op| op & 0xf == CIGAR_HARD_CLIP)
                .map(|&op| (op >> 4) as usize)
                .sum::<usize>();

        if read_len >= self.nbases {
            self.flush_cycles();
            self.realloc_buffers(read_len, info);
        }

        self.max_len = self.max_len.max(read_len);

        if order == READ_ORDER_FIRST {
            self.max_len_1st = self.max_len_1st.max(read_len);
        }

        if order == READ_ORDER_LAST {
            self.max_len_2nd = self.max_len_2nd.max(read_len);
        }

        let excluded = Flags::UNMAPPED
            | Flags::SECONDARY
            | Flags::SUPPLEMENTARY
            | Flags::QC_FAIL
            | Flags::DUPLICATE;

        if bam.flag & excluded.bits() == 0 {
            self.mapping_qualities[usize::from(bam.mapq)] += 1;
        }

        let mut gc_count = 0;

        if bam.is_original() {
            self.read_lengths[read_len] += 1;

            if order == READ_ORDER_FIRST {
                self.read_lengths_1st[read_len] += 1;
            }

            if order == READ_ORDER_LAST {
                self.read_lengths_2nd[read_len] += 1;
            }

            gc_count = self.collect_orig_read_stats(bam, info);
        }

        if bam.has(Flags::UNMAPPED) {
            return Ok(());
        }

        self.count_indels(bam)?;

        if bam.is_paired_and_mapped() && bam.is_original() {
            let mut isize = bam.isize.wrapping_abs();

            if info.nisize > 0 && isize > info.nisize {
                isize = info.nisize;
            }

            if isize > 0 || bam.tid == bam.mtid {
                let pos_fst = bam.mpos - bam.pos;
                let is_fst: i64 = if bam.has(Flags::FIRST_SEGMENT) { 1 } else { -1 };
                let is_fwd: i64 = if bam.has(Flags::REVERSE_COMPLEMENTED) {
                    -1
                } else {
                    1
                };
                let is_mfwd: i64 = if bam.has(Flags::MATE_REVERSE_COMPLEMENTED) {
                    -1
                } else {
                    1
                };

                let orientation = if is_fwd * is_mfwd > 0 {
                    OTHER
                } else if is_fst * pos_fst > 0 {
                    if is_fst * is_fwd > 0 { INWARD } else { OUTWARD }
                } else if is_fst * pos_fst < 0 {
                    if is_fst * is_fwd > 0 { OUTWARD } else { INWARD }
                } else {
                    // Exactly overlapping reads are taken as inward.
                    INWARD
                };

                self.isize.increment(isize, orientation);
            }
        }

        if let Some(nm) = integer_field(bam.raw, *b"NM")? {
            self.nmismatches = self.nmismatches.wrapping_add_signed(nm);
        }

        if bam.cigar.is_empty() {
            bail!("FIXME: mapped read with no cigar?");
        }

        let mut readlen = seq_len as i64;

        if self.regions.is_some() {
            // Count only on-target bases
            let mut iref = bam.pos + 1;

            for &op in &bam.cigar {
                let (cig, ncig_init) = (op & 0xf, i64::from(op >> 4));
                let mut ncig = ncig_init;

                if ncig == 0 {
                    continue;
                }

                match cig {
                    CIGAR_DELETION => readlen += ncig,
                    CIGAR_MATCH | CIGAR_SEQUENCE_MATCH | CIGAR_SEQUENCE_MISMATCH => {
                        if iref < self.reg_from {
                            ncig -= self.reg_from - iref;
                        }

                        if iref + ncig_init - 1 > self.reg_to {
                            ncig -= iref + ncig_init - 1 - self.reg_to;
                        }

                        self.nbases_mapped_cigar =
                            self.nbases_mapped_cigar.wrapping_add_signed(ncig.max(0));
                        iref += ncig_init;
                    }
                    CIGAR_INSERTION => {
                        iref += ncig;

                        if iref >= self.reg_from && iref <= self.reg_to {
                            self.nbases_mapped_cigar += ncig as u64;
                        }
                    }
                    _ => {}
                }
            }
        } else {
            for &op in &bam.cigar {
                let (cig, len) = (op & 0xf, u64::from(op >> 4));

                if matches!(
                    cig,
                    CIGAR_MATCH | CIGAR_INSERTION | CIGAR_SEQUENCE_MATCH | CIGAR_SEQUENCE_MISMATCH
                ) {
                    self.nbases_mapped_cigar += len;
                }

                if cig == CIGAR_DELETION {
                    readlen += len as i64;
                }
            }
        }

        if self.tid == bam.tid && bam.pos < self.pos {
            self.is_sorted = false;
        }

        self.pos = bam.pos;

        if !self.is_sorted {
            return Ok(());
        }

        if self.tid == -1 || self.tid != bam.tid {
            self.round_buffer_flush(-1, info)?;
        }

        // Drop pairs that can no longer overlap, to bound memory.
        self.last_read_flush += 1;

        if self.pair_count > PAIR_MAX && self.last_read_flush > PAIR_MAX {
            let removed = cleanup_overlaps(pairs, bam.pos);
            self.pair_count = self.pair_count.wrapping_sub(removed);
            self.last_read_flush = 0;
        }

        if self.last_pair_tid != bam.tid {
            let removed = cleanup_overlaps(pairs, i64::MAX - 1);
            self.pair_count = self.pair_count.wrapping_sub(removed);
            self.last_pair_tid = bam.tid;
            self.last_read_flush = 0;
        }

        let gcd_bin_size = i64::from(info.gcd_bin_size);

        if info.reference.is_some() {
            let mut inc_ref = 0;
            let mut inc_gcd = false;

            if self.rseq_pos == -1
                || self.tid != bam.tid
                || self.rseq_pos + self.nrseq_buf < bam.pos + readlen
            {
                inc_ref = bam.pos + readlen;
                inc_gcd = true;
            } else if self.gcd_pos + gcd_bin_size < bam.pos + readlen {
                inc_gcd = true;

                if self.rseq_pos + self.nrseq_buf < bam.pos + gcd_bin_size {
                    inc_ref = bam.pos + gcd_bin_size;
                }
            }

            if inc_gcd {
                self.igcd += 1;
                self.ensure_gcd(info);

                if inc_ref != 0 {
                    self.read_ref_seq(info, bam.tid, bam.pos, inc_ref)?;
                }

                self.gcd_pos = bam.pos;
                self.gcd[self.igcd].gc = self.gc_content(self.gcd_pos, info.gcd_bin_size);
            }

            self.count_mismatches_per_cycle(bam, read_len)?;
        } else if self.gcd_pos == -1 || self.tid != bam.tid || bam.pos - self.gcd_pos > gcd_bin_size
        {
            self.tid = bam.tid;
            self.gcd_pos = bam.pos;
            self.igcd += 1;
            self.ensure_gcd(info);
        }

        self.ensure_gcd(info);
        self.gcd[self.igcd].depth = self.gcd[self.igcd].depth.wrapping_add(1);

        if info.reference.is_none() {
            self.gcd[self.igcd].gc += gc_count as f32 / seq_len as f32;
        }

        // Coverage distribution
        self.round_buffer_flush(bam.pos, info)?;

        if self.regions.is_some() {
            let mut p = bam.pos;
            let (mut j, mut i) = (0, 0);

            while j < bam.cigar.len() && i < self.nchunks {
                let (op, oplen) = (bam.cigar[j] & 0xf, i64::from(bam.cigar[j] >> 4));

                if matches!(
                    op,
                    CIGAR_MATCH | CIGAR_SEQUENCE_MATCH | CIGAR_SEQUENCE_MISMATCH
                ) {
                    let pmin = p.max(self.chunks[i].0 - 1);
                    let pmax = (p + oplen).min(self.chunks[i].1);

                    if pmax > pmin {
                        if info.remove_overlaps {
                            self.remove_overlaps(bam, pairs, pmin, pmax)?;
                        } else {
                            self.insert_read(pmin, pmax)?;
                        }
                    }
                }

                let pnew = p + if consumes_reference(op) { oplen } else { 0 };

                if pnew >= self.chunks[i].1 {
                    i += 1;
                } else {
                    j += 1;
                    p = pnew;
                }
            }
        } else {
            let mut p = bam.pos;

            for &op in &bam.cigar {
                let (op, oplen) = (op & 0xf, i64::from(op >> 4));

                if matches!(
                    op,
                    CIGAR_MATCH | CIGAR_SEQUENCE_MATCH | CIGAR_SEQUENCE_MISMATCH
                ) {
                    if info.remove_overlaps {
                        self.remove_overlaps(bam, pairs, p, p + oplen)?;
                    } else {
                        self.insert_read(p, p + oplen)?;
                    }
                }

                if consumes_reference(op) {
                    p += oplen;
                }
            }
        }

        if info.remove_overlaps {
            self.remove_overlaps(bam, pairs, -1, -1)?;
        }

        Ok(())
    }

    /// `realloc_gcd_buffer`: makes `gcd[igcd]` valid, growing as `hts_expand0` does.
    fn ensure_gcd(&mut self, info: &Info) {
        if self.igcd >= self.gcd.len() {
            let n = (self.igcd + 1).next_power_of_two();
            self.gcd.resize(n, GcDepth::default());
            self.realloc_rseq_buffer(info);
        }
    }

    /// `collect_refstats`
    fn collect_reference_stats(&mut self, info: &Info) -> anyhow::Result<()> {
        let total_count = info.reference_names.len();
        let all_targets = self.regions.as_ref().is_none_or(Vec::is_empty);

        let regions: Vec<Regions> = if all_targets {
            (0..total_count)
                .map(|_| Regions {
                    intervals: vec![(1, i64::MAX)],
                    cpos: 0,
                })
                .collect()
        } else {
            self.regions.clone().unwrap_or_default()
        };

        let mut stats = ReferenceStats {
            total_count,
            count: 0,
            combined_length: 0,
            min_length: 0,
            max_length: 0,
            average_length: 0.0,
            average_gc: 0.0,
            rows: Vec::new(),
        };

        let mut gc_sum: f32 = 0.0;
        let mut warned = false;
        let mut buf = Vec::new();

        for (i, reg) in regions.iter().enumerate() {
            if reg.intervals.is_empty() {
                continue;
            }

            let Some(name) = info.reference_names.get(i) else {
                continue;
            };

            let header_length = i64::try_from(info.reference_lengths[i])?;

            for &(beg, original_end) in &reg.intervals {
                let start = beg;
                let end = if original_end == i64::MAX {
                    header_length
                } else {
                    original_end
                };

                if end < start {
                    stats.rows.push(None);
                    continue;
                }

                stats.count += 1;

                let (row_name, length) = if all_targets {
                    (String::from_utf8_lossy(name).into_owned(), header_length)
                } else {
                    let name = String::from_utf8_lossy(name);

                    // A whole sequence is named alone, a range as `name:start-end`.
                    let row_name = if beg == 1 && original_end == i64::MAX {
                        name.into_owned()
                    } else {
                        format!("{name}:{start}-{end}")
                    };

                    (row_name, (end - start + 1).min(header_length))
                };

                stats.combined_length += length;

                if stats.min_length == 0 || stats.min_length > length {
                    stats.min_length = length;
                }

                stats.max_length = stats.max_length.max(length);

                let Some(reference) = &info.reference else {
                    gc_sum = -1.0;
                    stats.rows.push(Some((row_name, length, -1.0, -1)));
                    continue;
                };

                let (mut gc, mut at, mut n) = (0i64, 0i64, 0i64);
                let mut start0 = start - 1;
                let mut rem = end - start0;

                while rem > 0 {
                    let chunk = rem.min(info.ref_chunk_size);

                    if !reference.fetch(name, start0, start0 + chunk - 1, &mut buf)? {
                        if !warned {
                            eprintln!(
                                "Warning: Some sequences not present in the reference, e.g. \"{}\".                             This message is printed only once.",
                                String::from_utf8_lossy(name)
                            );
                            warned = true;
                        }

                        break;
                    }

                    for &b in &buf {
                        match b {
                            b'G' | b'g' | b'C' | b'c' => gc += 1,
                            b'A' | b'a' | b'T' | b't' => at += 1,
                            b'N' | b'n' => n += 1,
                            _ => {}
                        }
                    }

                    let got = i64::try_from(buf.len())?;

                    if got < chunk {
                        break;
                    }

                    start0 += got;
                    rem -= got;
                }

                let gc_content = if gc + at > 0 {
                    gc as f32 / (gc + at) as f32
                } else {
                    0.0
                };

                gc_sum += gc_content;
                stats.rows.push(Some((row_name, length, gc_content, n)));
            }
        }

        if stats.count > 0 {
            stats.average_length = stats.combined_length as f32 / stats.count as f32;
            stats.average_gc = if gc_sum < 0.0 {
                -1.0
            } else {
                gc_sum / stats.count as f32
            };
        } else {
            stats.average_length = -1.0;
            stats.average_gc = -1.0;
        }

        self.reference_stats = Some(stats);

        Ok(())
    }

    /// `output_stats`
    fn output<W: Write>(
        &mut self,
        out: &mut W,
        info: &Info,
        command_line: &str,
        split_tag: Option<&str>,
        sparse: bool,
    ) -> io::Result<()> {
        self.flush_cycles();

        // Insert size average and standard deviation, from the main bulk only.
        self.isize.halve();

        let nitems = self.isize.len();
        let (mut nisize, mut nisize_inward, mut nisize_outward, mut nisize_other) =
            (0u64, 0u64, 0u64, 0u64);

        for isize in 0..nitems {
            let [i, o, x] = self.isize.get(isize);
            nisize_inward += i;
            nisize_outward += o;
            nisize_other += x;
            nisize += i + o + x;
        }

        let mut ibulk = 0;
        let mut bulk = 0.0f64;
        let mut avg_isize = 0.0f64;

        for isize in 0..nitems {
            let [i, o, x] = self.isize.get(isize);
            let num = i + o + x;

            if num > 0 {
                ibulk = isize + 1;
            }

            bulk += num as f64;
            avg_isize += (isize as u64).wrapping_mul(num) as f64;

            if bulk / nisize as f64 > f64::from(info.isize_main_bulk) {
                ibulk = isize + 1;
                nisize = bulk as u64;
                break;
            }
        }

        let divisor = if nisize > 0 { nisize } else { 1 };
        avg_isize /= divisor as f64;

        let mut sd_isize = 0.0f64;

        for isize in 1..ibulk {
            let [i, o, x] = self.isize.get(isize);
            let d = f64::from(isize) - avg_isize;
            sd_isize += (i + o + x) as f64 * d * d / divisor as f64;
        }

        let sd_isize = sd_isize.sqrt();

        writeln!(
            out,
            "# This file was produced by ntools stats ({}) and can be plotted using plot-bamstats",
            env!("CARGO_PKG_VERSION")
        )?;

        match (&self.split_name, split_tag) {
            (Some(name), Some(tag)) => writeln!(
                out,
                "# This file contains statistics only for reads with tag: {tag}={}",
                String::from_utf8_lossy(name)
            )?,
            _ => writeln!(out, "# This file contains statistics for all reads.")?,
        }

        writeln!(out, "# The command line was:  {command_line}")?;
        writeln!(
            out,
            "# CHK, Checksum\t[2]Read Names\t[3]Sequences\t[4]Qualities"
        )?;
        writeln!(
            out,
            "# CHK, CRC32 of reads which passed filtering followed by addition (32bit overflow)"
        )?;
        writeln!(
            out,
            "CHK\t{:08x}\t{:08x}\t{:08x}",
            self.checksum_names, self.checksum_reads, self.checksum_quals
        )?;

        let nreads = self.nreads_1st + self.nreads_2nd + self.nreads_other;

        writeln!(
            out,
            "# Summary Numbers. Use `grep ^SN | cut -f 2-` to extract this part."
        )?;
        writeln!(
            out,
            "SN\traw total sequences:\t{}\t# excluding supplementary and secondary reads",
            self.nreads_filtered + nreads
        )?;
        writeln!(out, "SN\tfiltered sequences:\t{}", self.nreads_filtered)?;
        writeln!(out, "SN\tsequences:\t{nreads}")?;
        writeln!(
            out,
            "SN\tis sorted:\t{}\t# {} by coordinate",
            u8::from(self.is_sorted),
            if self.is_sorted {
                "sorted"
            } else {
                "not sorted"
            }
        )?;
        writeln!(out, "SN\t1st fragments:\t{}", self.nreads_1st)?;
        writeln!(out, "SN\tlast fragments:\t{}", self.nreads_2nd)?;
        writeln!(
            out,
            "SN\treads mapped:\t{}",
            self.nreads_paired_and_mapped + self.nreads_single_mapped
        )?;
        writeln!(
            out,
            "SN\treads mapped and paired:\t{}\t# paired-end technology bit set + both mates mapped",
            self.nreads_paired_and_mapped
        )?;
        writeln!(out, "SN\treads unmapped:\t{}", self.nreads_unmapped)?;
        writeln!(
            out,
            "SN\treads properly paired:\t{}\t# proper-pair bit set",
            self.nreads_properly_paired
        )?;
        writeln!(
            out,
            "SN\treads paired:\t{}\t# paired-end technology bit set",
            self.nreads_paired_tech
        )?;
        writeln!(
            out,
            "SN\treads duplicated:\t{}\t# PCR or optical duplicate bit set",
            self.nreads_dup
        )?;
        writeln!(
            out,
            "SN\treads MQ0:\t{}\t# mapped and MQ=0",
            self.nreads_mq0
        )?;
        writeln!(out, "SN\treads QC failed:\t{}", self.nreads_qc_failed)?;
        writeln!(
            out,
            "SN\tnon-primary alignments:\t{}",
            self.nreads_secondary
        )?;
        writeln!(
            out,
            "SN\tsupplementary alignments:\t{}",
            self.nreads_supplementary
        )?;
        writeln!(
            out,
            "SN\ttotal length:\t{}\t# ignores clipping",
            self.total_len
        )?;
        writeln!(
            out,
            "SN\ttotal first fragment length:\t{}\t# ignores clipping",
            self.total_len_1st
        )?;
        writeln!(
            out,
            "SN\ttotal last fragment length:\t{}\t# ignores clipping",
            self.total_len_2nd
        )?;
        writeln!(
            out,
            "SN\tbases mapped:\t{}\t# ignores clipping",
            self.nbases_mapped
        )?;
        writeln!(
            out,
            "SN\tbases mapped (cigar):\t{}\t# more accurate",
            self.nbases_mapped_cigar.cast_signed()
        )?;
        writeln!(out, "SN\tbases trimmed:\t{}", self.nbases_trimmed)?;
        writeln!(out, "SN\tbases duplicated:\t{}", self.total_len_dup)?;
        writeln!(
            out,
            "SN\tmismatches:\t{}\t# from NM fields",
            self.nmismatches.cast_signed()
        )?;

        let error_rate = if self.nbases_mapped_cigar > 0 {
            f64::from(self.nmismatches as f32 / self.nbases_mapped_cigar as f32)
        } else {
            0.0
        };

        writeln!(
            out,
            "SN\terror rate:\t{}\t# mismatches / bases mapped (cigar)",
            format_e(error_rate)
        )?;

        let avg_read_length: f32 = if nreads > 0 {
            self.total_len as f32 / nreads as f32
        } else {
            0.0
        };

        let ratio = |a: u64, b: u64| -> f64 {
            if b > 0 {
                f64::from(a as f32 / b as f32)
            } else {
                0.0
            }
        };

        writeln!(
            out,
            "SN\taverage length:\t{:.0}",
            f64::from(avg_read_length)
        )?;
        writeln!(
            out,
            "SN\taverage first fragment length:\t{:.0}",
            ratio(self.total_len_1st, self.nreads_1st)
        )?;
        writeln!(
            out,
            "SN\taverage last fragment length:\t{:.0}",
            ratio(self.total_len_2nd, self.nreads_2nd)
        )?;
        writeln!(out, "SN\tmaximum length:\t{}", self.max_len)?;
        writeln!(
            out,
            "SN\tmaximum first fragment length:\t{}",
            self.max_len_1st
        )?;
        writeln!(
            out,
            "SN\tmaximum last fragment length:\t{}",
            self.max_len_2nd
        )?;
        writeln!(
            out,
            "SN\taverage quality:\t{:.1}",
            if self.total_len > 0 {
                self.sum_qual as f64 / self.total_len as f64
            } else {
                0.0
            }
        )?;
        writeln!(out, "SN\tinsert size average:\t{avg_isize:.1}")?;
        writeln!(out, "SN\tinsert size standard deviation:\t{sd_isize:.1}")?;
        writeln!(out, "SN\tinward oriented pairs:\t{nisize_inward}")?;
        writeln!(out, "SN\toutward oriented pairs:\t{nisize_outward}")?;
        writeln!(out, "SN\tpairs with other orientation:\t{nisize_other}")?;
        writeln!(
            out,
            "SN\tpairs on different chromosomes:\t{}",
            self.nreads_anomalous / 2
        )?;
        writeln!(
            out,
            "SN\tpercentage of properly paired reads (%):\t{:.1}",
            ratio(100 * self.nreads_properly_paired, nreads)
        )?;

        if self.target_count > 0 {
            writeln!(out, "SN\tbases inside the target:\t{}", self.target_count)?;

            let threshold = usize::try_from(info.cov_threshold).unwrap_or(0);
            let cov_sum: u64 = self.cov.iter().skip(threshold + 1).sum();

            writeln!(
                out,
                "SN\tpercentage of target genome with coverage > {} (%):\t{:.2}",
                info.cov_threshold,
                f64::from((100 * cov_sum) as f32 / self.target_count as f32)
            )?;
        }

        if self.max_len < self.nbases {
            self.max_len += 1;
        }

        if self.max_qual + 1 < QUALITY_BINS && self.max_qual < 255 {
            self.max_qual += 1;
        }

        writeln!(
            out,
            "# First Fragment Qualities. Use `grep ^FFQ | cut -f 2-` to extract this part."
        )?;
        writeln!(
            out,
            "# Columns correspond to qualities and rows to cycles. First column is the cycle number."
        )?;
        write_quality_rows(
            out,
            "FFQ",
            &self.quals_1st,
            self.quals_rows,
            self.max_len_1st,
            self.max_qual,
        )?;

        writeln!(
            out,
            "# Last Fragment Qualities. Use `grep ^LFQ | cut -f 2-` to extract this part."
        )?;
        writeln!(
            out,
            "# Columns correspond to qualities and rows to cycles. First column is the cycle number."
        )?;
        write_quality_rows(
            out,
            "LFQ",
            &self.quals_2nd,
            self.quals_rows,
            self.max_len_2nd,
            self.max_qual,
        )?;

        if let Some(mpc) = &self.mpc_buf {
            writeln!(
                out,
                "# Mismatches per cycle and quality. Use `grep ^MPC | cut -f 2-` to extract this part."
            )?;
            writeln!(
                out,
                "# Columns correspond to qualities, rows to cycles. First column is the cycle number, second"
            )?;
            writeln!(
                out,
                "# is the number of N's and the rest is the number of mismatches"
            )?;

            write_quality_rows(out, "MPC", mpc, QUALITY_BINS, self.max_len, self.max_qual)?;
        }

        writeln!(
            out,
            "# GC Content of first fragments. Use `grep ^GCF | cut -f 2-` to extract this part."
        )?;
        write_gc_rows(out, "GCF", &self.gc_1st, self.ngc)?;
        writeln!(
            out,
            "# GC Content of last fragments. Use `grep ^GCL | cut -f 2-` to extract this part."
        )?;
        write_gc_rows(out, "GCL", &self.gc_2nd, self.ngc)?;

        writeln!(
            out,
            "# ACGT content per cycle. Use `grep ^GCC | cut -f 2-` to extract this part. The columns are: cycle; A,C,G,T base counts as a percentage of all A/C/G/T bases [%]; and N and O counts as a percentage of all A/C/G/T bases [%]"
        )?;

        for ibase in 0..self.max_len {
            let (c1, c2) = (
                BaseCounts::from(&self.acgtno_cycles_1st[ibase]),
                BaseCounts::from(&self.acgtno_cycles_2nd[ibase]),
            );
            let sum = c1.a + c1.c + c1.g + c1.t + c2.a + c2.c + c2.g + c2.t;

            if sum == 0 {
                continue;
            }

            write_percent_row(
                out,
                "GCC",
                ibase + 1,
                &[
                    c1.a + c2.a,
                    c1.c + c2.c,
                    c1.g + c2.g,
                    c1.t + c2.t,
                    c1.n + c2.n,
                    c1.other + c2.other,
                ],
                sum,
            )?;
        }

        writeln!(
            out,
            "# ACGT content per cycle, read oriented. Use `grep ^GCT | cut -f 2-` to extract this part. The columns are: cycle; A,C,G,T base counts as a percentage of all A/C/G/T bases [%]"
        )?;

        for ibase in 0..self.max_len {
            let c = BaseCounts::from(&self.acgtno_revcomp[ibase]);
            let sum = c.a + c.c + c.g + c.t;

            if sum == 0 {
                continue;
            }

            write_percent_row(out, "GCT", ibase + 1, &[c.a, c.c, c.g, c.t], sum)?;
        }

        for (label, total_label, ordinal, cycles) in [
            ("FBC", "FTC", "first", &self.acgtno_cycles_1st),
            ("LBC", "LTC", "last", &self.acgtno_cycles_2nd),
        ] {
            writeln!(
                out,
                "# ACGT content per cycle for {ordinal} fragments. Use `grep ^{label} | cut -f 2-` to extract this part. The columns are: cycle; A,C,G,T base counts as a percentage of all A/C/G/T bases [%]; and N and O counts as a percentage of all A/C/G/T bases [%]"
            )?;

            let mut totals = BaseCounts::default();

            for (ibase, c) in cycles.iter().take(self.max_len).enumerate() {
                let c = BaseCounts::from(c);
                let sum = c.a + c.c + c.g + c.t;
                totals.a += c.a;
                totals.c += c.c;
                totals.g += c.g;
                totals.t += c.t;
                totals.n += c.n;

                if sum > 0 {
                    write_percent_row(
                        out,
                        label,
                        ibase + 1,
                        &[c.a, c.c, c.g, c.t, c.n, c.other],
                        sum,
                    )?;
                }
            }

            writeln!(
                out,
                "# ACGT raw counters for {ordinal} fragments. Use `grep ^{total_label} | cut -f 2-` to extract this part. The columns are: A,C,G,T,N base counters"
            )?;
            writeln!(
                out,
                "{total_label}\t{}\t{}\t{}\t{}\t{}",
                totals.a, totals.c, totals.g, totals.t, totals.n
            )?;
        }

        writeln!(
            out,
            "# Insert sizes. Use `grep ^IS | cut -f 2-` to extract this part. The columns are: insert size, pairs total, inward oriented pairs, outward oriented pairs, other pairs"
        )?;

        for isize in 0..ibulk {
            let [i, o, x] = self.isize.get(isize);

            if !sparse || i + o + x > 0 {
                writeln!(out, "IS\t{isize}\t{}\t{i}\t{o}\t{x}", i + o + x)?;
            }
        }

        for (label, description, lengths, max) in [
            ("RL", "", &self.read_lengths, self.max_len),
            (
                "FRL",
                " - first fragments",
                &self.read_lengths_1st,
                self.max_len_1st,
            ),
            (
                "LRL",
                " - last fragments",
                &self.read_lengths_2nd,
                self.max_len_2nd,
            ),
        ] {
            writeln!(
                out,
                "# Read lengths{description}. Use `grep ^{label} | cut -f 2-` to extract this part. The columns are: read length, count"
            )?;

            for ilen in 0..max {
                if lengths[ilen + 1] > 0 {
                    writeln!(out, "{label}\t{}\t{}", ilen + 1, lengths[ilen + 1])?;
                }
            }
        }

        writeln!(
            out,
            "# Mapping qualities for reads !(UNMAP|SECOND|SUPPL|QCFAIL|DUP). Use `grep ^MAPQ | cut -f 2-` to extract this part. The columns are: mapq, count"
        )?;

        for (mapq, &n) in self.mapping_qualities.iter().enumerate() {
            if n > 0 {
                writeln!(out, "MAPQ\t{mapq}\t{n}")?;
            }
        }

        writeln!(
            out,
            "# Indel distribution. Use `grep ^ID | cut -f 2-` to extract this part. The columns are: length, number of insertions, number of deletions"
        )?;

        for ilen in 0..self.nindels {
            if self.insertions[ilen] > 0 || self.deletions[ilen] > 0 {
                writeln!(
                    out,
                    "ID\t{}\t{}\t{}",
                    ilen + 1,
                    self.insertions[ilen],
                    self.deletions[ilen]
                )?;
            }
        }

        writeln!(
            out,
            "# Indels per cycle. Use `grep ^IC | cut -f 2-` to extract this part. The columns are: cycle, number of insertions (fwd), .. (rev) , number of deletions (fwd), .. (rev)"
        )?;

        for ilen in 0..=self.nbases {
            let counts = [
                self.ins_cycles_1st[ilen],
                self.ins_cycles_2nd[ilen],
                self.del_cycles_1st[ilen],
                self.del_cycles_2nd[ilen],
            ];

            if counts.iter().any(|&n| n > 0) {
                writeln!(
                    out,
                    "IC\t{}\t{}\t{}\t{}\t{}",
                    ilen + 1,
                    counts[0],
                    counts[1],
                    counts[2],
                    counts[3]
                )?;
            }
        }

        if self.is_sorted {
            self.write_coverage(out, info)?;
            self.write_gc_depth(out, info, avg_read_length)?;
        }

        if let Some(rs) = &self.reference_stats {
            writeln!(
                out,
                "# Reference statistics. Use `grep ^RFS | cut -f 2-` to extract this part."
            )?;
            writeln!(
                out,
                "# Total count, Output count, Average GC, Min length, Max length, Average length, Total length in first row."
            )?;
            writeln!(
                out,
                "# Sequence name, Length, GC content, Unknown count in following rows."
            )?;
            writeln!(
                out,
                "RFS\t{}\t{}\t{:.2}\t{}\t{}\t{:.2}\t{}",
                rs.total_count,
                rs.count,
                f64::from(rs.average_gc),
                rs.min_length,
                rs.max_length,
                f64::from(rs.average_length),
                rs.combined_length
            )?;

            for (name, length, gc, n) in rs.rows.iter().flatten() {
                writeln!(out, "RFS\t{name}\t{length}\t{:.2}\t{n}", f64::from(*gc))?;
            }
        }

        Ok(())
    }

    fn write_coverage<W: Write>(&self, out: &mut W, info: &Info) -> io::Result<()> {
        writeln!(
            out,
            "# Coverage distribution. Use `grep ^COV | cut -f 2-` to extract this part."
        )?;

        let (min, step) = (info.cov_min, info.cov_step);

        if self.cov[0] > 0 {
            writeln!(out, "COV\t[<{min}]\t{}\t{}", min - 1, self.cov[0])?;
        }

        for icov in 1..self.ncov - 1 {
            if self.cov[icov] > 0 {
                let icov_i32 = i32::try_from(icov).unwrap_or(i32::MAX);
                let lo = min + (icov_i32 - 1) * step;
                let hi = min + icov_i32 * step - 1;
                writeln!(out, "COV\t[{lo}-{hi}]\t{hi}\t{}", self.cov[icov])?;
            }
        }

        if self.cov[self.ncov - 1] > 0 {
            let last = min + (i32::try_from(self.ncov).unwrap_or(i32::MAX) - 2) * step - 1;
            writeln!(out, "COV\t[{last}<]\t{last}\t{}", self.cov[self.ncov - 1])?;
        }

        Ok(())
    }

    fn write_gc_depth<W: Write>(
        &mut self,
        out: &mut W,
        info: &Info,
        avg_read_length: f32,
    ) -> io::Result<()> {
        writeln!(
            out,
            "# GC-depth. Use `grep ^GCD | cut -f 2-` to extract this part. The columns are: GC%, unique sequence percentiles, 10th, 25th, 50th, 75th and 90th depth percentile"
        )?;

        // As in samtools, the last bin keeps its unrounded GC, and is sorted but not reported.
        let reported = self.igcd.min(self.gcd.len());

        for bin in &mut self.gcd[..reported] {
            if info.reference.is_some() {
                bin.gc = (100.0 * f64::from(bin.gc)).round_ties_even() as f32;
            } else if bin.depth > 0 {
                bin.gc =
                    (100.0 * f64::from(bin.gc) / f64::from(bin.depth)).round_ties_even() as f32;
            }
        }

        if !self.gcd.is_empty() {
            self.gcd[..=self.igcd].sort_by(|a, b| {
                a.gc.partial_cmp(&b.gc)
                    .unwrap_or(std::cmp::Ordering::Equal)
                    .then(a.depth.cmp(&b.depth))
            });
        }

        let bin_size = info.gcd_bin_size as f32;
        let mut igcd = 0;

        while igcd < self.igcd {
            let gc = self.gcd[igcd].gc;
            let mut nbins = 0;

            while igcd + nbins < self.igcd && f64::from(self.gcd[igcd + nbins].gc - gc).abs() < 0.1
            {
                nbins += 1;
            }

            let bins = &self.gcd[igcd..];
            let p = |percentile: u32| {
                f64::from(gcd_percentile(bins, nbins, percentile) * avg_read_length / bin_size)
            };

            writeln!(
                out,
                "GCD\t{:.1}\t{:.3}\t{:.3}\t{:.3}\t{:.3}\t{:.3}\t{:.3}",
                f64::from(gc),
                ((igcd + nbins + 1) as u32) as f64 * 100.0 / (self.igcd + 1) as f64,
                p(10),
                p(25),
                p(50),
                p(75),
                p(90)
            )?;

            igcd += nbins;
        }

        Ok(())
    }
}

/// `lidx2ridx`
fn lidx2ridx(offset: i32, size: i32, refpos: i64, pos: i64) -> i32 {
    let size64 = i64::from(size);
    i32::try_from((i64::from(offset) + (pos - refpos) % size64) % size64).unwrap_or(0)
}

/// `cleanup_overlaps`: drops pairs whose chunks all end before `max`, returning how many.
fn cleanup_overlaps(pairs: &mut HashMap<Vec<u8>, Pair>, max: i64) -> u32 {
    let before = pairs.len();
    pairs.retain(|_, pair| pair.chunks.last().is_some_and(|&(_, end)| end >= max));
    u32::try_from(before - pairs.len()).unwrap_or(u32::MAX)
}

/// `bwa_trim_read`: the number of bases BWA would trim from the 3' end.
fn bwa_trim_read(trim_qual: i32, quals: &[u8], len: usize, reverse: bool) -> usize {
    if len < BWA_MIN_READ_LENGTH {
        return 0;
    }

    let max_trimmed = len - BWA_MIN_READ_LENGTH + 1;
    let (mut sum, mut max_sum, mut max_l) = (0i32, 0i32, 0usize);

    for l in 0..max_trimmed {
        sum += trim_qual - i32::from(quals[if reverse { l } else { len - 1 - l }]);

        if sum < 0 {
            break;
        }

        if sum > max_sum {
            max_sum = sum;
            // BWA clips one base less than the maximum.
            max_l = l;
        }
    }

    max_l
}

/// `gcd_percentile`, in single precision.
fn gcd_percentile(gcd: &[GcDepth], n: usize, p: u32) -> f32 {
    let x = p as f32 * (n as u32 + 1) as f32 / 100.0;
    #[expect(
        clippy::cast_possible_truncation,
        reason = "samtools truncates the float to an int"
    )]
    let k = x as i32;

    if k <= 0 {
        return gcd[0].depth as f32;
    }

    let k = k as usize;

    if k >= n {
        return gcd[n - 1].depth as f32;
    }

    let d = x - k as f32;
    gcd[k - 1].depth as f32 + d * gcd[k].depth.wrapping_sub(gcd[k - 1].depth) as f32
}

fn write_quality_rows<W: Write>(
    out: &mut W,
    label: &str,
    quals: &[u64],
    rows: usize,
    max_len: usize,
    max_qual: usize,
) -> io::Result<()> {
    let mut line = Vec::new();

    for (ibase, counts) in quals.chunks_exact(rows).take(max_len).enumerate() {
        line.clear();
        line.extend_from_slice(label.as_bytes());
        line.push(b'\t');
        push_u64(&mut line, ibase as u64 + 1);

        for &n in &counts[..=max_qual] {
            line.push(b'\t');
            push_u64(&mut line, n);
        }

        line.push(b'\n');
        out.write_all(&line)?;
    }

    Ok(())
}

/// Writes `label`, the cycle, and each count as a percentage of `sum` (`%.2f`).
fn write_percent_row<W: Write>(
    out: &mut W,
    label: &str,
    cycle: usize,
    counts: &[u64],
    sum: u64,
) -> io::Result<()> {
    let mut line = Vec::with_capacity(64);
    line.extend_from_slice(label.as_bytes());
    line.push(b'\t');
    push_u64(&mut line, cycle as u64);

    for &n in counts {
        line.push(b'\t');
        push_fixed2(&mut line, 100.0 * n as f64 / sum as f64);
    }

    line.push(b'\n');
    out.write_all(&line)
}

fn push_u64(dst: &mut Vec<u8>, n: u64) {
    let mut buf = [0; lexical_core::BUFFER_SIZE];
    dst.extend_from_slice(lexical_core::write(n, &mut buf));
}

/// Appends `printf("%.2f", value)` for a nonnegative value.
fn push_fixed2(dst: &mut Vec<u8>, value: f64) {
    // `value * 100` is within an ulp of the exact product, so unless it lies near a rounding
    // boundary, rounding it gives the digits printf would; otherwise defer to exact formatting.
    let scaled = value * 100.0;
    let rounded = scaled.round();

    if value.is_finite()
        && (0.0..1e15).contains(&value)
        && ((scaled - rounded).abs() - 0.5).abs() > 1e-6
    {
        let n = rounded as u64;
        push_u64(dst, n / 100);
        dst.push(b'.');
        dst.push(b'0' + (n % 100 / 10) as u8);
        dst.push(b'0' + (n % 10) as u8);
    } else {
        dst.extend_from_slice(format!("{value:.2}").as_bytes());
    }
}

fn write_gc_rows<W: Write>(out: &mut W, label: &str, counts: &[u64], ngc: usize) -> io::Result<()> {
    let mut prev = 0;

    for ibase in 0..ngc {
        if counts[ibase] == counts[prev] {
            continue;
        }

        writeln!(
            out,
            "{label}\t{:.2}\t{}",
            (ibase + prev) as f64 * 0.5 * 100.0 / (ngc - 1) as f64,
            counts[prev]
        )?;
        prev = ibase;
    }

    Ok(())
}

/// The read group IDs to include for `-I`: those with this ID or sample name.
fn read_groups_for(header: &sam::Header, id: &str) -> Vec<Vec<u8>> {
    use sam::header::record::value::map::read_group::tag::SAMPLE;

    header
        .read_groups()
        .iter()
        .filter(|(key, map)| {
            key.as_slice() == id.as_bytes()
                || map
                    .other_fields()
                    .get(&SAMPLE)
                    .is_some_and(|sm| sm.as_slice() == id.as_bytes())
        })
        .map(|(key, _)| key.to_vec())
        .collect()
}

/// Formats like C's `printf("%e", value)`.
fn format_e(value: f64) -> String {
    let s = format!("{value:.6e}");
    let (mantissa, exponent) = s.split_once('e').expect("scientific notation");
    let exponent: i32 = exponent.parse().expect("valid exponent");
    let sign = if exponent < 0 { '-' } else { '+' };
    format!("{mantissa}e{sign}{:02}", exponent.abs())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_format_e() {
        assert_eq!(format_e(0.0), "0.000000e+00");
        assert_eq!(format_e(0.001_234_567_8), "1.234568e-03");
        assert_eq!(format_e(12345.0), "1.234500e+04");
    }

    #[test]
    fn test_bwa_trim_read() {
        let quals = [30; 50];
        assert_eq!(bwa_trim_read(20, &quals, 50, false), 0);

        let mut quals = [30; 50];
        quals[45..].fill(2);
        assert_eq!(bwa_trim_read(20, &quals, 50, false), 4);
        assert_eq!(bwa_trim_read(20, &quals, 30, false), 0);
    }

    #[test]
    fn test_parse_coverage() -> anyhow::Result<()> {
        assert_eq!(parse_coverage("1,1000,1")?, (1, 1000, 1));
        assert_eq!(parse_coverage("5,50,5x")?, (5, 50, 5));
        assert!(parse_coverage("1,2").is_err());
        Ok(())
    }
}
