//! `ntools stats`: alignment statistics for plot-bamstats (samtools stats).
//!
//! This is a port of `stats.c` from samtools 1.24, kept close to the original so that the
//! output matches it byte for byte: the same buffers and their growth rules, single-precision
//! (`float`) arithmetic where samtools uses it, and its quirks where they show in the output
//! (for example, the GC-depth table leaves out its last bin, and the coverage buffer keeps only
//! some of its contents when a long read enlarges it).

mod collect;
mod coverage;
mod cycles;
mod record;
mod regions;
mod report;

use std::{
    collections::HashMap,
    fs::File,
    io::{BufWriter, Write},
    path::{Path, PathBuf},
};

use anyhow::{Context, bail};
use noodles_sam::{self as sam, alignment::record::Flags};

use crate::{
    alignment::{AlignmentReader, Queries, RawEncoder, ReadOptions, Rec},
    bam_raw, flags, output,
    reference::IndexedFasta as Reference,
    regions::RegionSet,
};

use self::{
    cycles::{BaseCounts, PendingCycles},
    record::Bam,
    report::read_groups_for,
};
use crate::global_args::{IgnoredWriteIndex, InputFormatOptions, Verbosity};

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

    #[command(flatten)]
    input_format_options: InputFormatOptions,

    #[command(flatten)]
    write_index: IgnoredWriteIndex,

    #[command(flatten)]
    verbosity: Verbosity,
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
    insert_size_rows: InsertSizeRows,
}

/// Which insert sizes the report lists.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
enum InsertSizeRows {
    All,
    /// Only insert sizes with pairs (`-x`).
    WithPairs,
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
    barcodes: Vec<Barcode>,
    /// Barcode errors reported so far (`error_number`).
    barcode_errors: u32,
}

/// `barcode_info_t` and its counts: a barcode tag and its quality tag.
struct Barcode {
    tag: [u8; 2],
    quality_tag: [u8; 2],
    /// The barcode length, from the first barcode seen; 0 until then.
    len: usize,
    /// The position of a separator (a character other than A, C, G, T, or N), if any.
    separator: Option<usize>,
    max_qual: Option<usize>,
    counts: Vec<BaseCounts>,
    /// `[position * QUALITY_BINS + quality]`
    quals: Vec<u64>,
}

impl Barcode {
    fn new(tag: [u8; 2], quality_tag: [u8; 2]) -> Self {
        Self {
            tag,
            quality_tag,
            len: 0,
            separator: None,
            max_qual: None,
            counts: Vec::new(),
            quals: Vec::new(),
        }
    }

    /// The (segment, cycle) a position is reported as: positions after a separator form a
    /// second segment.
    fn segment(&self, i: usize) -> (usize, usize) {
        match self.separator {
            Some(sep) if i >= sep => (2, i - sep),
            _ => (1, i + 1),
        }
    }
}

/// The maximum barcode errors reported (`ERROR_LIMIT`).
const BARCODE_ERROR_LIMIT: u32 = 200;

pub fn run(mut args: Args) -> anyhow::Result<()> {
    let index = if args.customized_index {
        if args.regions.is_empty() {
            bail!("No index file provided");
        }

        Some(PathBuf::from(args.regions.remove(0)))
    } else {
        None
    };

    let mut read_options = ReadOptions {
        threads: args.threads,
        reference: args.cram_reference.clone(),
    };
    args.input_format_options.configure(&mut read_options)?;

    let (mut reader, path) = AlignmentReader::open(&args.input, &read_options)?;
    args.input_format_options.apply(&mut reader)?;
    let header = reader.read_header()?;

    let mut info = Info::from_args(&args, &header)?;

    let mut all = Stats::new(
        &mut info,
        &header,
        args.group_id.as_deref(),
        args.targets.as_deref(),
    )?;
    let mut pairs: HashMap<Vec<u8>, Pair> = HashMap::new();
    let mut splits = Splits::new(args.split_tag.as_deref())?;
    let mut encoder = RawEncoder::new();

    let mut process = |rec: Rec<'_>,
                       all: &mut Stats,
                       info: &Info,
                       pairs: &mut HashMap<Vec<u8>, Pair>|
     -> anyhow::Result<()> {
        let raw: &[u8] = rec.encode_raw(&header, &mut encoder)?;
        let bam = Bam::parse(raw)?;
        splits.collect(&bam, info, pairs, args.targets.as_deref())?;
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

    write_reports(&args, &info, all, splits)?;

    Ok(())
}

impl Info {
    fn from_args(args: &Args, header: &sam::Header) -> anyhow::Result<Self> {
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

        Ok(Info {
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
            insert_size_rows: if args.sparse {
                InsertSizeRows::WithPairs
            } else {
                InsertSizeRows::All
            },
        })
    }
}

/// Statistics by the value of a tag (`-S`), each written to its own file.
struct Splits {
    tag: Option<[u8; 2]>,
    stats: Vec<(Vec<u8>, Stats)>,
    index: HashMap<Vec<u8>, usize>,
}

impl Splits {
    fn new(tag: Option<&str>) -> anyhow::Result<Self> {
        let tag = tag
            .map(|tag| -> anyhow::Result<[u8; 2]> {
                tag.as_bytes()
                    .try_into()
                    .with_context(|| format!("invalid tag: {tag}"))
            })
            .transpose()?;

        Ok(Self {
            tag,
            stats: Vec::new(),
            index: HashMap::new(),
        })
    }

    /// Collects a record into the statistics of its tag value, if splitting.
    fn collect(
        &mut self,
        bam: &Bam<'_>,
        info: &Info,
        pairs: &mut HashMap<Vec<u8>, Pair>,
        targets: Option<&Path>,
    ) -> anyhow::Result<()> {
        let Some(tag) = self.tag else {
            return Ok(());
        };

        let raw = bam.raw;

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

        let i = match self.index.get(&value) {
            Some(&i) => i,
            None => {
                let mut stats = Stats::new_with_info(info, targets)?;
                stats.split_name = Some(value.clone());
                self.stats.push((value.clone(), stats));
                self.index.insert(value, self.stats.len() - 1);
                self.stats.len() - 1
            }
        };

        self.stats[i].1.collect(bam, info, pairs)
    }
}

/// Writes the report for all reads to stdout, and one per split to `PREFIX_VALUE.bamstat`.
fn write_reports(args: &Args, info: &Info, mut all: Stats, splits: Splits) -> anyhow::Result<()> {
    all.round_buffer_flush(-1, info)?;

    let command_line = std::iter::once("stats".to_owned())
        .chain(std::env::args().skip(2))
        .collect::<Vec<_>>()
        .join(" ");

    let mut out = BufWriter::new(output::raw_stdout()?);
    all.output(&mut out, info, &command_line, None)?;
    out.flush()?;

    for (name, mut stats) in splits.stats {
        stats.round_buffer_flush(-1, info)?;

        let prefix = args.split_prefix.as_deref().unwrap_or(&args.input);
        let path = format!("{prefix}_{}.bamstat", String::from_utf8_lossy(&name));
        let file =
            File::create(&path).with_context(|| format!("Could not open '{path}' for writing."))?;
        let mut out = BufWriter::new(file);
        stats.output(&mut out, info, &command_line, args.split_tag.as_deref())?;
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
            barcodes: vec![
                Barcode::new(*b"BC", *b"QT"),
                Barcode::new(*b"CR", *b"CY"),
                Barcode::new(*b"OX", *b"BZ"),
                Barcode::new(*b"RX", *b"QX"),
            ],
            barcode_errors: 0,
        };

        stats.realloc_rseq_buffer(info);

        if let Some(path) = targets {
            stats.init_regions(path, info)?;
        }

        Ok(stats)
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
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_parse_coverage() -> anyhow::Result<()> {
        assert_eq!(parse_coverage("1,1000,1")?, (1, 1000, 1));
        assert_eq!(parse_coverage("5,50,5x")?, (5, 50, 5));
        assert!(parse_coverage("1,2").is_err());
        Ok(())
    }
}
