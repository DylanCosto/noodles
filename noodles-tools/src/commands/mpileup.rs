//! `ntools mpileup`: text pileup (samtools mpileup).
//!
//! A port of `bam_plcmd.c` from samtools 1.24 on top of the htslib pileup engine port in
//! [`crate::pileup`] and the BAQ port in [`crate::baq`]. Each input is read on its own thread,
//! and reads are filtered and realigned (BAQ) in batches across the `-@` threads, keeping their
//! order; samtools does this work on one thread.
//!
//! Not supported: `--output-mods` (`-M`).

use std::{
    collections::HashSet,
    fs::{self, File},
    io::{self, BufWriter, Write},
    path::{Path, PathBuf},
    sync::{Arc, Mutex},
    thread,
};

use anyhow::{Context, bail};
use noodles_core::Region;
use noodles_pileup::{
    self as pileup, Entry, Overlaps, Pileup, Read,
    baq::{self, SEQ_NT16_TABLE},
    mods::BaseMods,
};
use noodles_sam as sam;

use crate::{
    alignment::{AlignmentReader, Queries, ReadOptions, Rec},
    flags, output,
    reference::IndexedFasta,
};

/// Reads per batch sent from a reader thread.
const BATCH_SIZE: usize = 4096;
/// Batches queued per input.
const QUEUE_DEPTH: usize = 4;

// Extra column flags, in samtools' output order (`MPLP_PRINT_*`).
const PRINT_MAPQ_CHAR: u32 = 1 << 11;
const PRINT_QPOS: u32 = 1 << 12;
const PRINT_QNAME: u32 = 1 << 13;
const PRINT_FLAG: u32 = 1 << 14;
const PRINT_RNAME: u32 = 1 << 15;
const PRINT_POS: u32 = 1 << 16;
const PRINT_MAPQ: u32 = 1 << 17;
const PRINT_RNEXT: u32 = 1 << 19;
const PRINT_PNEXT: u32 = 1 << 20;
const PRINT_RLEN: u32 = 1 << 24;
const PRINT_MODS: u32 = 1 << 25;
const PRINT_QPOS5: u32 = 1 << 26;
const PRINT_LAST: u32 = 1 << 27;

#[derive(Debug, clap::Args)]
#[command(disable_help_flag = true)]
pub struct Args {
    /// Input SAM, BAM, or CRAM files; with -X, followed by their indexes.
    inputs: Vec<String>,

    /// Count anomalous read pairs (not properly paired).
    #[arg(short = 'A', long = "count-orphans")]
    count_orphans: bool,

    /// A file listing the input files, one per line.
    #[arg(short = 'b', long = "bam-list")]
    bam_list: Option<PathBuf>,

    /// Disable BAQ (per-base alignment quality).
    #[arg(short = 'B', long = "no-BAQ", visible_alias = "no-baq")]
    no_baq: bool,

    /// Adjust mapping quality; 50 is recommended for BWA (0 disables).
    #[arg(
        short = 'C',
        long = "adjust-MQ",
        visible_alias = "adjust-mq",
        default_value_t = 0
    )]
    adjust_mq: i32,

    /// Maximum reads per file at a position (0 for no limit).
    #[arg(short = 'd', long = "max-depth", default_value_t = 8000)]
    max_depth: i64,

    /// Recalculate BAQ on the fly, ignoring existing BQ tags.
    #[arg(short = 'E', long = "redo-BAQ", visible_alias = "redo-baq")]
    redo_baq: bool,

    /// Indexed reference FASTA.
    #[arg(short = 'f', long = "fasta-ref")]
    fasta_ref: Option<PathBuf>,

    /// Exclude read groups listed in this file.
    #[arg(short = 'G', long = "exclude-RG", visible_alias = "exclude-rg")]
    exclude_rg: Option<PathBuf>,

    /// Positions (chr pos) or regions (BED) to output.
    #[arg(short = 'l', long = "positions")]
    positions: Option<PathBuf>,

    /// Region in which pileup is generated.
    #[arg(short = 'r', long = "region")]
    region: Option<String>,

    /// Ignore read groups (one sample per file).
    #[arg(short = 'R', long = "ignore-RG", visible_alias = "ignore-rg")]
    ignore_rg: bool,

    /// Skip alignments with a mapping quality below this.
    #[arg(
        short = 'q',
        long = "min-MQ",
        visible_alias = "min-mq",
        default_value_t = 0
    )]
    min_mq: i32,

    /// Skip bases with a base quality below this.
    #[arg(
        short = 'Q',
        long = "min-BQ",
        visible_alias = "min-bq",
        default_value_t = 13
    )]
    min_bq: i32,

    /// Required flags: skip reads with none of these set.
    #[arg(long = "rf", visible_alias = "incl-flags", value_parser = flags::parse)]
    required_flags: Option<sam::alignment::record::Flags>,

    /// Filter flags: skip reads with any of these set.
    #[arg(long = "ff", visible_alias = "excl-flags", value_parser = flags::parse)]
    filter_flags: Option<sam::alignment::record::Flags>,

    /// Disable read-pair overlap detection.
    #[arg(
        short = 'x',
        long = "ignore-overlaps-removal",
        visible_alias = "disable-overlap-removal"
    )]
    no_overlaps: bool,

    /// Use the index files given after the inputs.
    #[arg(short = 'X', long = "customized-index")]
    customized_index: bool,

    /// Output all positions (twice: including unused reference sequences).
    #[arg(short = 'a', action = clap::ArgAction::Count)]
    all: u8,

    /// Write output to this file.
    #[arg(short = 'o', long = "output")]
    output: Option<PathBuf>,

    /// Output base positions on reads.
    #[arg(short = 'O', long = "output-BP", visible_alias = "output-bp")]
    output_bp: bool,

    /// Output base positions on reads, 5' to 3'.
    #[arg(long = "output-BP-5", visible_alias = "output-bp-5")]
    output_bp5: bool,

    /// Output mapping quality.
    #[arg(short = 's', long = "output-MQ", visible_alias = "output-mq")]
    output_mq: bool,

    /// Output read names.
    #[arg(long = "output-QNAME", visible_alias = "output-qname")]
    output_qname: bool,

    /// Output extra columns: QNAME, FLAG, RNAME, POS, MAPQ, RNEXT, PNEXT, RLEN, or tags.
    #[arg(long = "output-extra")]
    output_extra: Option<String>,

    /// Separator for --output-extra tag values.
    #[arg(long = "output-sep", default_value = ",")]
    output_sep: String,

    /// Placeholder for absent --output-extra tags.
    #[arg(long = "output-empty", default_value = "*")]
    output_empty: String,

    /// Mark deletions on the reverse strand with `#`.
    #[arg(long = "reverse-del")]
    reverse_del: bool,

    /// Do not output inserted bases (twice: nor their lengths).
    #[arg(long = "no-output-ins", action = clap::ArgAction::Count)]
    no_output_ins: u8,

    /// Do not output base modifications in insertions.
    #[arg(long = "no-output-ins-mods")]
    no_output_ins_mods: bool,

    /// Do not output deleted bases (twice: nor their lengths).
    #[arg(long = "no-output-del", action = clap::ArgAction::Count)]
    no_output_del: u8,

    /// Do not output read start and end marks.
    #[arg(long = "no-output-ends")]
    no_output_ends: bool,

    /// Quality is in the Illumina 1.3+ encoding.
    #[arg(short = '6', long = "illumina1.3+")]
    illumina13: bool,

    /// Output base modifications (not supported).
    #[arg(short = 'M', long = "output-mods")]
    output_mods: bool,

    /// Reference FASTA (indexed) for CRAM input [default: -f].
    #[arg(long = "reference")]
    reference: Option<PathBuf>,

    /// Additional threads for decompression and BAQ.
    #[arg(short = '@', long = "threads", default_value_t = 0)]
    threads: usize,

    /// Print help.
    #[arg(long, action = clap::ArgAction::Help)]
    help: Option<bool>,
}

/// `mplp_conf_t`, as the reader threads need it.
struct Conf {
    min_mq: i32,
    capq_thres: i32,
    rflag_require: u16,
    rflag_filter: u16,
    realn: bool,
    redo_baq: bool,
    no_orphan: bool,
    illumina13: bool,
    /// `-l`, applied to reads only without `-a`.
    bed: Option<Arc<Bed>>,
    all: u8,
    exclude_rg: Option<HashSet<Vec<u8>>>,
}

/// `-l` intervals by sequence name: 0-based, half-open, sorted by start (`bed_read`), with the
/// running maximum of their ends.
#[derive(Debug, Default)]
struct Bed {
    by_name: std::collections::HashMap<Vec<u8>, Intervals>,
}

/// Intervals sorted by start, and the running maximum of their ends.
type Intervals = (Vec<(i64, i64)>, Vec<i64>);

/// A cached reference sequence by name; `None` if the FASTA does not have it.
type CachedSequence = (Vec<u8>, Option<Arc<Vec<u8>>>);

impl Bed {
    fn read(path: &Path) -> anyhow::Result<Self> {
        let text = fs::read(path)
            .with_context(|| format!("Could not read file \"{}\"", path.display()))?;
        let mut bed = Self::default();

        for (n, line) in text.split(|&b| b == b'\n').enumerate() {
            let line = line.trim_ascii_start();

            if line.is_empty() || line[0] == b'#' {
                continue;
            }

            let name_end = line
                .iter()
                .position(u8::is_ascii_whitespace)
                .unwrap_or(line.len());
            let name = &line[..name_end];
            let rest = String::from_utf8_lossy(&line[name_end..]);
            let mut fields = rest.split_ascii_whitespace().map(str::parse::<u64>);

            let (num, mut beg, mut end) = match (fields.next(), fields.next()) {
                (Some(Ok(beg)), Some(Ok(end))) => (2, beg, end),
                (Some(Ok(beg)), _) => (1, beg, 0),
                _ => (0, 0, 0),
            };

            if num == 1 {
                // A position, counting from 1.
                end = beg;
                beg = beg.wrapping_sub(1);
            }

            if num < 1 || end < beg {
                if name == b"browser" || name == b"track" {
                    continue;
                }

                bail!(
                    "[bed_read] Parse error reading \"{}\" at line {}",
                    path.display(),
                    n + 1
                );
            }

            bed.by_name
                .entry(name.to_vec())
                .or_default()
                .0
                .push((beg as i64, end as i64));
        }

        bed.index();
        Ok(bed)
    }

    /// Sorts the intervals and computes the running maximum of their ends.
    fn index(&mut self) {
        for (intervals, max_ends) in self.by_name.values_mut() {
            intervals.sort_unstable();

            *max_ends = intervals
                .iter()
                .scan(i64::MIN, |max, &(_, end)| {
                    *max = (*max).max(end);
                    Some(*max)
                })
                .collect();
        }
    }

    /// `bed_overlap`: whether `[beg, end)` overlaps an interval on `name`.
    fn overlaps(&self, name: &[u8], beg: i64, end: i64) -> bool {
        let Some((intervals, max_ends)) = self.by_name.get(name) else {
            return false;
        };

        // Of the intervals starting before `end`, whether any ends after `beg`.
        let n = intervals.partition_point(|&(b, _)| b < end);
        n > 0 && max_ends[n - 1] > beg
    }
}

/// Whole reference sequences by name, the last few used (`mplp_get_ref` keeps three per input;
/// these are shared by all inputs and threads).
struct ReferenceCache {
    fasta: Option<IndexedFasta>,
    cache: Mutex<Vec<CachedSequence>>,
}

impl ReferenceCache {
    const SLOTS: usize = 4;

    fn get(&self, name: &[u8]) -> anyhow::Result<Option<Arc<Vec<u8>>>> {
        let Some(fasta) = &self.fasta else {
            return Ok(None);
        };

        let mut cache = self.cache.lock().expect("lock poisoned");

        if let Some(i) = cache.iter().position(|(n, _)| n == name) {
            let entry = cache.remove(i);
            let sequence = entry.1.clone();
            cache.insert(0, entry);
            return Ok(sequence);
        }

        let mut buf = Vec::new();
        let sequence = fasta
            .fetch(name, 0, i64::MAX, &mut buf)?
            .then(|| Arc::new(buf));

        cache.insert(0, (name.to_vec(), sequence.clone()));
        cache.truncate(Self::SLOTS);

        Ok(sequence)
    }

    /// The sequence of a reference sequence ID of `names`, if any.
    fn get_id(&self, names: &[Vec<u8>], tid: i32) -> anyhow::Result<Option<Arc<Vec<u8>>>> {
        match usize::try_from(tid).ok().and_then(|i| names.get(i)) {
            Some(name) => self.get(name),
            None => Ok(None),
        }
    }
}

/// `mplp_func`: filters a read and adjusts its qualities, or returns `None` to skip it.
fn prepare(
    raw: Vec<u8>,
    conf: &Conf,
    names: &[Vec<u8>],
    references: &ReferenceCache,
) -> anyhow::Result<Option<Read>> {
    let mut read = Read::new(raw)?;

    if read.reference_sequence_id() < 0 || read.flags() & 0x4 != 0 {
        return Ok(None);
    }

    if conf.rflag_require != 0 && conf.rflag_require & read.flags() == 0 {
        return Ok(None);
    }

    if conf.rflag_filter != 0 && conf.rflag_filter & read.flags() != 0 {
        return Ok(None);
    }

    if let Some(bed) = &conf.bed
        && conf.all == 0
    {
        let name = names
            .get(read.reference_sequence_id() as usize)
            .map_or(&[][..], Vec::as_slice);

        if !bed.overlaps(name, read.position(), end_pos(&read)) {
            return Ok(None);
        }
    }

    if let Some(exclude) = &conf.exclude_rg
        && let Some((_, value)) = read.tag(*b"RG")
    {
        let value = &value[..memchr::memchr(0, value).unwrap_or(value.len())];

        if exclude.contains(value) {
            return Ok(None);
        }
    }

    if conf.illumina13 {
        for q in read.quality_scores_mut() {
            *q = q.saturating_sub(31);
        }
    }

    let reference = references.get_id(names, read.reference_sequence_id())?;

    if let Some(reference) = &reference {
        if (reference.len() as i64) <= read.position() {
            eprintln!(
                "[mplp_func] Skipping because {} is outside of {} [ref:{}]",
                read.position(),
                reference.len(),
                read.reference_sequence_id()
            );
            return Ok(None);
        }

        if conf.realn {
            let mode = baq::BaqMode {
                apply: true,
                extend: true,
                redo: conf.redo_baq,
            };
            let _ = baq::prob_realn(&mut read, reference, mode);
        }

        if conf.capq_thres > 10 {
            match baq::cap_mapq(&read, reference, conf.capq_thres) {
                None => return Ok(None),
                Some(q) if q < 0 => return Ok(None),
                Some(q) => {
                    if i32::from(read.mapping_quality()) > q {
                        read.set_mapping_quality(q as u8);
                    }
                }
            }
        }
    }

    if i32::from(read.mapping_quality()) < conf.min_mq {
        return Ok(None);
    }

    if conf.no_orphan && read.flags() & 0x1 != 0 && read.flags() & 0x2 == 0 {
        return Ok(None);
    }

    Ok(Some(read))
}

/// htslib's `bam_endpos`.
fn end_pos(read: &Read) -> i64 {
    let span: i64 = read
        .cigar()
        .iter()
        .filter(|&&op| matches!(op & 0xf, 0 | 2 | 3 | 7 | 8))
        .map(|&op| i64::from(op >> 4))
        .sum();

    read.position() + span.max(1)
}

/// Prepares a batch of records, on `workers` threads when there are several, keeping their
/// order.
fn prepare_batch(
    batch: Vec<Vec<u8>>,
    conf: &Conf,
    names: &[Vec<u8>],
    references: &ReferenceCache,
    workers: usize,
) -> anyhow::Result<Vec<Read>> {
    if workers <= 1 || batch.len() < 64 {
        let mut reads = Vec::with_capacity(batch.len());

        for raw in batch {
            if let Some(read) = prepare(raw, conf, names, references)? {
                reads.push(read);
            }
        }

        return Ok(reads);
    }

    let chunk_size = batch.len().div_ceil(workers);
    let mut chunks: Vec<Vec<Vec<u8>>> = Vec::with_capacity(workers);
    let mut batch = batch.into_iter();

    loop {
        let chunk: Vec<Vec<u8>> = batch.by_ref().take(chunk_size).collect();

        if chunk.is_empty() {
            break;
        }

        chunks.push(chunk);
    }

    let results: Vec<anyhow::Result<Vec<Read>>> = thread::scope(|scope| {
        let handles: Vec<_> = chunks
            .into_iter()
            .map(|chunk| {
                scope.spawn(move || {
                    let mut reads = Vec::with_capacity(chunk.len());

                    for raw in chunk {
                        if let Some(read) = prepare(raw, conf, names, references)? {
                            reads.push(read);
                        }
                    }

                    Ok(reads)
                })
            })
            .collect();

        handles
            .into_iter()
            .map(|h| h.join().expect("worker panicked"))
            .collect()
    });

    let mut reads = Vec::new();

    for result in results {
        reads.extend(result?);
    }

    Ok(reads)
}

/// A pileup over one input, fed from its reader thread.
struct Source {
    pileup: Pileup,
    batches: crossbeam_channel::Receiver<anyhow::Result<Vec<Read>>>,
    pending: std::vec::IntoIter<Read>,
}

impl Source {
    /// The next prepared read, or `None` at the end of the input.
    fn next_read(&mut self) -> anyhow::Result<Option<Read>> {
        loop {
            if let Some(read) = self.pending.next() {
                return Ok(Some(read));
            }

            match self.batches.recv() {
                Ok(batch) => self.pending = batch?.into_iter(),
                Err(_) => return Ok(None),
            }
        }
    }

    /// `bam_plp64_auto`: returns whether a position is ready.
    fn auto(&mut self) -> anyhow::Result<bool> {
        if self.pileup.next_position()? {
            return Ok(true);
        }

        if self.pileup.is_eof() {
            return Ok(false);
        }

        loop {
            match self.next_read()? {
                Some(read) => {
                    self.pileup.push(Some(read))?;

                    if self.pileup.next_position()? {
                        return Ok(true);
                    }
                }
                None => {
                    self.pileup.push(None)?;
                    return Ok(self.pileup.next_position()?);
                }
            }
        }
    }
}

pub fn run(args: Args) -> anyhow::Result<()> {
    if args.no_baq && args.redo_baq {
        bail!("Error: The -B option cannot be combined with -E");
    }

    let (inputs, indexes) = input_files(&args)?;

    if inputs.is_empty() {
        bail!("[mpileup] no input file/data given");
    }

    let mut print_flags = 0;

    if args.output_mq {
        print_flags |= PRINT_MAPQ_CHAR;
    }

    if args.output_bp {
        print_flags |= PRINT_QPOS;
    }

    if args.output_bp5 {
        print_flags |= PRINT_QPOS5;
    }

    if args.output_qname {
        print_flags |= PRINT_QNAME;
    }

    let mut aux_tags: Vec<[u8; 2]> = Vec::new();

    if let Some(extra) = &args.output_extra {
        for name in extra.split(',').filter(|s| !s.is_empty()) {
            let flag = match name {
                "QNAME" => Some(PRINT_QNAME),
                "FLAG" => Some(PRINT_FLAG),
                "RNAME" => Some(PRINT_RNAME),
                "POS" => Some(PRINT_POS),
                "MAPQ" => Some(PRINT_MAPQ),
                "RNEXT" => Some(PRINT_RNEXT),
                "PNEXT" => Some(PRINT_PNEXT),
                "RLEN" => Some(PRINT_RLEN),
                _ => None,
            };

            match flag {
                Some(flag) => print_flags |= flag,
                None if name.len() == 2 => {
                    let b = name.as_bytes();
                    aux_tags.push([b[0], b[1]]);
                }
                None => eprintln!(
                    "[build_auxlist] tag '{name}' has more than two characters or not supported"
                ),
            }
        }
    }

    let bed = args
        .positions
        .as_deref()
        .map(Bed::read)
        .transpose()?
        .map(Arc::new);

    let exclude_rg = match &args.exclude_rg {
        Some(path) => match fs::read(path) {
            Ok(text) => Some(
                text.split(|b| b.is_ascii_whitespace())
                    .filter(|s| !s.is_empty())
                    .map(<[u8]>::to_vec)
                    .collect(),
            ),
            Err(_) => {
                eprintln!(
                    "[bam_mpileup] Fail to open file {}. Continue anyway.",
                    path.display()
                );
                Some(HashSet::new())
            }
        },
        None => None,
    };

    let fasta = args
        .fasta_ref
        .as_deref()
        .map(IndexedFasta::open)
        .transpose()?;
    let cram_reference = args.reference.clone().or_else(|| args.fasta_ref.clone());

    let read_options = ReadOptions {
        threads: args.threads,
        reference: cram_reference,
    };

    let max_depth = if args.max_depth <= 0 {
        eprintln!("[mpileup] Max depth set to maximum value ({})", i32::MAX);
        i32::MAX as usize
    } else {
        if args.max_depth * inputs.len() as i64 > 1 << 20 {
            eprintln!("[mpileup] Combined max depth is above 1M. Potential memory hog!");
        }

        args.max_depth as usize
    };

    let conf = Conf {
        min_mq: args.min_mq,
        capq_thres: args.adjust_mq,
        rflag_require: args.required_flags.map_or(0, |f| f.bits()),
        rflag_filter: args
            .filter_flags
            .map_or(0x4 | 0x100 | 0x200 | 0x400, |f| f.bits()),
        realn: !args.no_baq,
        redo_baq: args.redo_baq,
        no_orphan: !args.count_orphans,
        illumina13: args.illumina13,
        bed: bed.clone(),
        all: args.all,
        exclude_rg,
    };

    let references = ReferenceCache {
        fasta,
        cache: Mutex::new(Vec::new()),
    };

    let sink: Box<dyn Write> = match &args.output {
        Some(path) => Box::new(
            File::create(path).with_context(|| format!("failed to write to {}", path.display()))?,
        ),
        None => Box::new(output::raw_stdout()?),
    };
    let mut out = BufWriter::with_capacity(1 << 20, sink);

    let options = OutputOptions {
        min_bq: args.min_bq,
        rev_del: args.reverse_del,
        no_ins: args.no_output_ins,
        mods: if args.output_mods {
            if args.no_output_ins_mods || args.no_output_ins > 0 {
                Mods::Bases
            } else {
                Mods::BasesAndInsertions
            }
        } else {
            Mods::None
        },
        no_del: args.no_output_del,
        no_ends: args.no_output_ends,
        flags: print_flags,
        aux_tags,
        sep: args.output_sep.bytes().next().unwrap_or(b','),
        empty: args.output_empty.bytes().next().unwrap_or(b'*'),
        all: args.all,
    };

    let context = PileupContext {
        inputs: &inputs,
        indexes: indexes.as_deref(),
        read_options: &read_options,
        conf: &conf,
        references: &references,
        options: &options,
        bed: bed.as_deref(),
        max_depth,
        overlaps: if args.no_overlaps {
            Overlaps::Ignore
        } else {
            Overlaps::Adjust
        },
        workers: args.threads,
    };

    match whole_chromosomes(&args, &inputs, indexes.as_deref(), &read_options)? {
        Some((headers, names)) => {
            report_samples(&headers, args.ignore_rg);
            pileup_by_chromosome(&context, &names, args.threads, &mut out)?;
        }
        None => pileup_region(
            &context,
            args.region.as_deref(),
            Some(args.ignore_rg),
            &mut out,
        )?,
    }

    out.flush()?;

    Ok(())
}

/// What the pileup of one region needs.
struct PileupContext<'a> {
    inputs: &'a [String],
    indexes: Option<&'a [String]>,
    read_options: &'a ReadOptions,
    conf: &'a Conf,
    references: &'a ReferenceCache,
    options: &'a OutputOptions,
    bed: Option<&'a Bed>,
    max_depth: usize,
    overlaps: Overlaps,
    /// Threads for BAQ in each reader.
    workers: usize,
}

/// Runs the pileup of `region` (a samtools region, or every read) and writes it to `out`.
/// `report` reports the number of samples, ignoring read groups if set.
fn pileup_region<W: Write>(
    context: &PileupContext<'_>,
    region: Option<&str>,
    report: Option<bool>,
    out: &mut W,
) -> anyhow::Result<()> {
    thread::scope(|scope| -> anyhow::Result<()> {
        let mut sources = Vec::with_capacity(context.inputs.len());
        let mut header_rxs = Vec::with_capacity(context.inputs.len());

        for (i, input) in context.inputs.iter().enumerate() {
            let (tx, rx) = crossbeam_channel::bounded(QUEUE_DEPTH);
            let (header_tx, header_rx) = crossbeam_channel::bounded(1);
            let index = context.indexes.map(|indexes| PathBuf::from(&indexes[i]));

            scope.spawn(move || {
                let result = read_input(
                    input,
                    context.read_options,
                    region,
                    index.as_deref(),
                    context.conf,
                    context.references,
                    context.workers,
                    &header_tx,
                    &tx,
                );

                if let Err(e) = result {
                    let _ = header_tx.try_send(None);
                    let _ = tx.send(Err(e));
                }
            });

            header_rxs.push(header_rx);
            let mut pileup = Pileup::new(context.max_depth, context.overlaps);

            if context.options.mods != Mods::None {
                pileup.enable_base_modifications();
            }

            sources.push(Source {
                pileup,
                batches: rx,
                pending: Vec::new().into_iter(),
            });
        }

        // Names and lengths come from the first file's header, as in samtools.
        let mut headers = Vec::with_capacity(context.inputs.len());

        for (rx, source) in header_rxs.iter().zip(&mut sources) {
            match rx.recv().ok().flatten() {
                Some(header) => headers.push(header),
                None => {
                    // The reader failed; report its error.
                    source.next_read()?;
                    bail!("failed to read a header");
                }
            }
        }

        let header: &sam::Header = &headers[0];
        let names: Vec<Vec<u8>> = header
            .reference_sequences()
            .keys()
            .map(|name| name.to_vec())
            .collect();
        let lengths: Vec<i64> = header
            .reference_sequences()
            .values()
            .map(|rs| rs.length().get() as i64)
            .collect();

        if let Some(ignore_rg) = report {
            report_samples(&headers, ignore_rg);
        }

        // The requested region, 0-based and half-open, of the first file.
        let region = match region {
            Some(raw) => Some(parse_region(raw, &names)?),
            None => None,
        };

        let result = mpileup(
            &mut sources,
            context.options,
            &names,
            &lengths,
            context.references,
            context.bed,
            region,
            out,
        );

        // Stop the readers if the pileup ended early.
        drop(sources);
        result
    })
}

/// The output a chromosome keeps in memory before spilling the rest to a temporary file.
const SPILL_THRESHOLD: usize = 64 << 20;

/// The inputs' headers and the reference sequence names, for [`pileup_by_chromosome`].
type Chromosomes = (Vec<sam::Header>, Vec<Vec<u8>>);

/// Returns the headers and reference sequence names if the pileup can run one chromosome at a
/// time on several threads: every input is indexed, no region is given, and there are several
/// sequences, all with names that are valid regions.
///
/// Chromosomes are independent in the pileup (reads, overlapping mates, and the depth limit
/// never span two), and a whole-chromosome region gives exactly that chromosome's part of the
/// output, also with `-a`/`-aa`. So the output is the same as reading the inputs through.
fn whole_chromosomes(
    args: &Args,
    inputs: &[String],
    indexes: Option<&[String]>,
    read_options: &ReadOptions,
) -> anyhow::Result<Option<Chromosomes>> {
    use crate::alignment::AlignmentIndex;

    if args.threads == 0 || args.region.is_some() {
        return Ok(None);
    }

    let mut headers = Vec::with_capacity(inputs.len());

    for (i, input) in inputs.iter().enumerate() {
        if input == "-" {
            return Ok(None);
        }

        let (mut reader, path) = AlignmentReader::open(input, read_options)?;
        let header = reader.read_header()?;

        let Some(path) = path else {
            return Ok(None);
        };

        let explicit = indexes.map(|indexes| Path::new(&indexes[i]));

        if AlignmentIndex::read(reader.format(), &path, explicit).is_err() {
            return Ok(None);
        }

        headers.push(header);
    }

    let names: Vec<Vec<u8>> = headers[0]
        .reference_sequences()
        .keys()
        .map(|name| name.to_vec())
        .collect();

    if names.len() < 2
        || names
            .iter()
            .any(|name| name.contains(&b':') || name.is_empty())
    {
        return Ok(None);
    }

    Ok(Some((headers, names)))
}

/// Runs the pileup of each chromosome as a whole-chromosome region, `jobs` at a time, and writes
/// their output in order. Chromosomes ahead of the one being written keep their output in
/// memory, then in a temporary file.
fn pileup_by_chromosome<W: Write>(
    context: &PileupContext<'_>,
    names: &[Vec<u8>],
    jobs: usize,
    out: &mut W,
) -> anyhow::Result<()> {
    let regions: Vec<String> = names
        .iter()
        .map(|name| String::from_utf8_lossy(name).into_owned())
        .collect();

    let jobs = jobs.max(1);

    // With fewer chromosomes than threads, each one realigns on several.
    let active = jobs.min(regions.len());
    let workers = if jobs / active > 1 { jobs / active } else { 0 };

    let read_options = ReadOptions {
        threads: 0,
        reference: context.read_options.reference.clone(),
    };

    let job_context = PileupContext {
        read_options: &read_options,
        workers,
        ..*context
    };
    let job_context = &job_context;

    thread::scope(|scope| -> anyhow::Result<()> {
        let mut pending = std::collections::VecDeque::new();
        let mut spawned = 0;

        for i in 0..regions.len() {
            while spawned < regions.len() && spawned < i + jobs {
                let region = &regions[spawned];

                pending.push_back(scope.spawn(move || -> anyhow::Result<SpillWriter> {
                    let mut writer = SpillWriter::default();
                    pileup_region(job_context, Some(region), None, &mut writer)?;
                    Ok(writer)
                }));

                spawned += 1;
            }

            let handle = pending.pop_front().expect("spawned ahead");
            let output = handle.join().expect("pileup thread panicked")?;
            output.copy_to(out)?;
        }

        Ok(())
    })
}

/// Output kept in memory up to [`SPILL_THRESHOLD`], then in a temporary file.
#[derive(Default)]
struct SpillWriter {
    memory: Vec<u8>,
    file: Option<BufWriter<File>>,
}

impl SpillWriter {
    /// Writes everything to `out`, in order.
    fn copy_to<W: Write>(self, out: &mut W) -> anyhow::Result<()> {
        use std::io::Seek as _;

        out.write_all(&self.memory)?;

        if let Some(file) = self.file {
            let mut file = file.into_inner().map_err(io::IntoInnerError::into_error)?;
            file.rewind()?;
            io::copy(&mut file, out)?;
        }

        Ok(())
    }
}

impl Write for SpillWriter {
    fn write(&mut self, buf: &[u8]) -> io::Result<usize> {
        if let Some(file) = &mut self.file {
            return file.write(buf);
        }

        if self.memory.len() + buf.len() > SPILL_THRESHOLD {
            let mut file = BufWriter::with_capacity(1 << 20, tempfile::tempfile()?);
            file.write_all(buf)?;
            self.file = Some(file);
        } else {
            self.memory.extend_from_slice(buf);
        }

        Ok(buf.len())
    }

    fn flush(&mut self) -> io::Result<()> {
        match &mut self.file {
            Some(file) => file.flush(),
            None => Ok(()),
        }
    }
}

/// Reports the number of samples, as samtools does (read group samples, or one per file).
fn report_samples(headers: &[sam::Header], ignore_rg: bool) {
    use sam::header::record::value::map::read_group::tag::SAMPLE;

    let sample_count: usize = headers
        .iter()
        .map(|h| {
            let samples: HashSet<_> = h
                .read_groups()
                .values()
                .filter_map(|rg| rg.other_fields().get(&SAMPLE))
                .collect();

            if ignore_rg || samples.is_empty() {
                1
            } else {
                samples.len()
            }
        })
        .sum();

    eprintln!(
        "[mpileup] {sample_count} samples in {} input files",
        headers.len()
    );
}

/// Parses `-r` as `sam_itr_querys` does: 0-based start and exclusive end.
fn parse_region(raw: &str, names: &[Vec<u8>]) -> anyhow::Result<(i32, i64, i64)> {
    let parsed: Region = raw
        .parse()
        .with_context(|| format!("[E::mpileup] fail to parse region '{raw}'"))?;
    let name: &[u8] = parsed.name().as_ref();
    let tid = names
        .iter()
        .position(|n| n.as_slice() == name)
        .with_context(|| format!("[E::mpileup] fail to parse region '{raw}'"))?;
    let interval = parsed.interval();
    let beg = interval.start().map_or(0, |p| usize::from(p) as i64 - 1);
    let end = interval.end().map_or(i64::MAX, |p| usize::from(p) as i64);
    Ok((tid as i32, beg, end))
}

/// The input files and their index files (`-b`, `-X`).
fn input_files(args: &Args) -> anyhow::Result<(Vec<String>, Option<Vec<String>>)> {
    if let Some(list) = &args.bam_list {
        if args.customized_index {
            bail!("Error: The -b option cannot be combined with -X");
        }

        let text = fs::read_to_string(list).with_context(|| format!("{}", list.display()))?;
        let files: Vec<String> = text
            .lines()
            .map(str::trim_end)
            .filter(|l| !l.is_empty())
            .map(String::from)
            .collect();

        for file in &files {
            if !file.contains(':') && !Path::new(file).exists() {
                bail!(
                    "The file list \"{}\" appears broken, could not locate: {file}",
                    list.display()
                );
            }
        }

        if files.is_empty() {
            bail!("No files read from {}", list.display());
        }

        return Ok((files, None));
    }

    if args.customized_index {
        if !args.inputs.len().is_multiple_of(2) {
            bail!("Odd number of filenames detected! Each BAM file should have an index file");
        }

        let n = args.inputs.len() / 2;
        return Ok((args.inputs[..n].to_vec(), Some(args.inputs[n..].to_vec())));
    }

    Ok((args.inputs.clone(), None))
}

/// Reads one input, sending batches of prepared reads.
#[expect(clippy::too_many_arguments, reason = "a thread's context")]
fn read_input(
    input: &str,
    read_options: &ReadOptions,
    region: Option<&str>,
    index: Option<&Path>,
    conf: &Conf,
    references: &ReferenceCache,
    workers: usize,
    header_tx: &crossbeam_channel::Sender<Option<sam::Header>>,
    tx: &crossbeam_channel::Sender<anyhow::Result<Vec<Read>>>,
) -> anyhow::Result<()> {
    let (mut reader, path) = AlignmentReader::open(input, read_options)?;
    let header = reader.read_header()?;
    let names: Vec<Vec<u8>> = header
        .reference_sequences()
        .keys()
        .map(|name| name.to_vec())
        .collect();

    let _ = header_tx.send(Some(header.clone()));
    let mut batch: Vec<Vec<u8>> = Vec::with_capacity(BATCH_SIZE);
    let mut encoder = noodles_bam::io::Writer::from(Vec::new());

    // An `Err` with no cause signals that the consumer hung up.
    let send = |batch: &mut Vec<Vec<u8>>| -> anyhow::Result<()> {
        let reads = prepare_batch(std::mem::take(batch), conf, &names, references, workers)?;

        if tx.send(Ok(reads)).is_err() {
            bail!(Disconnected);
        }

        Ok(())
    };

    let mut visit = |rec: Rec<'_>| -> anyhow::Result<()> {
        use noodles_sam::alignment::io::Write as _;

        let raw = match &rec {
            Rec::Bam(r) => r.as_ref().to_vec(),
            _ => {
                encoder.get_mut().clear();
                encoder.write_alignment_record(&header, rec.as_dyn())?;
                encoder.get_ref()[size_of::<u32>()..].to_vec()
            }
        };

        batch.push(raw);

        if batch.len() >= BATCH_SIZE {
            send(&mut batch)?;
        }

        Ok(())
    };

    let result = match region {
        Some(region) => {
            let path = path.with_context(|| format!("[mpileup] fail to load index for {input}"))?;
            let regions = [region.to_owned()];
            reader.for_each_in_queries(
                &path,
                index,
                &header,
                &Queries::Regions(&regions),
                &mut visit,
            )
        }
        None => reader.for_each(&header, &mut visit),
    };

    match result {
        Err(e) if e.is::<Disconnected>() => return Ok(()),
        result => result?,
    }

    if !batch.is_empty() {
        match send(&mut batch) {
            Err(e) if e.is::<Disconnected>() => {}
            result => result?,
        }
    }

    Ok(())
}

#[derive(Debug)]
struct Disconnected;

impl std::fmt::Display for Disconnected {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.write_str("the pileup stopped reading")
    }
}

impl std::error::Error for Disconnected {}

/// Which base modifications (MM and ML tags) to print, after their bases.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
enum Mods {
    None,
    Bases,
    BasesAndInsertions,
}

struct OutputOptions {
    min_bq: i32,
    rev_del: bool,
    no_ins: u8,
    mods: Mods,
    no_del: u8,
    no_ends: bool,
    flags: u32,
    aux_tags: Vec<[u8; 2]>,
    sep: u8,
    empty: u8,
    all: u8,
}

/// `mpileup`: the pileup loop and its output.
#[expect(clippy::too_many_arguments, reason = "the pileup's context")]
fn mpileup<W: Write>(
    sources: &mut [Source],
    options: &OutputOptions,
    names: &[Vec<u8>],
    lengths: &[i64],
    references: &ReferenceCache,
    bed: Option<&Bed>,
    region: Option<(i32, i64, i64)>,
    out: &mut W,
) -> anyhow::Result<()> {
    let n = sources.len();
    let (tid0, beg0, end0) = region.unwrap_or((0, 0, i64::MAX));
    let has_region = region.is_some();

    // `bam_mplp_t` state
    let mut min_tid: i32 = -1;
    let mut min_pos: i64 = i64::MAX;
    let mut tids = vec![-1i32; n];
    let mut poss = vec![i64::MAX; n];
    let mut ready = vec![false; n];
    let mut active = vec![false; n];

    let mut last_tid: i32 = -1;
    let mut last_pos: i64 = -1;
    let mut got_ref = false;
    let mut one_seq = false;
    let mut reference: Option<Arc<Vec<u8>>> = None;

    let mut line = Vec::with_capacity(1 << 16);
    let mut seq = Vec::new();
    let mut qual = Vec::new();
    let mut ins = Vec::new();

    let ref_at = |reference: &Option<Arc<Vec<u8>>>, pos: i64| -> u8 {
        match reference {
            Some(r) if pos < r.len() as i64 => r[pos as usize],
            _ => b'N',
        }
    };

    let name_of = |tid: i32| -> &[u8] { &names[tid as usize] };

    loop {
        // `bam_mplp64_auto`
        let mut new_min_tid = u32::MAX;
        let mut new_min_pos = i64::MAX;

        for i in 0..n {
            if poss[i] == min_pos && tids[i] == min_tid {
                ready[i] = sources[i].auto()?;

                if ready[i] {
                    let (tid, pos) = sources[i].pileup.position();
                    tids[i] = tid;
                    poss[i] = pos;
                } else {
                    tids[i] = 0;
                    poss[i] = 0;
                }
            }

            if ready[i] {
                let tid = tids[i] as u32;

                if tid < new_min_tid {
                    new_min_tid = tid;
                    new_min_pos = poss[i];
                } else if tid == new_min_tid && poss[i] < new_min_pos {
                    new_min_pos = poss[i];
                }
            }
        }

        min_pos = new_min_pos;
        min_tid = new_min_tid as i32;

        if new_min_pos == i64::MAX {
            break;
        }

        let (tid, pos) = (min_tid, min_pos);

        for i in 0..n {
            active[i] = ready[i] && poss[i] == min_pos && tids[i] == min_tid;
        }

        one_seq = true;

        if has_region && (pos < beg0 || pos >= end0) {
            continue;
        }

        if options.all > 0 {
            // Positions of earlier sequences not yet written
            while tid > last_tid {
                if last_tid >= 0 && !has_region {
                    loop {
                        last_pos += 1;

                        if last_pos >= lengths[last_tid as usize] {
                            break;
                        }

                        if let Some(bed) = bed
                            && !bed.overlaps(name_of(last_tid), last_pos, last_pos + 1)
                        {
                            continue;
                        }

                        write_empty(
                            out,
                            &mut line,
                            options,
                            name_of(last_tid),
                            last_pos,
                            n,
                            &reference,
                        )?;
                    }
                }

                last_tid += 1;
                got_ref = false;
                last_pos = -1;

                if options.all < 2 {
                    break;
                }

                if tid > last_tid {
                    reference = references.get_id(names, last_tid)?;
                    got_ref = reference.is_some();
                }
            }
        }

        if !got_ref || last_tid != tid {
            reference = references.get_id(names, tid)?;
            got_ref = reference.is_some();
            last_tid = tid;
        }

        if options.all > 0 {
            // Positions of this sequence not yet written
            loop {
                last_pos += 1;

                if last_pos >= pos {
                    break;
                }

                if has_region && last_pos < beg0 {
                    continue;
                }

                if let Some(bed) = bed
                    && !bed.overlaps(name_of(tid), last_pos, last_pos + 1)
                {
                    continue;
                }

                write_empty(
                    out,
                    &mut line,
                    options,
                    name_of(tid),
                    last_pos,
                    n,
                    &reference,
                )?;
            }

            last_pos = pos;
        }

        if let Some(bed) = bed
            && tid >= 0
            && !bed.overlaps(name_of(tid), pos, pos + 1)
        {
            continue;
        }

        line.clear();
        line.extend_from_slice(name_of(tid));
        line.push(b'\t');
        push_int(&mut line, pos + 1);
        line.push(b'\t');
        line.push(ref_at(&reference, pos));

        for (i, source) in sources.iter_mut().enumerate() {
            write_sample(
                &mut line,
                &mut source.pileup,
                active[i],
                options,
                names,
                pos,
                &reference,
                &mut seq,
                &mut qual,
                &mut ins,
            )?;
        }

        line.push(b'\n');
        out.write_all(&line)?;
    }

    if options.all > 0 {
        // The end of the last sequence (and, with -aa, those after it)
        if last_tid < 0 && has_region && options.all > 1 {
            last_tid = tid0;
            last_pos = beg0 - 1;
        } else if last_tid < 0 && !one_seq && options.all > 1 {
            last_tid = 0;
        }

        while last_tid >= 0 && (last_tid as usize) < names.len() {
            reference = references.get_id(names, last_tid)?;

            loop {
                last_pos += 1;

                if last_pos >= lengths[last_tid as usize] || last_pos >= end0 {
                    break;
                }

                if let Some(bed) = bed
                    && !bed.overlaps(name_of(last_tid), last_pos, last_pos + 1)
                {
                    continue;
                }

                write_empty(
                    out,
                    &mut line,
                    options,
                    name_of(last_tid),
                    last_pos,
                    n,
                    &reference,
                )?;
            }

            last_tid += 1;
            last_pos = -1;

            if options.all < 2 || has_region {
                break;
            }
        }
    }

    Ok(())
}

/// `print_empty_pileup`
fn write_empty<W: Write>(
    out: &mut W,
    line: &mut Vec<u8>,
    options: &OutputOptions,
    name: &[u8],
    pos: i64,
    n: usize,
    reference: &Option<Arc<Vec<u8>>>,
) -> io::Result<()> {
    line.clear();
    line.extend_from_slice(name);
    line.push(b'\t');
    push_int(line, pos + 1);
    line.push(b'\t');
    line.push(match reference {
        Some(r) if pos < r.len() as i64 => r[pos as usize],
        _ => b'N',
    });

    for _ in 0..n {
        line.extend_from_slice(b"\t0\t*\t*");
        push_empty_columns(line, options);
    }

    line.push(b'\n');
    out.write_all(line)
}

fn push_empty_columns(line: &mut Vec<u8>, options: &OutputOptions) {
    if options.flags == 0 && options.aux_tags.is_empty() {
        return;
    }

    let mut flag = PRINT_MAPQ_CHAR;

    while flag < PRINT_LAST {
        if flag != PRINT_MODS && options.flags & flag != 0 {
            line.extend_from_slice(b"\t*");
        }

        flag <<= 1;
    }

    for _ in &options.aux_tags {
        line.extend_from_slice(b"\t*");
    }
}

/// The base quality samtools filters and prints: that of the base at `qpos`, or 0 past the end.
fn base_quality(read: &Read, entry: &Entry) -> i32 {
    if (entry.qpos as usize) < read.sequence_len() {
        i32::from(read.quality_scores()[entry.qpos as usize])
    } else {
        0
    }
}

/// The columns of one input at one position.
#[expect(clippy::too_many_arguments, reason = "the line's context")]
fn write_sample(
    line: &mut Vec<u8>,
    pileup: &mut Pileup,
    active: bool,
    options: &OutputOptions,
    names: &[Vec<u8>],
    pos: i64,
    reference: &Option<Arc<Vec<u8>>>,
    seq: &mut Vec<u8>,
    qual: &mut Vec<u8>,
    ins: &mut Vec<u8>,
) -> anyhow::Result<()> {
    seq.clear();
    qual.clear();
    let mut cnt = 0i64;
    let n_entries = if active { pileup.entries().len() } else { 0 };

    for k in 0..n_entries {
        let entry = pileup.entries()[k];
        let (read, mods) = pileup.read_and_base_mods(&entry);
        let c = base_quality(read, &entry);

        if c >= options.min_bq {
            push_base(seq, read, mods, &entry, pos, reference, options, ins)?;
            qual.push(if c + 33 < 126 { (c + 33) as u8 } else { 126 });
            cnt += 1;
        }
    }

    line.push(b'\t');
    push_int(line, cnt);
    line.push(b'\t');

    if n_entries == 0 {
        line.extend_from_slice(b"*\t*");
        push_empty_columns(line, options);
        return Ok(());
    }

    if seq.is_empty() {
        line.push(b'*');
    } else {
        line.extend_from_slice(seq);
    }

    line.push(b'\t');

    if qual.is_empty() {
        line.push(b'*');
    } else {
        line.extend_from_slice(qual);
    }

    let mut flag = if options.flags == 0 {
        PRINT_LAST
    } else {
        PRINT_MAPQ_CHAR
    };

    while flag < PRINT_LAST {
        if flag != PRINT_MODS && options.flags & flag != 0 {
            line.push(b'\t');
            let mut n = 0;

            for entry in &pileup.entries()[..n_entries] {
                let read = pileup.read(entry);

                if base_quality(read, entry) < options.min_bq {
                    continue;
                }

                if n > 0 && flag != PRINT_MAPQ_CHAR {
                    line.push(b',');
                }

                n += 1;

                match flag {
                    PRINT_MAPQ_CHAR => {
                        line.push((read.mapping_quality() as u32 + 33).min(126) as u8)
                    }
                    PRINT_QPOS => push_int(line, entry.qpos + 1),
                    PRINT_QPOS5 => {
                        let pos5 = if read.is_reverse_complemented() {
                            read.sequence_len() as i64 - entry.qpos + i64::from(entry.is_del)
                        } else {
                            entry.qpos + 1
                        };
                        push_int(line, pos5);
                    }
                    PRINT_QNAME => line.extend_from_slice(read.name()),
                    PRINT_FLAG => push_int(line, i64::from(read.flags())),
                    PRINT_RNAME => {
                        if read.reference_sequence_id() >= 0 {
                            line.extend_from_slice(&names[read.reference_sequence_id() as usize]);
                        } else {
                            line.push(b'*');
                        }
                    }
                    PRINT_POS => push_int(line, read.position() + 1),
                    PRINT_MAPQ => push_int(line, i64::from(read.mapping_quality())),
                    PRINT_RNEXT => {
                        if read.mate_reference_sequence_id() >= 0 {
                            line.extend_from_slice(
                                &names[read.mate_reference_sequence_id() as usize],
                            );
                        } else {
                            line.push(b'*');
                        }
                    }
                    PRINT_PNEXT => push_int(line, read.mate_position() + 1),
                    PRINT_RLEN => push_int(line, read.sequence_len() as i64),
                    _ => {}
                }
            }

            if n == 0 {
                line.push(b'*');
            }
        }

        flag <<= 1;
    }

    for &tag in &options.aux_tags {
        line.push(b'\t');
        let mut n = 0;

        for entry in &pileup.entries()[..n_entries] {
            let read = pileup.read(entry);

            if base_quality(read, entry) < options.min_bq {
                continue;
            }

            if n > 0 {
                line.push(options.sep);
            }

            n += 1;
            push_tag_value(line, read, tag, options.empty);
        }

        if n == 0 {
            line.push(b'*');
        }
    }

    Ok(())
}

/// An --output-extra tag value, as samtools prints it.
fn push_tag_value(line: &mut Vec<u8>, read: &Read, tag: [u8; 2], empty: u8) {
    let Some((ty, value)) = read.tag(tag) else {
        line.push(empty);
        return;
    };

    match ty {
        b'Z' | b'H' => line.extend_from_slice(&value[..value.len() - 1]),
        b'c' => push_int(line, i64::from(value[0].cast_signed())),
        b'C' => push_int(line, i64::from(value[0])),
        b's' => push_int(line, i64::from(i16::from_le_bytes([value[0], value[1]]))),
        b'S' => push_int(line, i64::from(u16::from_le_bytes([value[0], value[1]]))),
        b'i' => push_int(
            line,
            i64::from(i32::from_le_bytes([value[0], value[1], value[2], value[3]])),
        ),
        b'I' => push_int(
            line,
            i64::from(u32::from_le_bytes([value[0], value[1], value[2], value[3]])),
        ),
        b'f' => push_kputd(
            line,
            f64::from(f32::from_le_bytes([value[0], value[1], value[2], value[3]])),
        ),
        b'A' => line.push(value[0]),
        _ => line.push(b'*'),
    }
}

/// `pileup_seq`: one read's base column entry.
#[expect(clippy::too_many_arguments, reason = "the read's context")]
fn push_base(
    seq: &mut Vec<u8>,
    read: &Read,
    mut mods: Option<&mut BaseMods>,
    entry: &Entry,
    pos: i64,
    reference: &Option<Arc<Vec<u8>>>,
    options: &OutputOptions,
    ins: &mut Vec<u8>,
) -> anyhow::Result<()> {
    const LOWER: &[u8; 16] = b",acmgrsvtwyhkdbn";
    const UPPER: &[u8; 16] = b".ACMGRSVTWYHKDBN";

    let reverse = read.is_reverse_complemented();

    if !options.no_ends && entry.is_head {
        seq.push(b'^');
        let mapq = read.mapping_quality();
        seq.push(if mapq > 93 { 126 } else { mapq + 33 });
    }

    if entry.is_del {
        seq.push(if entry.is_refskip {
            if reverse { b'<' } else { b'>' }
        } else if reverse && options.rev_del {
            b'#'
        } else {
            b'*'
        });
    } else {
        let mut c = if (entry.qpos as usize) < read.sequence_len() {
            read.base(entry.qpos as usize)
        } else {
            15
        };

        if let Some(r) = reference {
            let rb = if pos < r.len() as i64 {
                SEQ_NT16_TABLE[usize::from(r[pos as usize])]
            } else {
                15
            };

            if c == rb {
                c = 0;
            }
        }

        seq.push(if reverse { LOWER } else { UPPER }[usize::from(c)]);

        if let Some(state) = mods.as_deref_mut() {
            let mut calls = Vec::new();
            state.at_query_position(read, entry.qpos as usize, &mut calls)?;
            pileup::mods::write_calls(seq, &calls);
        }
    }

    let mut del_len = if entry.indel < 0 {
        (-entry.indel) as usize
    } else {
        0
    };

    if entry.indel > 0 {
        let mods = mods.filter(|_| options.mods == Mods::BasesAndInsertions);
        let (len, after) = pileup::insertion_with_mods(read, entry, mods, ins)?;
        del_len = after;

        if options.no_ins < 2 {
            seq.push(b'+');
            push_int(seq, len as i64);
        }

        if options.no_ins == 0 {
            // Modifications (in brackets) keep their case.
            let mut in_mod = false;
            let pad = if reverse && options.rev_del {
                b'#'
            } else {
                b'*'
            };

            seq.extend(ins.iter().map(|&b| {
                if b == b'[' {
                    in_mod = true;
                } else if b == b']' {
                    in_mod = false;
                }

                if in_mod {
                    b
                } else if !reverse {
                    b.to_ascii_uppercase()
                } else if b == b'*' {
                    pad
                } else {
                    b.to_ascii_lowercase()
                }
            }));
        }
    }

    if del_len > 0 {
        if options.no_del < 2 {
            push_int(seq, -(del_len as i64));
        }

        if options.no_del == 0 {
            for j in 1..=del_len as i64 {
                let c = match reference {
                    Some(r) if ((pos as i32 as i64) + j) < r.len() as i64 => r[(pos + j) as usize],
                    _ => b'N',
                };

                seq.push(if reverse {
                    c.to_ascii_lowercase()
                } else {
                    c.to_ascii_uppercase()
                });
            }
        }
    }

    if !options.no_ends && entry.is_tail {
        seq.push(b'$');
    }

    Ok(())
}

fn push_int(dst: &mut Vec<u8>, n: i64) {
    let mut buf = [0; lexical_core::BUFFER_SIZE];
    dst.extend_from_slice(lexical_core::write(n, &mut buf));
}

/// htslib's `kputd`: up to six significant digits, trailing zeros dropped, and `%g` outside
/// [0.0001, 999999].
fn push_kputd(dst: &mut Vec<u8>, d: f64) {
    if d == 0.0 {
        dst.extend_from_slice(if d.is_sign_negative() { b"-0" } else { b"0" });
        return;
    }

    let mut d = d;

    if d < 0.0 {
        dst.push(b'-');
        d = -d;
    }

    if !(0.0001..=999_999.0).contains(&d) {
        dst.extend_from_slice(crate::depth::format_g(d, 6).as_bytes());
        return;
    }

    // Six significant digits as an integer, and how many of them follow the decimal point.
    let (scale, decimals) = if d < 0.001 {
        (1e9, 9)
    } else if d < 0.01 {
        (1e8, 8)
    } else if d < 0.1 {
        (1e7, 7)
    } else if d < 1.0 {
        (1e6, 6)
    } else if d < 10.0 {
        (1e5, 5)
    } else if d < 100.0 {
        (1e4, 4)
    } else if d < 1000.0 {
        (1e3, 3)
    } else if d < 10_000.0 {
        (1e2, 2)
    } else if d < 100_000.0 {
        (1e1, 1)
    } else {
        (1.0, 0)
    };

    #[expect(
        clippy::cast_possible_truncation,
        clippy::cast_sign_loss,
        reason = "kputd stores rint() of the scaled value in a uint32_t"
    )]
    let i = (d * scale).round_ties_even() as u32;

    // The fixed-point text of i / scale with trailing zeros culled, as kputd's digit shuffling
    // produces.
    let digits = format!("{i:0width$}", width = decimals + 1);
    let (int_part, frac_part) = digits.split_at(digits.len() - decimals);
    let frac = frac_part.trim_end_matches('0');

    dst.extend_from_slice(int_part.as_bytes());

    if !frac.is_empty() {
        dst.push(b'.');
        // kputd keeps at most six significant digits after the leading zeros.
        dst.extend_from_slice(frac.as_bytes());
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_push_kputd() {
        let t = |d: f64| {
            let mut buf = Vec::new();
            push_kputd(&mut buf, d);
            String::from_utf8(buf).expect("ASCII")
        };

        assert_eq!(t(0.0), "0");
        assert_eq!(t(-0.0), "-0");
        assert_eq!(t(1.5), "1.5");
        assert_eq!(t(-2.25), "-2.25");
        assert_eq!(t(100.0), "100");
        assert_eq!(t(0.000_123_456_7), "0.000123457");
        assert_eq!(t(0.999_999_9), "1");
        assert_eq!(t(123_456.7), "123457");
        assert_eq!(t(1e-5), "1e-05");
        assert_eq!(t(1e7), "1e+07");
    }

    #[test]
    fn test_bed_overlaps() {
        let mut bed = Bed::default();
        bed.by_name.insert(
            b"chr1".to_vec(),
            (vec![(200, 300), (10, 20), (0, 100_000)], Vec::new()),
        );
        bed.index();

        assert!(bed.overlaps(b"chr1", 15, 16));
        assert!(bed.overlaps(b"chr1", 99_999, 100_000));
        assert!(!bed.overlaps(b"chr1", 100_000, 100_001));
        assert!(!bed.overlaps(b"chr2", 15, 16));
    }
}
