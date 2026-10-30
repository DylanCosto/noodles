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
    io::{BufWriter, Write},
    path::{Path, PathBuf},
    sync::{Arc, Mutex},
    thread,
};

mod chromosomes;
mod format;
mod input;
mod reference;
mod regions;

pub(crate) use format::push_kputd;

use anyhow::{Context, bail};
use noodles_pileup::{Overlaps, Pileup};
use noodles_sam as sam;

use crate::{alignment::ReadOptions, flags, output, reference::IndexedFasta};

use self::{
    chromosomes::{pileup_by_chromosome, whole_chromosomes},
    format::{push_int, write_empty, write_sample},
    input::{Source, input_files, read_input},
    reference::ReferenceCache,
    regions::{Bed, parse_region},
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

pub fn run(args: Args) -> anyhow::Result<()> {
    if args.no_baq && args.redo_baq {
        bail!("Error: The -B option cannot be combined with -E");
    }

    let (inputs, indexes) = input_files(&args)?;

    if inputs.is_empty() {
        bail!("[mpileup] no input file/data given");
    }

    let bed = args
        .positions
        .as_deref()
        .map(Bed::read)
        .transpose()?
        .map(Arc::new);

    let exclude_rg = args.exclude_rg.as_deref().map(read_excluded_read_groups);

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

    let max_depth = max_depth(&args, inputs.len());

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

    let options = OutputOptions::from_args(&args);

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

    let samples = if args.ignore_rg {
        Samples::OnePerFile
    } else {
        Samples::FromReadGroups
    };

    match whole_chromosomes(&args, &inputs, indexes.as_deref(), &read_options)? {
        Some((headers, names)) => {
            report_samples(&headers, samples);
            pileup_by_chromosome(&context, &names, args.threads, &mut out)?;
        }
        None => pileup_region(&context, args.region.as_deref(), Some(samples), &mut out)?,
    }

    out.flush()?;

    Ok(())
}

impl OutputOptions {
    fn from_args(args: &Args) -> Self {
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

        Self {
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
        }
    }
}

/// Read groups to leave out (`-G`), whitespace-separated. An unreadable file excludes none, as
/// in samtools.
fn read_excluded_read_groups(path: &Path) -> HashSet<Vec<u8>> {
    match fs::read(path) {
        Ok(text) => text
            .split(|b| b.is_ascii_whitespace())
            .filter(|s| !s.is_empty())
            .map(<[u8]>::to_vec)
            .collect(),
        Err(_) => {
            eprintln!(
                "[bam_mpileup] Fail to open file {}. Continue anyway.",
                path.display()
            );
            HashSet::new()
        }
    }
}

/// The per-position read limit (`-d`), warning as samtools does.
fn max_depth(args: &Args, input_count: usize) -> usize {
    if args.max_depth <= 0 {
        eprintln!("[mpileup] Max depth set to maximum value ({})", i32::MAX);
        i32::MAX as usize
    } else {
        if args.max_depth * input_count as i64 > 1 << 20 {
            eprintln!("[mpileup] Combined max depth is above 1M. Potential memory hog!");
        }

        args.max_depth as usize
    }
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
/// With `report`, the number of samples is reported first.
fn pileup_region<W: Write>(
    context: &PileupContext<'_>,
    region: Option<&str>,
    report: Option<Samples>,
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

        check_reference_sequences(&headers, context.inputs)?;

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

        if let Some(samples) = report {
            report_samples(&headers, samples);
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

/// Checks that each input's reference sequences are those of the first, in the same order,
/// as positions are compared by reference sequence index. (samtools mixes up inputs whose
/// sequences differ, or fails reading one with more.)
pub(super) fn check_reference_sequences(
    headers: &[sam::Header],
    inputs: &[String],
) -> anyhow::Result<()> {
    let first = headers[0].reference_sequences();

    for (header, input) in headers.iter().zip(inputs).skip(1) {
        let matches = header.reference_sequences().len() <= first.len()
            && header
                .reference_sequences()
                .keys()
                .zip(first.keys())
                .all(|(a, b)| a == b);

        if !matches {
            bail!(
                "[mpileup] the reference sequences of {input} differ from those of {}",
                inputs[0]
            );
        }
    }

    Ok(())
}

/// How samples are counted for the report.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
enum Samples {
    /// The read groups' samples, or one for a file without any.
    FromReadGroups,
    /// One per file (`--ignore-RG`).
    OnePerFile,
}

/// Reports the number of samples, as samtools does.
fn report_samples(headers: &[sam::Header], samples: Samples) {
    use sam::header::record::value::map::read_group::tag::SAMPLE;

    let sample_count: usize = headers
        .iter()
        .map(|h| {
            let read_group_samples: HashSet<_> = h
                .read_groups()
                .values()
                .filter_map(|rg| rg.other_fields().get(&SAMPLE))
                .collect();

            if samples == Samples::OnePerFile || read_group_samples.is_empty() {
                1
            } else {
                read_group_samples.len()
            }
        })
        .sum();

    eprintln!(
        "[mpileup] {sample_count} samples in {} input files",
        headers.len()
    );
}

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

    let mut inputs = MultiPileup::new(n);
    let mut empty = EmptyLines {
        options,
        names,
        input_count: n,
        bed,
        line: Vec::new(),
    };

    let mut last_tid: i32 = -1;
    let mut last_pos: i64 = -1;
    let mut got_ref = false;
    let mut one_seq = false;
    let mut reference: Option<Arc<Vec<u8>>> = None;

    let mut line = Vec::with_capacity(1 << 16);
    let mut seq = Vec::new();
    let mut qual = Vec::new();
    let mut ins = Vec::new();

    let name_of = |tid: i32| -> &[u8] { &names[tid as usize] };

    while let Some((tid, pos)) = inputs.next(sources)? {
        one_seq = true;

        if has_region && (pos < beg0 || pos >= end0) {
            continue;
        }

        if options.all > 0 {
            // Positions of earlier sequences not yet written
            while tid > last_tid {
                if last_tid >= 0 && !has_region {
                    let end = lengths[last_tid as usize];
                    let reference = reference.as_deref().map(Vec::as_slice);
                    empty.fill(out, last_tid, last_pos + 1..end, i64::MIN, reference)?;
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
            let skip_before = if has_region { beg0 } else { i64::MIN };
            let reference = reference.as_deref().map(Vec::as_slice);
            empty.fill(out, tid, last_pos + 1..pos, skip_before, reference)?;
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
        line.push(match &reference {
            Some(r) if pos < r.len() as i64 => r[pos as usize],
            _ => b'N',
        });

        for (source, &active) in sources.iter_mut().zip(&inputs.active) {
            write_sample(
                &mut line,
                active.then_some(&mut source.pileup),
                options,
                names,
                pos,
                reference.as_deref().map(Vec::as_slice),
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

            let end = lengths[last_tid as usize].min(end0);
            let reference_bases = reference.as_deref().map(Vec::as_slice);
            empty.fill(out, last_tid, last_pos + 1..end, i64::MIN, reference_bases)?;

            last_tid += 1;
            last_pos = -1;

            if options.all < 2 || has_region {
                break;
            }
        }
    }

    Ok(())
}

/// `bam_mplp_t`: the inputs' pileups, kept at the same position.
struct MultiPileup {
    min_tid: i32,
    min_pos: i64,
    tids: Vec<i32>,
    positions: Vec<i64>,
    ready: Vec<bool>,
    /// Whether each input has reads at the current position.
    active: Vec<bool>,
}

impl MultiPileup {
    fn new(input_count: usize) -> Self {
        Self {
            min_tid: -1,
            min_pos: i64::MAX,
            tids: vec![-1; input_count],
            positions: vec![i64::MAX; input_count],
            ready: vec![false; input_count],
            active: vec![false; input_count],
        }
    }

    /// `bam_mplp64_auto`: advances the inputs at the current position and returns the next
    /// position with reads in any input, or `None` at the end.
    fn next(&mut self, sources: &mut [Source]) -> anyhow::Result<Option<(i32, i64)>> {
        let mut new_min_tid = u32::MAX;
        let mut new_min_pos = i64::MAX;

        for (i, source) in sources.iter_mut().enumerate() {
            if self.positions[i] == self.min_pos && self.tids[i] == self.min_tid {
                self.ready[i] = source.auto()?;

                if self.ready[i] {
                    let (tid, pos) = source.pileup.position();
                    self.tids[i] = tid;
                    self.positions[i] = pos;
                } else {
                    self.tids[i] = 0;
                    self.positions[i] = 0;
                }
            }

            if self.ready[i] {
                let tid = self.tids[i] as u32;

                if tid < new_min_tid {
                    new_min_tid = tid;
                    new_min_pos = self.positions[i];
                } else if tid == new_min_tid && self.positions[i] < new_min_pos {
                    new_min_pos = self.positions[i];
                }
            }
        }

        self.min_pos = new_min_pos;
        self.min_tid = new_min_tid as i32;

        if new_min_pos == i64::MAX {
            return Ok(None);
        }

        for i in 0..self.active.len() {
            self.active[i] =
                self.ready[i] && self.positions[i] == self.min_pos && self.tids[i] == self.min_tid;
        }

        Ok(Some((self.min_tid, self.min_pos)))
    }
}

/// Writes `-a` lines for positions without reads.
struct EmptyLines<'a> {
    options: &'a OutputOptions,
    names: &'a [Vec<u8>],
    input_count: usize,
    bed: Option<&'a Bed>,
    line: Vec<u8>,
}

impl EmptyLines<'_> {
    /// Writes a line for each position of `tid` in `positions`, except those before
    /// `skip_before` and those outside the BED file.
    fn fill<W: Write>(
        &mut self,
        out: &mut W,
        tid: i32,
        positions: std::ops::Range<i64>,
        skip_before: i64,
        reference: Option<&[u8]>,
    ) -> anyhow::Result<()> {
        let name = &self.names[tid as usize];

        for pos in positions {
            if pos < skip_before {
                continue;
            }

            if let Some(bed) = self.bed
                && !bed.overlaps(name, pos, pos + 1)
            {
                continue;
            }

            write_empty(
                out,
                &mut self.line,
                self.options,
                name,
                pos,
                self.input_count,
                reference,
            )?;
        }

        Ok(())
    }
}
