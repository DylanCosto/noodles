//! `ntools fixmate`: fill in mate information from name-grouped alignments (samtools fixmate).
//!
//! This follows `bam_mating_core` and `bam_sanitize` in samtools 1.24 step by step, including
//! the order in which tags are replaced.

use std::{collections::BTreeMap, path::PathBuf, thread};

use anyhow::{Context, bail};
use bstr::BString;
use noodles_core::Position;
use noodles_sam::{
    self as sam,
    alignment::{
        RecordBuf,
        io::Write as _,
        record::{
            Flags, MappingQuality,
            cigar::{Op, op::Kind},
            data::field::Tag,
        },
        record_buf::data::field::Value,
    },
    header::record::value::map::header::{sort_order::COORDINATE, tag::SORT_ORDER},
};

use crate::{
    alignment::{AlignmentReader, ReadOptions},
    output::{self, AlignmentWriter, OutputFormat, WriteOptions},
    program,
    sanitize::Sanitize,
};

/// Base qualities below this are not counted in mate scores (`MD_MIN_QUALITY`).
const MIN_SCORED_QUALITY: u8 = 15;

/// Largest CIGAR operation length (28 bits).
const MAX_OP_LEN: usize = (1 << 28) - 1;

const MATE_CIGAR: Tag = Tag::new(b'M', b'C');
const MATE_MAPPING_QUALITY: Tag = Tag::new(b'M', b'Q');
const MATE_SCORE: Tag = Tag::new(b'm', b's');
const TEMPLATE_CIGAR: Tag = Tag::new(b'c', b't');

/// Tags removed from unmapped records when sanitizing auxiliary data.
const UNMAPPED_REMOVED_TAGS: [Tag; 4] = [
    Tag::EDIT_DISTANCE,
    Tag::MISMATCHED_POSITIONS,
    Tag::new(b'C', b'G'),
    Tag::new(b'S', b'M'),
];

#[derive(Debug, clap::Args)]
pub struct Args {
    /// Name-grouped input SAM, BAM, or CRAM file.
    input: String,

    /// Output file.
    output: PathBuf,

    /// Remove unmapped reads and secondary alignments.
    #[arg(short = 'r')]
    remove_reads: bool,

    /// Disable the FR proper pair check.
    #[arg(short = 'p')]
    no_proper_pair_check: bool,

    /// Add template CIGAR (ct) tags.
    #[arg(short = 'c')]
    template_cigar: bool,

    /// Add mate score (ms) tags.
    #[arg(short = 'm')]
    mate_score: bool,

    /// Output uncompressed BAM.
    #[arg(short = 'u')]
    uncompressed: bool,

    /// Sanitize alignment fields [all].
    #[arg(short = 'z', long = "sanitize", default_value = "all")]
    sanitize: String,

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

#[derive(Clone, Copy, Debug)]
struct Options {
    remove_reads: bool,
    proper_pair_check: bool,
    template_cigar: bool,
    mate_score: bool,
}

pub fn run(args: Args) -> anyhow::Result<()> {
    let sanitize = Sanitize::parse(&args.sanitize)?;

    let read_options = ReadOptions {
        threads: args.threads,
        reference: args.reference.clone(),
    };

    let (mut reader, _) = AlignmentReader::open(&args.input, &read_options)?;
    let mut header = reader.read_header()?;

    let is_coordinate_sorted = header
        .header()
        .and_then(|hd| hd.other_fields().get(&SORT_ORDER))
        .is_some_and(|so| so == COORDINATE);

    if is_coordinate_sorted {
        bail!("[bam_mating_core] ERROR: Coordinate sorted, require grouped/sorted by queryname.");
    }

    if !args.no_pg {
        program::add(&mut header)?;
    }

    let header = header;

    let format = args
        .output_format
        .or_else(|| OutputFormat::from_path(&args.output))
        .unwrap_or(OutputFormat::Bam);

    let write_options = WriteOptions {
        format,
        level: output::parse_level(if args.uncompressed { 0 } else { 6 })?,
        threads: args.threads,
        reference: args.reference,
    };

    let mut writer = AlignmentWriter::create(Some(&args.output), &write_options, &header)?;
    writer.write_header(&header)?;

    let options = Options {
        remove_reads: args.remove_reads,
        proper_pair_check: !args.no_proper_pair_check,
        template_cigar: args.template_cigar,
        mate_score: args.mate_score,
    };

    let lengths: Vec<usize> = header
        .reference_sequences()
        .values()
        .map(|rs| rs.length().get())
        .collect();

    if args.threads > 0
        && matches!(reader, AlignmentReader::Bam(_))
        && writer.preformatted_kind() == Some(OutputFormat::Bam)
    {
        let config = BatchConfig {
            header: &header,
            lengths: &lengths,
            sanitize,
            options,
        };

        fix_in_parallel(&mut reader, &config, args.threads, &mut writer)?;
        writer.finish(&header)?;
        return Ok(());
    }

    let mut template: Vec<RecordBuf> = Vec::new();
    let mut framed = Vec::new();

    let flush =
        |template: &mut Vec<RecordBuf>, writer: &mut AlignmentWriter| -> anyhow::Result<()> {
            fix_template(template, options)?;

            for record in template.iter() {
                let flags = record.flags();

                if options.remove_reads && flags.intersects(Flags::SECONDARY | Flags::UNMAPPED) {
                    continue;
                }

                writer.write(&header, &crate::alignment::Rec::Buf(record))?;
            }

            template.clear();
            Ok(())
        };

    reader.for_each(&header, |record| {
        let mut record = match &record {
            // The BAM decoder builds record buffers much faster than the generic conversion.
            crate::alignment::Rec::Bam(r) => {
                framed.clear();
                crate::runs::write_raw_record(&mut framed, r.as_ref())?;
                let mut buf = RecordBuf::default();
                noodles_bam::io::Reader::from(&framed[..]).read_record_buf(&header, &mut buf)?;
                buf
            }
            _ => RecordBuf::try_from_alignment_record(&header, record.as_dyn())?,
        };
        sanitize_record(&mut record, &lengths, sanitize)?;

        if template
            .first()
            .is_some_and(|first| first.name() != record.name())
        {
            flush(&mut template, &mut writer)?;
        }

        template.push(record);
        Ok(())
    })?;

    if !template.is_empty() {
        flush(&mut template, &mut writer)?;
    }

    writer.finish(&header)?;

    Ok(())
}

/// Raw record bytes per batch sent to a worker.
const BATCH_CAPACITY: usize = 1 << 20;

/// What a worker needs to fix a batch of templates.
struct BatchConfig<'a> {
    header: &'a sam::Header,
    lengths: &'a [usize],
    sanitize: Sanitize,
    options: Options,
}

/// Fixes templates on `threads` workers: this thread batches raw records at template
/// boundaries, workers decode, fix, and re-encode them, and a writer thread writes batches in
/// input order.
fn fix_in_parallel(
    reader: &mut AlignmentReader,
    config: &BatchConfig<'_>,
    threads: usize,
    writer: &mut AlignmentWriter,
) -> anyhow::Result<()> {
    let header = config.header;

    thread::scope(|scope| {
        let (work_tx, work_rx) = crossbeam_channel::bounded::<(u64, Vec<u8>)>(threads * 2);
        let (done_tx, done_rx) =
            crossbeam_channel::bounded::<(u64, anyhow::Result<Vec<u8>>)>(threads * 2);

        for _ in 0..threads {
            let work_rx = work_rx.clone();
            let done_tx = done_tx.clone();

            scope.spawn(move || {
                for (id, batch) in work_rx {
                    let result = fix_batch(&batch, config);

                    if done_tx.send((id, result)).is_err() {
                        break;
                    }
                }
            });
        }

        drop(work_rx);
        drop(done_tx);

        let write = scope.spawn(move || -> anyhow::Result<()> {
            let mut pending = BTreeMap::new();
            let mut next = 0;

            for (id, result) in done_rx {
                pending.insert(id, result?);

                while let Some(records) = pending.remove(&next) {
                    writer.write_bam_records(&records)?;
                    next += 1;
                }
            }

            Ok(())
        });

        let mut batch = Vec::with_capacity(BATCH_CAPACITY);
        let mut last_name = Vec::new();
        let mut id = 0;

        let produced = reader.for_each(header, |record| {
            let crate::alignment::Rec::Bam(record) = record else {
                bail!("parallel fixmate needs BAM input");
            };

            let raw = record.as_ref();
            let name = crate::sort_order::name(raw);

            // Batches end only between templates.
            if batch.len() >= BATCH_CAPACITY && name != last_name.as_slice() {
                let full = std::mem::replace(&mut batch, Vec::with_capacity(BATCH_CAPACITY));
                work_tx
                    .send((id, full))
                    .map_err(|_| anyhow::anyhow!("fixmate writer stopped"))?;
                id += 1;
            }

            last_name.clear();
            last_name.extend_from_slice(name);
            crate::runs::write_raw_record(&mut batch, raw)?;

            Ok(())
        });

        let sent = match produced {
            Ok(()) if !batch.is_empty() => work_tx
                .send((id, batch))
                .map_err(|_| anyhow::anyhow!("fixmate writer stopped")),
            other => other,
        };

        drop(work_tx);

        let written = write
            .join()
            .map_err(|_| anyhow::anyhow!("fixmate writer thread panicked"))?;

        // A writer error explains a failed send, so report it first.
        written.and(sent)
    })
}

/// Decodes, fixes, and re-encodes the templates of one batch.
fn fix_batch(batch: &[u8], config: &BatchConfig<'_>) -> anyhow::Result<Vec<u8>> {
    let mut reader = noodles_bam::io::Reader::from(batch);
    let mut encoder = noodles_bam::io::Writer::from(Vec::with_capacity(batch.len()));
    let mut template: Vec<RecordBuf> = Vec::new();

    let mut flush = |template: &mut Vec<RecordBuf>| -> anyhow::Result<()> {
        fix_template(template, config.options)?;

        for record in template.iter() {
            let removed = config.options.remove_reads
                && record
                    .flags()
                    .intersects(Flags::SECONDARY | Flags::UNMAPPED);

            if !removed {
                encoder.write_alignment_record(config.header, record)?;
            }
        }

        template.clear();
        Ok(())
    };

    loop {
        let mut record = RecordBuf::default();

        if reader.read_record_buf(config.header, &mut record)? == 0 {
            break;
        }

        sanitize_record(&mut record, config.lengths, config.sanitize)?;

        if template
            .first()
            .is_some_and(|first| first.name() != record.name())
        {
            flush(&mut template)?;
        }

        template.push(record);
    }

    if !template.is_empty() {
        flush(&mut template)?;
    }

    Ok(encoder.into_inner())
}

/// The 0-based position, or -1.
fn pos(record: &RecordBuf) -> i64 {
    record.alignment_start().map_or(-1, |p| {
        i64::try_from(usize::from(p)).unwrap_or(i64::MAX) - 1
    })
}

fn set_pos(record: &mut RecordBuf, pos: i64) {
    *record.alignment_start_mut() = usize::try_from(pos + 1).ok().and_then(Position::new);
}

fn tid(record: &RecordBuf) -> Option<usize> {
    record.reference_sequence_id()
}

fn reference_span(ops: &[Op]) -> usize {
    ops.iter()
        .filter(|op| op.kind().consumes_reference())
        .map(|op| op.len())
        .sum()
}

/// HTSlib's `bam_endpos`: the 0-based exclusive end, covering at least one base.
fn end_pos(record: &RecordBuf) -> i64 {
    let span = if record.flags().is_unmapped() {
        0
    } else {
        reference_span(record.cigar().as_ref())
    };

    pos(record) + i64::try_from(span.max(1)).unwrap_or(i64::MAX)
}

/// HTSlib's `bam_sanitize`.
fn sanitize_record(
    record: &mut RecordBuf,
    lengths: &[usize],
    sanitize: Sanitize,
) -> anyhow::Result<()> {
    if sanitize.has(Sanitize::POS) && tid(record).is_none() {
        set_pos(record, -1);

        if sanitize.has(Sanitize::UNMAP) {
            *record.flags_mut() |= Flags::UNMAPPED;
        }
    }

    if sanitize.has(Sanitize::CIGAR) && !record.flags().is_unmapped() {
        if pos(record) < 0 && sanitize.has(Sanitize::UNMAP) {
            *record.flags_mut() |= Flags::UNMAPPED;
        } else {
            let length = tid(record)
                .and_then(|id| lengths.get(id))
                .map_or(0, |&n| i64::try_from(n).unwrap_or(i64::MAX));

            if pos(record) >= length && sanitize.has(Sanitize::UNMAP) {
                *record.flags_mut() |= Flags::UNMAPPED;

                if sanitize.has(Sanitize::POS) {
                    *record.reference_sequence_id_mut() = None;
                    set_pos(record, -1);
                }
            } else if end_pos(record) > length {
                trim(record, length);
            }
        }
    }

    if record.flags().is_unmapped() {
        if sanitize.has(Sanitize::CIGAR) && !record.cigar().as_ref().is_empty() {
            record.cigar_mut().as_mut().clear();
        }

        if sanitize.has(Sanitize::MQUAL) {
            *record.mapping_quality_mut() = MappingQuality::new(0);
        }

        if sanitize.has(Sanitize::AUX) {
            for tag in UNMAPPED_REMOVED_TAGS {
                while record.data_mut().shift_remove(&tag).is_some() {}
            }
        }
    }

    if sanitize.has(Sanitize::CIGARX) && !record.flags().is_unmapped() {
        for op in record.cigar_mut().as_mut() {
            if matches!(op.kind(), Kind::SequenceMatch | Kind::SequenceMismatch) {
                *op = Op::new(Kind::Match, op.len());
            }
        }
    }

    if sanitize.has(Sanitize::CIGDUP) && !record.flags().is_unmapped() {
        merge_cigar_ops(record.cigar_mut().as_mut());
    }

    Ok(())
}

/// Merges adjacent operations of the same kind and drops empty ones (`FIX_CIGDUP`).
fn merge_cigar_ops(ops: &mut Vec<Op>) {
    let mut merged: Vec<Op> = Vec::with_capacity(ops.len());

    for &op in ops.iter() {
        if op.is_empty() {
            continue;
        }

        match merged.last_mut() {
            Some(last) if last.kind() == op.kind() => {
                let total = last.len() + op.len();

                if total > MAX_OP_LEN {
                    *last = Op::new(op.kind(), MAX_OP_LEN);
                    merged.push(Op::new(op.kind(), total - MAX_OP_LEN));
                } else {
                    *last = Op::new(op.kind(), total);
                }
            }
            _ => merged.push(op),
        }
    }

    *ops = merged;
}

/// HTSlib's `bam_trim`: soft clips the alignment past reference position `end`.
fn trim(record: &mut RecordBuf, end: i64) {
    let ops: Vec<Op> = record.cigar().as_ref().to_vec();
    let mut position = pos(record);
    let mut i = 0;
    let mut kind = Kind::Match;
    let mut len = 0;

    while i < ops.len() {
        kind = ops[i].kind();
        len = i64::try_from(ops[i].len()).unwrap_or(i64::MAX);

        if kind.consumes_reference() {
            position += len;

            if position > end {
                break;
            }
        }

        i += 1;
    }

    if i == ops.len() {
        return;
    }

    let mut kept: Vec<Op> = ops[..i].to_vec();
    let mut trailing: Vec<Op> = Vec::new();

    let to_len = |n: i64| usize::try_from(n).unwrap_or(0);

    if position - len < end {
        // The operation straddles the end: split it.
        kept.push(Op::new(kind, to_len(end - (position - len))));
        trailing.push(Op::new(Kind::SoftClip, to_len(position - end)));
    } else if position > end {
        *record.flags_mut() |= Flags::UNMAPPED;
        *record.flags_mut() &= !Flags::PROPERLY_SEGMENTED;
    } else {
        trailing.push(Op::new(Kind::SoftClip, ops[i].len()));
    }

    for op in &ops[i + 1..] {
        if op.kind() == Kind::HardClip {
            trailing.push(*op);
        } else if let Some(last) = trailing.last_mut() {
            *last = Op::new(Kind::SoftClip, last.len() + op.len());
        } else {
            trailing.push(Op::new(Kind::SoftClip, op.len()));
        }
    }

    kept.extend(trailing);
    *record.cigar_mut() = kept.into_iter().collect();
}

fn cigar_string(record: &RecordBuf) -> String {
    let ops = record.cigar().as_ref();

    if ops.is_empty() {
        return "*".to_owned();
    }

    ops.iter()
        .map(|op| format!("{}{}", op.len(), char::from(kind_char(op.kind()))))
        .collect()
}

fn kind_char(kind: Kind) -> u8 {
    match kind {
        Kind::Match => b'M',
        Kind::Insertion => b'I',
        Kind::Deletion => b'D',
        Kind::Skip => b'N',
        Kind::SoftClip => b'S',
        Kind::HardClip => b'H',
        Kind::Pad => b'P',
        Kind::SequenceMatch => b'=',
        Kind::SequenceMismatch => b'X',
    }
}

/// Removes a tag (keeping the order of the others) and appends a new value, as samtools'
/// `bam_aux_del` followed by `bam_aux_append` does.
fn replace_tag(record: &mut RecordBuf, tag: Tag, value: Value) {
    let data = record.data_mut();
    data.shift_remove(&tag);
    data.insert(tag, value);
}

/// `sync_mate_inner`
fn sync_mate_inner(src: &RecordBuf, dest: &mut RecordBuf) {
    *dest.mate_reference_sequence_id_mut() = tid(src);
    *dest.mate_alignment_start_mut() = src.alignment_start();

    let flags = dest.flags_mut();
    flags.set(
        Flags::MATE_REVERSE_COMPLEMENTED,
        src.flags().is_reverse_complemented(),
    );

    if src.flags().is_unmapped() {
        flags.insert(Flags::MATE_UNMAPPED);
    }
}

/// `sync_mq_mc`
fn sync_mq_mc(src: &RecordBuf, dest: &mut RecordBuf) {
    if !src.flags().is_unmapped() {
        let mq = src.mapping_quality().map_or(255, u8::from);
        replace_tag(dest, MATE_MAPPING_QUALITY, Value::Int32(i32::from(mq)));
    }

    if !src.flags().is_unmapped() || !dest.flags().is_unmapped() {
        replace_tag(
            dest,
            MATE_CIGAR,
            Value::String(BString::from(cigar_string(src))),
        );
    }
}

/// `sync_unmapped_pos_inner`
fn sync_unmapped_pos_inner(src: &RecordBuf, dest: &mut RecordBuf) {
    if dest.flags().is_unmapped() && !src.flags().is_unmapped() {
        *dest.reference_sequence_id_mut() = tid(src);
        *dest.alignment_start_mut() = src.alignment_start();
    }
}

/// `calc_mate_score`: the sum of base qualities of at least 15. Missing qualities count as 255,
/// as they are stored in BAM.
fn mate_score(record: &RecordBuf) -> u32 {
    let scores = record.quality_scores().as_ref();

    if scores.is_empty() {
        let n = u32::try_from(record.sequence().len()).unwrap_or(u32::MAX);
        return n.saturating_mul(255);
    }

    scores
        .iter()
        .filter(|&&q| q >= MIN_SCORED_QUALITY)
        .map(|&q| u32::from(q))
        .sum()
}

/// `plausibly_properly_paired`
fn plausibly_properly_paired(a: &RecordBuf, b: &RecordBuf) -> bool {
    if a.flags().is_unmapped() || b.flags().is_unmapped() || tid(a) != tid(b) {
        return false;
    }

    let five_prime = |r: &RecordBuf| {
        if r.flags().is_reverse_complemented() {
            end_pos(r)
        } else {
            pos(r)
        }
    };

    let (first, second) = if five_prime(a) > five_prime(b) {
        (b, a)
    } else {
        (a, b)
    };

    !first.flags().is_reverse_complemented() && second.flags().is_reverse_complemented()
}

/// `bam_template_cigar`: sets `ct` on the lower-positioned record.
fn template_cigar(a: &mut RecordBuf, b: &mut RecordBuf) {
    let unusable = |r: &RecordBuf| tid(r).is_none() || pos(r) < 0 || r.flags().is_unmapped();

    if tid(a) != tid(b) || unusable(a) || unusable(b) {
        return;
    }

    let swapped = pos(a) > pos(b);
    let (b1, b2): (&RecordBuf, &RecordBuf) = if swapped { (b, a) } else { (a, b) };

    let segment = |r: &RecordBuf| {
        if r.flags().is_first_segment() {
            '1'
        } else {
            '2'
        }
    };
    let strand = |r: &RecordBuf| {
        if r.flags().is_reverse_complemented() {
            'R'
        } else {
            'F'
        }
    };

    let mut s = String::new();
    s.push(segment(b1));
    s.push(strand(b1));
    for op in b1.cigar().as_ref() {
        s.push_str(&format!("{}{}", op.len(), char::from(kind_char(op.kind()))));
    }
    s.push_str(&(pos(b2) - end_pos(b1)).to_string());
    s.push('T');
    s.push(segment(b2));
    s.push(strand(b2));
    for op in b2.cigar().as_ref() {
        s.push_str(&format!("{}{}", op.len(), char::from(kind_char(op.kind()))));
    }

    a.data_mut().shift_remove(&TEMPLATE_CIGAR);
    b.data_mut().shift_remove(&TEMPLATE_CIGAR);

    let target = if swapped { b } else { a };
    target
        .data_mut()
        .insert(TEMPLATE_CIGAR, Value::String(BString::from(s)));
}

/// Returns mutable references to two distinct elements.
fn pair_mut(records: &mut [RecordBuf], i: usize, j: usize) -> (&mut RecordBuf, &mut RecordBuf) {
    assert_ne!(i, j, "a record cannot be its own mate");

    if i < j {
        let (left, right) = records.split_at_mut(j);
        (&mut left[i], &mut right[0])
    } else {
        let (left, right) = records.split_at_mut(i);
        (&mut right[0], &mut left[j])
    }
}

/// Fixes one template (records sharing a name), as `bam_mating_core` does.
fn fix_template(records: &mut [RecordBuf], options: Options) -> anyhow::Result<()> {
    let not_primary = Flags::SECONDARY | Flags::SUPPLEMENTARY;

    let mut pre: Option<usize> = None;
    let mut cur: Option<usize> = None;
    let mut by_segment: [Option<usize>; 2] = [None, None];
    let mut pre_end = 0;

    for n in 0..records.len() {
        let flags = records[n].flags();

        if flags.intersects(not_primary) {
            continue;
        }

        let is_r2 = usize::from(flags.is_last_segment());

        let Some(p) = pre else {
            pre = Some(n);
            by_segment[is_r2] = Some(n);
            pre_end = if flags.is_unmapped() {
                0
            } else {
                end_pos(&records[n])
            };
            continue;
        };

        cur = Some(n);
        by_segment[is_r2] = Some(n);
        let cur_end = if flags.is_unmapped() {
            0
        } else {
            end_pos(&records[n])
        };

        let (a, b) = pair_mut(records, p, n);

        *a.flags_mut() |= Flags::SEGMENTED;
        *b.flags_mut() |= Flags::SEGMENTED;

        // sync_mate
        sync_unmapped_pos_inner(a, b);
        sync_unmapped_pos_inner(b, a);
        sync_mate_inner(a, b);
        sync_mate_inner(b, a);
        sync_mq_mc(a, b);
        sync_mq_mc(b, a);

        let unmapped = Flags::UNMAPPED | Flags::MATE_UNMAPPED;

        if tid(a) == tid(b) && !b.flags().intersects(unmapped) && !a.flags().intersects(unmapped) {
            let cur5 = if b.flags().is_reverse_complemented() {
                cur_end
            } else {
                pos(b)
            };
            let pre5 = if a.flags().is_reverse_complemented() {
                pre_end
            } else {
                pos(a)
            };

            *b.template_length_mut() = i32::try_from(pre5 - cur5).context("TLEN overflows")?;
            *a.template_length_mut() = i32::try_from(cur5 - pre5).context("TLEN overflows")?;
        } else {
            *b.template_length_mut() = 0;
            *a.template_length_mut() = 0;
        }

        if options.template_cigar {
            template_cigar(a, b);
        }

        if options.proper_pair_check && !plausibly_properly_paired(a, b) {
            *a.flags_mut() &= !Flags::PROPERLY_SEGMENTED;
            *b.flags_mut() &= !Flags::PROPERLY_SEGMENTED;
        }

        if options.mate_score {
            let (score_a, score_b) = (mate_score(a), mate_score(b));
            replace_tag(b, MATE_SCORE, score_value(score_a));
            replace_tag(a, MATE_SCORE, score_value(score_b));
        }

        if options.remove_reads {
            let orphaned = Flags::MATE_REVERSE_COMPLEMENTED | Flags::PROPERLY_SEGMENTED;

            if a.flags().is_unmapped() {
                *b.flags_mut() &= !orphaned;
            }

            if b.flags().is_unmapped() {
                *a.flags_mut() &= !orphaned;
            }
        }
    }

    if let (None, Some(p)) = (cur, pre) {
        let record = &mut records[p];
        *record.mate_reference_sequence_id_mut() = None;
        *record.mate_alignment_start_mut() = None;
        *record.template_length_mut() = 0;
        *record.flags_mut() &= !(Flags::MATE_REVERSE_COMPLEMENTED | Flags::PROPERLY_SEGMENTED);
    }

    for n in 0..records.len() {
        let flags = records[n].flags();

        if !flags.intersects(not_primary) {
            continue;
        }

        // Supplementary alignments get their mate's primary information.
        if flags.is_supplementary() && flags.is_segmented() {
            let mate = by_segment[usize::from(!flags.is_last_segment())];

            if let Some(m) = mate {
                let (mate, record) = pair_mut(records, m, n);
                sync_mate_inner(mate, record);
                sync_mq_mc(mate, record);
            }
        }
    }

    Ok(())
}

/// samtools appends mate scores as `i` (32-bit) values.
fn score_value(score: u32) -> Value {
    i32::try_from(score).map_or(Value::UInt32(score), Value::Int32)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_merge_cigar_ops() {
        let mut ops = vec![
            Op::new(Kind::Match, 5),
            Op::new(Kind::Match, 3),
            Op::new(Kind::Insertion, 0),
            Op::new(Kind::Deletion, 2),
            Op::new(Kind::Match, 4),
        ];
        merge_cigar_ops(&mut ops);
        assert_eq!(
            ops,
            [
                Op::new(Kind::Match, 8),
                Op::new(Kind::Deletion, 2),
                Op::new(Kind::Match, 4)
            ]
        );
    }
}
