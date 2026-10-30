//! `ntools calmd`: recalculate MD and NM tags (samtools calmd).
//!
//! A port of `bam_md.c` from samtools 1.24, with BAQ (`-r`) from [`crate::baq`].

use std::{path::PathBuf, sync::Arc};

use anyhow::bail;

use noodles_pileup::{
    Read,
    baq::{self, BaqMode},
    md::{self, FillMdOptions},
};

use crate::global_args::{IgnoredWriteIndex, InputFormatOptions, OutputFormatOptions, Verbosity};
use crate::{
    alignment::{AlignmentReader, RawEncoder, ReadOptions},
    output::{self, AlignmentWriter, OutputFormat, WriteOptions},
    program,
    reference::IndexedFasta,
};

#[derive(Debug, clap::Args)]
#[command(disable_help_flag = true)]
pub struct Args {
    /// Input SAM, BAM, or CRAM file.
    input: String,

    /// Indexed reference FASTA [default: --reference].
    reference_fasta: Option<PathBuf>,

    /// Change bases identical to the reference to `=`.
    #[arg(short = 'e')]
    use_equal: bool,

    /// Output uncompressed BAM.
    #[arg(short = 'u')]
    uncompressed: bool,

    /// Output compressed BAM.
    #[arg(short = 'b')]
    bam: bool,

    /// Output format [default: BAM with -b or -u, else SAM].
    #[arg(long = "output-fmt", ignore_case = true)]
    output_format: Option<OutputFormat>,

    /// Ignored (the input format is detected).
    #[arg(short = 'S')]
    sam: bool,

    /// Modify the qualities (with -r, cap them by BAQ).
    #[arg(short = 'A')]
    apply_baq: bool,

    /// Quiet mode: do not report changed tags.
    #[arg(short = 'Q')]
    quiet: bool,

    /// Compute BQ tags (without -A) or cap base qualities by BAQ (with -A).
    #[arg(short = 'r')]
    realign: bool,

    /// Extended BAQ.
    #[arg(short = 'E')]
    extended_baq: bool,

    /// Drop all tags but RG.
    #[arg(short = 'd')]
    drop_tags: bool,

    /// Bin base qualities.
    #[arg(short = 'q')]
    bin_qualities: bool,

    /// Ignored (samtools accepts it).
    #[arg(short = 'h')]
    hash_qnm: bool,

    /// Do not update MD or NM.
    #[arg(short = 'N')]
    no_update: bool,

    /// Mask matching bases of reads with at least this many mismatches (0 disables).
    #[arg(short = 'n', default_value_t = 0)]
    max_nm: i32,

    /// Cap mapping qualities; 50 is recommended for BWA (0 disables).
    #[arg(short = 'C', default_value_t = 0)]
    cap_q: i32,

    /// Reference FASTA (indexed).
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

    #[command(flatten)]
    input_format_options: InputFormatOptions,

    #[command(flatten)]
    output_format_options: OutputFormatOptions,

    #[command(flatten)]
    write_index: IgnoredWriteIndex,

    #[command(flatten)]
    verbosity: Verbosity,
}

pub fn run(args: Args) -> anyhow::Result<()> {
    let Some(reference_path) = args
        .reference_fasta
        .clone()
        .or_else(|| args.reference.clone())
    else {
        bail!("Usage: ntools calmd [-eubrAESQ] <aln.bam> <ref.fasta>");
    };

    let mut read_options = ReadOptions {
        threads: args.threads,
        reference: Some(reference_path.clone()),
    };
    args.input_format_options.configure(&mut read_options)?;

    let (mut reader, _) = AlignmentReader::open(&args.input, &read_options)?;
    args.input_format_options.apply(&mut reader)?;
    let mut header = reader.read_header()?;
    let input_header = header.clone();

    if header.reference_sequences().is_empty() {
        eprintln!("[bam_fillmd] warning: input SAM does not have header, performing a no-op.");
    }

    if !args.no_pg {
        program::add(&mut header)?;
    }

    let fasta = IndexedFasta::open(&reference_path).map_err(|e| {
        anyhow::anyhow!(
            "Failed to open reference file '{}': {e}",
            reference_path.display()
        )
    })?;

    let write_options = WriteOptions {
        format: args
            .output_format
            .unwrap_or(if args.bam || args.uncompressed {
                OutputFormat::Bam
            } else {
                OutputFormat::Sam
            }),
        level: output::parse_level(if args.uncompressed {
            0
        } else {
            args.output_format_options.level(6)?
        })?,
        threads: args.threads,
        reference: Some(reference_path.clone()),
    };

    let mut writer = AlignmentWriter::create(None, &write_options, &header)?;
    writer.write_header(&header)?;

    let options = FillMdOptions {
        update_tags: !args.no_update,
        use_equal: args.use_equal,
        max_nm: args.max_nm,
        drop_tags: args.drop_tags,
        bin_qualities: args.bin_qualities,
    };

    let mode = BaqMode {
        apply: args.apply_baq,
        extend: args.extended_baq,
        redo: false,
    };

    let names: Vec<Vec<u8>> = input_header
        .reference_sequences()
        .keys()
        .map(|name| name.to_vec())
        .collect();

    let mut last_tid = -2;
    let mut reference: Option<Arc<Vec<u8>>> = None;
    let mut cache: Vec<Option<Arc<Vec<u8>>>> = Vec::new();
    let mut skipped = 0u32;
    let mut encoder = RawEncoder::new();
    let mut batch: Vec<Job> = Vec::with_capacity(BATCH_SIZE);

    let context = Context {
        names: &names,
        options,
        quiet: args.quiet,
        mode,
        realign: args.realign,
        cap_q: args.cap_q,
    };

    let mut flush = |batch: &mut Vec<Job>, writer: &mut AlignmentWriter| -> anyhow::Result<()> {
        let results = process_batch(std::mem::take(batch), &context, args.threads)?;

        for (raw, filled) in results {
            if !filled {
                skipped += 1;
            }

            writer.write_raw(&header, &raw)?;
        }

        Ok(())
    };

    reader.for_each(&input_header, |rec| {
        let raw = rec.encode_raw(&input_header, &mut encoder)?.to_vec();

        let tid = crate::bam_raw::le_u32(&raw[0..4]).cast_signed();

        if tid >= 0 && tid != last_tid {
            reference = fetch(&fasta, &names, &mut cache, tid, last_tid)?;
            last_tid = tid;

            if reference.is_none() {
                let name = names
                    .get(tid as usize)
                    .map_or_else(|| "(unknown)".into(), |n| String::from_utf8_lossy(n));
                eprintln!("[bam_fillmd] fail to find sequence '{name}' in the reference.");

                if args.realign || args.cap_q > 10 {
                    bail!("fail to find sequence '{name}' in the reference");
                }
            }
        }

        let reference = if tid >= 0 { reference.clone() } else { None };
        batch.push(Job { raw, reference });

        if batch.len() >= BATCH_SIZE {
            flush(&mut batch, &mut writer)?;
        }

        Ok(())
    })?;

    flush(&mut batch, &mut writer)?;
    writer.finish(&header)?;

    if skipped > 0 {
        eprintln!("[calmd] Warning: {skipped} records skipped due to no query sequence");
    }

    Ok(())
}

/// A processed record, and whether it had a sequence.
type Processed = (Vec<u8>, bool);

/// Records per batch processed across threads.
const BATCH_SIZE: usize = 4096;

/// A record and the reference sequence of its reference sequence ID, if any.
struct Job {
    raw: Vec<u8>,
    reference: Option<Arc<Vec<u8>>>,
}

/// What processing a record needs.
struct Context<'a> {
    names: &'a [Vec<u8>],
    options: FillMdOptions,
    quiet: bool,
    mode: BaqMode,
    realign: bool,
    cap_q: i32,
}

/// Realigns, caps, and fills MD/NM in one record. Returns the record and whether it had a
/// sequence.
fn process(job: Job, context: &Context<'_>) -> anyhow::Result<Processed> {
    let Some(reference) = job.reference else {
        return Ok((job.raw, true));
    };

    let mut read = Read::new(job.raw)?;

    if context.realign
        && baq::prob_realn(&mut read, &reference, context.mode) == Err(baq::BaqSkipped::Invalid)
    {
        bail!("BAQ alignment failed");
    }

    if context.cap_q > 10 {
        let q = baq::cap_mapq(&read, &reference, context.cap_q).unwrap_or(-1);

        // As in samtools, -1 (skip) becomes a mapping quality of 255.
        if i32::from(read.mapping_quality()) > q {
            read.set_mapping_quality(q as u8);
        }
    }

    let filled = match md::fill_md(&mut read, &reference, &context.options) {
        Some(report) => {
            if !context.quiet {
                let name = String::from_utf8_lossy(read.name());

                if let Some((old, new)) = report.nm {
                    eprintln!(
                        "[bam_fillmd1] different NM for read '{name}': {} -> {new}",
                        old as i32
                    );
                }

                if let Some((old, new)) = &report.md {
                    eprintln!(
                        "[bam_fillmd1] different MD for read '{name}': '{}' -> '{}'",
                        String::from_utf8_lossy(old),
                        String::from_utf8_lossy(new)
                    );
                }
            }

            true
        }
        None => {
            if !context.quiet {
                let ref_name = &context.names[read.reference_sequence_id() as usize];
                eprintln!(
                    "[bam_fillmd1] no sequence in alignment record for '{}' at {}:{}, skipped",
                    String::from_utf8_lossy(read.name()),
                    String::from_utf8_lossy(ref_name),
                    read.position() + 1
                );
            }

            false
        }
    };

    Ok((read.into_bytes(), filled))
}

/// Processes a batch, on `threads` threads when there are several, keeping its order.
fn process_batch(
    batch: Vec<Job>,
    context: &Context<'_>,
    threads: usize,
) -> anyhow::Result<Vec<Processed>> {
    if threads <= 1 || batch.len() < 64 {
        return batch.into_iter().map(|job| process(job, context)).collect();
    }

    let chunk_size = batch.len().div_ceil(threads);
    let mut chunks: Vec<Vec<Job>> = Vec::with_capacity(threads);
    let mut jobs = batch.into_iter();

    loop {
        let chunk: Vec<Job> = jobs.by_ref().take(chunk_size).collect();

        if chunk.is_empty() {
            break;
        }

        chunks.push(chunk);
    }

    let results: Vec<anyhow::Result<Vec<Processed>>> = std::thread::scope(|scope| {
        let handles: Vec<_> = chunks
            .into_iter()
            .map(|chunk| {
                scope.spawn(move || chunk.into_iter().map(|job| process(job, context)).collect())
            })
            .collect();

        handles
            .into_iter()
            .map(|h| h.join().expect("calmd worker panicked"))
            .collect()
    });

    let mut out = Vec::new();

    for result in results {
        out.extend(result?);
    }

    Ok(out)
}

/// `get_ref`: the whole reference sequence for `tid`, keeping every one fetched once the input
/// goes back to an earlier sequence (it is then not sorted by position).
fn fetch(
    fasta: &IndexedFasta,
    names: &[Vec<u8>],
    cache: &mut Vec<Option<Arc<Vec<u8>>>>,
    tid: i32,
    last_tid: i32,
) -> anyhow::Result<Option<Arc<Vec<u8>>>> {
    let i = tid as usize;

    if let Some(Some(sequence)) = cache.get(i) {
        return Ok(Some(Arc::clone(sequence)));
    }

    let Some(name) = names.get(i) else {
        return Ok(None);
    };

    let mut buf = Vec::new();

    if !fasta.fetch(name, 0, i64::MAX, &mut buf)? {
        return Ok(None);
    }

    let sequence = Arc::new(buf);

    if !cache.is_empty() || last_tid > tid {
        cache.resize(names.len(), None);
        cache[i] = Some(Arc::clone(&sequence));
    }

    Ok(Some(sequence))
}
