//! `ntools calmd`: recalculate MD and NM tags (samtools calmd).
//!
//! A port of `bam_md.c` from samtools 1.24, with BAQ (`-r`) from [`crate::baq`].

use std::{path::PathBuf, sync::Arc};

use anyhow::bail;

use crate::{
    alignment::{AlignmentReader, ReadOptions, Rec},
    baq::{self, BaqMode, SEQ_NT16_TABLE},
    output::{self, AlignmentWriter, OutputFormat, WriteOptions},
    pileup::Read,
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
}

/// What `bam_fillmd1_core` changes.
#[derive(Clone, Copy, Debug)]
struct Options {
    use_equal: bool,
    drop_tags: bool,
    bin_qualities: bool,
    update: bool,
    max_nm: i32,
    quiet: bool,
}

pub fn run(args: Args) -> anyhow::Result<()> {
    let Some(reference_path) = args
        .reference_fasta
        .clone()
        .or_else(|| args.reference.clone())
    else {
        bail!("Usage: ntools calmd [-eubrAESQ] <aln.bam> <ref.fasta>");
    };

    let read_options = ReadOptions {
        threads: args.threads,
        reference: Some(reference_path.clone()),
    };

    let (mut reader, _) = AlignmentReader::open(&args.input, &read_options)?;
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
        format: if args.bam || args.uncompressed {
            OutputFormat::Bam
        } else {
            OutputFormat::Sam
        },
        level: output::parse_level(if args.uncompressed { 0 } else { 6 })?,
        threads: args.threads,
        reference: None,
    };

    let mut writer = AlignmentWriter::create(None, &write_options, &header)?;
    writer.write_header(&header)?;

    let options = Options {
        use_equal: args.use_equal,
        drop_tags: args.drop_tags,
        bin_qualities: args.bin_qualities,
        update: !args.no_update,
        max_nm: args.max_nm,
        quiet: args.quiet,
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
    let mut encoder = noodles_bam::io::Writer::from(Vec::new());
    let mut batch: Vec<Job> = Vec::with_capacity(BATCH_SIZE);

    let context = Context {
        names: &names,
        options,
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
        use noodles_sam::alignment::io::Write as _;

        let raw = match &rec {
            Rec::Bam(r) => r.as_ref().to_vec(),
            _ => {
                encoder.get_mut().clear();
                encoder.write_alignment_record(&input_header, rec.as_dyn())?;
                encoder.get_ref()[size_of::<u32>()..].to_vec()
            }
        };

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
    options: Options,
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
        if i32::from(read.mapq()) > q {
            read.set_mapq(q as u8);
        }
    }

    let name = &context.names[read.tid as usize];
    let filled = fill_md(&mut read, &reference, name, &context.options);

    Ok((read.raw, filled))
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

/// `bam_fillmd1_core`. Returns `false` if the read has no sequence.
fn fill_md(read: &mut Read, reference: &[u8], ref_name: &[u8], options: &Options) -> bool {
    // The CIGAR is read while the sequence changes.
    let cigar = std::mem::take(&mut read.cigar);
    let filled = fill_md_with(read, &cigar, reference, ref_name, options);
    read.cigar = cigar;
    filled
}

fn fill_md_with(
    read: &mut Read,
    cigar: &[u32],
    reference: &[u8],
    ref_name: &[u8],
    options: &Options,
) -> bool {
    const UNMAPPED: u16 = 0x4;

    let ref_len = reference.len() as i64;
    let l_qseq = read.l_seq;

    if l_qseq == 0 {
        if !options.quiet {
            eprintln!(
                "[bam_fillmd1] no sequence in alignment record for '{}' at {}:{}, skipped",
                String::from_utf8_lossy(read.name()),
                String::from_utf8_lossy(ref_name),
                read.pos + 1
            );
        }

        return false;
    }

    let in_bounds = |rpos: i64| rpos < ref_len && reference[rpos as usize] != 0;

    let mut md: Vec<u8> = Vec::new();
    let mut matched = 0i32;
    let mut nm = 0i32;
    let (mut qpos, mut rpos) = (0usize, read.pos);

    'ops: for &op in cigar {
        let (kind, oplen) = (op & 0xf, (op >> 4) as usize);

        match kind {
            0 | 7 | 8 => {
                for j in 0..oplen {
                    let z = qpos + j;
                    let r = rpos + j as i64;

                    if !in_bounds(r) || z >= l_qseq {
                        break 'ops;
                    }

                    let c1 = read.base(z);
                    let c2 = SEQ_NT16_TABLE[usize::from(reference[r as usize])];

                    if (c1 == c2 && c1 != 15 && c2 != 15) || c1 == 0 {
                        if options.use_equal {
                            read.set_base(z, 0);
                        }

                        matched += 1;
                    } else {
                        push_int(&mut md, matched);
                        md.push(reference[r as usize].to_ascii_uppercase());
                        matched = 0;
                        nm += 1;
                    }
                }

                rpos += oplen as i64;
                qpos += oplen;
            }
            2 => {
                push_int(&mut md, matched);
                md.push(b'^');

                let mut j = 0;

                while j < oplen && in_bounds(rpos + j as i64) {
                    md.push(reference[(rpos + j as i64) as usize].to_ascii_uppercase());
                    j += 1;
                }

                matched = 0;
                rpos += j as i64;
                nm += j as i32;

                if j < oplen {
                    break 'ops;
                }
            }
            1 | 4 => {
                qpos += oplen;

                if kind == 1 {
                    nm += oplen as i32;
                }
            }
            3 => rpos += oplen as i64,
            _ => {}
        }
    }

    push_int(&mut md, matched);

    // -n: mask the matching bases of reads with too many mismatches.
    if options.max_nm > 0 && nm >= options.max_nm {
        let (mut qpos, mut rpos) = (0usize, read.pos);

        'mask: for &op in cigar {
            let (kind, oplen) = (op & 0xf, (op >> 4) as usize);

            match kind {
                0 | 7 | 8 => {
                    for j in 0..oplen {
                        let z = qpos + j;
                        let r = rpos + j as i64;

                        if !in_bounds(r) || z >= l_qseq {
                            break 'mask;
                        }

                        let c1 = read.base(z);
                        let c2 = SEQ_NT16_TABLE[usize::from(reference[r as usize])];

                        if (c1 == c2 && c1 != 15 && c2 != 15) || c1 == 0 {
                            read.set_base(z, 15);
                            read.qual_mut()[z] = 0;
                        }
                    }

                    rpos += oplen as i64;
                    qpos += oplen;
                }
                2 | 3 => rpos += oplen as i64,
                1 | 4 => qpos += oplen,
                _ => {}
            }
        }
    }

    let is_mapped = read.flag & UNMAPPED == 0;
    let name = String::from_utf8_lossy(read.name()).into_owned();

    if options.update && is_mapped {
        match read.find_tag(*b"NM") {
            None => append_nm(read, nm),
            Some((start, end)) => {
                let old = integer_value(&read.raw[start + 2..end]);

                if i64::from(nm) != old {
                    if !options.quiet {
                        eprintln!(
                            "[bam_fillmd1] different NM for read '{name}': {} -> {nm}",
                            old as i32
                        );
                    }

                    read.remove_tag(*b"NM");
                    append_nm(read, nm);
                }
            }
        }

        match read.find_tag(*b"MD") {
            None => read.append_string_tag(*b"MD", &md),
            Some((start, end)) => {
                let value = &read.raw[start + 3..end];
                let old = &value[..memchr::memchr(0, value).unwrap_or(value.len())];

                if !old.eq_ignore_ascii_case(&md) {
                    if !options.quiet {
                        eprintln!(
                            "[bam_fillmd1] different MD for read '{name}': '{}' -> '{}'",
                            String::from_utf8_lossy(old),
                            String::from_utf8_lossy(&md)
                        );
                    }

                    read.remove_tag(*b"MD");
                    read.append_string_tag(*b"MD", &md);
                }
            }
        }
    }

    if options.drop_tags {
        // bam_aux_drop_other: keep only the RG field.
        let rg = read.find_tag(*b"RG").map(|(s, e)| read.raw[s..e].to_vec());
        let data_start = crate::bam_raw::data_range(&read.raw).map_or(read.raw.len(), |r| r.start);
        read.raw.truncate(data_start);

        if let Some(rg) = rg {
            read.raw.extend_from_slice(&rg);
        }
    }

    if options.bin_qualities {
        for q in read.qual_mut() {
            if *q >= 3 {
                *q = *q / 10 * 10 + 7;
            }
        }
    }

    true
}

/// `bam_aux_append(b, "NM", 'i', 4, &nm)`
fn append_nm(read: &mut Read, nm: i32) {
    read.raw.extend_from_slice(b"NMi");
    read.raw.extend_from_slice(&nm.to_le_bytes());
}

/// `bam_aux2i` of a field's type and value: integers, or 0 for other types.
fn integer_value(field: &[u8]) -> i64 {
    let v = &field[1..];

    match field[0] {
        b'c' => i64::from(v[0].cast_signed()),
        b'C' => i64::from(v[0]),
        b's' => i64::from(i16::from_le_bytes([v[0], v[1]])),
        b'S' => i64::from(u16::from_le_bytes([v[0], v[1]])),
        b'i' => i64::from(i32::from_le_bytes([v[0], v[1], v[2], v[3]])),
        b'I' => i64::from(u32::from_le_bytes([v[0], v[1], v[2], v[3]])),
        _ => 0,
    }
}

fn push_int(dst: &mut Vec<u8>, n: i32) {
    let mut buf = [0; lexical_core::BUFFER_SIZE];
    dst.extend_from_slice(lexical_core::write(n, &mut buf));
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_integer_value() {
        assert_eq!(integer_value(&[b'C', 7]), 7);
        assert_eq!(integer_value(&[b'c', 0xff]), -1);
        assert_eq!(integer_value(&[b'i', 1, 1, 0, 0]), 257);
        assert_eq!(integer_value(&[b'Z', b'x', 0]), 0);
    }
}
