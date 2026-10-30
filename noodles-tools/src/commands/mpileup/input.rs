//! Reading inputs: records are filtered and prepared (BAQ, MAPQ capping) in batches on worker threads.

use std::{fs, path::Path, thread};

use anyhow::{Context, bail};
use noodles_pileup::{Pileup, Read, baq};
use noodles_sam as sam;

use super::{Args, BATCH_SIZE, Conf, reference::ReferenceCache};
use crate::alignment::{AlignmentReader, Queries, RawEncoder, ReadOptions, Rec};

/// `mplp_func`: filters a read and adjusts its qualities, or returns `None` to skip it.
pub(super) fn prepare(
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
pub(super) fn end_pos(read: &Read) -> i64 {
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
pub(super) fn prepare_batch(
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
pub(super) struct Source {
    pub(super) pileup: Pileup,
    pub(super) batches: crossbeam_channel::Receiver<anyhow::Result<Vec<Read>>>,
    pub(super) pending: std::vec::IntoIter<Read>,
}

impl Source {
    /// The next prepared read, or `None` at the end of the input.
    pub(super) fn next_read(&mut self) -> anyhow::Result<Option<Read>> {
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
    pub(super) fn auto(&mut self) -> anyhow::Result<bool> {
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

/// The input files and their index files (`-b`, `-X`).
pub(super) fn input_files(args: &Args) -> anyhow::Result<(Vec<String>, Option<Vec<String>>)> {
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
pub(super) fn read_input(
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
    let mut encoder = RawEncoder::new();

    // An `Err` with no cause signals that the consumer hung up.
    let send = |batch: &mut Vec<Vec<u8>>| -> anyhow::Result<()> {
        let reads = prepare_batch(std::mem::take(batch), conf, &names, references, workers)?;

        if tx.send(Ok(reads)).is_err() {
            bail!(Disconnected);
        }

        Ok(())
    };

    let mut visit = |rec: Rec<'_>| -> anyhow::Result<()> {
        let raw = rec.encode_raw(&header, &mut encoder)?.to_vec();

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
pub(super) struct Disconnected;

impl std::fmt::Display for Disconnected {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.write_str("the pileup stopped reading")
    }
}

impl std::error::Error for Disconnected {}
