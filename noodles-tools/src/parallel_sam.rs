//! SAM output formatted on a worker pool.
//!
//! Records arrive as raw BAM on the calling thread and are batched; workers format batches into
//! SAM text, and a writer thread writes the text in batch order, so output is identical to
//! single-threaded formatting.

use std::{
    collections::BTreeMap,
    io::{self, Write},
    sync::Arc,
    thread::{self, JoinHandle},
};

use crossbeam_channel::{Receiver, Sender};
use noodles_bam as bam;
use noodles_sam::{self as sam, alignment::io::Write as _};

use crate::runs::write_raw_record;

/// Raw record bytes per batch. Large enough to amortize channel traffic, small enough that
/// batches spread evenly across workers.
const BATCH_CAPACITY: usize = 1 << 20;

/// Batches in flight per worker before the reader blocks, which bounds memory.
const QUEUE_DEPTH_PER_WORKER: usize = 4;

/// A unit of work: raw BAM records to format, or text (e.g., the header) to pass through.
enum Job {
    Records(Vec<u8>),
    Text(Vec<u8>),
}

type Batch = (u64, Job);
type Formatted = (u64, io::Result<Vec<u8>>);

pub struct ParallelSamWriter {
    batch: Vec<u8>,
    next_id: u64,
    work: Option<Sender<Batch>>,
    workers: Vec<JoinHandle<()>>,
    writer: Option<JoinHandle<io::Result<()>>>,
}

impl ParallelSamWriter {
    /// Starts `workers` formatting threads and a writer thread that writes to `sink`.
    pub fn new<W>(sink: W, header: &sam::Header, workers: usize) -> Self
    where
        W: Write + Send + 'static,
    {
        let workers = workers.max(1);
        let header = Arc::new(header.clone());

        let (work_tx, work_rx) =
            crossbeam_channel::bounded::<Batch>(workers * QUEUE_DEPTH_PER_WORKER);
        let (done_tx, done_rx) = crossbeam_channel::unbounded::<Formatted>();

        let worker_handles = (0..workers)
            .map(|_| {
                let work_rx = work_rx.clone();
                let done_tx = done_tx.clone();
                let header = Arc::clone(&header);
                thread::spawn(move || format_batches(&header, &work_rx, &done_tx))
            })
            .collect();

        drop(done_tx);

        let writer = thread::spawn(move || write_in_order(sink, &done_rx));

        Self {
            batch: Vec::with_capacity(BATCH_CAPACITY),
            next_id: 0,
            work: Some(work_tx),
            workers: worker_handles,
            writer: Some(writer),
        }
    }

    /// Queues already formatted text, such as the SAM header, in output order.
    pub fn write_text(&mut self, text: Vec<u8>) -> io::Result<()> {
        self.send_batch()?;
        self.send(Job::Text(text))
    }

    /// Queues an encoded BAM record (without `block_size`).
    pub fn push(&mut self, record: &[u8]) -> io::Result<()> {
        write_raw_record(&mut self.batch, record)?;

        if self.batch.len() >= BATCH_CAPACITY {
            self.send_batch()?;
        }

        Ok(())
    }

    fn send_batch(&mut self) -> io::Result<()> {
        if self.batch.is_empty() {
            return Ok(());
        }

        let batch = std::mem::replace(&mut self.batch, Vec::with_capacity(BATCH_CAPACITY));
        self.send(Job::Records(batch))
    }

    fn send(&mut self, job: Job) -> io::Result<()> {
        let id = self.next_id;
        self.next_id += 1;

        let work = self
            .work
            .as_ref()
            .expect("work channel is open until finish");

        if work.send((id, job)).is_ok() {
            return Ok(());
        }

        // Every worker stopped, which happens only after the writer thread failed. Report its
        // error rather than the closed channel.
        match self.writer.take().map(JoinHandle::join) {
            Some(Ok(Err(e))) => Err(e),
            _ => Err(io::Error::other("SAM writer stopped")),
        }
    }

    /// Formats remaining records, waits for all output to be written, and reports any error.
    pub fn finish(mut self) -> io::Result<()> {
        let sent = self.send_batch();

        // Closing the work channel ends the workers, which closes the output channel.
        drop(self.work.take());

        for handle in self.workers.drain(..) {
            handle
                .join()
                .map_err(|_| io::Error::other("SAM formatting thread panicked"))?;
        }

        let written = match self.writer.take() {
            Some(handle) => handle
                .join()
                .map_err(|_| io::Error::other("SAM writer thread panicked"))?,
            // `send` already joined the writer and returned its error.
            None => Ok(()),
        };

        written.and(sent)
    }
}

fn format_batches(header: &sam::Header, work: &Receiver<Batch>, done: &Sender<Formatted>) {
    let mut record = bam::Record::default();

    for (id, job) in work {
        let result = match job {
            Job::Records(batch) => format_batch(header, &batch, &mut record),
            Job::Text(text) => Ok(text),
        };

        if done.send((id, result)).is_err() {
            return;
        }
    }
}

fn format_batch(
    header: &sam::Header,
    batch: &[u8],
    record: &mut bam::Record,
) -> io::Result<Vec<u8>> {
    // SAM text is usually about twice the size of the raw BAM record.
    let mut writer = sam::io::Writer::new(Vec::with_capacity(batch.len() * 2));
    let mut reader = bam::io::Reader::from(batch);

    while reader.read_record(record)? != 0 {
        writer.write_alignment_record(header, &*record)?;
    }

    Ok(writer.into_inner())
}

fn write_in_order<W: Write>(mut sink: W, done: &Receiver<Formatted>) -> io::Result<()> {
    let mut pending = BTreeMap::new();
    let mut next = 0;

    for (id, result) in done {
        pending.insert(id, result?);

        while let Some(text) = pending.remove(&next) {
            sink.write_all(&text)?;
            next += 1;
        }
    }

    sink.flush()
}
