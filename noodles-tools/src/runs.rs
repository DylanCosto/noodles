//! Sorted runs of raw BAM records: in-memory chunks, spill files, and k-way merging.
//!
//! Every merge breaks ties by source order, and sources are numbered in input order, so sorting
//! is stable, as in samtools.

use std::{
    cmp::Ordering,
    collections::BinaryHeap,
    fs::File,
    io::{self, BufReader, Read, Write},
    path::{Path, PathBuf},
    thread,
};

use noodles_bgzf::io::writer::CompressionLevel;

use crate::{
    bgzf_stream::{BgzfReader, BgzfWriter},
    sort_order::{Order, coordinate_key},
};

/// The fixed-size prefix of every encoded BAM record, up to the read name.
const MIN_RECORD_LEN: usize = 32;

/// Reads one `block_size`-prefixed raw record into `buf`. Returns `false` at end of stream.
pub fn read_raw_record<R: Read>(reader: &mut R, buf: &mut Vec<u8>) -> io::Result<bool> {
    let mut len = [0; 4];

    match reader.read(&mut len[..1])? {
        0 => return Ok(false),
        _ => reader.read_exact(&mut len[1..])?,
    }

    let len = usize::try_from(u32::from_le_bytes(len))
        .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))?;

    buf.resize(len, 0);
    reader.read_exact(buf)?;

    // Keys are read from fixed offsets and the name, so those must be in bounds.
    if len < MIN_RECORD_LEN || MIN_RECORD_LEN + usize::from(buf[8]) > len {
        return Err(io::Error::new(
            io::ErrorKind::InvalidData,
            "truncated BAM record",
        ));
    }

    Ok(true)
}

pub fn write_raw_record<W: Write>(writer: &mut W, record: &[u8]) -> io::Result<()> {
    let len =
        u32::try_from(record.len()).map_err(|e| io::Error::new(io::ErrorKind::InvalidInput, e))?;
    writer.write_all(&len.to_le_bytes())?;
    writer.write_all(record)
}

#[derive(Clone, Copy, Debug)]
struct Entry {
    /// The coordinate key, precomputed so coordinate sorting compares integers.
    key: u64,
    start: usize,
    end: usize,
}

/// An in-memory batch of raw records in one contiguous buffer.
#[derive(Debug, Default)]
pub struct Chunk {
    data: Vec<u8>,
    entries: Vec<Entry>,
}

impl Chunk {
    pub fn push(&mut self, record: &[u8]) {
        let start = self.data.len();
        self.data.extend_from_slice(record);
        self.entries.push(Entry {
            key: coordinate_key(record),
            start,
            end: self.data.len(),
        });
    }

    pub fn is_empty(&self) -> bool {
        self.entries.is_empty()
    }

    /// Approximate memory use, for comparison against the sort memory limit.
    pub fn size(&self) -> usize {
        self.data.len() + self.entries.len() * size_of::<Entry>()
    }

    pub fn clear(&mut self) {
        self.data.clear();
        self.entries.clear();
    }

    /// Sorts the chunk with `threads` additional threads and calls `f` with each record in order.
    pub fn drain_sorted<F>(&mut self, order: Order, threads: usize, mut f: F) -> io::Result<()>
    where
        F: FnMut(&[u8]) -> io::Result<()>,
    {
        let data = &self.data;
        let parts = threads + 1;
        let part_len = self.entries.len().div_ceil(parts).max(1);

        let compare = |a: &Entry, b: &Entry| match order {
            Order::Coordinate => a.key.cmp(&b.key),
            _ => order.compare(&data[a.start..a.end], &data[b.start..b.end]),
        };

        // Stable sorts of contiguous parts, merged with ties going to the earlier part, give a
        // stable sort of the whole chunk.
        thread::scope(|scope| {
            for part in self.entries.chunks_mut(part_len) {
                scope.spawn(move || part.sort_by(compare));
            }
        });

        let parts: Vec<&[Entry]> = self.entries.chunks(part_len).collect();
        let mut heap = BinaryHeap::with_capacity(parts.len());

        for (source, part) in parts.iter().enumerate() {
            if let Some(&entry) = part.first() {
                heap.push(HeapEntry {
                    entry,
                    source,
                    index: 0,
                    compare: &compare,
                });
            }
        }

        while let Some(top) = heap.pop() {
            f(&data[top.entry.start..top.entry.end])?;

            let index = top.index + 1;

            if let Some(&entry) = parts[top.source].get(index) {
                heap.push(HeapEntry {
                    entry,
                    index,
                    ..top
                });
            }
        }

        self.clear();

        Ok(())
    }
}

struct HeapEntry<'c, C> {
    entry: Entry,
    source: usize,
    index: usize,
    compare: &'c C,
}

impl<C: Fn(&Entry, &Entry) -> Ordering> Ord for HeapEntry<'_, C> {
    fn cmp(&self, other: &Self) -> Ordering {
        // `BinaryHeap` is a max-heap, so the smallest record must compare greatest.
        (self.compare)(&self.entry, &other.entry)
            .then(self.source.cmp(&other.source))
            .reverse()
    }
}

impl<C: Fn(&Entry, &Entry) -> Ordering> PartialOrd for HeapEntry<'_, C> {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}

impl<C: Fn(&Entry, &Entry) -> Ordering> PartialEq for HeapEntry<'_, C> {
    fn eq(&self, other: &Self) -> bool {
        self.cmp(other) == Ordering::Equal
    }
}

impl<C: Fn(&Entry, &Entry) -> Ordering> Eq for HeapEntry<'_, C> {}

/// Writes a chunk as a sorted spill file of raw records, and returns its path.
pub fn spill(
    chunk: &mut Chunk,
    order: Order,
    dir: &Path,
    index: usize,
    threads: usize,
) -> io::Result<PathBuf> {
    let path = dir.join(format!("{index:04}.bam"));
    let file = File::create(&path)?;
    let mut writer = BgzfWriter::new(file, CompressionLevel::FAST, threads);

    chunk.drain_sorted(order, threads, |record| {
        write_raw_record(&mut writer, record)
    })?;

    writer.finish()?;

    Ok(path)
}

/// A source of sorted raw records for a k-way merge.
pub trait RecordSource {
    /// Reads the next record into `buf`, returning `false` when exhausted.
    fn next_record(&mut self, buf: &mut Vec<u8>) -> io::Result<bool>;
}

/// A spill file written by [`spill`].
pub struct SpillFile(BgzfReader<BufReader<File>>);

impl SpillFile {
    pub fn open(path: &Path, threads: usize) -> io::Result<Self> {
        let file = File::open(path)?;
        Ok(Self(BgzfReader::new(BufReader::new(file), threads)))
    }
}

impl RecordSource for SpillFile {
    fn next_record(&mut self, buf: &mut Vec<u8>) -> io::Result<bool> {
        read_raw_record(&mut self.0, buf)
    }
}

struct Head {
    record: Vec<u8>,
    source: usize,
    order: Order,
}

impl Ord for Head {
    fn cmp(&self, other: &Self) -> Ordering {
        self.order
            .compare(&self.record, &other.record)
            .then(self.source.cmp(&other.source))
            .reverse()
    }
}

impl PartialOrd for Head {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}

impl PartialEq for Head {
    fn eq(&self, other: &Self) -> bool {
        self.cmp(other) == Ordering::Equal
    }
}

impl Eq for Head {}

/// Merges sorted sources, calling `f` with each record in order. Ties go to the earlier source.
pub fn merge<S, F>(sources: &mut [S], order: Order, mut f: F) -> anyhow::Result<()>
where
    S: RecordSource,
    F: FnMut(&[u8]) -> anyhow::Result<()>,
{
    let mut heap = BinaryHeap::with_capacity(sources.len());

    for (source, s) in sources.iter_mut().enumerate() {
        let mut record = Vec::new();

        if s.next_record(&mut record)? {
            heap.push(Head {
                record,
                source,
                order,
            });
        }
    }

    while let Some(mut head) = heap.pop() {
        f(&head.record)?;

        if sources[head.source].next_record(&mut head.record)? {
            heap.push(head);
        }
    }

    Ok(())
}
