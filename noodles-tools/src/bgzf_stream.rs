//! Single- and multithreaded BGZF streams behind one type.

use std::{
    io::{self, BufRead, Read, Seek, SeekFrom, Write},
    num::NonZero,
};

use noodles_bgzf::{self as bgzf, VirtualPosition, gzi, io::writer::CompressionLevel};

/// A BGZF reader that decompresses on the calling thread or on a worker pool.
pub enum BgzfReader<R> {
    Single(bgzf::io::Reader<R>),
    Multi(bgzf::io::MultithreadedReader<R>),
}

impl<R> BgzfReader<R>
where
    R: Read + Send + 'static,
{
    /// Creates a reader using `threads` additional decompression workers (`0` is single-threaded).
    pub fn new(inner: R, threads: usize) -> Self {
        match NonZero::new(threads) {
            Some(n) => Self::Multi(bgzf::io::MultithreadedReader::with_worker_count(n, inner)),
            None => Self::Single(bgzf::io::Reader::new(inner)),
        }
    }
}

impl<R> Read for BgzfReader<R>
where
    R: Read + Send + 'static,
{
    fn read(&mut self, buf: &mut [u8]) -> io::Result<usize> {
        match self {
            Self::Single(r) => r.read(buf),
            Self::Multi(r) => r.read(buf),
        }
    }

    fn read_exact(&mut self, buf: &mut [u8]) -> io::Result<()> {
        match self {
            Self::Single(r) => r.read_exact(buf),
            Self::Multi(r) => r.read_exact(buf),
        }
    }
}

impl<R> BufRead for BgzfReader<R>
where
    R: Read + Send + 'static,
{
    fn fill_buf(&mut self) -> io::Result<&[u8]> {
        match self {
            Self::Single(r) => r.fill_buf(),
            Self::Multi(r) => r.fill_buf(),
        }
    }

    fn consume(&mut self, amt: usize) {
        match self {
            Self::Single(r) => r.consume(amt),
            Self::Multi(r) => r.consume(amt),
        }
    }
}

impl<R> bgzf::io::Read for BgzfReader<R>
where
    R: Read + Send + 'static,
{
    fn virtual_position(&self) -> VirtualPosition {
        match self {
            Self::Single(r) => r.virtual_position(),
            Self::Multi(r) => r.virtual_position(),
        }
    }
}

impl<R> bgzf::io::BufRead for BgzfReader<R> where R: Read + Send + 'static {}

impl<R> bgzf::io::Seek for BgzfReader<R>
where
    R: Read + Seek + Send + 'static,
{
    fn seek_to_virtual_position(&mut self, pos: VirtualPosition) -> io::Result<VirtualPosition> {
        match self {
            Self::Single(r) => r.seek_to_virtual_position(pos),
            Self::Multi(r) => r.seek_to_virtual_position(pos),
        }
    }

    fn seek_with_index(&mut self, index: &gzi::Index, pos: SeekFrom) -> io::Result<u64> {
        match self {
            Self::Single(r) => r.seek_with_index(index, pos),
            Self::Multi(r) => r.seek_with_index(index, pos),
        }
    }
}

/// A BGZF writer that compresses on the calling thread or on a worker pool.
pub enum BgzfWriter<W>
where
    W: Write + Send + 'static,
{
    Single(bgzf::io::Writer<W>),
    Multi(bgzf::io::MultithreadedWriter<W>),
}

impl<W> BgzfWriter<W>
where
    W: Write + Send + 'static,
{
    /// Creates a writer using `threads` additional compression workers (`0` is single-threaded).
    pub fn new(inner: W, level: CompressionLevel, threads: usize) -> Self {
        match NonZero::new(threads) {
            Some(n) => Self::Multi(
                bgzf::io::multithreaded_writer::Builder::default()
                    .set_compression_level(level)
                    .set_worker_count(n)
                    .build_from_writer(inner),
            ),
            None => Self::Single(
                bgzf::io::writer::Builder::default()
                    .set_compression_level(level)
                    .build_from_writer(inner),
            ),
        }
    }

    /// Flushes remaining blocks, writes the EOF marker, and returns the inner writer.
    pub fn finish(self) -> io::Result<W> {
        match self {
            Self::Single(w) => w.finish(),
            Self::Multi(mut w) => w.finish(),
        }
    }
}

impl<W> Write for BgzfWriter<W>
where
    W: Write + Send + 'static,
{
    fn write(&mut self, buf: &[u8]) -> io::Result<usize> {
        match self {
            Self::Single(w) => w.write(buf),
            Self::Multi(w) => w.write(buf),
        }
    }

    fn write_all(&mut self, buf: &[u8]) -> io::Result<()> {
        match self {
            Self::Single(w) => w.write_all(buf),
            Self::Multi(w) => w.write_all(buf),
        }
    }

    fn flush(&mut self) -> io::Result<()> {
        match self {
            Self::Single(w) => w.flush(),
            Self::Multi(w) => w.flush(),
        }
    }
}
