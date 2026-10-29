use std::{
    io::{self, BufRead, Read},
    vec,
};

use noodles_bgzf as bgzf;

use super::IndexedRecords;
use crate::binning_index::index::{Header, reference_sequence::bin::Chunk};

enum State {
    Seek,
    Read(bgzf::VirtualPosition),
    Done,
}

/// How a format reader classifies a record read from a query.
#[doc(hidden)]
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub enum RecordMatch {
    /// The record does not intersect the region.
    Skip,
    /// The record intersects the region.
    Include,
    /// The record starts past the region, as does every later record in its chunk.
    NextChunk,
}

/// A query reader.
///
/// This reader returns the uncompressed data between all the given chunks.
pub struct Query<'r, R> {
    reader: &'r mut R,
    chunks: vec::IntoIter<Chunk>,
    chunks_are_sorted: bool,
    state: State,
}

impl<'r, R> Query<'r, R> {
    /// Creates a query reader.
    pub fn new(reader: &'r mut R, chunks: Vec<Chunk>) -> Self {
        Self {
            reader,
            chunks_are_sorted: chunks_are_sorted(&chunks),
            chunks: chunks.into_iter(),
            state: State::Seek,
        }
    }

    /// Applies the match of a complete record just read, returning whether to yield it.
    ///
    /// A record past the region ends the query only when every later chunk follows the current
    /// one in the file. Otherwise, only the rest of the current chunk is skipped, preserving
    /// custom chunk order and overlapping chunks, including their duplicate records.
    #[doc(hidden)]
    pub fn accept(&mut self, record_match: RecordMatch) -> bool {
        match record_match {
            RecordMatch::Include => return true,
            RecordMatch::NextChunk if self.chunks_are_sorted => self.state = State::Done,
            RecordMatch::NextChunk => self.skip_current_chunk(),
            RecordMatch::Skip => {}
        }
        false
    }

    /// Discards unread bytes in the active chunk, retaining all later chunks.
    ///
    /// Call this only after consuming a complete record. If no chunk is active, this is a no-op.
    #[doc(hidden)]
    pub fn skip_current_chunk(&mut self) {
        if matches!(self.state, State::Read(_)) {
            self.state = State::Seek;
        }
    }
}

pub(crate) fn chunks_are_sorted(chunks: &[Chunk]) -> bool {
    chunks
        .windows(2)
        .all(|pair| pair[0].end() <= pair[1].start())
}

impl<R> Query<'_, R>
where
    R: bgzf::io::BufRead + bgzf::io::Seek,
{
    /// Creates an iterator that parses indexed records.
    pub fn indexed_records(self, header: &Header) -> IndexedRecords<Self> {
        IndexedRecords::new(self, header)
    }
}

impl<R> Read for Query<'_, R>
where
    R: bgzf::io::BufRead + bgzf::io::Seek,
{
    fn read(&mut self, buf: &mut [u8]) -> io::Result<usize> {
        let mut src = self.fill_buf()?;
        let amt = src.read(buf)?;
        self.consume(amt);
        Ok(amt)
    }
}

impl<R> BufRead for Query<'_, R>
where
    R: bgzf::io::BufRead + bgzf::io::Seek,
{
    fn fill_buf(&mut self) -> io::Result<&[u8]> {
        loop {
            match self.state {
                State::Seek => {
                    self.state = match self.chunks.next() {
                        Some(chunk) => {
                            self.reader.seek_to_virtual_position(chunk.start())?;
                            State::Read(chunk.end())
                        }
                        None => State::Done,
                    }
                }
                State::Read(chunk_end) => {
                    if self.reader.virtual_position() < chunk_end {
                        return self.reader.fill_buf();
                    } else {
                        self.state = State::Seek;
                    }
                }
                State::Done => return Ok(&[]),
            }
        }
    }

    fn consume(&mut self, amt: usize) {
        self.reader.consume(amt);
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn chunk(start: u64, end: u64) -> Chunk {
        Chunk::new(
            bgzf::VirtualPosition::from(start),
            bgzf::VirtualPosition::from(end),
        )
    }

    #[test]
    fn test_accept_ends_query_for_chunks_in_file_order() {
        let mut reader = bgzf::io::Reader::new(io::Cursor::new(Vec::new()));
        let mut query = Query::new(&mut reader, vec![chunk(0, 1), chunk(2, 3)]);

        assert!(query.accept(RecordMatch::Include));
        assert!(!query.accept(RecordMatch::Skip));
        assert!(!query.accept(RecordMatch::NextChunk));
        assert!(matches!(query.state, State::Done));
    }

    #[test]
    fn test_accept_skips_only_the_current_chunk_for_unordered_chunks() {
        let mut reader = bgzf::io::Reader::new(io::Cursor::new(Vec::new()));
        let mut query = Query::new(&mut reader, vec![chunk(2, 3), chunk(0, 1)]);
        query.state = State::Read(bgzf::VirtualPosition::from(3));

        assert!(!query.accept(RecordMatch::NextChunk));
        assert!(matches!(query.state, State::Seek));
    }
}
