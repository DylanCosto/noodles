use std::io::{self, BufRead, Seek};

use noodles_core::{Position, Region, region::Interval};

use crate::{Record, record::Sequence, repository::Adapter};

/// An indexed reader adapter.
pub struct IndexedReader<R> {
    reader: crate::io::IndexedReader<R>,
}

impl<R> IndexedReader<R>
where
    R: BufRead + Seek,
{
    /// Creates an indexed reader adapter.
    ///
    /// # Examples
    ///
    /// ```
    /// # use std::io;
    /// use noodles_fasta::{self as fasta, fai, repository::adapters::IndexedReader};
    /// let index = fai::Index::default();
    /// let reader = fasta::io::IndexedReader::new(io::empty(), index);
    /// let adapter = IndexedReader::new(reader);
    /// ```
    pub fn new(reader: crate::io::IndexedReader<R>) -> Self {
        Self { reader }
    }
}

impl<R> IndexedReader<R>
where
    R: BufRead,
{
    fn index_record(&self, name: &[u8]) -> Option<&crate::fai::Record> {
        self.reader
            .index()
            .as_ref()
            .iter()
            .find(|record| record.name() == name)
    }
}

impl<R> Adapter for IndexedReader<R>
where
    R: BufRead + Seek + Send + Sync,
{
    /// Returns `None` when the index has no sequence named `name`.
    fn get(&mut self, name: &[u8]) -> Option<io::Result<Record>> {
        self.index_record(name)?;
        let region = Region::new(name, ..);
        Some(self.reader.query(&region))
    }

    /// Returns `None` when the index has no sequence named `name`.
    fn get_range(&mut self, name: &[u8], interval: Interval) -> Option<io::Result<Sequence>> {
        let length = self.index_record(name)?.length();

        let Ok(length) = usize::try_from(length) else {
            let region = Region::new(name, interval);
            return Some(self.reader.query(&region).map(Record::into_sequence));
        };

        // Stop at the indexed length. A query that starts past it would seek into the next record,
        // e.g., `sq0:13-14` of a 10-base `sq0` would return bases of the next sequence.
        let start = interval.start().unwrap_or(Position::MIN);
        let end = interval
            .end()
            .map(usize::from)
            .unwrap_or(length)
            .min(length);

        let Some(end) = Position::new(end).filter(|&end| end >= start) else {
            return Some(Ok(Sequence::default()));
        };

        let region = Region::new(name, start..=end);
        Some(self.reader.query(&region).map(Record::into_sequence))
    }
}
