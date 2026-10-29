use crate::container::block;

const LOW_READER_COUNT: usize = 64;
const LOW_READER_COUNT_I32: i32 = LOW_READER_COUNT as i32;

pub struct ExternalDataReaders<'c> {
    low_readers: [Option<&'c [u8]>; LOW_READER_COUNT],
    // Sorted by content ID. Tag blocks use large content IDs (e.g., 5063770 for `MD:Z`), and a
    // binary search over a slice's few such blocks is much cheaper than hashing on every read.
    high_readers: Vec<(block::ContentId, &'c [u8])>,
}

impl<'c> ExternalDataReaders<'c> {
    pub fn new() -> Self {
        Self {
            low_readers: [None; LOW_READER_COUNT],
            high_readers: Vec::new(),
        }
    }

    pub fn insert(&mut self, id: block::ContentId, reader: &'c [u8]) {
        match id {
            i @ 0..LOW_READER_COUNT_I32 => {
                self.low_readers[i as usize] = Some(reader);
            }
            _ => match self.high_readers.binary_search_by_key(&id, |(i, _)| *i) {
                Ok(i) => self.high_readers[i].1 = reader,
                Err(i) => self.high_readers.insert(i, (id, reader)),
            },
        }
    }

    pub fn get_mut(&mut self, id: &block::ContentId) -> Option<&mut &'c [u8]> {
        match *id {
            i @ 0..LOW_READER_COUNT_I32 => self.low_readers[i as usize].as_mut(),
            _ => self
                .high_readers
                .binary_search_by_key(id, |(i, _)| *i)
                .ok()
                .map(|i| &mut self.high_readers[i].1),
        }
    }
}

impl<'c> ExternalDataReaders<'c> {
    /// Returns a stable index of the reader for `id`, for [`Self::get_mut_by_index`].
    pub(crate) fn index_of(&self, id: &block::ContentId) -> Option<usize> {
        match *id {
            i @ 0..LOW_READER_COUNT_I32 => self.low_readers[i as usize].map(|_| i as usize),
            _ => self
                .high_readers
                .binary_search_by_key(id, |(i, _)| *i)
                .ok()
                .map(|i| LOW_READER_COUNT + i),
        }
    }

    /// Returns the reader at an index from [`Self::index_of`].
    ///
    /// # Panics
    ///
    /// Panics if the index is not one `index_of` returned for these readers.
    pub(crate) fn get_mut_by_index(&mut self, i: usize) -> &mut &'c [u8] {
        match i.checked_sub(LOW_READER_COUNT) {
            None => self.low_readers[i]
                .as_mut()
                .expect("index_of returns present readers"),
            Some(j) => &mut self.high_readers[j].1,
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_new() {
        let readers = ExternalDataReaders::new();
        assert!(readers.low_readers.iter().all(|reader| reader.is_none()));
        assert!(readers.high_readers.is_empty());
    }

    #[test]
    fn test_insert() {
        let mut readers = ExternalDataReaders::new();

        readers.insert(0, b"z");
        readers.insert(0, b"a");
        readers.insert(LOW_READER_COUNT_I32 - 1, b"b");
        readers.insert(LOW_READER_COUNT_I32, b"c");

        assert_eq!(readers.low_readers[0], Some(&b"a"[..]));
        assert_eq!(readers.low_readers[63], Some(&b"b"[..]));
        assert_eq!(readers.high_readers, [(64, &b"c"[..])]);
    }

    #[test]
    fn test_get_mut() {
        let mut readers = ExternalDataReaders::new();

        readers.insert(0, b"a");
        readers.insert(LOW_READER_COUNT_I32 - 1, b"b");
        readers.insert(LOW_READER_COUNT_I32, b"c");

        assert_eq!(readers.get_mut(&0), Some(&mut &b"a"[..]));
        assert_eq!(readers.get_mut(&63), Some(&mut &b"b"[..]));
        assert_eq!(readers.get_mut(&64), Some(&mut &b"c"[..]));
    }

    #[test]
    fn test_high_ids_are_kept_sorted() {
        let mut readers = ExternalDataReaders::new();
        let ids = [5063770, i32::MAX, -1, 64, i32::MIN, 4281155, 64];
        for (i, id) in ids.into_iter().enumerate() {
            readers.insert(id, &b"abcdefg"[i..=i]);
        }
        assert_eq!(readers.get_mut(&64), Some(&mut &b"g"[..]));
        assert_eq!(readers.get_mut(&i32::MIN), Some(&mut &b"e"[..]));
        assert_eq!(readers.get_mut(&5063770), Some(&mut &b"a"[..]));
        assert_eq!(readers.get_mut(&65), None);
    }
}
