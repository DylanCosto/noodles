//! A string set with HTSlib's khash layout (`KHASH_SET_INIT_STR`/`KHASH_MAP_INIT_STR`), for
//! output in samtools' iteration order.

/// khash's `__ac_HASH_UPPER`.
const HASH_UPPER: f64 = 0.77;

/// Bucket states (khash's two flag bits per bucket).
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
enum Bucket {
    Empty,
    Deleted,
    Occupied,
}

/// An insert-only khash string set: keys are placed, and the table grows and rehashes in place,
/// exactly as `kh_put` and `kh_resize` do, so [`Self::iter`] visits keys in khash's order.
#[derive(Debug, Default)]
pub struct StringSet {
    keys: Vec<Vec<u8>>,
    buckets: Vec<Bucket>,
    size: usize,
    occupied: usize,
    upper_bound: usize,
}

impl StringSet {
    pub fn new() -> Self {
        Self::default()
    }

    /// Adds `key` if it is absent (`kh_get`, then `kh_put`), returning whether it was new. A key
    /// already present does not grow the table, as it would through `kh_put` alone.
    pub fn insert(&mut self, key: &[u8]) -> bool {
        if self.contains(key) {
            return false;
        }

        self.put(key);
        true
    }

    /// `kh_get`
    fn contains(&self, key: &[u8]) -> bool {
        if self.buckets.is_empty() {
            return false;
        }

        let mask = self.buckets.len() - 1;
        let mut i = hash(key) as usize & mask;
        let mut step = 0;

        while self.buckets[i] != Bucket::Empty {
            if self.keys[i] == key {
                return true;
            }

            step += 1;
            i = (i + step) & mask;
        }

        false
    }

    /// `kh_put` of a key that is absent.
    fn put(&mut self, key: &[u8]) {
        if self.occupied >= self.upper_bound {
            let n = self.buckets.len();

            if n > self.size << 1 {
                self.resize(n - 1);
            } else {
                self.resize(n + 1);
            }
        }

        let mask = self.buckets.len() - 1;
        let mut i = hash(key) as usize & mask;
        let mut step = 0;

        // Nothing is deleted outside a resize, so the key goes in the first empty bucket.
        while self.buckets[i] != Bucket::Empty {
            step += 1;
            i = (i + step) & mask;
        }

        self.keys[i] = key.to_vec();
        self.buckets[i] = Bucket::Occupied;
        self.size += 1;
        self.occupied += 1;
    }

    /// The keys in bucket order (`kh_begin` to `kh_end`).
    pub fn iter(&self) -> impl Iterator<Item = &[u8]> {
        self.buckets
            .iter()
            .zip(&self.keys)
            .filter(|(bucket, _)| **bucket == Bucket::Occupied)
            .map(|(_, key)| key.as_slice())
    }

    pub fn is_empty(&self) -> bool {
        self.size == 0
    }

    /// `kh_resize`: rehashes in place, kicking existing keys out of their new buckets.
    fn resize(&mut self, new_n: usize) {
        let new_n = new_n.next_power_of_two().max(4);

        if self.size >= (new_n as f64 * HASH_UPPER + 0.5) as usize {
            return;
        }

        let old_n = self.buckets.len();
        let mut new_occupied = vec![false; new_n];

        if old_n < new_n {
            self.keys.resize_with(new_n, Vec::new);
        }

        let new_mask = new_n - 1;

        for j in 0..old_n {
            if self.buckets[j] != Bucket::Occupied {
                continue;
            }

            let mut key = std::mem::take(&mut self.keys[j]);
            self.buckets[j] = Bucket::Deleted;

            loop {
                let mut i = hash(&key) as usize & new_mask;
                let mut step = 0;

                while new_occupied[i] {
                    step += 1;
                    i = (i + step) & new_mask;
                }

                new_occupied[i] = true;

                if i < old_n && self.buckets[i] == Bucket::Occupied {
                    std::mem::swap(&mut key, &mut self.keys[i]);
                    self.buckets[i] = Bucket::Deleted;
                } else {
                    self.keys[i] = key;
                    break;
                }
            }
        }

        self.keys.truncate(new_n);
        self.buckets = new_occupied
            .into_iter()
            .map(|occupied| {
                if occupied {
                    Bucket::Occupied
                } else {
                    Bucket::Empty
                }
            })
            .collect();
        self.occupied = self.size;
        self.upper_bound = (new_n as f64 * HASH_UPPER + 0.5) as usize;
    }
}

/// khash's `kh_str_hash_func` (FNV-1a in HTSlib 1.24).
fn hash(s: &[u8]) -> u32 {
    const OFFSET_BASIS: u32 = 2_166_136_261;
    const PRIME: u32 = 16_777_619;

    s.iter()
        .fold(OFFSET_BASIS, |h, &b| (h ^ u32::from(b)).wrapping_mul(PRIME))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_insert() {
        let mut set = StringSet::new();
        assert!(set.is_empty());
        assert!(set.insert(b"a"));
        assert!(!set.insert(b"a"));
        assert!(set.insert(b"b"));

        let mut keys: Vec<_> = set.iter().collect();
        keys.sort_unstable();
        assert_eq!(keys, [b"a", b"b"]);
    }

    #[test]
    fn test_grows() {
        let mut set = StringSet::new();

        for i in 0..100 {
            assert!(set.insert(format!("s{i}").as_bytes()));
        }

        assert_eq!(set.iter().count(), 100);
        assert!(!set.insert(b"s42"));
    }
}
