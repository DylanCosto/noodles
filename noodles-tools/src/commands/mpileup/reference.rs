//! The reference sequence cache (`-f`).

use std::sync::{Arc, Mutex};

use crate::reference::IndexedFasta;

/// A cached reference sequence by name; `None` if the FASTA does not have it.
pub(super) type CachedSequence = (Vec<u8>, Option<Arc<Vec<u8>>>);

/// Whole reference sequences by name, the last few used (`mplp_get_ref` keeps three per input;
/// these are shared by all inputs and threads).
pub(super) struct ReferenceCache {
    pub(super) fasta: Option<IndexedFasta>,
    pub(super) cache: Mutex<Vec<CachedSequence>>,
}

impl ReferenceCache {
    pub(super) const SLOTS: usize = 4;

    pub(super) fn get(&self, name: &[u8]) -> anyhow::Result<Option<Arc<Vec<u8>>>> {
        let Some(fasta) = &self.fasta else {
            return Ok(None);
        };

        let mut cache = self.cache.lock().expect("lock poisoned");

        if let Some(i) = cache.iter().position(|(n, _)| n == name) {
            let entry = cache.remove(i);
            let sequence = entry.1.clone();
            cache.insert(0, entry);
            return Ok(sequence);
        }

        let mut buf = Vec::new();
        let sequence = fasta
            .fetch(name, 0, i64::MAX, &mut buf)?
            .then(|| Arc::new(buf));

        cache.insert(0, (name.to_vec(), sequence.clone()));
        cache.truncate(Self::SLOTS);

        Ok(sequence)
    }

    /// The sequence of a reference sequence ID of `names`, if any.
    pub(super) fn get_id(
        &self,
        names: &[Vec<u8>],
        tid: i32,
    ) -> anyhow::Result<Option<Arc<Vec<u8>>>> {
        match usize::try_from(tid).ok().and_then(|i| names.get(i)) {
            Some(name) => self.get(name),
            None => Ok(None),
        }
    }
}
