//! khash's string hash.

/// khash's `__ac_X31_hash_string`.
pub(crate) fn x31_hash(s: &[u8]) -> u32 {
    let Some((&first, rest)) = s.split_first() else {
        return 0;
    };

    rest.iter().fold(u32::from(first), |h, &c| {
        (h << 5).wrapping_sub(h).wrapping_add(u32::from(c))
    })
}

/// khash's `__ac_Wang_hash`.
pub(crate) fn wang_hash(mut key: u32) -> u32 {
    key = key.wrapping_add(!(key << 15));
    key ^= key >> 10;
    key = key.wrapping_add(key << 3);
    key ^= key >> 6;
    key = key.wrapping_add(!(key << 11));
    key ^= key >> 16;
    key
}
