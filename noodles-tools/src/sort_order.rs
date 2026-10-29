//! samtools-compatible record orderings over raw (encoded) BAM records.
//!
//! Keys are read directly from the fixed-layout BAM record fields, so sorting and merging never
//! decode records.

use std::cmp::Ordering;

use crate::bam_raw::{field, le_u16, le_u32};

const REVERSE_COMPLEMENTED: u16 = 0x10;
/// READ1 and READ2 bits; samtools orders same-name records by these.
const SEGMENT_BITS: u16 = 0xc0;

/// A record order.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub enum Order {
    /// By reference sequence, position, then strand; unmapped records with no reference last.
    Coordinate,
    /// By read name, comparing digit runs numerically (`samtools sort -n`).
    NameNatural,
    /// By read name, bytewise (`samtools sort -N`).
    NameAscii,
    /// Grouped by read name in hash order (`samtools collate`): records with the same name are
    /// adjacent, in input order.
    Collate,
}

impl Order {
    /// The `@HD` `SO` and `SS` values samtools writes for this order.
    pub fn header_values(self) -> (&'static str, Option<&'static str>) {
        match self {
            Self::Coordinate => ("coordinate", None),
            Self::NameNatural => ("queryname", Some("queryname:natural")),
            Self::NameAscii => ("queryname", Some("queryname:lexicographical")),
            Self::Collate => ("unsorted", None),
        }
    }

    pub fn compare(self, a: &[u8], b: &[u8]) -> Ordering {
        match self {
            Self::Coordinate => coordinate_key(a).cmp(&coordinate_key(b)),
            Self::NameNatural => strnum_cmp(name(a), name(b))
                .then_with(|| (flags(a) & SEGMENT_BITS).cmp(&(flags(b) & SEGMENT_BITS))),
            Self::NameAscii => name(a)
                .cmp(name(b))
                .then_with(|| (flags(a) & SEGMENT_BITS).cmp(&(flags(b) & SEGMENT_BITS))),
            Self::Collate => name_hash(a)
                .cmp(&name_hash(b))
                .then_with(|| name(a).cmp(name(b))),
        }
    }
}

/// The samtools coordinate key: `tid` (as unsigned, so `-1` sorts last), `pos + 1`, then strand.
pub fn coordinate_key(record: &[u8]) -> u64 {
    let tid = u64::from(le_u32(&record[field::REFERENCE_SEQUENCE_ID]));
    // `pos` is at least -1, so `pos + 1` is nonnegative and fits in 31 bits.
    let pos = u64::from(le_u32(&record[field::POSITION]).wrapping_add(1));
    let is_reverse = u64::from(flags(record) & REVERSE_COMPLEMENTED != 0);
    (tid << 32) | (pos << 1) | is_reverse
}

/// FNV-1a over the read name, the grouping key for [`Order::Collate`].
pub fn name_hash(record: &[u8]) -> u64 {
    const OFFSET_BASIS: u64 = 0xcbf2_9ce4_8422_2325;
    const PRIME: u64 = 0x0100_0000_01b3;

    name(record)
        .iter()
        .fold(OFFSET_BASIS, |h, &b| (h ^ u64::from(b)).wrapping_mul(PRIME))
}

/// Returns the read name without its NUL terminator.
pub fn name(record: &[u8]) -> &[u8] {
    let len = usize::from(record[field::NAME_LENGTH]);
    let end = field::NAME + len.saturating_sub(1);
    &record[field::NAME..end]
}

pub fn flags(record: &[u8]) -> u16 {
    le_u16(&record[field::FLAGS])
}

/// Rewrites the reference and mate reference sequence IDs through `map` (index: input ID).
/// IDs of `-1` (none) are unchanged.
pub fn remap_reference_sequence_ids(record: &mut [u8], map: &[i32]) {
    for range in [
        field::REFERENCE_SEQUENCE_ID,
        field::MATE_REFERENCE_SEQUENCE_ID,
    ] {
        let id = le_u32(&record[range.clone()]).cast_signed();

        if let Ok(i) = usize::try_from(id) {
            record[range].copy_from_slice(&map[i].to_le_bytes());
        }
    }
}

/// samtools' natural string comparison (`strnum_cmp`): digit runs compare numerically, ignoring
/// leading zeros, and when two names are otherwise equal, the one with fewer leading zeros sorts
/// first.
pub fn strnum_cmp(a: &[u8], b: &[u8]) -> Ordering {
    // A translation of the C, where reading past the end yields the NUL terminator.
    let at = |s: &[u8], i: usize| s.get(i).copied().unwrap_or(0);
    let (mut i, mut j) = (0, 0);

    while at(a, i) != 0 && at(b, j) != 0 {
        if !at(a, i).is_ascii_digit() || !at(b, j).is_ascii_digit() {
            if at(a, i) != at(b, j) {
                return at(a, i).cmp(&at(b, j));
            }

            i += 1;
            j += 1;
        } else {
            while at(a, i) == b'0' {
                i += 1;
            }

            while at(b, j) == b'0' {
                j += 1;
            }

            while at(a, i).is_ascii_digit() && at(a, i) == at(b, j) {
                i += 1;
                j += 1;
            }

            let difference = at(a, i).cmp(&at(b, j));

            while at(a, i).is_ascii_digit() && at(b, j).is_ascii_digit() {
                i += 1;
                j += 1;
            }

            if at(a, i).is_ascii_digit() {
                return Ordering::Greater;
            } else if at(b, j).is_ascii_digit() {
                return Ordering::Less;
            } else if difference != Ordering::Equal {
                return difference;
            }
        }
    }

    if at(a, i) == at(b, j) {
        i.cmp(&j)
    } else if at(a, i) != 0 {
        Ordering::Greater
    } else {
        Ordering::Less
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_strnum_cmp() {
        use Ordering::{Equal, Greater, Less};

        let cases: [(&str, &str, Ordering); 10] = [
            ("a", "b", Less),
            ("r2", "r10", Less),
            ("r10", "r2", Greater),
            ("r10", "r10", Equal),
            ("r01", "r1", Greater),
            ("r1", "r01", Less),
            ("r1a", "r1b", Less),
            ("r12x", "r13a", Less),
            ("abc", "ab", Greater),
            ("x9:1", "x10:0", Less),
        ];

        for (a, b, expected) in cases {
            assert_eq!(
                strnum_cmp(a.as_bytes(), b.as_bytes()),
                expected,
                "{a:?} vs {b:?}"
            );
        }
    }
}
