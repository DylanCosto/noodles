//! samtools-compatible record orderings over raw (encoded) BAM records.
//!
//! Keys are read directly from the fixed-layout BAM record fields, so sorting and merging never
//! decode records.

use std::cmp::Ordering;

use crate::{
    bam_raw::{self, field, le_u16, le_u32},
    template_coordinate,
};

const REVERSE_COMPLEMENTED: u16 = 0x10;

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
    /// By the value of an auxiliary field, with records lacking it first, then by `then`
    /// (`samtools sort -t`).
    Tag { tag: [u8; 2], then: TagTieBreak },
    /// By the coordinates of both ends of the template, library, `CB`, `MI`, name, and end
    /// (`samtools sort --template-coordinate`). Records must be keyed (see
    /// [`crate::template_coordinate`]).
    TemplateCoordinate,
}

/// The order of records with equal [`Order::Tag`] values.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub enum TagTieBreak {
    Coordinate,
    NameNatural,
    NameAscii,
}

impl TagTieBreak {
    fn order(self) -> Order {
        match self {
            Self::Coordinate => Order::Coordinate,
            Self::NameNatural => Order::NameNatural,
            Self::NameAscii => Order::NameAscii,
        }
    }
}

impl Order {
    /// The `@HD` `SO` and `SS` values samtools writes for this order.
    pub fn header_values(self) -> (&'static str, Option<String>) {
        match self {
            Self::Coordinate => ("coordinate", None),
            Self::NameNatural => ("queryname", Some("queryname:natural".into())),
            Self::NameAscii => ("queryname", Some("queryname:lexicographical".into())),
            Self::Collate => ("unsorted", None),
            Self::TemplateCoordinate => ("unsorted", Some("unsorted:template-coordinate".into())),
            Self::Tag { tag, then } => {
                let then = match then {
                    TagTieBreak::Coordinate => "coordinate",
                    TagTieBreak::NameNatural => "queryname:natural",
                    TagTieBreak::NameAscii => "queryname:lexicographical",
                };

                let tag = String::from_utf8_lossy(&tag);
                ("unsorted", Some(format!("unsorted:{tag}:{then}")))
            }
        }
    }

    pub fn compare(self, a: &[u8], b: &[u8]) -> Ordering {
        match self {
            Self::Coordinate => coordinate_key(a).cmp(&coordinate_key(b)),
            Self::NameNatural => {
                strnum_cmp(name(a), name(b)).then_with(|| name_flags_key(a).cmp(&name_flags_key(b)))
            }
            Self::NameAscii => name(a)
                .cmp(name(b))
                .then_with(|| name_flags_key(a).cmp(&name_flags_key(b))),
            Self::Collate => name_hash(a)
                .cmp(&name_hash(b))
                .then_with(|| name(a).cmp(name(b))),
            Self::Tag { tag, .. } => {
                self.compare_tagged(a, tag_value(a, tag), b, tag_value(b, tag))
            }
            Self::TemplateCoordinate => template_coordinate::compare(a, b),
        }
    }

    /// Compares records by [`Order::Tag`] given their tag values (see [`tag_value`]), which
    /// callers may have looked up ahead of time. Other orders ignore the values.
    pub fn compare_tagged(
        self,
        a: &[u8],
        a_value: Option<&[u8]>,
        b: &[u8],
        b_value: Option<&[u8]>,
    ) -> Ordering {
        match self {
            Self::Tag { then, .. } => {
                compare_tag_values(a_value, b_value).then_with(|| then.order().compare(a, b))
            }
            _ => self.compare(a, b),
        }
    }
}

/// samtools' tie-break for records with the same name: READ1, READ2, then secondary, then
/// supplementary (`bam1_cmp_core`).
fn name_flags_key(record: &[u8]) -> u16 {
    let f = flags(record);
    ((f & 0xc0) << 8) | ((f & 0x100) << 3) | ((f & 0x800) >> 3)
}

/// Returns the first auxiliary field with `tag`, from its type byte to the end of the record,
/// or `None` if the record has no such field (or malformed auxiliary data).
pub fn tag_value(record: &[u8], tag: [u8; 2]) -> Option<&[u8]> {
    tag_offset(record, tag).map(|offset| &record[offset..])
}

/// The offset of the type byte of the first auxiliary field with `tag`.
pub fn tag_offset(record: &[u8], tag: [u8; 2]) -> Option<usize> {
    bam_raw::find_field(record, tag)
        .ok()
        .flatten()
        .map(|range| range.start + 2)
}

/// A tag value's class, as the type letter samtools normalizes it to (`normalize_type`).
fn normalized_type(ty: u8) -> u8 {
    match ty {
        b'c' | b'C' | b's' | b'S' | b'i' | b'I' => b'c',
        b'f' | b'd' => b'f',
        b'H' | b'Z' => b'H',
        _ => ty,
    }
}

fn integer_value(value: &[u8]) -> i64 {
    let b = &value[1..];

    match value[0] {
        b'c' => i64::from(b[0].cast_signed()),
        b'C' => i64::from(b[0]),
        b's' => i64::from(i16::from_le_bytes([b[0], b[1]])),
        b'S' => i64::from(u16::from_le_bytes([b[0], b[1]])),
        b'i' => i64::from(i32::from_le_bytes([b[0], b[1], b[2], b[3]])),
        b'I' => i64::from(u32::from_le_bytes([b[0], b[1], b[2], b[3]])),
        _ => 0,
    }
}

#[expect(
    clippy::cast_precision_loss,
    reason = "HTSlib's bam_aux2f converts integers to double"
)]
fn float_value(value: &[u8]) -> f64 {
    let b = &value[1..];

    match value[0] {
        b'f' => f64::from(f32::from_le_bytes([b[0], b[1], b[2], b[3]])),
        b'd' => f64::from_le_bytes([b[0], b[1], b[2], b[3], b[4], b[5], b[6], b[7]]),
        _ => integer_value(value) as f64,
    }
}

fn string_value(value: &[u8]) -> &[u8] {
    let s = &value[1..];
    memchr::memchr(0, s).map_or(s, |end| &s[..end])
}

/// samtools' `bam1_cmp_by_tag` without its tie-break: records without the tag first, then by
/// type class, then by value. `Equal` leaves the order to the tie-break.
pub fn compare_tag_values(a: Option<&[u8]>, b: Option<&[u8]>) -> Ordering {
    let (a, b) = match (a, b) {
        (None, None) => return Ordering::Equal,
        (None, Some(_)) => return Ordering::Less,
        (Some(_), None) => return Ordering::Greater,
        (Some(a), Some(b)) => (a, b),
    };

    let class = match (normalized_type(a[0]), normalized_type(b[0])) {
        (a_type, b_type) if a_type == b_type => a_type,
        // Integers compare with floats as floats.
        (b'c', b'f') | (b'f', b'c') => b'f',
        (a_type, b_type) => return a_type.cmp(&b_type),
    };

    match class {
        b'c' => integer_value(a).cmp(&integer_value(b)),
        b'f' => {
            let (x, y) = (float_value(a), float_value(b));

            // NaN, which samtools orders inconsistently, sorts last here so the order is total.
            x.partial_cmp(&y)
                .unwrap_or_else(|| x.is_nan().cmp(&y.is_nan()))
        }
        b'A' => a[1].cmp(&b[1]),
        b'H' => string_value(a).cmp(string_value(b)),
        _ => Ordering::Equal,
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
