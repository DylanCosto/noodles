//! samtools' template-coordinate order (`sort --template-coordinate`).
//!
//! Records are ordered by the 5' unclipped coordinates of both ends of their template (the mate's
//! from `MNEXT`, `PNEXT`, and the `MC` tag), the library of their read group, the `CB` and `MI`
//! tags, the read name, and which end of the template the record is.
//!
//! The key depends on the header (for libraries), so it is computed once per record and appended
//! to the raw record: a keyed record is `record ++ key ++ key length (u32 LE)`. Keeping the record
//! first leaves its fixed fields where the spill files and other orders expect them.

use std::{cmp::Ordering, collections::HashMap, io};

use bstr::BString;
use noodles_sam::{self as sam, header::record::value::map::read_group::tag::LIBRARY};

use crate::{
    bam_raw::{self, le_u32},
    sort_order::name,
};

const PAIRED: u16 = 0x01;
const UNMAPPED: u16 = 0x04;
const MATE_UNMAPPED: u16 = 0x08;
const REVERSE_COMPLEMENTED: u16 = 0x10;
const MATE_REVERSE_COMPLEMENTED: u16 = 0x20;

/// The length of the fixed part of a key: two reference sequence IDs, two positions, two strands,
/// and whether the record is the upper end of its template.
const FIXED_LEN: usize = 4 + 4 + 8 + 8 + 1 + 1 + 1;
/// The offset of the upper-end byte, compared after the strings.
const UPPER: usize = FIXED_LEN - 1;

/// Libraries (`LB`) by read group ID.
#[derive(Debug, Default)]
pub struct Libraries(HashMap<BString, BString>);

impl Libraries {
    /// The libraries of a header's read groups (`lookup_libraries`).
    pub fn from_header(header: &sam::Header) -> Self {
        Self(
            header
                .read_groups()
                .iter()
                .filter_map(|(id, read_group)| {
                    read_group
                        .other_fields()
                        .get(&LIBRARY)
                        .map(|lb| (id.clone(), lb.clone()))
                })
                .collect(),
        )
    }

    /// Libraries by the read group IDs of an input, given the merged read group ID of each
    /// (`samtools merge` translates `RG` before computing keys).
    pub fn translated<'a, I>(&self, read_groups: I) -> Self
    where
        I: IntoIterator<Item = (&'a BString, &'a BString)>,
    {
        Self(
            read_groups
                .into_iter()
                .filter_map(|(input_id, merged_id)| {
                    self.0
                        .get(merged_id)
                        .map(|lb| (input_id.clone(), lb.clone()))
                })
                .collect(),
        )
    }

    /// The library of a record's read group, or the empty string. `RG` values of type `H` count
    /// only when `hex_ids` is [`HexIds::Accept`].
    pub fn get(&self, record: &[u8], hex_ids: HexIds) -> io::Result<&[u8]> {
        if self.0.is_empty() {
            return Ok(b"");
        }

        let value = match bam_raw::find_field(record, *b"RG")? {
            Some(range) => match record[range.start + 2] {
                b'Z' => string_value(&record[range]),
                b'H' if hex_ids == HexIds::Accept => string_value(&record[range]),
                _ => None,
            },
            None => None,
        };

        Ok(value
            .and_then(|id| self.0.get(id))
            .map_or(b"".as_slice(), |lb| lb.as_slice()))
    }
}

/// Whether `RG` values of type `H` name read groups.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub enum HexIds {
    /// As in `samtools merge`, which converts them to `Z` when translating read groups.
    Accept,
    /// As in `samtools sort`, which only looks up `Z` values.
    Ignore,
}

/// The value of a `Z` or `H` field (`bam_aux2Z`), given the field from its tag.
fn string_value(field: &[u8]) -> Option<&[u8]> {
    match field[2] {
        b'Z' | b'H' => Some(&field[3..field.len() - 1]),
        _ => None,
    }
}

/// Appends `record`, followed by its template-coordinate key with `library`, to `dst`
/// (`template_coordinate_key`).
///
/// # Errors
///
/// Returns an error, with samtools' message, if a record with a mapped mate has no `MC` tag, or
/// if `MC`, `CB`, or `MI` is not a string.
pub fn push_keyed(dst: &mut Vec<u8>, record: &[u8], library: &[u8]) -> io::Result<()> {
    let raw = bam_raw::Record::new(record)?;
    let flags = raw.flags();

    let mut tid1 = i32::MAX;
    let mut tid2 = i32::MAX;
    let mut pos1 = i64::MAX;
    let mut pos2 = i64::MAX;
    let mut neg1 = false;
    let mut neg2 = false;

    if flags & UNMAPPED == 0 {
        let cigar = match raw.long_cigar_field()? {
            Some(range) => &record[range.start + 8..range.end],
            None => raw.raw_cigar(),
        };

        tid1 = raw.reference_sequence_id();
        neg1 = flags & REVERSE_COMPLEMENTED != 0;
        let pos = i64::from(raw.position());
        pos1 = if neg1 {
            unclipped_end(pos, cigar)
        } else {
            unclipped_start(pos, cigar)
        };
    }

    let [mc, cb, mi] = bam_raw::find_fields(record, [*b"MC", *b"CB", *b"MI"])?;

    if flags & PAIRED != 0 && flags & MATE_UNMAPPED == 0 {
        let Some(range) = mc else {
            return Err(key_error(
                "[bam_sort] error: no MC tag. Please run samtools fixmate on file first.",
            ));
        };

        let cigar = string_value(&record[range]).ok_or_else(|| {
            key_error(
                "[bam_sort] error: MC tag wrong type. Please use the MC tag provided by samtools \
                 fixmate.",
            )
        })?;

        tid2 = raw.mate_reference_sequence_id();
        neg2 = flags & MATE_REVERSE_COMPLEMENTED != 0;
        let pos = i64::from(raw.mate_position());
        pos2 = if neg2 {
            unclipped_other_end(pos, cigar)
        } else {
            unclipped_other_start(pos, cigar)
        };
    }

    let tag_string = |range: Option<std::ops::Range<usize>>, tag: &str| match range {
        Some(range) => string_value(&record[range]).ok_or_else(|| {
            key_error(&format!(
                "[bam_sort] error: {tag} tag wrong type (not a string)."
            ))
        }),
        None => Ok(b"".as_slice()),
    };

    let cell = tag_string(cb, "CB")?;
    let molecule = strip_mid_suffix(tag_string(mi, "MI")?);

    let is_upper =
        !(tid1 < tid2 || (tid1 == tid2 && pos1 < pos2) || (tid1 == tid2 && pos1 == pos2 && !neg1));

    if is_upper {
        (tid1, tid2) = (tid2, tid1);
        (pos1, pos2) = (pos2, pos1);
        (neg1, neg2) = (neg2, neg1);
    }

    dst.extend_from_slice(record);
    let key_start = dst.len();

    // The fixed part compares bytewise: big-endian with the sign bit flipped, and reverse strands
    // (0) before forward ones (1).
    dst.extend_from_slice(&(tid1.cast_unsigned() ^ (1 << 31)).to_be_bytes());
    dst.extend_from_slice(&(tid2.cast_unsigned() ^ (1 << 31)).to_be_bytes());
    dst.extend_from_slice(&(pos1.cast_unsigned() ^ (1 << 63)).to_be_bytes());
    dst.extend_from_slice(&(pos2.cast_unsigned() ^ (1 << 63)).to_be_bytes());
    dst.push(u8::from(!neg1));
    dst.push(u8::from(!neg2));
    dst.push(u8::from(is_upper));

    for s in [library, cell, molecule] {
        let len =
            u32::try_from(s.len()).map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))?;
        dst.extend_from_slice(&len.to_le_bytes());
        dst.extend_from_slice(s);
    }

    let key_len = u32::try_from(dst.len() - key_start)
        .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))?;
    dst.extend_from_slice(&key_len.to_le_bytes());

    Ok(())
}

fn key_error(message: &str) -> io::Error {
    io::Error::new(io::ErrorKind::InvalidData, message.to_owned())
}

/// Splits a keyed record into the record and its key.
fn split(keyed: &[u8]) -> (&[u8], &[u8]) {
    let (rest, len) = keyed.split_at(keyed.len() - 4);
    let len = usize::try_from(le_u32(len)).unwrap_or(usize::MAX);
    rest.split_at(rest.len() - len)
}

/// The record of a keyed record.
pub fn record(keyed: &[u8]) -> &[u8] {
    split(keyed).0
}

/// Compares keyed records (`bam1_cmp_template_coordinate`).
pub fn compare(a: &[u8], b: &[u8]) -> Ordering {
    let (a_record, a_key) = split(a);
    let (b_record, b_key) = split(b);

    a_key[..UPPER]
        .cmp(&b_key[..UPPER])
        .then_with(|| {
            let (mut a_rest, mut b_rest) = (&a_key[FIXED_LEN..], &b_key[FIXED_LEN..]);
            let mut next = || (next_string(&mut a_rest), next_string(&mut b_rest));

            let (a_library, b_library) = next();
            let (a_cell, b_cell) = next();
            let (a_molecule, b_molecule) = next();

            a_library
                .cmp(b_library)
                .then_with(|| a_cell.cmp(b_cell))
                .then_with(|| compare_mid(a_molecule, b_molecule))
        })
        .then_with(|| name(a_record).cmp(name(b_record)))
        .then_with(|| a_key[UPPER].cmp(&b_key[UPPER]))
}

fn next_string<'a>(src: &mut &'a [u8]) -> &'a [u8] {
    let (len, rest) = src.split_at(4);
    let len = usize::try_from(le_u32(len)).unwrap_or(usize::MAX);
    let (s, rest) = rest.split_at(len);
    *src = rest;
    s
}

/// Removes a trailing slash and single character from a molecular identifier (`1/A` -> `1`).
fn strip_mid_suffix(mid: &[u8]) -> &[u8] {
    match mid {
        [rest @ .., b'/', _] => rest,
        _ => mid,
    }
}

/// `template_coordinate_key_compare_mid` on suffix-stripped identifiers: a prefix sorts first,
/// and characters compare as C `char`s, which are signed on the platforms samtools is usually
/// built for.
fn compare_mid(a: &[u8], b: &[u8]) -> Ordering {
    a.iter()
        .map(|&c| c.cast_signed())
        .cmp(b.iter().map(|&c| c.cast_signed()))
}

/// The encoded CIGAR operations as (length, operation code) pairs.
fn ops(cigar: &[u8]) -> impl DoubleEndedIterator<Item = (i64, u32)> + '_ {
    cigar.as_chunks::<4>().0.iter().map(|op| {
        let op = u32::from_le_bytes(*op);
        (i64::from(op >> 4), op & 0x0f)
    })
}

const MATCH: u32 = 0;
const DELETION: u32 = 2;
const SKIP: u32 = 3;
const SOFT_CLIP: u32 = 4;
const HARD_CLIP: u32 = 5;
const SEQUENCE_MATCH: u32 = 7;
const SEQUENCE_MISMATCH: u32 = 8;

/// The 1-based start, including leading soft clips (`unclipped_start`).
fn unclipped_start(pos: i64, cigar: &[u8]) -> i64 {
    let mut clipped = 0;

    for (len, op) in ops(cigar) {
        match op {
            SOFT_CLIP => clipped += len,
            HARD_CLIP => {}
            _ => break,
        }
    }

    pos - clipped + 1
}

/// The 1-based end, including trailing soft clips (`unclipped_end`).
fn unclipped_end(pos: i64, cigar: &[u8]) -> i64 {
    let reference_len: i64 = ops(cigar)
        .filter(|&(_, op)| {
            matches!(
                op,
                MATCH | DELETION | SKIP | SEQUENCE_MATCH | SEQUENCE_MISMATCH
            )
        })
        .map(|(len, _)| len)
        .sum();

    // `bam_endpos`: a record without CIGAR operations covers one base.
    let end = if cigar.is_empty() {
        pos + 1
    } else {
        pos + reference_len
    };

    let mut clipped = 0;

    for (len, op) in ops(cigar).rev() {
        match op {
            SOFT_CLIP => clipped += len,
            HARD_CLIP => {}
            _ => break,
        }
    }

    end + clipped
}

/// Parses the next CIGAR operation of a textual CIGAR (as samtools does: a missing length is 1),
/// returning the length and the operation character, or `None` at the end of the string or at
/// `*`.
fn next_text_op(cigar: &mut &[u8]) -> Option<(i64, u8)> {
    match cigar.first() {
        None | Some(b'*') => None,
        Some(c) if c.is_ascii_digit() => {
            let digits = cigar.iter().take_while(|c| c.is_ascii_digit()).count();
            let len = cigar[..digits].iter().fold(0i64, |n, &d| {
                n.saturating_mul(10).saturating_add(i64::from(d - b'0'))
            });
            let op = cigar.get(digits).copied().unwrap_or(0);
            *cigar = cigar.get(digits + 1..).unwrap_or_default();
            Some((len, op))
        }
        Some(&op) => {
            *cigar = &cigar[1..];
            Some((1, op))
        }
    }
}

/// The mate's 1-based unclipped start from its 0-based position and `MC` CIGAR
/// (`unclipped_other_start`).
fn unclipped_other_start(pos: i64, mut cigar: &[u8]) -> i64 {
    let mut clipped = 0i64;

    while let Some((len, op)) = next_text_op(&mut cigar) {
        match op {
            b'S' => clipped = clipped.saturating_add(len),
            b'H' => {}
            _ => break,
        }
    }

    pos - clipped + 1
}

/// The mate's 1-based unclipped end from its 0-based position and `MC` CIGAR
/// (`unclipped_other_end`). As in samtools, leading soft clips are skipped but later ones count.
fn unclipped_other_end(pos: i64, mut cigar: &[u8]) -> i64 {
    let mut reference_len = 0i64;
    let mut skip = true;

    while let Some((len, op)) = next_text_op(&mut cigar) {
        match op {
            b'M' | b'D' | b'N' | b'=' | b'X' => {
                reference_len = reference_len.saturating_add(len);
                skip = false;
            }
            b'S' if !skip => reference_len = reference_len.saturating_add(len),
            0 => break,
            _ => {}
        }
    }

    pos.saturating_add(reference_len)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_unclipped_other() {
        assert_eq!(unclipped_other_start(99, b"2S8M"), 98);
        assert_eq!(unclipped_other_start(99, b"2H3S8M"), 97);
        assert_eq!(unclipped_other_start(99, b"*"), 100);
        assert_eq!(unclipped_other_end(99, b"10M2S"), 111);
        assert_eq!(unclipped_other_end(99, b"3S4M1D2N3M2S4H"), 111);
        assert_eq!(unclipped_other_end(99, b"*"), 99);
        assert_eq!(unclipped_other_end(99, b"10"), 99);
    }

    #[test]
    fn test_compare_mid() {
        assert_eq!(strip_mid_suffix(b"1/A"), b"1");
        assert_eq!(strip_mid_suffix(b"/A"), b"");
        assert_eq!(strip_mid_suffix(b"10"), b"10");
        assert_eq!(compare_mid(b"1", b"10"), Ordering::Less);
        assert_eq!(compare_mid(b"\x80", b"a"), Ordering::Less);
    }
}
