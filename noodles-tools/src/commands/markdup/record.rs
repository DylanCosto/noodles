//! Fields of encoded records that duplicate marking reads and writes.

use std::{collections::HashMap, io};

use anyhow::{Context, bail};
use noodles_bam as bam;
use noodles_sam as sam;
use noodles_sam::alignment::record::{
    Flags,
    cigar::op::Kind,
    data::field::{Tag, Value},
};

use super::{MATE_SCORE, MIN_SCORED_QUALITY};
use crate::bam_raw::{self, field, le_u16, le_u32};

/// Whether a duplicate may have records that `-S` marks: it has an `SA` or `XA` tag, or an
/// unmapped mate.
pub(super) fn has_other_records(record: &[u8]) -> io::Result<bool> {
    Ok(flags_of(record).is_mate_unmapped()
        || bam_raw::find_field(record, *b"SA")?.is_some()
        || bam_raw::find_field(record, *b"XA")?.is_some())
}

pub(super) fn i32_field(record: &[u8], range: std::ops::Range<usize>) -> i32 {
    le_u32(&record[range]).cast_signed()
}

pub(super) fn flags_of(record: &[u8]) -> Flags {
    Flags::from_bits_retain(le_u16(&record[field::FLAGS]))
}

pub(super) fn set_flags(record: &mut [u8], flags: Flags) {
    record[field::FLAGS].copy_from_slice(&flags.bits().to_le_bytes());
}

pub(super) fn name(record: &[u8]) -> &[u8] {
    crate::sort_order::name(record)
}

/// `calc_score`: the sum of base qualities of at least 15 (missing qualities count as 255).
pub(super) fn read_score(record: &[u8]) -> i64 {
    let Ok(data) = bam_raw::data_range(record) else {
        return 0;
    };

    let base_count = usize::try_from(le_u32(&record[field::SEQUENCE_LENGTH])).unwrap_or(0);
    let quality_scores = &record[data.start - base_count..data.start];

    quality_scores
        .iter()
        .filter(|&&q| q >= MIN_SCORED_QUALITY)
        .map(|&q| i64::from(q))
        .sum()
}

/// `get_mate_score`, from the `ms` tag.
pub(super) fn mate_score(record: &bam::Record) -> anyhow::Result<i64> {
    let data = record.data();

    match sam::alignment::record::Data::get(&data, &MATE_SCORE).transpose()? {
        Some(value) => value
            .as_int()
            .context("[markdup] error, ms tag is not an integer"),
        None => {
            bail!("[markdup] error, no ms score tag. Please run samtools fixmate on file first.")
        }
    }
}

pub(super) fn read_group_of(
    record: &bam::Record,
    read_groups: &HashMap<Vec<u8>, usize>,
) -> io::Result<usize> {
    let data = record.data();

    Ok(
        match sam::alignment::record::Data::get(&data, &Tag::READ_GROUP).transpose()? {
            Some(Value::String(id)) => read_groups.get(id.as_ref() as &[u8]).copied().unwrap_or(0),
            _ => 0,
        },
    )
}

/// `has_mate`: the read's mate is present and mapped.
pub(super) fn has_mate(record: &bam::Record) -> bool {
    let raw = record.as_ref();
    let flags = flags_of(raw);
    let mtid = i32_field(raw, field::MATE_REFERENCE_SEQUENCE_ID);
    let mpos = i32_field(raw, 24..28);

    flags.is_segmented() && !flags.is_mate_unmapped() && !(mtid == -1 && mpos == -1)
}

/// The CIGAR as `(kind, len)` pairs.
pub(super) fn cigar_ops(record: &bam::Record) -> io::Result<Vec<(Kind, i64)>> {
    record
        .cigar()
        .iter()
        .map(|result| result.map(|op| (op.kind(), i64::try_from(op.len()).unwrap_or(i64::MAX))))
        .collect()
}

/// HTSlib's `bam_endpos`.
pub(super) fn end_pos(record: &bam::Record, ops: &[(Kind, i64)]) -> i64 {
    let raw = record.as_ref();
    let pos = i64::from(i32_field(raw, field::POSITION));

    let span: i64 = if flags_of(raw).is_unmapped() {
        0
    } else {
        ops.iter()
            .filter(|(kind, _)| kind.consumes_reference())
            .map(|(_, len)| len)
            .sum()
    };

    pos + span.max(1)
}

pub(super) fn is_clip(kind: Kind) -> bool {
    matches!(kind, Kind::SoftClip | Kind::HardClip)
}

/// `unclipped_start(b, 1)`
pub(super) fn unclipped_start(record: &bam::Record, ops: &[(Kind, i64)]) -> i64 {
    let pos = i64::from(i32_field(record.as_ref(), field::POSITION));
    let clipped: i64 = ops
        .iter()
        .take_while(|(kind, _)| is_clip(*kind))
        .map(|(_, len)| len)
        .sum();

    pos - clipped + 1
}

/// `unclipped_end(b, 1)`
pub(super) fn unclipped_end(record: &bam::Record, ops: &[(Kind, i64)]) -> i64 {
    let clipped: i64 = ops
        .iter()
        .rev()
        .take_while(|(kind, _)| is_clip(*kind))
        .map(|(_, len)| len)
        .sum();

    end_pos(record, ops) + clipped
}

/// Parses a CIGAR string as samtools' `unclipped_other_*` do: numbers via `strtol`, and an
/// operation without a number counts 1.
pub(super) fn parse_cigar_string(s: &[u8]) -> Vec<(u8, i64)> {
    let mut ops = Vec::new();
    let mut i = 0;

    while i < s.len() && s[i] != b'*' {
        let start = i;

        while i < s.len() && s[i].is_ascii_digit() {
            i += 1;
        }

        let n = if i > start {
            std::str::from_utf8(&s[start..i])
                .ok()
                .and_then(|t| t.parse().ok())
                .unwrap_or(i64::MAX)
        } else {
            1
        };

        if i < s.len() {
            ops.push((s[i], n));
            i += 1;
        }
    }

    ops
}

/// `unclipped_other_start(mpos, MC, 1)`
pub(super) fn unclipped_other_start(mpos: i64, cigar: &[(u8, i64)]) -> i64 {
    let mut clipped = 0;

    for &(op, n) in cigar {
        match op {
            b'S' | b'H' => clipped += n,
            _ => break,
        }
    }

    mpos - clipped + 1
}

/// `unclipped_other_end(mpos, MC, 1)`
pub(super) fn unclipped_other_end(mpos: i64, cigar: &[(u8, i64)]) -> i64 {
    let mut reference_position = 0;
    let mut skip = true;

    for &(op, n) in cigar {
        match op {
            b'M' | b'D' | b'N' | b'=' | b'X' => {
                reference_position += n;
                skip = false;
            }
            b'S' | b'H' if !skip => reference_position += n,
            _ => {}
        }
    }

    mpos + reference_position
}

/// Appends an integer tag with the smallest type that holds it, as HTSlib does.
pub(super) fn append_integer_tag(record: &mut Vec<u8>, tag: Tag, n: i64) {
    let tag: [u8; 2] = tag.into();
    record.extend_from_slice(&tag);

    if let Ok(n) = u8::try_from(n) {
        record.push(b'C');
        record.push(n);
    } else if let Ok(n) = i8::try_from(n) {
        record.push(b'c');
        record.extend_from_slice(&n.to_le_bytes());
    } else if let Ok(n) = u16::try_from(n) {
        record.push(b'S');
        record.extend_from_slice(&n.to_le_bytes());
    } else if let Ok(n) = i16::try_from(n) {
        record.push(b's');
        record.extend_from_slice(&n.to_le_bytes());
    } else if let Ok(n) = u32::try_from(n) {
        record.push(b'I');
        record.extend_from_slice(&n.to_le_bytes());
    } else {
        record.push(b'i');
        let n = i32::try_from(n).unwrap_or(i32::MAX);
        record.extend_from_slice(&n.to_le_bytes());
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_unclipped_other_positions() {
        let cigar = parse_cigar_string(b"5S10M2D3M4S2H");
        assert_eq!(unclipped_other_start(100, &cigar), 96);
        assert_eq!(unclipped_other_end(100, &cigar), 100 + 10 + 2 + 3 + 4 + 2);
        assert!(parse_cigar_string(b"*").is_empty());
    }
}
