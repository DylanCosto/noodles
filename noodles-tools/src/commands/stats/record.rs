//! The view of a record that statistics are collected from.

use noodles_sam::alignment::record::Flags;

use super::{
    CIGAR_DELETION, CIGAR_MATCH, CIGAR_SEQUENCE_MATCH, CIGAR_SEQUENCE_MISMATCH, CIGAR_SKIP,
};
use crate::bam_raw::{self, Record, le_u16, le_u32};

/// The fixed fields and parts of an encoded BAM record.
pub(super) struct Bam<'a> {
    pub(super) tid: i32,
    pub(super) pos: i64,
    pub(super) mapq: u8,
    pub(super) flag: u16,
    pub(super) mtid: i32,
    pub(super) mpos: i64,
    pub(super) isize: i32,
    pub(super) name: &'a [u8],
    pub(super) cigar: Vec<u32>,
    pub(super) seq: &'a [u8],
    pub(super) qual: &'a [u8],
    pub(super) l_seq: usize,
    pub(super) raw: &'a [u8],
}

impl<'a> Bam<'a> {
    pub(super) fn parse(raw: &'a [u8]) -> anyhow::Result<Self> {
        let record = Record::new(raw)?;
        let mut cigar = Vec::new();
        record.read_cigar(&mut cigar)?;

        Ok(Self {
            tid: record.reference_sequence_id(),
            pos: i64::from(record.position()),
            mapq: record.mapping_quality(),
            flag: record.flags(),
            mtid: record.mate_reference_sequence_id(),
            mpos: i64::from(record.mate_position()),
            isize: record.template_length(),
            name: record.name(),
            cigar,
            seq: record.packed_sequence(),
            qual: record.quality_scores(),
            l_seq: record.sequence_len(),
            raw,
        })
    }

    pub(super) fn has(&self, flag: Flags) -> bool {
        self.flag & flag.bits() != 0
    }

    pub(super) fn base(&self, i: usize) -> u8 {
        (self.seq[i / 2] >> ((!i & 1) << 2)) & 0xf
    }

    pub(super) fn is_paired_and_mapped(&self) -> bool {
        self.has(Flags::SEGMENTED) && !self.has(Flags::UNMAPPED) && !self.has(Flags::MATE_UNMAPPED)
    }

    pub(super) fn is_original(&self) -> bool {
        !self.has(Flags::SECONDARY) && !self.has(Flags::SUPPLEMENTARY)
    }

    /// `READ_ORDER_*`: 1 for a first (or single) read, 2 for a last read, 0 or 3 otherwise.
    pub(super) fn order(&self) -> u32 {
        if self.has(Flags::SEGMENTED) {
            u32::from(self.has(Flags::FIRST_SEGMENT)) + 2 * u32::from(self.has(Flags::LAST_SEGMENT))
        } else {
            READ_ORDER_FIRST
        }
    }

    /// htslib's `bam_endpos`.
    pub(super) fn end_pos(&self) -> i64 {
        let span: i64 = if self.has(Flags::UNMAPPED) {
            0
        } else {
            self.cigar
                .iter()
                .filter(|&&op| consumes_reference(op & 0xf))
                .map(|&op| i64::from(op >> 4))
                .sum()
        };

        self.pos + span.max(1)
    }
}

pub(super) const READ_ORDER_FIRST: u32 = 1;
pub(super) const READ_ORDER_LAST: u32 = 2;

pub(super) fn consumes_reference(op: u32) -> bool {
    matches!(
        op,
        CIGAR_MATCH | CIGAR_DELETION | CIGAR_SKIP | CIGAR_SEQUENCE_MATCH | CIGAR_SEQUENCE_MISMATCH
    )
}

/// htslib's `bam_aux2i`: the value of an integer field, or 0.
pub(super) fn integer_field(raw: &[u8], tag: [u8; 2]) -> anyhow::Result<Option<i64>> {
    let Some(range) = bam_raw::find_field(raw, tag)? else {
        return Ok(None);
    };

    let value = &raw[range.start + 3..range.end];

    Ok(Some(match raw[range.start + 2] {
        b'c' => i64::from(value[0].cast_signed()),
        b'C' => i64::from(value[0]),
        b's' => i64::from(i16::from_le_bytes([value[0], value[1]])),
        b'S' => i64::from(le_u16(value)),
        b'i' => i64::from(le_u32(value).cast_signed()),
        b'I' => i64::from(le_u32(value)),
        _ => 0,
    }))
}

/// The value bytes of a field, as samtools passes `bam_aux_get(...) + 1` to string functions:
/// up to the first NUL.
pub(super) fn string_value(raw: &[u8], tag: [u8; 2]) -> anyhow::Result<Option<&[u8]>> {
    let Some(range) = bam_raw::find_field(raw, tag)? else {
        return Ok(None);
    };

    let value = &raw[range.start + 3..range.end];
    let end = memchr::memchr(0, value).unwrap_or(value.len());
    Ok(Some(&value[..end]))
}
