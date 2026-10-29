//! Direct encoding of CRAM records as BAM records.

use std::{io, mem};

use noodles_bam as bam;
use noodles_core::Position;
use noodles_sam::{
    self as sam,
    alignment::record::{
        cigar::{Op, op::Kind},
        data::field::Tag,
    },
};

use super::{
    Record,
    cigar::iter::Iter as CigarOps,
    data::field::{Value, value::Array},
};

// § 4.2.1 "BIN field calculation" (2021-06-03): "unmapped reads with `POS` 0 (which becomes -1 in
// BAM) therefore use `reg2bin(-1, 0)` which is computed as 4680."
const UNMAPPED_BIN: u16 = 4680;

// § 4.2 "The BAM format" (2023-05-24): `mapq` is 255 when unavailable.
const MISSING_MAPPING_QUALITY: u8 = 255;

// § 4.2.3 "SEQ and QUAL encoding" (2023-05-24): missing quality scores are 0xff, and stored ones
// are in [0, 93].
const MISSING_QUALITY_SCORE: u8 = 0xff;
const MAX_QUALITY_SCORE: u8 = 93;

// § 4.2.2 "N_CIGAR_OP field" (2023-11-16): op lengths are stored in 28 bits.
const MAX_OP_LENGTH: usize = (1 << 28) - 1;

const MISSING_NAME: &[u8] = b"*";
const MAX_NAME_LENGTH: usize = 254;
const NUL: u8 = 0x00;

impl Record<'_> {
    /// Encodes this record as a BAM record, without its leading `block_size`, into `dst`.
    ///
    /// The bytes are the same as those [`bam::io::Writer`] writes for this record after
    /// `block_size`, but they are written from the decoded CRAM fields directly instead of
    /// through the generic [`sam::alignment::Record`] accessors. The sequence is reconstructed in
    /// bulk. `dst` is cleared first.
    ///
    /// # Errors
    ///
    /// Returns an error if the record cannot be represented as a BAM record, e.g., when a
    /// reference sequence ID has no entry in `header`.
    pub fn encode_bam_record(&self, header: &sam::Header, dst: &mut Vec<u8>) -> io::Result<()> {
        dst.clear();

        // Rare or invalid records (overflowing CIGARs, quality scores stored as read features,
        // out-of-range values) are left to the generic encoder, which also reports the errors.
        if self.try_encode_bam_record(header, dst).is_none() {
            self.encode_bam_record_generically(header, dst)?;
        }

        Ok(())
    }

    fn try_encode_bam_record(&self, header: &sam::Header, dst: &mut Vec<u8>) -> Option<()> {
        const CIGAR_OP_COUNT_RANGE: std::ops::Range<usize> = 12..14;
        const READ_LENGTH_RANGE: std::ops::Range<usize> = 16..20;

        let reference_sequence_count = header.reference_sequences().len();

        let reference_sequence_id =
            encode_reference_sequence_id(reference_sequence_count, self.reference_sequence_id)?;
        let alignment_start = encode_position(self.alignment_start)?;
        let mate_reference_sequence_id = encode_reference_sequence_id(
            reference_sequence_count,
            self.mate_reference_sequence_id,
        )?;
        let mate_alignment_start = encode_position(self.mate_alignment_start)?;

        let name = match self.name.as_deref() {
            Some(name) if is_valid_name(name) => name,
            Some(_) => return None,
            None => MISSING_NAME,
        };

        let mapping_quality = self
            .mapping_quality
            .map(u8::from)
            .unwrap_or(MISSING_MAPPING_QUALITY);

        let bin = match (self.alignment_start, self.alignment_end()) {
            (Some(start), Some(end)) => region_to_bin(start, end),
            _ => UNMAPPED_BIN,
        };

        dst.extend(reference_sequence_id.to_le_bytes());
        dst.extend(alignment_start.to_le_bytes());
        // SAFETY: `is_valid_name` bounds the length to 254.
        dst.push((name.len() + 1) as u8);
        dst.push(mapping_quality);
        dst.extend(bin.to_le_bytes());
        dst.extend([0; 2]); // n_cigar_op, set below
        dst.extend(u16::from(self.bam_flags).to_le_bytes());
        dst.extend([0; 4]); // l_seq, set below
        dst.extend(mate_reference_sequence_id.to_le_bytes());
        dst.extend(mate_alignment_start.to_le_bytes());
        dst.extend(self.template_length.to_le_bytes());

        dst.extend_from_slice(name);
        dst.push(NUL);

        let (op_count, cigar_read_length) = self.write_bam_cigar(dst)?;
        let op_count = u16::try_from(op_count).ok()?;
        dst[CIGAR_OP_COUNT_RANGE].copy_from_slice(&op_count.to_le_bytes());

        let base_count = self.write_bam_sequence(dst);

        // § 1.4.10 "`SEQ`" (2022-08-22): the sequence length must equal the CIGAR read length.
        if base_count > 0 && cigar_read_length > 0 && base_count != cigar_read_length {
            return None;
        }

        let n = u32::try_from(base_count).ok()?;
        dst[READ_LENGTH_RANGE].copy_from_slice(&n.to_le_bytes());

        self.write_bam_quality_scores(base_count, dst)?;
        self.write_bam_data(dst)?;

        Some(())
    }

    /// Writes the simplified CIGAR ops and returns their count and read length.
    fn write_bam_cigar(&self, dst: &mut Vec<u8>) -> Option<(usize, usize)> {
        fn write_op(dst: &mut Vec<u8>, op: Op) -> Option<()> {
            if op.len() > MAX_OP_LENGTH {
                return None;
            }

            // SAFETY: `op.len() <= MAX_OP_LENGTH < 2^28`.
            let n = ((op.len() as u32) << 4) | encode_kind(op.kind());
            dst.extend(n.to_le_bytes());

            Some(())
        }

        if self.bam_flags.is_unmapped() {
            return Some((0, 0));
        }

        let mut op_count = 0;
        let mut read_length = 0;
        let mut prev_op: Option<Op> = None;

        // Adjacent ops of the same kind are merged, as `TrySimplify` does for the generic CIGAR.
        for op in CigarOps::new(&self.features, self.read_length) {
            if op.kind().consumes_read() {
                read_length += op.len();
            }

            match prev_op {
                Some(prev) if prev.kind() == op.kind() => {
                    prev_op = Some(Op::new(prev.kind(), prev.len() + op.len()));
                }
                Some(prev) => {
                    write_op(dst, prev)?;
                    op_count += 1;
                    prev_op = Some(op);
                }
                None => prev_op = Some(op),
            }
        }

        if let Some(op) = prev_op {
            write_op(dst, op)?;
            op_count += 1;
        }

        Some((op_count, read_length))
    }

    /// Writes the 4-bit packed sequence and returns its length.
    fn write_bam_sequence(&self, dst: &mut Vec<u8>) -> usize {
        let start = dst.len();

        if self.bam_flags.is_unmapped() || self.cram_flags.sequence_is_missing() {
            dst.extend_from_slice(&self.sequence);
        } else {
            self.sequence_decoder().extend_to(dst);
        }

        let base_count = dst.len() - start;
        pack_bases_in_place(&mut dst[start..]);
        dst.truncate(start + base_count.div_ceil(2));

        base_count
    }

    fn write_bam_quality_scores(&self, base_count: usize, dst: &mut Vec<u8>) -> Option<()> {
        if !(self.bam_flags.is_unmapped() || self.cram_flags.quality_scores_are_stored_as_array()) {
            return None;
        }

        let scores = &self.quality_scores[..];

        if scores.len() == base_count {
            // `max` vectorizes, where a short-circuiting `all` does not.
            if scores.iter().copied().max().unwrap_or_default() > MAX_QUALITY_SCORE {
                return None;
            }

            dst.extend_from_slice(scores);
        } else if scores.is_empty() {
            dst.resize(dst.len() + base_count, MISSING_QUALITY_SCORE);
        } else {
            return None;
        }

        Some(())
    }

    fn write_bam_data(&self, dst: &mut Vec<u8>) -> Option<()> {
        // Without a header, the generic record has no data, not even a read group.
        let Some(header) = self.header else {
            return Some(());
        };

        for (tag, value) in &self.data {
            if *tag == Tag::CIGAR {
                continue;
            }

            dst.extend(tag.as_ref());
            write_bam_value(dst, value)?;
        }

        if let Some(id) = self.read_group_id {
            let (name, _) = header.read_groups().get_index(id)?;
            dst.extend(Tag::READ_GROUP.as_ref());
            dst.push(b'Z');
            write_string(dst, name)?;
        }

        Some(())
    }

    fn encode_bam_record_generically(
        &self,
        header: &sam::Header,
        dst: &mut Vec<u8>,
    ) -> io::Result<()> {
        use sam::alignment::io::Write as _;

        let mut sequence = Vec::new();
        let record = self.with_sequence(&mut sequence);

        dst.clear();
        let mut writer = bam::io::Writer::from(mem::take(dst));
        let result = writer.write_alignment_record(header, &record);
        *dst = writer.into_inner();
        result?;

        dst.drain(..mem::size_of::<u32>());

        Ok(())
    }
}

fn encode_reference_sequence_id(reference_sequence_count: usize, id: Option<usize>) -> Option<i32> {
    const UNMAPPED: i32 = -1;

    match id {
        Some(id) if id < reference_sequence_count => i32::try_from(id).ok(),
        Some(_) => None,
        None => Some(UNMAPPED),
    }
}

fn encode_position(position: Option<Position>) -> Option<i32> {
    const MISSING: i32 = -1;

    match position {
        Some(position) => i32::try_from(usize::from(position) - 1).ok(),
        None => Some(MISSING),
    }
}

fn is_valid_name(name: &[u8]) -> bool {
    (1..=MAX_NAME_LENGTH).contains(&name.len())
        && name != MISSING_NAME
        && name.iter().all(|&b| b.is_ascii_graphic() && b != b'@')
}

// § 5.3 "C source code for computing bin number and overlapping bins" (2021-06-03)
#[allow(clippy::eq_op)]
fn region_to_bin(alignment_start: Position, alignment_end: Position) -> u16 {
    let start = usize::from(alignment_start) - 1;
    let end = usize::from(alignment_end) - 1;

    let bin = if start >> 14 == end >> 14 {
        ((1 << 15) - 1) / 7 + (start >> 14)
    } else if start >> 17 == end >> 17 {
        ((1 << 12) - 1) / 7 + (start >> 17)
    } else if start >> 20 == end >> 20 {
        ((1 << 9) - 1) / 7 + (start >> 20)
    } else if start >> 23 == end >> 23 {
        ((1 << 6) - 1) / 7 + (start >> 23)
    } else if start >> 26 == end >> 26 {
        ((1 << 3) - 1) / 7 + (start >> 26)
    } else {
        0
    };

    // Overflowing bin IDs are truncated, as in the BAM encoder.
    bin as u16
}

fn encode_kind(kind: Kind) -> u32 {
    match kind {
        Kind::Match => 0,
        Kind::Insertion => 1,
        Kind::Deletion => 2,
        Kind::Skip => 3,
        Kind::SoftClip => 4,
        Kind::HardClip => 5,
        Kind::Pad => 6,
        Kind::SequenceMatch => 7,
        Kind::SequenceMismatch => 8,
    }
}

/// Packs the bases in `buf` two per byte into its first half, in place.
fn pack_bases_in_place(buf: &mut [u8]) {
    // § 4.2.3 "SEQ and QUAL encoding" (2023-11-16): "The case-insensitive base codes [...] are
    // mapped to [0, 15] respectively with all other characters mapping to 'N' (value 15)".
    const CODES: [u8; 256] = build_codes();
    // § 4.2.3 (2021-06-03): "When `l_seq` is odd the bottom 4 bits of the last byte are
    // undefined, but we recommend writing these as zero."
    const PAD: u8 = 0;

    let len = buf.len();

    // Byte i reads bases 2i and 2i + 1, which are at or after i, so it never reads a packed byte.
    for i in 0..len / 2 {
        let l = CODES[usize::from(buf[2 * i])];
        let r = CODES[usize::from(buf[2 * i + 1])];
        buf[i] = (l << 4) | r;
    }

    if len % 2 == 1 {
        buf[len / 2] = (CODES[usize::from(buf[len - 1])] << 4) | PAD;
    }
}

const fn build_codes() -> [u8; 256] {
    const BASES: [u8; 16] = *b"=ACMGRSVTWYHKDBN";
    const N: u8 = 0x0f;

    let mut codes = [N; 256];
    let mut i = 0;

    while i < BASES.len() {
        let base = BASES[i];
        // SAFETY: `i < BASES.len() <= u8::MAX`.
        codes[base as usize] = i as u8;
        codes[base.to_ascii_lowercase() as usize] = i as u8;
        i += 1;
    }

    codes
}

fn write_bam_value(dst: &mut Vec<u8>, value: &Value<'_>) -> Option<()> {
    match value {
        Value::Character(c) => dst.extend([b'A', *c]),
        Value::Int8(n) => dst.extend([b'c', *n as u8]),
        Value::UInt8(n) => dst.extend([b'C', *n]),
        Value::Int16(n) => {
            dst.push(b's');
            dst.extend(n.to_le_bytes());
        }
        Value::UInt16(n) => {
            dst.push(b'S');
            dst.extend(n.to_le_bytes());
        }
        Value::Int32(n) => {
            dst.push(b'i');
            dst.extend(n.to_le_bytes());
        }
        Value::UInt32(n) => {
            dst.push(b'I');
            dst.extend(n.to_le_bytes());
        }
        Value::Float(n) => {
            dst.push(b'f');
            dst.extend(n.to_le_bytes());
        }
        Value::String(s) => {
            dst.push(b'Z');
            write_string(dst, s)?;
        }
        Value::Hex(s) => {
            // § 1.5 "The alignment section: optional fields" (2024-11-06): "`([0-9A-F][0-9A-F])*`".
            if !s.len().is_multiple_of(2)
                || !s.iter().all(|b| matches!(b, b'0'..=b'9' | b'A'..=b'F'))
            {
                return None;
            }

            dst.push(b'H');
            dst.extend_from_slice(s);
            dst.push(NUL);
        }
        Value::Array(array) => write_array(dst, array)?,
    }

    Some(())
}

fn write_string(dst: &mut Vec<u8>, s: &[u8]) -> Option<()> {
    // § 1.5 "The alignment section: optional fields" (2024-11-06): "`[ !-~]*`".
    if !s.iter().all(|b| matches!(b, b' '..=b'~')) {
        return None;
    }

    dst.extend_from_slice(s);
    dst.push(NUL);

    Some(())
}

fn write_array(dst: &mut Vec<u8>, array: &Array<'_>) -> Option<()> {
    let src = array.as_bam_bytes()?;
    dst.push(b'B');
    dst.extend_from_slice(src);
    Some(())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_pack_bases_in_place() {
        let mut buf = b"ACGTn".to_vec();
        pack_bases_in_place(&mut buf);
        assert_eq!(&buf[..3], [0x12, 0x48, 0xf0]);

        let mut buf = b"=X".to_vec();
        pack_bases_in_place(&mut buf);
        assert_eq!(&buf[..1], [0x0f]);
    }

    #[test]
    fn test_region_to_bin() -> Result<(), noodles_core::position::TryFromIntError> {
        assert_eq!(
            region_to_bin(Position::try_from(8)?, Position::try_from(13)?),
            4681
        );
        assert_eq!(
            region_to_bin(Position::try_from(63245986)?, Position::try_from(63245986)?),
            8541
        );
        Ok(())
    }
}
