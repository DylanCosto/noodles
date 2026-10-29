//! An encoded BAM record with mutable qualities and tags.

use std::{io, ops::Range};

use noodles_bam as bam;

const NAME_LENGTH: usize = 8;
const NAME: usize = 32;

pub(crate) const CIGAR_SOFT_CLIP: u32 = 4;
pub(crate) const CIGAR_SKIP: u32 = 3;

/// An alignment record as htslib holds it while piling up: its encoded BAM bytes, whose
/// qualities, sequence, mapping quality, and tags may be changed, and its decoded fixed fields.
///
/// A record whose CIGAR is a `kSmN` placeholder takes its CIGAR from the `CG` tag, as htslib
/// reads it.
#[derive(Clone, Debug)]
pub struct Read {
    pub(crate) raw: Vec<u8>,
    pub(crate) tid: i32,
    pub(crate) pos: i64,
    pub(crate) flag: u16,
    pub(crate) mtid: i32,
    pub(crate) mpos: i64,
    pub(crate) isize: i32,
    pub(crate) l_seq: usize,
    pub(crate) cigar: Vec<u32>,
    pub(crate) name_end: usize,
    pub(crate) seq_start: usize,
    pub(crate) qual_start: usize,
}

impl Read {
    /// Decodes the fixed fields of an encoded BAM record (without its `block_size`).
    ///
    /// # Errors
    ///
    /// Returns an error if the record is truncated.
    pub fn new(raw: Vec<u8>) -> io::Result<Self> {
        if raw.len() < NAME {
            return Err(truncated());
        }

        let u32_at =
            |offset: usize| u32::from_le_bytes(raw[offset..offset + 4].try_into().unwrap());
        let i32_at = |offset: usize| u32_at(offset).cast_signed();

        let name_len = usize::from(raw[NAME_LENGTH]);
        let n_cigar = usize::from(u16::from_le_bytes([raw[12], raw[13]]));
        let flag = u16::from_le_bytes([raw[14], raw[15]]);
        let l_seq = usize::try_from(u32_at(16))
            .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))?;

        let cigar_start = NAME + name_len;
        let seq_start = cigar_start + 4 * n_cigar;
        let qual_start = seq_start + l_seq.div_ceil(2);

        if name_len == 0 || qual_start + l_seq > raw.len() {
            return Err(truncated());
        }

        let mut cigar = decode_ops(&raw[cigar_start..seq_start]);

        if cigar.len() == 2
            && cigar[0] & 0xf == CIGAR_SOFT_CLIP
            && (cigar[0] >> 4) as usize == l_seq
            && cigar[1] & 0xf == CIGAR_SKIP
            && let Some(range) = find_field(&raw, qual_start + l_seq, *b"CG")?
            && raw[range.start + 2] == b'B'
            && raw[range.start + 3] == b'I'
        {
            cigar = decode_ops(&raw[range.start + 8..range.end]);
        }

        Ok(Self {
            tid: i32_at(0),
            pos: i64::from(i32_at(4)),
            flag,
            mtid: i32_at(20),
            mpos: i64::from(i32_at(24)),
            isize: i32_at(28),
            l_seq,
            cigar,
            name_end: cigar_start - 1,
            seq_start,
            qual_start,
            raw,
        })
    }

    /// The encoded record.
    pub fn as_bytes(&self) -> &[u8] {
        &self.raw
    }

    /// Returns the encoded record.
    pub fn into_bytes(self) -> Vec<u8> {
        self.raw
    }

    /// The reference sequence ID, or -1.
    pub fn reference_sequence_id(&self) -> i32 {
        self.tid
    }

    /// The 0-based position, or -1.
    pub fn position(&self) -> i64 {
        self.pos
    }

    /// The SAM flags.
    pub fn flags(&self) -> u16 {
        self.flag
    }

    /// The mate's reference sequence ID, or -1.
    pub fn mate_reference_sequence_id(&self) -> i32 {
        self.mtid
    }

    /// The mate's 0-based position, or -1.
    pub fn mate_position(&self) -> i64 {
        self.mpos
    }

    /// The template length.
    pub fn template_length(&self) -> i32 {
        self.isize
    }

    /// The number of bases.
    pub fn sequence_len(&self) -> usize {
        self.l_seq
    }

    /// The CIGAR operations, encoded as in BAM (`len << 4 | op`).
    pub fn cigar(&self) -> &[u32] {
        &self.cigar
    }

    /// The read name.
    pub fn name(&self) -> &[u8] {
        &self.raw[NAME..self.name_end]
    }

    /// The mapping quality.
    pub fn mapping_quality(&self) -> u8 {
        self.raw[9]
    }

    /// Sets the mapping quality.
    pub fn set_mapping_quality(&mut self, mapq: u8) {
        self.raw[9] = mapq;
    }

    /// Whether the read is on the reverse strand.
    pub fn is_reverse_complemented(&self) -> bool {
        self.flag & 0x10 != 0
    }

    /// The 4-bit code (`=ACMGRSVTWYHKDBN`) of base `i`.
    ///
    /// # Panics
    ///
    /// Panics if `i` is out of bounds.
    pub fn base(&self, i: usize) -> u8 {
        assert!(i < self.l_seq, "base {i} of {}", self.l_seq);
        (self.raw[self.seq_start + i / 2] >> ((!i & 1) << 2)) & 0xf
    }

    /// Sets the 4-bit code of base `i`.
    ///
    /// # Panics
    ///
    /// Panics if `i` is out of bounds.
    pub fn set_base(&mut self, i: usize, code: u8) {
        assert!(i < self.l_seq, "base {i} of {}", self.l_seq);
        let byte = &mut self.raw[self.seq_start + i / 2];

        if i & 1 == 0 {
            *byte = (*byte & 0x0f) | (code << 4);
        } else {
            *byte = (*byte & 0xf0) | (code & 0x0f);
        }
    }

    /// The base qualities.
    pub fn quality_scores(&self) -> &[u8] {
        &self.raw[self.qual_start..self.qual_start + self.l_seq]
    }

    /// The base qualities, to change.
    pub fn quality_scores_mut(&mut self) -> &mut [u8] {
        &mut self.raw[self.qual_start..self.qual_start + self.l_seq]
    }

    /// The type and value bytes of the first auxiliary field with `tag`.
    pub fn tag(&self, tag: [u8; 2]) -> Option<(u8, &[u8])> {
        let range = self.tag_range(tag)?;
        Some((
            self.raw[range.start + 2],
            &self.raw[range.start + 3..range.end],
        ))
    }

    /// The value of the first string (`Z` or `H`) field with `tag`, without its NUL.
    pub fn string_tag(&self, tag: [u8; 2]) -> Option<&[u8]> {
        match self.tag(tag)? {
            (b'Z' | b'H', value) => Some(&value[..value.len() - 1]),
            _ => None,
        }
    }

    /// The value of the first integer field with `tag` (htslib's `bam_aux2i`: 0 for fields of
    /// other types).
    pub fn integer_tag(&self, tag: [u8; 2]) -> Option<i64> {
        let (ty, v) = self.tag(tag)?;

        Some(match ty {
            b'c' => i64::from(v[0].cast_signed()),
            b'C' => i64::from(v[0]),
            b's' => i64::from(i16::from_le_bytes([v[0], v[1]])),
            b'S' => i64::from(u16::from_le_bytes([v[0], v[1]])),
            b'i' => i64::from(i32::from_le_bytes([v[0], v[1], v[2], v[3]])),
            b'I' => i64::from(u32::from_le_bytes([v[0], v[1], v[2], v[3]])),
            _ => 0,
        })
    }

    /// Removes the first auxiliary field with `tag`.
    pub fn remove_tag(&mut self, tag: [u8; 2]) {
        if let Some(range) = self.tag_range(tag) {
            self.raw.drain(range);
        }
    }

    /// Renames the first auxiliary field with `tag`.
    pub fn rename_tag(&mut self, tag: [u8; 2], new_tag: [u8; 2]) {
        if let Some(range) = self.tag_range(tag) {
            self.raw[range.start..range.start + 2].copy_from_slice(&new_tag);
        }
    }

    /// Appends a field of type `ty` with the encoded `value` (htslib's `bam_aux_append`).
    pub fn append_tag(&mut self, tag: [u8; 2], ty: u8, value: &[u8]) {
        self.raw.extend_from_slice(&tag);
        self.raw.push(ty);
        self.raw.extend_from_slice(value);
    }

    /// Appends a string (`Z`) field.
    pub fn append_string_tag(&mut self, tag: [u8; 2], value: &[u8]) {
        self.append_tag(tag, b'Z', value);
        self.raw.push(0);
    }

    /// Removes every auxiliary field but the first with `tag` (htslib's `bam_aux_drop_other`).
    pub fn retain_only_tag(&mut self, tag: [u8; 2]) {
        let kept = self.tag_range(tag).map(|range| self.raw[range].to_vec());
        self.raw.truncate(self.data_start());

        if let Some(field) = kept {
            self.raw.extend_from_slice(&field);
        }
    }

    /// `bam_cigar2rlen`: the reference length of the CIGAR.
    pub(crate) fn reference_length(&self) -> i64 {
        self.cigar
            .iter()
            .filter(|&&op| matches!(op & 0xf, 0 | 2 | 3 | 7 | 8))
            .map(|&op| i64::from(op >> 4))
            .sum()
    }

    fn data_start(&self) -> usize {
        self.qual_start + self.l_seq
    }

    pub(crate) fn tag_range(&self, tag: [u8; 2]) -> Option<Range<usize>> {
        find_field(&self.raw, self.data_start(), tag).ok().flatten()
    }
}

impl TryFrom<bam::Record> for Read {
    type Error = io::Error;

    fn try_from(record: bam::Record) -> Result<Self, Self::Error> {
        Self::new(Vec::from(record))
    }
}

fn decode_ops(src: &[u8]) -> Vec<u32> {
    src.as_chunks::<4>()
        .0
        .iter()
        .map(|b| u32::from_le_bytes(*b))
        .collect()
}

/// The byte range of the first auxiliary field with `tag`, searching from `start`.
fn find_field(raw: &[u8], start: usize, tag: [u8; 2]) -> io::Result<Option<Range<usize>>> {
    let mut i = start;

    while i < raw.len() {
        let src = &raw[i..];

        let &[t0, t1, ty, ref rest @ ..] = src else {
            return Err(truncated());
        };

        let value_len = match ty {
            b'A' | b'c' | b'C' => 1,
            b's' | b'S' => 2,
            b'i' | b'I' | b'f' => 4,
            b'Z' | b'H' => rest.iter().position(|&b| b == 0).ok_or_else(truncated)? + 1,
            b'B' => {
                let &[subtype, n0, n1, n2, n3, ..] = rest else {
                    return Err(truncated());
                };

                let width = match subtype {
                    b'c' | b'C' => 1,
                    b's' | b'S' => 2,
                    b'i' | b'I' | b'f' => 4,
                    _ => return Err(invalid("invalid array subtype")),
                };

                let count = usize::try_from(u32::from_le_bytes([n0, n1, n2, n3]))
                    .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))?;

                count
                    .checked_mul(width)
                    .and_then(|n| n.checked_add(5))
                    .ok_or_else(|| invalid("array too long"))?
            }
            _ => return Err(invalid("invalid field type")),
        };

        let len = 3 + value_len;

        if len > src.len() {
            return Err(truncated());
        }

        if [t0, t1] == tag {
            return Ok(Some(i..i + len));
        }

        i += len;
    }

    Ok(None)
}

fn truncated() -> io::Error {
    io::Error::new(io::ErrorKind::UnexpectedEof, "truncated BAM record")
}

fn invalid(msg: &str) -> io::Error {
    io::Error::new(io::ErrorKind::InvalidData, msg.to_owned())
}

#[cfg(test)]
mod tests {
    use super::*;

    /// A record named `r1` at 100 with CIGAR 4M, sequence ACGT, and qualities 30.
    pub(crate) fn record() -> Vec<u8> {
        let mut raw = Vec::new();
        raw.extend(0i32.to_le_bytes()); // ref_id
        raw.extend(99i32.to_le_bytes()); // pos
        raw.push(3); // l_read_name
        raw.push(60); // mapq
        raw.extend(0u16.to_le_bytes()); // bin
        raw.extend(1u16.to_le_bytes()); // n_cigar_op
        raw.extend(0u16.to_le_bytes()); // flag
        raw.extend(4u32.to_le_bytes()); // l_seq
        raw.extend((-1i32).to_le_bytes()); // next_ref_id
        raw.extend((-1i32).to_le_bytes()); // next_pos
        raw.extend(0i32.to_le_bytes()); // tlen
        raw.extend(b"r1\0");
        raw.extend((4u32 << 4).to_le_bytes()); // 4M
        raw.extend([0x12, 0x48]); // ACGT
        raw.extend([30; 4]);
        raw
    }

    #[test]
    fn test_tags() -> io::Result<()> {
        let mut read = Read::new(record())?;
        assert_eq!(read.name(), b"r1");
        assert_eq!(read.cigar(), [4 << 4]);
        assert_eq!(
            (0..4).map(|i| read.base(i)).collect::<Vec<_>>(),
            [1, 2, 4, 8]
        );

        read.append_tag(*b"NM", b'C', &[2]);
        read.append_string_tag(*b"MD", b"4");
        read.append_string_tag(*b"RG", b"a");
        assert_eq!(read.integer_tag(*b"NM"), Some(2));
        assert_eq!(read.string_tag(*b"MD"), Some(&b"4"[..]));

        read.rename_tag(*b"MD", *b"XM");
        assert_eq!(read.string_tag(*b"XM"), Some(&b"4"[..]));

        read.retain_only_tag(*b"RG");
        assert_eq!(read.tag(*b"NM"), None);
        assert_eq!(read.string_tag(*b"RG"), Some(&b"a"[..]));

        read.set_base(1, 0);
        assert_eq!(read.base(1), 0);
        assert_eq!(read.base(0), 1);

        Ok(())
    }
}
