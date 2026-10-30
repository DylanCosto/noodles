//! Field access and editing on encoded (raw) BAM records, excluding `block_size`.

use std::{io, ops::Range};

/// Offsets of the fixed-length fields.
pub mod field {
    use std::ops::Range;

    pub const REFERENCE_SEQUENCE_ID: Range<usize> = 0..4;
    pub const POSITION: Range<usize> = 4..8;
    pub const NAME_LENGTH: usize = 8;
    pub const MAPPING_QUALITY: usize = 9;
    pub const CIGAR_OP_COUNT: Range<usize> = 12..14;
    pub const FLAGS: Range<usize> = 14..16;
    pub const SEQUENCE_LENGTH: Range<usize> = 16..20;
    pub const MATE_REFERENCE_SEQUENCE_ID: Range<usize> = 20..24;
    pub const MATE_POSITION: Range<usize> = 24..28;
    pub const TEMPLATE_LENGTH: Range<usize> = 28..32;
    pub const NAME: usize = 32;
}

/// SAM's bases by their 4-bit BAM codes (HTSlib's `seq_nt16_str`).
pub const BASES: &[u8; 16] = b"=ACMGRSVTWYHKDBN";

/// Decodes a byte of two 4-bit bases to its two characters.
pub static BASE_PAIRS: [[u8; 2]; 256] = {
    let mut table = [[0; 2]; 256];
    let mut i = 0;
    while i < 256 {
        table[i] = [BASES[i >> 4], BASES[i & 0x0f]];
        i += 1;
    }
    table
};

/// A view of an encoded BAM record whose variable-length sections were checked to fit.
#[derive(Clone, Copy, Debug)]
pub struct Record<'a> {
    buf: &'a [u8],
    seq_start: usize,
    qual_start: usize,
    data_start: usize,
}

impl<'a> Record<'a> {
    /// Wraps an encoded record (without `block_size`).
    ///
    /// # Errors
    ///
    /// Returns an error if the record is shorter than its fixed fields, read name, CIGAR,
    /// sequence, and qualities.
    pub fn new(buf: &'a [u8]) -> io::Result<Self> {
        if buf.len() < field::NAME {
            return Err(truncated());
        }

        let name_len = usize::from(buf[field::NAME_LENGTH]);
        let n_cigar = usize::from(le_u16(&buf[field::CIGAR_OP_COUNT]));
        let l_seq = usize::try_from(le_u32(&buf[field::SEQUENCE_LENGTH]))
            .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))?;

        let seq_start = field::NAME + name_len + 4 * n_cigar;
        let qual_start = seq_start + l_seq.div_ceil(2);
        let data_start = qual_start + l_seq;

        if name_len == 0 || data_start > buf.len() {
            return Err(truncated());
        }

        Ok(Self {
            buf,
            seq_start,
            qual_start,
            data_start,
        })
    }

    /// The encoded record.
    #[inline]
    pub fn as_bytes(&self) -> &'a [u8] {
        self.buf
    }

    #[inline]
    fn i32_at(&self, offset: usize) -> i32 {
        i32::from_le_bytes([
            self.buf[offset],
            self.buf[offset + 1],
            self.buf[offset + 2],
            self.buf[offset + 3],
        ])
    }

    /// The reference sequence ID (-1 when unplaced).
    #[inline]
    pub fn reference_sequence_id(&self) -> i32 {
        self.i32_at(field::REFERENCE_SEQUENCE_ID.start)
    }

    /// The 0-based position (-1 when unplaced).
    #[inline]
    pub fn position(&self) -> i32 {
        self.i32_at(field::POSITION.start)
    }

    #[inline]
    pub fn mapping_quality(&self) -> u8 {
        self.buf[field::MAPPING_QUALITY]
    }

    #[inline]
    pub fn flags(&self) -> u16 {
        le_u16(&self.buf[field::FLAGS])
    }

    /// The mate's reference sequence ID (-1 when unplaced).
    #[inline]
    pub fn mate_reference_sequence_id(&self) -> i32 {
        self.i32_at(field::MATE_REFERENCE_SEQUENCE_ID.start)
    }

    /// The mate's 0-based position (-1 when unplaced).
    #[inline]
    pub fn mate_position(&self) -> i32 {
        self.i32_at(field::MATE_POSITION.start)
    }

    #[inline]
    pub fn template_length(&self) -> i32 {
        self.i32_at(field::TEMPLATE_LENGTH.start)
    }

    /// The read name, without its NUL terminator.
    #[inline]
    pub fn name(&self) -> &'a [u8] {
        let len = usize::from(self.buf[field::NAME_LENGTH]);
        &self.buf[field::NAME..field::NAME + len - 1]
    }

    /// The encoded CIGAR operations, as stored (see [`Self::read_cigar`] for the `CG` tag).
    #[inline]
    pub fn raw_cigar(&self) -> &'a [u8] {
        let name_len = usize::from(self.buf[field::NAME_LENGTH]);
        &self.buf[field::NAME + name_len..self.seq_start]
    }

    /// The byte range (in [`Self::as_bytes`]) of the `CG` field holding the CIGAR, when the
    /// stored CIGAR is the `kSmN` placeholder BAM uses for more than 65,535 operations.
    ///
    /// # Errors
    ///
    /// Returns an error if the auxiliary data is malformed.
    pub fn long_cigar_field(&self) -> io::Result<Option<Range<usize>>> {
        const SOFT_CLIP: u32 = 4;
        const SKIP: u32 = 3;

        let ops = self.raw_cigar().as_chunks::<4>().0;

        let [first, second] = ops else {
            return Ok(None);
        };

        let (first, second) = (u32::from_le_bytes(*first), u32::from_le_bytes(*second));

        if first & 0xf != SOFT_CLIP
            || first >> 4 != le_u32(&self.buf[field::SEQUENCE_LENGTH])
            || second & 0xf != SKIP
        {
            return Ok(None);
        }

        Ok(self
            .find_field(*b"CG")?
            .filter(|range| self.buf[range.start + 2] == b'B' && self.buf[range.start + 3] == b'I'))
    }

    /// Reads the CIGAR operations into `dst`, replacing a `kSmN` placeholder with the `CG` tag's
    /// operations as HTSlib does when it reads a record.
    ///
    /// # Errors
    ///
    /// Returns an error if the auxiliary data is malformed.
    pub fn read_cigar(&self, dst: &mut Vec<u32>) -> io::Result<()> {
        let ops = match self.long_cigar_field()? {
            Some(range) => &self.buf[range.start + 8..range.end],
            None => self.raw_cigar(),
        };

        dst.clear();
        dst.extend(
            ops.as_chunks::<4>()
                .0
                .iter()
                .map(|b| u32::from_le_bytes(*b)),
        );

        Ok(())
    }

    /// The number of bases.
    #[inline]
    pub fn sequence_len(&self) -> usize {
        self.data_start - self.qual_start
    }

    /// The bases, packed two per byte as 4-bit codes.
    #[inline]
    pub fn packed_sequence(&self) -> &'a [u8] {
        &self.buf[self.seq_start..self.qual_start]
    }

    /// The 4-bit code of base `i` (an index into [`BASES`]), for `i` less than
    /// [`Self::sequence_len`].
    ///
    /// # Panics
    ///
    /// Panics if `i` is past the packed sequence (an `i` equal to an odd sequence length reads
    /// the padding, and is caught only in debug builds).
    #[inline]
    pub fn base(&self, i: usize) -> u8 {
        debug_assert!(i < self.sequence_len(), "base {i} out of range");
        (self.packed_sequence()[i / 2] >> ((!i & 1) << 2)) & 0x0f
    }

    /// Appends the bases, decoded to characters, to `dst`.
    pub fn decode_sequence(&self, dst: &mut Vec<u8>) {
        let len = self.sequence_len();
        let packed = self.packed_sequence();

        // Writing into a presized slice avoids a capacity check per base pair.
        let start = dst.len();
        dst.resize(start + len, 0);
        let (pairs, rest) = dst[start..].as_chunks_mut::<2>();

        for (pair, &byte) in pairs.iter_mut().zip(packed) {
            *pair = BASE_PAIRS[usize::from(byte)];
        }

        if let [last] = rest {
            *last = BASE_PAIRS[usize::from(packed[len / 2])][0];
        }
    }

    #[inline]
    pub fn quality_scores(&self) -> &'a [u8] {
        &self.buf[self.qual_start..self.data_start]
    }

    /// The auxiliary data.
    #[inline]
    pub fn data(&self) -> &'a [u8] {
        &self.buf[self.data_start..]
    }

    /// The byte range (in [`Self::as_bytes`]) of the first auxiliary field with `tag`,
    /// including its tag and type.
    ///
    /// # Errors
    ///
    /// Returns an error if the auxiliary data is malformed.
    pub fn find_field(&self, tag: [u8; 2]) -> io::Result<Option<Range<usize>>> {
        find_field_from(self.buf, self.data_start, tag)
    }
}

#[inline]
pub fn le_u16(buf: &[u8]) -> u16 {
    u16::from_le_bytes([buf[0], buf[1]])
}

#[inline]
pub fn le_u32(buf: &[u8]) -> u32 {
    u32::from_le_bytes([buf[0], buf[1], buf[2], buf[3]])
}

/// Returns the byte range of the auxiliary data (optional fields).
pub fn data_range(record: &[u8]) -> io::Result<Range<usize>> {
    let name_len = usize::from(record[field::NAME_LENGTH]);
    let cigar_len = usize::from(le_u16(&record[field::CIGAR_OP_COUNT])) * 4;
    let base_count = usize::try_from(le_u32(&record[field::SEQUENCE_LENGTH]))
        .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))?;

    let start = field::NAME + name_len + cigar_len + base_count.div_ceil(2) + base_count;

    if start > record.len() {
        return Err(truncated());
    }

    Ok(start..record.len())
}

/// Returns the tag and byte length (including tag and type) of the field at the start of `src`.
fn next_field(src: &[u8]) -> io::Result<([u8; 2], usize)> {
    let &[t0, t1, ty, ref rest @ ..] = src else {
        return Err(truncated());
    };

    let value_len = match ty {
        b'A' | b'c' | b'C' => 1,
        b's' | b'S' => 2,
        b'i' | b'I' | b'f' => 4,
        // HTSlib's non-standard double.
        b'd' => 8,
        b'Z' | b'H' => memchr::memchr(0, rest).ok_or_else(truncated)? + 1,
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

    Ok(([t0, t1], len))
}

/// Writes `record` to `dst` without the auxiliary fields whose tags `remove` selects.
pub fn strip_fields<F>(record: &[u8], dst: &mut Vec<u8>, remove: F) -> io::Result<()>
where
    F: Fn([u8; 2]) -> bool,
{
    let data = data_range(record)?;

    dst.clear();
    dst.extend_from_slice(&record[..data.start]);

    let mut src = &record[data];

    while !src.is_empty() {
        let (tag, len) = next_field(src)?;

        if !remove(tag) {
            dst.extend_from_slice(&src[..len]);
        }

        src = &src[len..];
    }

    Ok(())
}

/// Returns the byte range of the first auxiliary field with `tag`, including its tag and type.
pub fn find_field(record: &[u8], tag: [u8; 2]) -> io::Result<Option<Range<usize>>> {
    let data = data_range(record)?;
    find_field_from(record, data.start, tag)
}

/// Returns the byte ranges of the first auxiliary fields with each of `tags`, in one pass.
pub fn find_fields<const N: usize>(
    record: &[u8],
    tags: [[u8; 2]; N],
) -> io::Result<[Option<Range<usize>>; N]> {
    let mut ranges = [const { None }; N];
    let mut start = data_range(record)?.start;

    while start < record.len() {
        let (field_tag, len) = next_field(&record[start..])?;

        for (tag, range) in tags.iter().zip(&mut ranges) {
            if *tag == field_tag && range.is_none() {
                *range = Some(start..start + len);
            }
        }

        start += len;
    }

    Ok(ranges)
}

/// [`find_field`] with the auxiliary data starting at `start`.
fn find_field_from(
    record: &[u8],
    mut start: usize,
    tag: [u8; 2],
) -> io::Result<Option<Range<usize>>> {
    while start < record.len() {
        let (field_tag, len) = next_field(&record[start..])?;

        if field_tag == tag {
            return Ok(Some(start..start + len));
        }

        start += len;
    }

    Ok(None)
}

/// Returns the value of the first field with `tag` if it is a string (`Z`) field, `Ok(None)` if
/// there is no such field, and an error if it has another type.
pub fn string_field(record: &[u8], tag: [u8; 2]) -> io::Result<Option<&[u8]>> {
    match find_field(record, tag)? {
        Some(range) if record[range.start + 2] == b'Z' => {
            Ok(Some(&record[range.start + 3..range.end - 1]))
        }
        Some(_) => Err(invalid("not a string field")),
        None => Ok(None),
    }
}

/// Sets a string (`Z`) field as HTSlib's `bam_aux_update_str` does: an existing field is replaced
/// in place, a new one is appended, and an existing field of another type is an error.
pub fn update_string_field(record: &mut Vec<u8>, tag: [u8; 2], value: &[u8]) -> io::Result<()> {
    const NUL: u8 = 0x00;

    match find_field(record, tag)? {
        Some(range) if record[range.start + 2] == b'Z' => {
            record.splice(range.start + 3..range.end - 1, value.iter().copied());
        }
        Some(_) => return Err(invalid("not a string field")),
        None => {
            record.extend_from_slice(&tag);
            record.push(b'Z');
            record.extend_from_slice(value);
            record.push(NUL);
        }
    }

    Ok(())
}

/// Unmaps a record as `samtools view -p` does for records that fail the filters: sets the
/// unmapped flag, and clears the mapping quality, template length, and CIGAR (a long CIGAR's `CG`
/// field too, which HTSlib would have moved into the CIGAR on reading).
///
/// # Errors
///
/// Returns an error if the record is malformed.
pub fn unmap(record: &mut Vec<u8>) -> io::Result<()> {
    const FLAG_UNMAPPED: u16 = 0x4;

    if let Some(range) = Record::new(record)?.long_cigar_field()? {
        record.drain(range);
    }

    let flags = le_u16(&record[field::FLAGS]) | FLAG_UNMAPPED;
    record[field::FLAGS].copy_from_slice(&flags.to_le_bytes());
    record[field::MAPPING_QUALITY] = 0;
    record[field::TEMPLATE_LENGTH].copy_from_slice(&0i32.to_le_bytes());

    let name_end = field::NAME + usize::from(record[field::NAME_LENGTH]);
    let cigar_len = 4 * usize::from(le_u16(&record[field::CIGAR_OP_COUNT]));
    record.drain(name_end..name_end + cigar_len);
    record[field::CIGAR_OP_COUNT].copy_from_slice(&0u16.to_le_bytes());

    Ok(())
}

/// Sets the `add` flag bits of a record, then clears the `remove` ones (samtools view's
/// `--add-flags` and `--remove-flags`).
///
/// # Errors
///
/// Returns an error if the record is too short to hold flags.
pub fn edit_flags(record: &mut [u8], add: u16, remove: u16) -> io::Result<()> {
    let field = record.get_mut(field::FLAGS).ok_or_else(truncated)?;
    let flags = (le_u16(field) | add) & !remove;
    field.copy_from_slice(&flags.to_le_bytes());
    Ok(())
}

/// `bam_remove_B` (samtools' `bam.c`, 1.24): collapses backward (`B`) CIGAR operations, merging
/// the overlapping bases into a consensus. A `B` that cannot be removed (first in the CIGAR, or
/// longer than the query before it) marks the record unmapped instead.
///
/// Records with a long CIGAR (in a `CG` field) are left as they are.
///
/// # Errors
///
/// Returns an error if the record is malformed.
pub fn remove_backward_ops(record: &mut Vec<u8>) -> io::Result<()> {
    const FLAG_UNMAPPED: u16 = 0x4;
    const BACK: u32 = 9;

    let view = Record::new(record)?;

    if view.flags() & FLAG_UNMAPPED != 0 || view.long_cigar_field()?.is_some() {
        return Ok(());
    }

    let mut cigar = Vec::new();
    view.read_cigar(&mut cigar)?;

    if !cigar.iter().any(|&op| op & 0xf == BACK) {
        return Ok(());
    }

    let set_unmapped = |record: &mut Vec<u8>| {
        let flags = le_u16(&record[field::FLAGS]) | FLAG_UNMAPPED;
        record[field::FLAGS].copy_from_slice(&flags.to_le_bytes());
    };

    if cigar[0] & 0xf == BACK {
        set_unmapped(record);
        return Ok(());
    }

    let consumes_query = |op: u32| (0x3c1a7 >> ((op & 0xf) << 1)) & 1 != 0;

    let mut seq: Vec<u8> = (0..view.sequence_len()).map(|i| view.base(i)).collect();
    let mut qual = view.quality_scores().to_vec();
    let no_qual = qual.first() == Some(&0xff);

    let mut new_cigar: Vec<u32> = Vec::with_capacity(cigar.len());
    let (mut i, mut j, mut end_j) = (0usize, 0usize, None::<usize>);

    for (k, &op) in cigar.iter().enumerate() {
        let len = (op >> 4) as usize;

        if op & 0xf == BACK {
            // A trailing `B` is ignored.
            if k == cigar.len() - 1 {
                break;
            }

            if len > j {
                set_unmapped(record);
                return Ok(());
            }

            // Look back through the query-consuming operations, shortening the one the `B`
            // ends in and dropping those after it.
            let mut u = 0;
            let mut t = new_cigar.len();

            while t > 0 {
                let prev = new_cigar[t - 1];
                let prev_len = (prev >> 4) as usize;

                if consumes_query(prev) {
                    if u + prev_len >= len {
                        new_cigar[t - 1] -= ((len - u) as u32) << 4;
                        break;
                    }

                    u += prev_len;
                }

                t -= 1;
            }

            // Squeeze out an operation shortened to nothing.
            if t > 0 && new_cigar[t - 1] >> 4 == 0 {
                t -= 1;
            }

            new_cigar.truncate(t);
            end_j = Some(j);
            j -= len;
        } else {
            new_cigar.push(op);

            if consumes_query(op) {
                if i != j {
                    for u in 0..len {
                        let c = seq[i + u];

                        if end_j.is_some_and(|end| j + u < end) {
                            // In the overlap: keep the better base.
                            if c == seq[j + u] {
                                qual[j + u] = qual[j + u].max(qual[i + u]);
                            } else if qual[j + u] < qual[i + u] {
                                seq[j + u] = c;
                                qual[j + u] = qual[i + u] - qual[j + u];
                            } else {
                                qual[j + u] -= qual[i + u];
                            }
                        } else {
                            seq[j + u] = c;
                            qual[j + u] = qual[i + u];
                        }
                    }
                }

                i += len;
                j += len;
            }
        }
    }

    // A record without qualities stays without them. HTSlib resets only the first byte (which
    // marks them missing), leaving the merged values of the others.
    if no_qual {
        qual.fill(0xff);
    }

    // Merge adjacent operations of the same kind, then drop empty ones.
    for k in 1..new_cigar.len() {
        if new_cigar[k] & 0xf == new_cigar[k - 1] & 0xf {
            new_cigar[k] += new_cigar[k - 1] & !0xf;
            new_cigar[k - 1] &= 0xf;
        }
    }

    new_cigar.retain(|op| op >> 4 != 0);

    let name_end = field::NAME + usize::from(record[field::NAME_LENGTH]);
    let data = record[view.data_start..].to_vec();

    record.truncate(name_end);
    record.extend(new_cigar.iter().flat_map(|op| op.to_le_bytes()));
    record.extend(seq[..j].chunks(2).map(|pair| {
        let (hi, lo) = (pair[0], pair.get(1).copied().unwrap_or(0));
        (hi << 4) | lo
    }));
    record.extend_from_slice(&qual[..j]);
    record.extend_from_slice(&data);

    let n_cigar = u16::try_from(new_cigar.len()).map_err(|_| invalid("CIGAR too long"))?;
    record[field::CIGAR_OP_COUNT].copy_from_slice(&n_cigar.to_le_bytes());
    let l_seq = u32::try_from(j).map_err(|_| invalid("sequence too long"))?;
    record[field::SEQUENCE_LENGTH].copy_from_slice(&l_seq.to_le_bytes());

    Ok(())
}

fn truncated() -> io::Error {
    io::Error::new(io::ErrorKind::UnexpectedEof, "truncated BAM record")
}

fn invalid(msg: &str) -> io::Error {
    io::Error::new(io::ErrorKind::InvalidData, msg.to_owned())
}

#[cfg(test)]
mod tests {
    use noodles_bam as bam;
    use noodles_sam::{
        self as sam,
        alignment::{
            RecordBuf,
            io::Write as _,
            record::data::field::Tag,
            record_buf::data::field::{Value, value::Array},
        },
    };

    use super::*;

    fn encode(record: &RecordBuf) -> io::Result<Vec<u8>> {
        let mut writer = bam::io::Writer::from(Vec::new());
        writer.write_alignment_record(&sam::Header::default(), record)?;
        Ok(writer.get_ref()[4..].to_vec())
    }

    fn sam_record(line: &str) -> io::Result<Vec<u8>> {
        let header: sam::Header = "@SQ\tSN:c\tLN:100\n".parse().map_err(io::Error::other)?;
        let record = sam::Record::try_from(line.as_bytes())?;
        let mut writer = bam::io::Writer::from(Vec::new());
        writer.write_alignment_record(&header, &record)?;
        Ok(writer.get_ref()[4..].to_vec())
    }

    #[test]
    fn test_remove_backward_ops() -> io::Result<()> {
        // noodles cannot parse `B` in SAM, so it is written as `D` and patched.
        let with_back = |line: &str| -> io::Result<Vec<u8>> {
            const DELETION: u8 = 2;
            const BACK: u8 = 9;

            let mut buf = sam_record(line)?;
            let start = field::NAME + usize::from(buf[field::NAME_LENGTH]);
            let n = usize::from(le_u16(&buf[field::CIGAR_OP_COUNT]));

            for op in buf[start..start + 4 * n].chunks_mut(4) {
                if op[0] & 0xf == DELETION {
                    op[0] = (op[0] & !0xf) | BACK;
                }
            }

            Ok(buf)
        };

        // `samtools view -B` collapses 5M2B5M to 8M, merging the two overlapping bases.
        let mut buf = with_back("r\t0\tc\t1\t60\t5M2D5M\t*\t0\t0\tACGTAGTCCA\t(((((5((((")?;
        remove_backward_ops(&mut buf)?;
        assert_eq!(
            buf,
            sam_record("r\t0\tc\t1\t60\t8M\t*\t0\t0\tACGGACCA\t(((.!(((")?
        );

        // A leading `B` cannot be removed; the record is unmapped instead.
        let mut buf = with_back("r\t0\tc\t1\t60\t2D5M\t*\t0\t0\tACGTA\t*")?;
        remove_backward_ops(&mut buf)?;
        let mut expected = with_back("r\t0\tc\t1\t60\t2D5M\t*\t0\t0\tACGTA\t*")?;
        expected[field::FLAGS.start] |= 0x4;
        assert_eq!(buf, expected);

        Ok(())
    }

    #[test]
    fn test_unmap() -> io::Result<()> {
        let mut buf = sam_record("r\t99\tc\t5\t60\t2S3M\t=\t9\t20\tACGTA\t*\tNM:i:0")?;
        unmap(&mut buf)?;
        assert_eq!(
            buf,
            sam_record("r\t103\tc\t5\t0\t*\t=\t9\t0\tACGTA\t*\tNM:i:0")?
        );
        Ok(())
    }

    #[test]
    fn test_strip_fields() -> io::Result<()> {
        let data = |tags: &[&[u8; 2]]| {
            tags.iter()
                .map(|&&t| {
                    let value = match &t {
                        b"ZZ" => Value::from("text"),
                        b"BB" => Value::Array(Array::Int16(vec![1, -2, 3])),
                        _ => Value::from(7),
                    };
                    (Tag::from(t), value)
                })
                .collect()
        };

        let record = RecordBuf::builder()
            .set_name("r0")
            .set_sequence(b"ACGT".to_vec().into())
            .set_quality_scores(vec![30; 4].into())
            .set_data(data(&[b"NM", b"ZZ", b"BB", b"AS"]))
            .build();

        let mut dst = Vec::new();
        strip_fields(&encode(&record)?, &mut dst, |tag| {
            tag == *b"ZZ" || tag == *b"AS"
        })?;

        let expected = RecordBuf::builder()
            .set_name("r0")
            .set_sequence(b"ACGT".to_vec().into())
            .set_quality_scores(vec![30; 4].into())
            .set_data(data(&[b"NM", b"BB"]))
            .build();

        assert_eq!(dst, encode(&expected)?);

        let src = encode(&record)?;
        assert!(strip_fields(&src[..src.len() - 1], &mut dst, |_| false).is_err());

        Ok(())
    }

    #[test]
    fn test_update_string_field() -> io::Result<()> {
        let record = |data: Vec<(Tag, Value)>| {
            RecordBuf::builder()
                .set_name("r0")
                .set_data(data.into_iter().collect())
                .build()
        };

        let nm = (Tag::from(*b"NM"), Value::from(1));
        let as_ = (Tag::from(*b"AS"), Value::from(9));

        // An existing field is replaced in place.
        let mut buf = encode(&record(vec![
            nm.clone(),
            (Tag::from(*b"dt"), Value::from("LB")),
            as_.clone(),
        ]))?;
        assert_eq!(string_field(&buf, *b"dt")?, Some(&b"LB"[..]));
        update_string_field(&mut buf, *b"dt", b"SQ!")?;
        let expected = record(vec![
            nm.clone(),
            (Tag::from(*b"dt"), Value::from("SQ!")),
            as_.clone(),
        ]);
        assert_eq!(buf, encode(&expected)?);

        // A new field is appended.
        let mut buf = encode(&record(vec![nm.clone()]))?;
        assert_eq!(string_field(&buf, *b"do")?, None);
        update_string_field(&mut buf, *b"do", b"r1")?;
        let expected = record(vec![nm.clone(), (Tag::from(*b"do"), Value::from("r1"))]);
        assert_eq!(buf, encode(&expected)?);

        // A field of another type is an error.
        let mut buf = encode(&record(vec![nm]))?;
        assert!(string_field(&buf, *b"NM").is_err());
        assert!(update_string_field(&mut buf, *b"NM", b"x").is_err());

        Ok(())
    }
}
