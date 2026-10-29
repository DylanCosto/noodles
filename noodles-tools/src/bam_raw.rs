//! Field access and editing on encoded (raw) BAM records, excluding `block_size`.

use std::{io, ops::Range};

/// Offsets of the fixed-length fields.
pub mod field {
    use std::ops::Range;

    pub const REFERENCE_SEQUENCE_ID: Range<usize> = 0..4;
    pub const POSITION: Range<usize> = 4..8;
    pub const NAME_LENGTH: usize = 8;
    pub const CIGAR_OP_COUNT: Range<usize> = 12..14;
    pub const FLAGS: Range<usize> = 14..16;
    pub const SEQUENCE_LENGTH: Range<usize> = 16..20;
    pub const MATE_REFERENCE_SEQUENCE_ID: Range<usize> = 20..24;
    pub const NAME: usize = 32;
}

pub fn le_u16(buf: &[u8]) -> u16 {
    u16::from_le_bytes([buf[0], buf[1]])
}

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
}
