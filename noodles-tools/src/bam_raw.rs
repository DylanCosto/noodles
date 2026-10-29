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

/// Reads a record's CIGAR operations into `dst`, replacing a `kSmN` placeholder with the `CG`
/// tag's operations as HTSlib does when it reads a record.
pub fn read_cigar(record: &[u8], dst: &mut Vec<u32>) -> io::Result<()> {
    const SOFT_CLIP: u32 = 4;
    const SKIP: u32 = 3;

    let name_len = usize::from(record[field::NAME_LENGTH]);
    let n_cigar = usize::from(le_u16(&record[field::CIGAR_OP_COUNT]));
    let start = field::NAME + name_len;
    let raw = record
        .get(start..start + 4 * n_cigar)
        .ok_or_else(truncated)?;

    dst.clear();
    dst.extend(
        raw.as_chunks::<4>()
            .0
            .iter()
            .map(|b| u32::from_le_bytes(*b)),
    );

    let l_seq = le_u32(&record[field::SEQUENCE_LENGTH]);

    if dst.len() == 2
        && dst[0] & 0xf == SOFT_CLIP
        && dst[0] >> 4 == l_seq
        && dst[1] & 0xf == SKIP
        && let Some(range) = find_field(record, *b"CG")?
        && record[range.start + 2] == b'B'
        && record[range.start + 3] == b'I'
    {
        dst.clear();
        dst.extend(
            record[range.start + 8..range.end]
                .as_chunks::<4>()
                .0
                .iter()
                .map(|b| u32::from_le_bytes(*b)),
        );
    }

    Ok(())
}

/// Returns the byte range of the first auxiliary field with `tag`, including its tag and type.
pub fn find_field(record: &[u8], tag: [u8; 2]) -> io::Result<Option<Range<usize>>> {
    let data = data_range(record)?;
    let mut start = data.start;

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
