//! Formats encoded BAM records directly as SAM text.
//!
//! This is the common-case fast path behind SAM output. Any record outside it (an invalid name or
//! string, out-of-range qualities, or a long CIGAR stored in a `CG` tag) is declined, and the
//! caller formats it with the noodles SAM writer instead, so output and errors are unchanged.

use crate::bam_raw::{BASE_PAIRS, field, le_u16, le_u32};

const CIGAR_OPS: &[u8; 9] = b"MIDNSHP=X";

/// `,N` for every byte value, with its length, for formatting `B:C` arrays (e.g., `ML`).
static U8_ITEMS: [([u8; 4], u8); 256] = {
    let mut table = [([0; 4], 0); 256];
    let mut n = 0;

    while n < 256 {
        let mut buf = [b',', 0, 0, 0];
        let len = if n >= 100 {
            buf[1] = b'0' + (n / 100) as u8;
            buf[2] = b'0' + (n / 10 % 10) as u8;
            buf[3] = b'0' + (n % 10) as u8;
            4
        } else if n >= 10 {
            buf[1] = b'0' + (n / 10) as u8;
            buf[2] = b'0' + (n % 10) as u8;
            3
        } else {
            buf[1] = b'0' + n as u8;
            2
        };

        table[n] = (buf, len);
        n += 1;
    }

    table
};

/// The largest quality score SAM can represent (`~` is 93 + 33).
const MAX_QUALITY_SCORE: u8 = 93;

/// Reference sequence names by ID.
#[derive(Clone, Debug)]
pub struct Names(Vec<Box<[u8]>>);

impl Names {
    pub fn new(header: &noodles_sam::Header) -> Self {
        Self(
            header
                .reference_sequences()
                .keys()
                .map(|name| name.to_vec().into_boxed_slice())
                .collect(),
        )
    }

    /// Returns `Some(None)` for no reference (`-1`), and `None` for an out-of-range ID.
    fn get(&self, id: i32) -> Option<Option<&[u8]>> {
        match usize::try_from(id) {
            Ok(i) => self.0.get(i).map(|name| Some(&name[..])),
            Err(_) => Some(None),
        }
    }
}

/// Appends `record` as a SAM line (with newline) to `dst`, or returns `false`, leaving `dst`
/// unchanged, if the record needs the generic writer.
pub fn format_record(names: &Names, record: &[u8], dst: &mut Vec<u8>) -> bool {
    let start = dst.len();

    if try_format(names, record, dst).is_none() {
        dst.truncate(start);
        return false;
    }

    true
}

fn try_format(names: &Names, record: &[u8], dst: &mut Vec<u8>) -> Option<()> {
    if record.len() < field::NAME {
        return None;
    }

    let id = le_u32(&record[field::REFERENCE_SEQUENCE_ID]).cast_signed();
    let position = le_u32(&record[field::POSITION]).cast_signed();
    let name_len = usize::from(record[field::NAME_LENGTH]);
    let mapping_quality = record[9];
    let cigar_len = usize::from(le_u16(&record[field::CIGAR_OP_COUNT]));
    let flags = le_u16(&record[field::FLAGS]);
    let base_count = usize::try_from(le_u32(&record[field::SEQUENCE_LENGTH])).ok()?;
    let mate_id = le_u32(&record[field::MATE_REFERENCE_SEQUENCE_ID]).cast_signed();
    let mate_position = le_u32(&record[24..28]).cast_signed();
    let template_length = le_u32(&record[28..32]).cast_signed();

    let name_end = field::NAME.checked_add(name_len)?;
    let cigar_end = name_end.checked_add(cigar_len.checked_mul(4)?)?;
    let sequence_end = cigar_end.checked_add(base_count.div_ceil(2))?;
    let quality_end = sequence_end.checked_add(base_count)?;

    if quality_end > record.len() || name_len == 0 {
        return None;
    }

    // QNAME, without its NUL terminator. noodles rejects names outside `!`..`~` or with `@`.
    let name = &record[field::NAME..name_end - 1];
    if name.is_empty() || !name.iter().all(|&c| c.is_ascii_graphic() && c != b'@') {
        return None;
    }
    dst.extend_from_slice(name);

    dst.push(b'\t');
    write_u64(dst, u64::from(flags));

    let reference_name = names.get(id)?;
    dst.push(b'\t');
    dst.extend_from_slice(reference_name.unwrap_or(b"*"));

    dst.push(b'\t');
    write_position(dst, position)?;

    dst.push(b'\t');
    write_u64(dst, u64::from(mapping_quality));

    // CIGAR. A two-op placeholder may stand for a long CIGAR in a `CG` tag; leave those to noodles.
    dst.push(b'\t');
    let cigar = &record[name_end..cigar_end];
    if cigar_len == 0 {
        dst.push(b'*');
    } else {
        // The placeholder is `<base count>S<span>N`.
        if cigar_len == 2 {
            let (first, second) = (le_u32(&cigar[0..4]), le_u32(&cigar[4..8]));
            let soft_clips_read =
                first & 0x0f == 4 && usize::try_from(first >> 4).ok()? == base_count;

            if soft_clips_read && second & 0x0f == 3 {
                return None;
            }
        }

        let (ops, _) = cigar.as_chunks::<4>();

        for op in ops {
            let n = u32::from_le_bytes(*op);
            let kind = *CIGAR_OPS.get(usize::try_from(n & 0x0f).ok()?)?;
            write_u64(dst, u64::from(n >> 4));
            dst.push(kind);
        }
    }

    dst.push(b'\t');
    let mate_name = names.get(mate_id)?;
    match mate_name {
        Some(_) if mate_id == id => dst.push(b'='),
        Some(name) => dst.extend_from_slice(name),
        None => dst.push(b'*'),
    }

    dst.push(b'\t');
    write_position(dst, mate_position)?;

    dst.push(b'\t');
    write_i64(dst, i64::from(template_length));

    // SEQ
    dst.push(b'\t');
    if base_count == 0 {
        dst.push(b'*');
    } else {
        let packed = &record[cigar_end..sequence_end];

        // Writing into a presized slice avoids a capacity check per base pair.
        let start = dst.len();
        dst.resize(start + base_count, 0);
        let out = &mut dst[start..];

        let (pairs, rest) = out.as_chunks_mut::<2>();
        for (pair, &byte) in pairs.iter_mut().zip(packed) {
            *pair = BASE_PAIRS[usize::from(byte)];
        }

        if let [last] = rest {
            *last = BASE_PAIRS[usize::from(packed[base_count / 2])][0];
        }
    }

    // QUAL
    dst.push(b'\t');
    let quality_scores = &record[sequence_end..quality_end];
    if base_count == 0 || quality_scores[0] == 0xff {
        if base_count > 0 && quality_scores.iter().any(|&q| q != 0xff) {
            return None;
        }
        dst.push(b'*');
    } else {
        // Folded rather than short-circuited, so that the check vectorizes.
        if quality_scores
            .iter()
            .fold(false, |any, &q| any | (q > MAX_QUALITY_SCORE))
        {
            return None;
        }
        let start = dst.len();
        dst.resize(start + base_count, 0);

        for (d, &q) in dst[start..].iter_mut().zip(quality_scores) {
            *d = q + b'!';
        }
    }

    write_data(dst, &record[quality_end..])?;

    dst.push(b'\n');

    Some(())
}

fn write_position(dst: &mut Vec<u8>, position: i32) -> Option<()> {
    let n = i64::from(position) + 1;

    if n < 0 {
        return None;
    }

    write_i64(dst, n);
    Some(())
}

fn write_u64(dst: &mut Vec<u8>, n: u64) {
    let mut buf = [0; 20];
    dst.extend_from_slice(lexical_core::write(n, &mut buf));
}

fn write_i64(dst: &mut Vec<u8>, n: i64) {
    let mut buf = [0; 20];
    dst.extend_from_slice(lexical_core::write(n, &mut buf));
}

fn write_f32(dst: &mut Vec<u8>, n: f32) -> Option<()> {
    // The same options as the noodles SAM writer.
    const FORMAT: u128 = lexical_core::format::STANDARD;
    const OPTIONS: lexical_core::WriteFloatOptions = lexical_core::WriteFloatOptionsBuilder::new()
        .trim_floats(true)
        .build_strict();

    if !n.is_finite() {
        return None;
    }

    let mut buf = [0; 64];
    dst.extend_from_slice(lexical_core::write_with_options::<_, FORMAT>(
        n, &mut buf, &OPTIONS,
    ));
    Some(())
}

fn write_data(dst: &mut Vec<u8>, mut src: &[u8]) -> Option<()> {
    while !src.is_empty() {
        let &[t0, t1, ty, ref rest @ ..] = src else {
            return None;
        };

        dst.extend_from_slice(&[b'\t', t0, t1, b':']);

        let consumed = match ty {
            b'A' => {
                let c = *rest.first()?;
                if !c.is_ascii_graphic() {
                    return None;
                }
                dst.extend_from_slice(&[b'A', b':', c]);
                1
            }
            b'c' | b'C' | b's' | b'S' | b'i' | b'I' => {
                let (n, len) = read_integer(ty, rest)?;
                dst.extend_from_slice(b"i:");
                write_i64(dst, n);
                len
            }
            b'f' => {
                let n = f32::from_le_bytes(rest.get(..4)?.try_into().ok()?);
                dst.extend_from_slice(b"f:");
                write_f32(dst, n)?;
                4
            }
            // HTSlib's non-standard double, printed with `%g`.
            b'd' => {
                let n = f64::from_le_bytes(rest.get(..8)?.try_into().ok()?);
                dst.extend_from_slice(b"d:");
                dst.extend_from_slice(crate::depth::format_g(n, 6).as_bytes());
                8
            }
            b'Z' | b'H' => {
                let len = memchr::memchr(0, rest)?;
                let s = &rest[..len];

                // Written as a range test so that the check over long strings (e.g., `MM`)
                // vectorizes.
                let valid = if ty == b'Z' {
                    s.iter()
                        .fold(true, |valid, &c| valid & (b' '..=b'~').contains(&c))
                } else {
                    s.len() % 2 == 0 && s.iter().all(|c| matches!(c, b'0'..=b'9' | b'A'..=b'F'))
                };

                if !valid {
                    return None;
                }

                dst.extend_from_slice(&[ty, b':']);
                dst.extend_from_slice(s);
                len + 1
            }
            b'B' => write_array(dst, rest)?,
            _ => return None,
        };

        src = rest.get(consumed..)?;
    }

    Some(())
}

/// Reads an integer of the given BAM type, returning it and its byte length.
fn read_integer(ty: u8, src: &[u8]) -> Option<(i64, usize)> {
    Some(match ty {
        b'c' => (i64::from((*src.first()?).cast_signed()), 1),
        b'C' => (i64::from(*src.first()?), 1),
        b's' => (
            i64::from(i16::from_le_bytes(src.get(..2)?.try_into().ok()?)),
            2,
        ),
        b'S' => (
            i64::from(u16::from_le_bytes(src.get(..2)?.try_into().ok()?)),
            2,
        ),
        b'i' => (
            i64::from(i32::from_le_bytes(src.get(..4)?.try_into().ok()?)),
            4,
        ),
        b'I' => (
            i64::from(u32::from_le_bytes(src.get(..4)?.try_into().ok()?)),
            4,
        ),
        _ => return None,
    })
}

/// Writes a `B` array value, returning the bytes consumed after the type.
fn write_array(dst: &mut Vec<u8>, src: &[u8]) -> Option<usize> {
    let &[subtype, n0, n1, n2, n3, ref values @ ..] = src else {
        return None;
    };

    let count = usize::try_from(u32::from_le_bytes([n0, n1, n2, n3])).ok()?;

    let width = match subtype {
        b'c' | b'C' => 1,
        b's' | b'S' => 2,
        b'i' | b'I' | b'f' => 4,
        _ => return None,
    };

    let len = count.checked_mul(width)?;
    let values = values.get(..len)?;

    dst.extend_from_slice(&[b'B', b':', subtype]);

    if subtype == b'C' {
        dst.reserve(values.len() * 4);

        for &n in values {
            let (buf, len) = &U8_ITEMS[usize::from(n)];
            dst.extend_from_slice(&buf[..usize::from(*len)]);
        }

        return Some(5 + len);
    }

    for value in values.chunks_exact(width) {
        dst.push(b',');

        if subtype == b'f' {
            write_f32(dst, f32::from_le_bytes(value.try_into().ok()?))?;
        } else {
            write_i64(dst, read_integer(subtype, value)?.0);
        }
    }

    Some(5 + len)
}

#[cfg(test)]
mod tests {
    use noodles_bam as bam;
    use noodles_sam::{
        self as sam,
        alignment::{
            RecordBuf,
            io::Write as _,
            record::{
                Flags, MappingQuality,
                cigar::{Op, op::Kind},
                data::field::Tag,
            },
            record_buf::data::field::{Value, value::Array},
        },
    };

    use super::*;

    fn header() -> sam::Header {
        "@SQ\tSN:sq0\tLN:1000\n@SQ\tSN:sq1\tLN:1000\n"
            .parse()
            .expect("valid header")
    }

    fn generic(header: &sam::Header, record: &RecordBuf) -> Vec<u8> {
        let mut writer = sam::io::Writer::new(Vec::new());
        writer
            .write_alignment_record(header, record)
            .expect("formats");
        writer.into_inner()
    }

    fn raw(header: &sam::Header, record: &RecordBuf) -> Vec<u8> {
        let mut writer = bam::io::Writer::from(Vec::new());
        writer
            .write_alignment_record(header, record)
            .expect("encodes");
        writer.get_ref()[4..].to_vec()
    }

    #[test]
    fn test_format_record_matches_generic_writer() {
        let header = header();
        let names = Names::new(&header);

        let data = [
            (Tag::from(*b"NM"), Value::from(3)),
            (Tag::from(*b"XN"), Value::Int8(-4)),
            (Tag::from(*b"XU"), Value::UInt32(4_000_000_000)),
            (Tag::from(*b"XA"), Value::Character(b'q')),
            (Tag::from(*b"XZ"), Value::from("a string with spaces")),
            (Tag::from(*b"XH"), Value::Hex(b"CAFE".into())),
            (Tag::from(*b"XF"), Value::Float(0.25)),
            (
                Tag::from(*b"XB"),
                Value::Array(Array::Int16(vec![-1, 2, 300])),
            ),
            (
                Tag::from(*b"XC"),
                Value::Array(Array::Float(vec![1.5, -0.125])),
            ),
        ]
        .into_iter()
        .collect();

        let records = [
            RecordBuf::default(),
            RecordBuf::builder()
                .set_name("r0")
                .set_flags(Flags::SEGMENTED | Flags::FIRST_SEGMENT)
                .set_reference_sequence_id(0)
                .set_alignment_start(noodles_core::Position::new(10).unwrap())
                .set_mapping_quality(MappingQuality::new(37).unwrap())
                .set_cigar(
                    [Op::new(Kind::SoftClip, 2), Op::new(Kind::Match, 3)]
                        .into_iter()
                        .collect(),
                )
                .set_mate_reference_sequence_id(0)
                .set_mate_alignment_start(noodles_core::Position::new(50).unwrap())
                .set_template_length(-45)
                .set_sequence(b"ACGTN".to_vec().into())
                .set_quality_scores(vec![0, 10, 20, 30, 40].into())
                .set_data(data)
                .build(),
            RecordBuf::builder()
                .set_name("r1")
                .set_reference_sequence_id(1)
                .set_alignment_start(noodles_core::Position::MIN)
                .set_mate_reference_sequence_id(0)
                .set_sequence(b"ACG".to_vec().into())
                .build(),
        ];

        for record in &records {
            let mut dst = Vec::new();
            assert!(format_record(&names, &raw(&header, record), &mut dst));
            assert_eq!(dst, generic(&header, record), "{record:?}");
        }
    }

    #[test]
    fn test_format_record_declines_invalid_strings() {
        let header = header();
        let names = Names::new(&header);

        let record = RecordBuf::builder()
            .set_name("r0")
            .set_data(
                [(Tag::from(*b"XZ"), Value::from("tabXhere"))]
                    .into_iter()
                    .collect(),
            )
            .build();

        // The BAM encoder rejects a tab in a string, so it is written into the raw bytes. The
        // record ends with the field's `X`, `here`, and NUL.
        let mut src = raw(&header, &record);
        assert!(src.ends_with(b"Xhere\0"));
        let x = src.len() - 6;
        src[x] = b'\t';

        let mut dst = b"kept".to_vec();
        assert!(!format_record(&names, &src, &mut dst));
        assert_eq!(dst, b"kept");
    }
}
