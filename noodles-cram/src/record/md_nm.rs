//! MD and NM tags regenerated from read features, as HTSlib does when decoding.
//!
//! HTSlib's CRAM writer drops MD and NM tags that match the alignment (unless asked to keep
//! them) and its reader rebuilds them for mapped records with a reference, appending MD and then
//! NM after the stored tags. This ports `cram_decode_seq`'s MD/NM generation (HTSlib 1.24), so
//! records decode with the same tags. A `cF` tag (CRAM flags) is removed; its bits 1 and 2 stop
//! MD and NM from being generated.

use std::{borrow::Cow, io::Write};

use noodles_sam::alignment::record::data::field::Tag;

use super::{Feature, Record, data::field::Value};
use crate::io::reader::container::slice::ReferenceSequence;

const CRAM_FLAGS: Tag = Tag::new(b'c', b'F');

/// Removes a `cF` tag and adds the MD and NM tags HTSlib generates for this record.
pub(crate) fn add_missing_md_nm(record: &mut Record<'_>) {
    let mut cram_flags = 0;

    if let Some(i) = record.data.iter().position(|(tag, _)| *tag == CRAM_FLAGS)
        && let Value::UInt8(flags) = record.data[i].1
    {
        cram_flags = flags;
        record.data.remove(i);
    }

    if record.bam_flags.is_unmapped() || record.cram_flags.sequence_is_missing() {
        return;
    }

    let has_md = cram_flags & 1 != 0
        || record
            .data
            .iter()
            .any(|(tag, _)| *tag == Tag::MISMATCHED_POSITIONS);
    let has_nm = cram_flags & 2 != 0
        || record
            .data
            .iter()
            .any(|(tag, _)| *tag == Tag::EDIT_DISTANCE);

    if has_md && has_nm {
        return;
    }

    let Some((md, nm)) = calculate(record) else {
        return;
    };

    if !has_md {
        record.data.push((
            Tag::MISMATCHED_POSITIONS,
            Value::String(Cow::Owned(md.into())),
        ));
    }

    if !has_nm {
        let value = if let Ok(n) = u8::try_from(nm) {
            Value::UInt8(n)
        } else if let Ok(n) = u16::try_from(nm) {
            Value::UInt16(n)
        } else {
            Value::UInt32(nm)
        };

        record.data.push((Tag::EDIT_DISTANCE, value));
    }
}

/// A record's reference bases: the loaded bases, the 0-based reference position of the first,
/// and the reference sequence's length.
struct Reference<'a> {
    bases: &'a [u8],
    offset: usize,
    len: usize,
}

impl Reference<'_> {
    /// The uppercased base at a 0-based reference position, or `None` outside the loaded bases.
    fn get(&self, pos: usize) -> Option<u8> {
        pos.checked_sub(self.offset)
            .and_then(|i| self.bases.get(i))
            .map(|b| b & !0x20)
    }
}

fn reference<'a>(record: &'a Record<'_>) -> Option<Reference<'a>> {
    let reference_sequence_id = record.reference_sequence_id?;
    let len = record
        .header?
        .reference_sequences()
        .get_index(reference_sequence_id)
        .map(|(_, rs)| usize::from(rs.length()))?;

    let (bases, offset) = match record.reference_sequence.as_ref()? {
        ReferenceSequence::Embedded {
            reference_start,
            sequence,
        } => (*sequence, usize::from(*reference_start) - 1),
        ReferenceSequence::ExternalRange {
            reference_start,
            sequence,
        } => ((**sequence).as_ref(), usize::from(*reference_start) - 1),
        ReferenceSequence::External { sequence } => ((**sequence).as_ref(), 0),
    };

    Some(Reference { bases, offset, len })
}

/// The MD string and NM count, or `None` if the record has no reference or reads outside the
/// loaded reference bases (where HTSlib fails to decode the slice).
fn calculate(record: &Record<'_>) -> Option<(Vec<u8>, u32)> {
    let reference = reference(record)?;
    let read_length = record.read_length;

    let mut md = Vec::new();
    let mut md_dist: i64 = 0;
    let mut nm: u32 = 0;
    let mut seq_pos = 1;
    let mut ref_pos = usize::from(record.alignment_start?) - 1;

    // `add_md_char`
    let push_mismatch = |md: &mut Vec<u8>, md_dist: &mut i64, base: u8| {
        write!(md, "{md_dist}").expect("writing to a Vec cannot fail");
        md.push(base);
        *md_dist = 0;
    };

    // Matching bases, where reference Ns are mismatches.
    let matches = |md: &mut Vec<u8>, md_dist: &mut i64, nm: &mut u32, start: usize, len: usize| {
        if start + len > reference.len {
            if *md_dist >= 0 {
                *md_dist += len as i64;
            }

            return Some(());
        }

        for i in 0..len {
            if reference.get(start + i)? == b'N' {
                push_mismatch(md, md_dist, b'N');
                *nm += 1;
            } else {
                *md_dist += 1;
            }
        }

        Some(())
    };

    let push_dist = |md: &mut Vec<u8>, md_dist: i64| {
        if md_dist >= 0 {
            write!(md, "{md_dist}").expect("writing to a Vec cannot fail");
        }
    };

    for feature in &record.features {
        let pos = usize::from(feature.position());

        if pos > seq_pos {
            let len = pos - seq_pos;
            matches(&mut md, &mut md_dist, &mut nm, ref_pos, len)?;
            ref_pos += len;
            seq_pos = pos;
        }

        match feature {
            Feature::SoftClip { bases, .. } => seq_pos += bases.len(),
            Feature::Substitution { .. } => {
                if ref_pos >= reference.len {
                    push_dist(&mut md, md_dist);
                    md_dist = -1;
                    nm = nm.wrapping_sub(1);
                } else {
                    let base = reference.get(ref_pos)?;
                    push_mismatch(&mut md, &mut md_dist, base);
                }

                nm = nm.wrapping_add(1);
                seq_pos += 1;
                ref_pos += 1;
            }
            Feature::Deletion { len, .. } => {
                let len = *len;
                push_dist(&mut md, md_dist);

                if ref_pos + len <= reference.len {
                    md.push(b'^');

                    for i in 0..len {
                        md.push(reference.get(ref_pos + i)?);
                    }

                    md_dist = 0;
                    nm += len as u32;
                } else {
                    if reference.len >= ref_pos {
                        md.push(b'^');

                        for i in ref_pos..reference.len {
                            md.push(reference.get(i)?);
                        }

                        md.push(b'0');
                        nm += (reference.len - ref_pos) as u32;
                    }

                    md_dist = -1;
                }

                ref_pos += len;
            }
            Feature::Insertion { bases, .. } => {
                seq_pos += bases.len();
                nm += bases.len() as u32;
            }
            Feature::InsertBase { .. } => {
                seq_pos += 1;
                nm += 1;
            }
            Feature::Bases { bases, .. } => {
                let len = bases.len();
                push_dist(&mut md, md_dist);

                let mut x = 0;

                while x < len {
                    if x > 0 {
                        md.push(b'0');
                    }

                    if ref_pos + x >= reference.len {
                        break;
                    }

                    md.push(reference.get(ref_pos + x)?);
                    x += 1;
                }

                nm += x as u32;
                md_dist = 0;
                seq_pos += len;
                ref_pos += len;
            }
            Feature::ReadBase { .. } => {
                push_dist(&mut md, md_dist);

                if ref_pos >= reference.len {
                    md_dist = -1;
                } else {
                    md.push(reference.get(ref_pos)?);
                    nm += 1;
                    md_dist = 0;
                }

                seq_pos += 1;
                ref_pos += 1;
            }
            Feature::ReferenceSkip { len, .. } => ref_pos += len,
            Feature::Scores { .. }
            | Feature::QualityScore { .. }
            | Feature::Padding { .. }
            | Feature::HardClip { .. } => {}
        }
    }

    if read_length >= seq_pos {
        let len = read_length - seq_pos + 1;
        matches(&mut md, &mut md_dist, &mut nm, ref_pos, len)?;
    }

    push_dist(&mut md, md_dist);

    Some((md, nm))
}

#[cfg(test)]
mod tests {
    use std::sync::Arc;

    use noodles_core::Position;
    use noodles_fasta as fasta;
    use noodles_sam::{self as sam, header::record::value::map::ReferenceSequence as Sq};

    use bstr::BStr;

    use super::*;

    fn header() -> sam::Header {
        sam::Header::builder()
            .add_reference_sequence(
                "sq0",
                sam::header::record::value::Map::<Sq>::new(
                    std::num::NonZero::new(12).expect("non-zero"),
                ),
            )
            .build()
    }

    fn record<'c>(header: &'c sam::Header, features: Vec<Feature<'c>>) -> Record<'c> {
        Record {
            header: Some(header),
            bam_flags: sam::alignment::record::Flags::empty(),
            reference_sequence_id: Some(0),
            alignment_start: Position::new(2),
            read_length: 8,
            reference_sequence: Some(ReferenceSequence::External {
                sequence: Arc::new(fasta::record::Sequence::from(b"ACGTNacgtACG".to_vec())),
            }),
            features,
            ..Record::default()
        }
    }

    fn md_nm(record: &Record<'_>) -> (String, u32) {
        let (md, nm) = calculate(record).expect("calculated");
        (String::from_utf8(md).expect("ASCII"), nm)
    }

    #[test]
    fn test_calculate() {
        let header = header();
        let position = |n| Position::new(n).expect("non-zero");

        // Reference CGTNacgt: the N is a mismatch.
        assert_eq!(md_nm(&record(&header, Vec::new())), ("3N4".into(), 1));

        let features = vec![
            Feature::Substitution {
                position: position(2),
                code: 0,
            },
            Feature::Deletion {
                position: position(4),
                len: 2,
            },
            Feature::InsertBase {
                position: position(6),
                base: b'A',
            },
        ];

        // C, G>?, T, ^NA, then A (inserted), CGTA...
        assert_eq!(md_nm(&record(&header, features)), ("1G1^NA4".into(), 4));
    }

    #[test]
    fn test_add_missing_md_nm() {
        let header = header();
        let mut r = record(&header, Vec::new());
        r.data.push((CRAM_FLAGS, Value::UInt8(2)));
        add_missing_md_nm(&mut r);

        assert_eq!(r.data.len(), 1);
        assert_eq!(r.data[0].0, Tag::MISMATCHED_POSITIONS);
        assert!(matches!(&r.data[0].1, Value::String(s) if s.as_ref() == BStr::new("3N4")));
    }
}
