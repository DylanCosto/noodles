use std::{io, sync::Arc};

use noodles_core::Position;
use noodles_fasta as fasta;
use noodles_sam::{
    self as sam,
    alignment::{
        record::data::field::{Tag, Value},
        record_buf::{QualityScores, Sequence, data::field::Value as ValueBuf},
    },
};

use super::{Feature, Record, ReferenceBases};
use crate::{
    container::compression_header::preservation_map::substitution_matrix::Base,
    record::{Flags, MateFlags},
};

impl Record {
    /// Converts an alignment record, comparing its bases with the sequence that
    /// `get_reference_sequence` returns for its reference sequence name. When that is `None`, the
    /// bases are stored verbatim.
    pub fn try_from_alignment_record<F>(
        header: &sam::Header,
        record: &dyn sam::alignment::Record,
        get_reference_sequence: F,
    ) -> io::Result<Self>
    where
        F: FnOnce(&[u8]) -> io::Result<Option<Arc<fasta::record::Sequence>>>,
    {
        let bam_flags = record.flags()?;
        let mut cram_flags = Flags::QUALITY_SCORES_ARE_STORED_AS_ARRAY;

        let sequence = if record.sequence().is_empty() {
            cram_flags.insert(Flags::SEQUENCE_IS_MISSING);
            Sequence::default()
        } else {
            Sequence::from(record.sequence().iter().collect::<Vec<_>>())
        };

        let reference_sequence_id = record.reference_sequence_id(header).transpose()?;
        let alignment_start = record.alignment_start().transpose()?;

        // As in HTSlib, a mapped record without a sequence (`*`) has the read length of its
        // CIGAR, so decoders can rebuild the CIGAR from the features.
        let read_length = if sequence.is_empty() && !bam_flags.is_unmapped() {
            record.cigar().read_length()?
        } else {
            sequence.len()
        };

        // Quality scores are stored as an array of the read length; missing ones (`*`) are
        // stored as 0xff, as in HTSlib.
        let quality_scores = if record.quality_scores().is_empty() {
            QualityScores::from(vec![0xff; read_length])
        } else {
            record.quality_scores().iter().collect::<io::Result<_>>()?
        };

        let (features, reference_bases) =
            if let (Some(id), Some(start)) = (reference_sequence_id, alignment_start) {
                let (reference_sequence_name, _) = header
                    .reference_sequences()
                    .get_index(id)
                    .expect("missing reference sequence ID");

                let reference_sequence = get_reference_sequence(reference_sequence_name)?;

                let reference_bases = if reference_sequence.is_some() {
                    ReferenceBases::External
                } else {
                    ReferenceBases::Verbatim
                };

                let features = cigar_to_features(
                    record.cigar().as_ref(),
                    reference_sequence.as_deref(),
                    cram_flags,
                    start,
                    &sequence,
                    &quality_scores,
                )?;

                (features, reference_bases)
            } else {
                (Vec::new(), ReferenceBases::Unused)
            };

        let data = record.data();
        let (data_buf, read_group_id) = get_filtered_data(header, data.as_ref())?;

        Ok(Self {
            bam_flags,
            cram_flags,
            reference_sequence_id,
            read_length,
            alignment_start,
            read_group_id,
            name: record.name().map(|s| s.into()),
            mate_flags: MateFlags::default(),
            mate_reference_sequence_id: record.mate_reference_sequence_id(header).transpose()?,
            mate_alignment_start: record.mate_alignment_start().transpose()?,
            template_length: record.template_length()?,
            mate_distance: None,
            data: data_buf,
            features,
            reference_bases,
            mapping_quality: record.mapping_quality().transpose()?,
            sequence: sequence.into(),
            quality_scores: quality_scores.into(),
        })
    }
}

/// Builds the features of a record's CIGAR. Without a reference sequence, aligned bases are stored
/// verbatim, as stretches of bases.
fn cigar_to_features(
    cigar: &dyn sam::alignment::record::Cigar,
    reference_sequence: Option<&fasta::record::Sequence>,
    flags: Flags,
    alignment_start: Position,
    sequence: &Sequence,
    quality_scores: &QualityScores,
) -> io::Result<Vec<Feature>> {
    use noodles_sam::alignment::record::cigar::op::Kind;

    let reference_bases: Option<&[u8]> = reference_sequence.map(|sequence| sequence.as_ref());
    let bases: &[u8] = sequence.as_ref();
    let scores: &[u8] = quality_scores.as_ref();
    let stores_scores = !flags.quality_scores_are_stored_as_array();

    let mut features = Vec::new();

    // 0-based offsets into the reference sequence and the read.
    let mut reference_offset = usize::from(alignment_start) - 1;
    let mut read_offset: usize = 0;

    for result in cigar.iter() {
        let op = result?;
        let len = op.len();
        let kind = op.kind();

        let read_end = read_offset
            .checked_add(len)
            .ok_or_else(|| io::Error::new(io::ErrorKind::InvalidInput, "CIGAR is too long"))?;

        // A missing sequence (`*`) has no bases to compare or store; HTSlib writes `N` for its
        // inserted and soft-clipped bases.
        if kind.consumes_read() && !bases.is_empty() && read_end > bases.len() {
            return Err(io::Error::new(
                io::ErrorKind::InvalidInput,
                "CIGAR and query sequence are of different length",
            ));
        }

        let position = read_position(read_offset)?;

        match kind {
            Kind::Match | Kind::SequenceMatch | Kind::SequenceMismatch => {
                if !bases.is_empty() {
                    if stores_scores {
                        push_scores(&mut features, position, scores, read_offset..read_end);
                    }

                    if let Some(reference_bases) = reference_bases {
                        push_mismatches(
                            &mut features,
                            reference_bases,
                            reference_offset,
                            &bases[read_offset..read_end],
                            scores,
                            read_offset,
                        )?;
                    } else if len > 0 {
                        features.push(Feature::Bases {
                            position,
                            bases: bases[read_offset..read_end].to_vec(),
                        });
                    }
                }
            }
            Kind::Insertion => {
                let inserted = read_bases(bases, read_offset..read_end);

                if let [base] = inserted[..] {
                    features.push(Feature::InsertBase { position, base });
                } else {
                    features.push(Feature::Insertion {
                        position,
                        bases: inserted,
                    });
                }

                if stores_scores {
                    push_scores(&mut features, position, scores, read_offset..read_end);
                }
            }
            Kind::Deletion => features.push(Feature::Deletion { position, len }),
            Kind::Skip => features.push(Feature::ReferenceSkip { position, len }),
            Kind::SoftClip => {
                features.push(Feature::SoftClip {
                    position,
                    bases: read_bases(bases, read_offset..read_end),
                });

                if stores_scores {
                    push_scores(&mut features, position, scores, read_offset..read_end);
                }
            }
            Kind::HardClip => features.push(Feature::HardClip { position, len }),
            Kind::Pad => features.push(Feature::Padding { position, len }),
        }

        if kind.consumes_reference() {
            reference_offset = reference_offset
                .checked_add(len)
                .ok_or_else(|| io::Error::new(io::ErrorKind::InvalidInput, "CIGAR is too long"))?;
        }

        if kind.consumes_read() {
            read_offset = read_end;
        }
    }

    Ok(features)
}

fn read_position(offset: usize) -> io::Result<Position> {
    offset
        .checked_add(1)
        .and_then(Position::new)
        .ok_or_else(|| io::Error::new(io::ErrorKind::InvalidInput, "invalid read position"))
}

/// Returns the read bases in `range`, or `N`s when the sequence is missing.
fn read_bases(bases: &[u8], range: std::ops::Range<usize>) -> Vec<u8> {
    if bases.is_empty() {
        vec![b'N'; range.len()]
    } else {
        bases[range].to_vec()
    }
}

/// Pushes a quality score feature for the scores in `range`, if the read has them.
fn push_scores(
    features: &mut Vec<Feature>,
    position: Position,
    scores: &[u8],
    range: std::ops::Range<usize>,
) {
    let Some(scores) = scores.get(range) else {
        return;
    };

    if let [quality_score] = *scores {
        features.push(Feature::QualityScore {
            position,
            quality_score,
        });
    } else {
        features.push(Feature::Scores {
            position,
            quality_scores: scores.to_vec(),
        });
    }
}

/// Pushes substitutions (or verbatim read bases) for the aligned bases `read_bases`, which start
/// at `read_offset` in the read and `reference_offset` in the reference.
///
/// As in HTSlib, bases past the end of the reference are stored verbatim.
pub(crate) fn push_mismatches(
    features: &mut Vec<Feature>,
    reference_bases: &[u8],
    reference_offset: usize,
    read_bases: &[u8],
    scores: &[u8],
    read_offset: usize,
) -> io::Result<()> {
    let reference_bases = reference_bases.get(reference_offset..).unwrap_or_default();
    let n = reference_bases.len().min(read_bases.len());

    // Most bases match: compare chunks first, and bases only in chunks that differ.
    const CHUNK_LEN: usize = 16;

    let chunks = reference_bases[..n]
        .chunks(CHUNK_LEN)
        .zip(read_bases[..n].chunks(CHUNK_LEN));

    for (chunk_index, (reference_chunk, read_chunk)) in chunks.enumerate() {
        if reference_chunk == read_chunk {
            continue;
        }

        for (j, (&reference_base, &read_base)) in reference_chunk.iter().zip(read_chunk).enumerate()
        {
            if !reference_base.eq_ignore_ascii_case(&read_base) {
                let i = chunk_index * CHUNK_LEN + j;
                let position = read_position(read_offset + i)?;
                features.push(mismatch(position, reference_base, read_base, || {
                    scores.get(read_offset + i).copied()
                }));
            }
        }
    }

    // Bases past the end of the reference.
    for (i, &base) in read_bases.iter().enumerate().skip(n) {
        features.push(Feature::ReadBase {
            position: read_position(read_offset + i)?,
            base,
            quality_score: scores.get(read_offset + i).copied().unwrap_or(0xff),
        });
    }

    Ok(())
}

/// A substitution, or the read base verbatim if either base cannot be substituted.
fn mismatch<F>(position: Position, reference_base: u8, read_base: u8, score: F) -> Feature
where
    F: FnOnce() -> Option<u8>,
{
    match (Base::try_from(reference_base), Base::try_from(read_base)) {
        (Ok(reference_base), Ok(read_base)) => Feature::Substitution {
            position,
            reference_base,
            read_base,
        },
        // HTSlib stores missing quality scores (`*`) as 0xff.
        _ => Feature::ReadBase {
            position,
            base: read_base,
            quality_score: score().unwrap_or(0xff),
        },
    }
}

/// Returns the data fields to store as tags, and the read group ID of the `RG` field.
///
/// As in HTSlib, an `RG` field naming a read group that the header lacks (or that is not a
/// string) is stored verbatim as a tag.
#[allow(clippy::type_complexity)]
fn get_filtered_data(
    header: &sam::Header,
    data: &dyn sam::alignment::record::Data<'_>,
) -> io::Result<(Vec<(Tag, ValueBuf)>, Option<usize>)> {
    let mut data_buf = Vec::new();
    let mut read_group_id = None;

    for result in data.iter() {
        let (tag, value) = result?;

        if tag == Tag::READ_GROUP
            && read_group_id.is_none()
            && let Value::String(name) = &value
            && let Some(id) = header.read_groups().get_index_of(*name)
        {
            read_group_id = Some(id);
            continue;
        }

        data_buf.push((tag, value.try_into()?));
    }

    Ok((data_buf, read_group_id))
}

#[cfg(test)]
mod tests {
    use noodles_sam::alignment::{
        record::cigar::{Op, op::Kind},
        record_buf::Cigar,
    };

    use super::*;

    fn mapped_record(
        cigar: &str,
        sequence: &[u8],
        quality_scores: &[u8],
    ) -> sam::alignment::RecordBuf {
        use sam::alignment::record_buf::Cigar;

        let cigar: Cigar = cigar
            .split_inclusive(|c: char| c.is_ascii_alphabetic() || c == '=')
            .map(|s| {
                let (n, op) = s.split_at(s.len() - 1);
                let kind = match op {
                    "M" => Kind::Match,
                    "I" => Kind::Insertion,
                    "D" => Kind::Deletion,
                    "S" => Kind::SoftClip,
                    _ => unreachable!(),
                };
                Op::new(kind, n.parse().unwrap())
            })
            .collect();

        sam::alignment::RecordBuf::builder()
            .set_flags(sam::alignment::record::Flags::empty())
            .set_reference_sequence_id(0)
            .set_alignment_start(Position::MIN)
            .set_cigar(cigar)
            .set_sequence(Sequence::from(sequence.to_vec()))
            .set_quality_scores(QualityScores::from(quality_scores.to_vec()))
            .build()
    }

    fn convert(record: &sam::alignment::RecordBuf) -> io::Result<Record> {
        let header =
            sam::Header::builder()
                .add_reference_sequence(
                    "sq0",
                    sam::header::record::value::Map::<
                        sam::header::record::value::map::ReferenceSequence,
                    >::new(std::num::NonZero::new(4).unwrap()),
                )
                .build();

        let repository = fasta::Repository::new(vec![fasta::Record::new(
            fasta::record::Definition::new("sq0", None),
            fasta::record::Sequence::from(b"ACGT".to_vec()),
        )]);

        Record::try_from_alignment_record(&header, record, |name| repository.get(name).transpose())
    }

    #[test]
    fn missing_sequence_takes_the_cigar_read_length() -> io::Result<()> {
        // SEQ and QUAL `*`, as HTSlib stores them: the read length is the CIGAR's, qualities
        // are 0xff, and inserted and soft-clipped bases are `N`.
        let record = convert(&mapped_record("1S2M2I", b"", b""))?;

        assert!(record.cram_flags.contains(Flags::SEQUENCE_IS_MISSING));
        assert_eq!(record.read_length, 5);
        assert_eq!(record.quality_scores, [0xff; 5]);
        assert_eq!(
            record.features,
            [
                Feature::SoftClip {
                    position: Position::MIN,
                    bases: b"N".to_vec(),
                },
                Feature::Insertion {
                    position: Position::try_from(4).unwrap(),
                    bases: b"NN".to_vec(),
                },
            ]
        );

        Ok(())
    }

    #[test]
    fn missing_quality_scores_are_stored_as_0xff() -> io::Result<()> {
        let record = convert(&mapped_record("2M", b"AC", b""))?;
        assert_eq!(record.read_length, 2);
        assert_eq!(record.quality_scores, [0xff; 2]);
        assert!(record.features.is_empty());
        Ok(())
    }

    #[test]
    fn bases_past_the_reference_end_are_stored_verbatim() -> io::Result<()> {
        // The reference is ACGT; the read covers positions 3 to 6.
        let mut record = mapped_record("4M", b"GTCA", b"ABCD");
        *record.alignment_start_mut() = Position::new(3);
        let record = convert(&record)?;

        assert_eq!(
            record.features,
            [
                Feature::ReadBase {
                    position: Position::try_from(3).unwrap(),
                    base: b'C',
                    quality_score: b'C',
                },
                Feature::ReadBase {
                    position: Position::try_from(4).unwrap(),
                    base: b'A',
                    quality_score: b'D',
                },
            ]
        );

        Ok(())
    }

    #[test]
    fn cigar_longer_than_the_sequence_is_an_error() {
        assert!(matches!(
            convert(&mapped_record("4M", b"AC", b"")),
            Err(e) if e.kind() == io::ErrorKind::InvalidInput
        ));
    }

    #[test]
    fn read_group_uses_dedicated_series_without_duplicate_auxiliary_tag() -> io::Result<()> {
        let header = sam::Header::builder()
            .add_read_group("rg0", Default::default())
            .build();
        let data = [(Tag::READ_GROUP, "rg0"), (Tag::COMMENT, "note")];
        let data = data
            .into_iter()
            .map(|(tag, s)| (tag, ValueBuf::from(s)))
            .collect();
        let record = sam::alignment::RecordBuf::builder().set_data(data).build();
        let actual = Record::try_from_alignment_record(&header, &record, |_| Ok(None))?;
        assert_eq!(actual.read_group_id, Some(0));
        assert_eq!(actual.data, [(Tag::COMMENT, ValueBuf::from("note"))]);
        Ok(())
    }

    #[test]
    fn read_group_missing_from_the_header_is_stored_as_a_tag() -> io::Result<()> {
        let header = sam::Header::builder()
            .add_read_group("rg0", Default::default())
            .build();
        let data = [(Tag::READ_GROUP, "nosuch"), (Tag::COMMENT, "note")];
        let data: sam::alignment::record_buf::Data = data
            .into_iter()
            .map(|(tag, s)| (tag, ValueBuf::from(s)))
            .collect();
        let record = sam::alignment::RecordBuf::builder()
            .set_data(data.clone())
            .build();
        let actual = Record::try_from_alignment_record(&header, &record, |_| Ok(None))?;
        assert_eq!(actual.read_group_id, None);
        assert_eq!(
            actual.data,
            data.iter()
                .map(|(tag, value)| (tag, value.clone()))
                .collect::<Vec<_>>()
        );
        Ok(())
    }

    #[test]
    fn test_cigar_to_features() -> Result<(), Box<dyn std::error::Error>> {
        let reference_sequence = fasta::record::Sequence::from(b"ACGT".to_vec());

        let flags = Flags::default();
        let alignment_start = Position::MIN;

        let cigar: Cigar = [Op::new(Kind::Match, 1)].into_iter().collect();
        let sequence = Sequence::from(b"A");
        let quality_scores = [45].into_iter().collect();
        let actual = cigar_to_features(
            &cigar,
            Some(&reference_sequence),
            flags,
            alignment_start,
            &sequence,
            &quality_scores,
        )?;
        let expected = vec![Feature::QualityScore {
            position: Position::try_from(1)?,
            quality_score: 45,
        }];
        assert_eq!(actual, expected);

        let cigar: Cigar = [Op::new(Kind::Match, 1)].into_iter().collect();
        let sequence = Sequence::from(b"C");
        let quality_scores = [45].into_iter().collect();
        let actual = cigar_to_features(
            &cigar,
            Some(&reference_sequence),
            flags,
            alignment_start,
            &sequence,
            &quality_scores,
        )?;
        let expected = vec![
            Feature::QualityScore {
                position: Position::try_from(1)?,
                quality_score: 45,
            },
            Feature::Substitution {
                position: Position::try_from(1)?,
                reference_base: Base::A,
                read_base: Base::C,
            },
        ];
        assert_eq!(actual, expected);

        let cigar: Cigar = [Op::new(Kind::Match, 2)].into_iter().collect();
        let sequence = Sequence::from(b"AC");
        let quality_scores = [45, 35].into_iter().collect();
        let actual = cigar_to_features(
            &cigar,
            Some(&reference_sequence),
            flags,
            alignment_start,
            &sequence,
            &quality_scores,
        )?;
        let expected = vec![Feature::Scores {
            position: Position::try_from(1)?,
            quality_scores: vec![45, 35],
        }];
        assert_eq!(actual, expected);

        let cigar: Cigar = [Op::new(Kind::Match, 2)].into_iter().collect();
        let sequence = Sequence::from(b"AT");
        let quality_scores = [45, 35].into_iter().collect();
        let actual = cigar_to_features(
            &cigar,
            Some(&reference_sequence),
            flags,
            alignment_start,
            &sequence,
            &quality_scores,
        )?;
        let expected = vec![
            Feature::Scores {
                position: Position::try_from(1)?,
                quality_scores: vec![45, 35],
            },
            Feature::Substitution {
                position: Position::try_from(2)?,
                reference_base: Base::C,
                read_base: Base::T,
            },
        ];
        assert_eq!(actual, expected);

        let cigar: Cigar = [Op::new(Kind::Insertion, 1), Op::new(Kind::Match, 1)]
            .into_iter()
            .collect();
        let sequence = Sequence::from(b"AC");
        let quality_scores = [45, 35].into_iter().collect();
        let actual = cigar_to_features(
            &cigar,
            Some(&reference_sequence),
            flags,
            alignment_start,
            &sequence,
            &quality_scores,
        )?;
        let expected = vec![
            Feature::InsertBase {
                position: Position::try_from(1)?,
                base: b'A',
            },
            Feature::QualityScore {
                position: Position::try_from(1)?,
                quality_score: 45,
            },
            Feature::QualityScore {
                position: Position::try_from(2)?,
                quality_score: 35,
            },
            Feature::Substitution {
                position: Position::try_from(2)?,
                reference_base: Base::A,
                read_base: Base::C,
            },
        ];
        assert_eq!(actual, expected);

        let cigar: Cigar = [Op::new(Kind::Insertion, 2), Op::new(Kind::Match, 1)]
            .into_iter()
            .collect();
        let sequence = Sequence::from(b"ACG");
        let quality_scores = [45, 35, 43].into_iter().collect();
        let actual = cigar_to_features(
            &cigar,
            Some(&reference_sequence),
            flags,
            alignment_start,
            &sequence,
            &quality_scores,
        )?;
        let expected = vec![
            Feature::Insertion {
                position: Position::try_from(1)?,
                bases: vec![b'A', b'C'],
            },
            Feature::Scores {
                position: Position::try_from(1)?,
                quality_scores: vec![45, 35],
            },
            Feature::QualityScore {
                position: Position::try_from(3)?,
                quality_score: 43,
            },
            Feature::Substitution {
                position: Position::try_from(3)?,
                reference_base: Base::A,
                read_base: Base::G,
            },
        ];
        assert_eq!(actual, expected);

        let cigar: Cigar = [Op::new(Kind::Deletion, 1), Op::new(Kind::Match, 2)]
            .into_iter()
            .collect();
        let sequence = Sequence::from(b"AC");
        let quality_scores = [45, 35].into_iter().collect();
        let actual = cigar_to_features(
            &cigar,
            Some(&reference_sequence),
            flags,
            alignment_start,
            &sequence,
            &quality_scores,
        )?;
        let expected = vec![
            Feature::Deletion {
                position: Position::try_from(1)?,
                len: 1,
            },
            Feature::Scores {
                position: Position::try_from(1)?,
                quality_scores: vec![45, 35],
            },
            Feature::Substitution {
                position: Position::try_from(1)?,
                reference_base: Base::C,
                read_base: Base::A,
            },
            Feature::Substitution {
                position: Position::try_from(2)?,
                reference_base: Base::G,
                read_base: Base::C,
            },
        ];
        assert_eq!(actual, expected);

        let cigar: Cigar = [Op::new(Kind::Skip, 1), Op::new(Kind::Match, 1)]
            .into_iter()
            .collect();
        let sequence = Sequence::from(b"A");
        let quality_scores = [45].into_iter().collect();
        let actual = cigar_to_features(
            &cigar,
            Some(&reference_sequence),
            flags,
            alignment_start,
            &sequence,
            &quality_scores,
        )?;
        let expected = vec![
            Feature::ReferenceSkip {
                position: Position::try_from(1)?,
                len: 1,
            },
            Feature::QualityScore {
                position: Position::try_from(1)?,
                quality_score: 45,
            },
            Feature::Substitution {
                position: Position::try_from(1)?,
                reference_base: Base::C,
                read_base: Base::A,
            },
        ];
        assert_eq!(actual, expected);

        let cigar: Cigar = [Op::new(Kind::SoftClip, 1), Op::new(Kind::Match, 1)]
            .into_iter()
            .collect();
        let sequence = Sequence::from(b"AC");
        let quality_scores = [45, 35].into_iter().collect();
        let actual = cigar_to_features(
            &cigar,
            Some(&reference_sequence),
            flags,
            alignment_start,
            &sequence,
            &quality_scores,
        )?;
        let expected = vec![
            Feature::SoftClip {
                position: Position::try_from(1)?,
                bases: vec![b'A'],
            },
            Feature::QualityScore {
                position: Position::try_from(1)?,
                quality_score: 45,
            },
            Feature::QualityScore {
                position: Position::try_from(2)?,
                quality_score: 35,
            },
            Feature::Substitution {
                position: Position::try_from(2)?,
                reference_base: Base::A,
                read_base: Base::C,
            },
        ];
        assert_eq!(actual, expected);

        let cigar: Cigar = [Op::new(Kind::SoftClip, 2), Op::new(Kind::Match, 1)]
            .into_iter()
            .collect();
        let sequence = Sequence::from(b"ACG");
        let quality_scores = [45, 35, 43].into_iter().collect();
        let actual = cigar_to_features(
            &cigar,
            Some(&reference_sequence),
            flags,
            alignment_start,
            &sequence,
            &quality_scores,
        )?;
        let expected = vec![
            Feature::SoftClip {
                position: Position::try_from(1)?,
                bases: vec![b'A', b'C'],
            },
            Feature::Scores {
                position: Position::try_from(1)?,
                quality_scores: vec![45, 35],
            },
            Feature::QualityScore {
                position: Position::try_from(3)?,
                quality_score: 43,
            },
            Feature::Substitution {
                position: Position::try_from(3)?,
                reference_base: Base::A,
                read_base: Base::G,
            },
        ];
        assert_eq!(actual, expected);

        let cigar: Cigar = [Op::new(Kind::HardClip, 1), Op::new(Kind::Match, 1)]
            .into_iter()
            .collect();
        let sequence = Sequence::from(b"A");
        let quality_scores = [45].into_iter().collect();
        let actual = cigar_to_features(
            &cigar,
            Some(&reference_sequence),
            flags,
            alignment_start,
            &sequence,
            &quality_scores,
        )?;
        let expected = vec![
            Feature::HardClip {
                position: Position::try_from(1)?,
                len: 1,
            },
            Feature::QualityScore {
                position: Position::try_from(1)?,
                quality_score: 45,
            },
        ];
        assert_eq!(actual, expected);

        let cigar: Cigar = [Op::new(Kind::Pad, 1), Op::new(Kind::Match, 1)]
            .into_iter()
            .collect();
        let sequence = Sequence::from(b"A");
        let quality_scores = [45].into_iter().collect();
        let actual = cigar_to_features(
            &cigar,
            Some(&reference_sequence),
            flags,
            alignment_start,
            &sequence,
            &quality_scores,
        )?;
        let expected = vec![
            Feature::Padding {
                position: Position::try_from(1)?,
                len: 1,
            },
            Feature::QualityScore {
                position: Position::try_from(1)?,
                quality_score: 45,
            },
        ];
        assert_eq!(actual, expected);

        Ok(())
    }

    #[test]
    fn test_cigar_to_features_with_quality_scores_stored_as_array()
    -> Result<(), Box<dyn std::error::Error>> {
        let reference_sequence = fasta::record::Sequence::from(b"ACGT".to_vec());

        let flags = Flags::QUALITY_SCORES_ARE_STORED_AS_ARRAY;
        let alignment_start = Position::MIN;

        let cigar: Cigar = [Op::new(Kind::Match, 1)].into_iter().collect();
        let sequence = Sequence::from(b"A");
        let quality_scores = [45].into_iter().collect();
        let actual = cigar_to_features(
            &cigar,
            Some(&reference_sequence),
            flags,
            alignment_start,
            &sequence,
            &quality_scores,
        )?;
        let expected = vec![];
        assert_eq!(actual, expected);

        let cigar: Cigar = [Op::new(Kind::Match, 1)].into_iter().collect();
        let sequence = Sequence::from(b"C");
        let quality_scores = [45].into_iter().collect();
        let actual = cigar_to_features(
            &cigar,
            Some(&reference_sequence),
            flags,
            alignment_start,
            &sequence,
            &quality_scores,
        )?;
        let expected = vec![Feature::Substitution {
            position: Position::try_from(1)?,
            reference_base: Base::A,
            read_base: Base::C,
        }];
        assert_eq!(actual, expected);

        let cigar: Cigar = [Op::new(Kind::Match, 2)].into_iter().collect();
        let sequence = Sequence::from(b"AC");
        let quality_scores = [45, 35].into_iter().collect();
        let actual = cigar_to_features(
            &cigar,
            Some(&reference_sequence),
            flags,
            alignment_start,
            &sequence,
            &quality_scores,
        )?;
        let expected = vec![];
        assert_eq!(actual, expected);

        let cigar: Cigar = [Op::new(Kind::Match, 2)].into_iter().collect();
        let sequence = Sequence::from(b"AT");
        let quality_scores = [45, 35].into_iter().collect();
        let actual = cigar_to_features(
            &cigar,
            Some(&reference_sequence),
            flags,
            alignment_start,
            &sequence,
            &quality_scores,
        )?;
        let expected = vec![Feature::Substitution {
            position: Position::try_from(2)?,
            reference_base: Base::C,
            read_base: Base::T,
        }];
        assert_eq!(actual, expected);

        let cigar: Cigar = [Op::new(Kind::Insertion, 1), Op::new(Kind::Match, 1)]
            .into_iter()
            .collect();
        let sequence = Sequence::from(b"AC");
        let quality_scores = [45, 35].into_iter().collect();
        let actual = cigar_to_features(
            &cigar,
            Some(&reference_sequence),
            flags,
            alignment_start,
            &sequence,
            &quality_scores,
        )?;
        let expected = vec![
            Feature::InsertBase {
                position: Position::try_from(1)?,
                base: b'A',
            },
            Feature::Substitution {
                position: Position::try_from(2)?,
                reference_base: Base::A,
                read_base: Base::C,
            },
        ];
        assert_eq!(actual, expected);

        let cigar: Cigar = [Op::new(Kind::Insertion, 2), Op::new(Kind::Match, 1)]
            .into_iter()
            .collect();
        let sequence = Sequence::from(b"ACG");
        let quality_scores = [45, 35, 43].into_iter().collect();
        let actual = cigar_to_features(
            &cigar,
            Some(&reference_sequence),
            flags,
            alignment_start,
            &sequence,
            &quality_scores,
        )?;
        let expected = vec![
            Feature::Insertion {
                position: Position::try_from(1)?,
                bases: vec![b'A', b'C'],
            },
            Feature::Substitution {
                position: Position::try_from(3)?,
                reference_base: Base::A,
                read_base: Base::G,
            },
        ];
        assert_eq!(actual, expected);

        let cigar: Cigar = [Op::new(Kind::Deletion, 1), Op::new(Kind::Match, 2)]
            .into_iter()
            .collect();
        let sequence = Sequence::from(b"AC");
        let quality_scores = [45, 35].into_iter().collect();
        let actual = cigar_to_features(
            &cigar,
            Some(&reference_sequence),
            flags,
            alignment_start,
            &sequence,
            &quality_scores,
        )?;
        let expected = vec![
            Feature::Deletion {
                position: Position::try_from(1)?,
                len: 1,
            },
            Feature::Substitution {
                position: Position::try_from(1)?,
                reference_base: Base::C,
                read_base: Base::A,
            },
            Feature::Substitution {
                position: Position::try_from(2)?,
                reference_base: Base::G,
                read_base: Base::C,
            },
        ];
        assert_eq!(actual, expected);

        let cigar: Cigar = [Op::new(Kind::Skip, 1), Op::new(Kind::Match, 1)]
            .into_iter()
            .collect();
        let sequence = Sequence::from(b"A");
        let quality_scores = [45].into_iter().collect();
        let actual = cigar_to_features(
            &cigar,
            Some(&reference_sequence),
            flags,
            alignment_start,
            &sequence,
            &quality_scores,
        )?;
        let expected = vec![
            Feature::ReferenceSkip {
                position: Position::try_from(1)?,
                len: 1,
            },
            Feature::Substitution {
                position: Position::try_from(1)?,
                reference_base: Base::C,
                read_base: Base::A,
            },
        ];
        assert_eq!(actual, expected);

        let cigar: Cigar = [Op::new(Kind::SoftClip, 1), Op::new(Kind::Match, 1)]
            .into_iter()
            .collect();
        let sequence = Sequence::from(b"AC");
        let quality_scores = [45, 35].into_iter().collect();
        let actual = cigar_to_features(
            &cigar,
            Some(&reference_sequence),
            flags,
            alignment_start,
            &sequence,
            &quality_scores,
        )?;
        let expected = vec![
            Feature::SoftClip {
                position: Position::try_from(1)?,
                bases: vec![b'A'],
            },
            Feature::Substitution {
                position: Position::try_from(2)?,
                reference_base: Base::A,
                read_base: Base::C,
            },
        ];
        assert_eq!(actual, expected);

        let cigar: Cigar = [Op::new(Kind::SoftClip, 2), Op::new(Kind::Match, 1)]
            .into_iter()
            .collect();
        let sequence = Sequence::from(b"ACG");
        let quality_scores = [45, 35, 43].into_iter().collect();
        let actual = cigar_to_features(
            &cigar,
            Some(&reference_sequence),
            flags,
            alignment_start,
            &sequence,
            &quality_scores,
        )?;
        let expected = vec![
            Feature::SoftClip {
                position: Position::try_from(1)?,
                bases: vec![b'A', b'C'],
            },
            Feature::Substitution {
                position: Position::try_from(3)?,
                reference_base: Base::A,
                read_base: Base::G,
            },
        ];
        assert_eq!(actual, expected);

        let cigar: Cigar = [Op::new(Kind::HardClip, 1), Op::new(Kind::Match, 1)]
            .into_iter()
            .collect();
        let sequence = Sequence::from(b"A");
        let quality_scores = [45].into_iter().collect();
        let actual = cigar_to_features(
            &cigar,
            Some(&reference_sequence),
            flags,
            alignment_start,
            &sequence,
            &quality_scores,
        )?;
        let expected = vec![Feature::HardClip {
            position: Position::try_from(1)?,
            len: 1,
        }];
        assert_eq!(actual, expected);

        let cigar: Cigar = [Op::new(Kind::Pad, 1), Op::new(Kind::Match, 1)]
            .into_iter()
            .collect();
        let sequence = Sequence::from(b"A");
        let quality_scores = [45].into_iter().collect();
        let actual = cigar_to_features(
            &cigar,
            Some(&reference_sequence),
            flags,
            alignment_start,
            &sequence,
            &quality_scores,
        )?;
        let expected = vec![Feature::Padding {
            position: Position::try_from(1)?,
            len: 1,
        }];
        assert_eq!(actual, expected);

        Ok(())
    }
}
