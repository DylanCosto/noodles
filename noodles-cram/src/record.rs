//! CRAM record.

mod cigar;
pub(crate) mod data;
pub(crate) mod feature;
mod flags;
mod mate_flags;
mod quality_scores;
mod sequence;

use std::{borrow::Cow, io, iter, marker::PhantomData};

use bstr::{BStr, ByteSlice};
use noodles_core::Position;
use noodles_sam::{
    self as sam,
    alignment::record::{MappingQuality, data::field::Tag},
};

use self::{
    cigar::Cigar,
    data::{Data, field::Value},
    quality_scores::QualityScores,
    sequence::Sequence,
};
pub(crate) use self::{feature::Feature, flags::Flags, mate_flags::MateFlags};
use crate::{
    container::compression_header::preservation_map::SubstitutionMatrix,
    io::reader::container::slice::ReferenceSequence,
};

/// A CRAM record.
#[derive(Clone, Debug, PartialEq)]
pub struct Record<'c> {
    pub(crate) id: u64,
    pub(crate) header: Option<&'c sam::Header>,
    pub(crate) reference_sequence: Option<ReferenceSequence<'c>>,
    pub(crate) substitution_matrix: SubstitutionMatrix,
    pub(crate) bam_flags: sam::alignment::record::Flags,
    pub(crate) cram_flags: Flags,
    pub(crate) reference_sequence_id: Option<usize>,
    pub(crate) read_length: usize,
    pub(crate) alignment_start: Option<Position>,
    pub(crate) read_group_id: Option<usize>,
    pub(crate) name: Option<Cow<'c, [u8]>>,
    pub(crate) mate_flags: MateFlags,
    pub(crate) mate_reference_sequence_id: Option<usize>,
    pub(crate) mate_alignment_start: Option<Position>,
    pub(crate) template_length: i32,
    pub(crate) mate_distance: Option<usize>,
    pub(crate) data: Vec<(Tag, Value<'c>)>,
    pub(crate) sequence: Cow<'c, [u8]>,
    pub(crate) features: Vec<Feature<'c>>,
    pub(crate) mapping_quality: Option<MappingQuality>,
    pub(crate) quality_scores: Cow<'c, [u8]>,
}

impl Record<'_> {
    /// The read bases as reconstructed from the reference sequence and read features, or `None`
    /// if they are stored as is (unmapped records and records with missing sequences).
    fn reconstructed_sequence(&self) -> Option<Sequence<'_, '_>> {
        if self.bam_flags.is_unmapped() || self.cram_flags.sequence_is_missing() {
            return None;
        }

        let (reference_sequence, alignment_start) = match self.reference_sequence.as_ref() {
            Some(ReferenceSequence::Embedded {
                reference_start,
                sequence,
            }) => {
                let alignment_start = usize::from(self.alignment_start.unwrap());
                let offset = usize::from(*reference_start);
                let offset_alignment_start = Position::new(alignment_start - offset + 1).unwrap();
                (Some(*sequence), offset_alignment_start)
            }
            Some(ReferenceSequence::External { sequence, .. }) => {
                (Some((**sequence).as_ref()), self.alignment_start.unwrap())
            }
            None => (None, Position::MIN),
        };

        Some(Sequence::new(
            reference_sequence,
            self.substitution_matrix.clone(),
            &self.features,
            alignment_start,
            self.read_length,
        ))
    }

    /// Converts this record to an alignment record buffer.
    ///
    /// This is equivalent to [`sam::alignment::RecordBuf::try_from_alignment_record`] but copies
    /// the sequence and quality scores in runs instead of base by base.
    pub(crate) fn try_into_alignment_record_buf(
        &self,
        header: &sam::Header,
    ) -> io::Result<sam::alignment::RecordBuf> {
        use sam::alignment::Record as _;

        let mut record_buf = sam::alignment::RecordBuf::default();

        if let Some(name) = self.name() {
            *record_buf.name_mut() = Some(name.into());
        }

        *record_buf.flags_mut() = self.flags()?;
        *record_buf.reference_sequence_id_mut() = self.reference_sequence_id(header).transpose()?;
        *record_buf.alignment_start_mut() = self.alignment_start;
        *record_buf.mapping_quality_mut() = self.mapping_quality;

        let cigar = record_buf.cigar_mut().as_mut();

        for result in self.cigar().iter() {
            cigar.push(result?);
        }

        *record_buf.mate_reference_sequence_id_mut() =
            self.mate_reference_sequence_id(header).transpose()?;
        *record_buf.mate_alignment_start_mut() = self.mate_alignment_start;
        *record_buf.template_length_mut() = self.template_length;

        let sequence = record_buf.sequence_mut().as_mut();

        match self.reconstructed_sequence() {
            Some(reconstructed) => reconstructed.extend_into(sequence),
            None => sequence.extend_from_slice(&self.sequence),
        }

        let quality_scores = record_buf.quality_scores_mut().as_mut();

        if self.bam_flags.is_unmapped() || self.cram_flags.quality_scores_are_stored_as_array() {
            quality_scores.extend_from_slice(&self.quality_scores);
        } else {
            for result in self.quality_scores().iter() {
                quality_scores.push(result?);
            }
        }

        let data = record_buf.data_mut();

        for result in self.data().iter() {
            let (tag, value) = result?;
            data.insert(tag, value.try_into()?);
        }

        Ok(record_buf)
    }

    fn alignment_span(&self) -> usize {
        calculate_alignment_span(self.read_length, &self.features)
    }

    pub(crate) fn alignment_end(&self) -> Option<Position> {
        self.alignment_start.and_then(|alignment_start| {
            let end = usize::from(alignment_start) + self.alignment_span() - 1;
            Position::new(end)
        })
    }
}

impl Default for Record<'_> {
    fn default() -> Self {
        Self {
            id: 0,
            header: None,
            reference_sequence: None,
            substitution_matrix: SubstitutionMatrix::default(),
            bam_flags: sam::alignment::record::Flags::UNMAPPED,
            cram_flags: Flags::default(),
            reference_sequence_id: None,
            read_length: 0,
            alignment_start: None,
            read_group_id: None,
            name: None,
            mate_flags: MateFlags::default(),
            mate_reference_sequence_id: None,
            mate_alignment_start: None,
            template_length: 0,
            mate_distance: None,
            data: Vec::new(),
            sequence: Cow::from(&[]),
            features: Vec::new(),
            mapping_quality: None,
            quality_scores: Cow::from(&[]),
        }
    }
}

impl sam::alignment::Record for Record<'_> {
    fn name(&self) -> Option<&BStr> {
        self.name.as_deref().map(|name| name.as_bstr())
    }

    fn flags(&self) -> io::Result<sam::alignment::record::Flags> {
        Ok(self.bam_flags)
    }

    fn reference_sequence_id<'r, 'h: 'r>(
        &'r self,
        _: &'h sam::Header,
    ) -> Option<io::Result<usize>> {
        self.reference_sequence_id.map(Ok)
    }

    fn alignment_start(&self) -> Option<io::Result<Position>> {
        self.alignment_start.map(Ok)
    }

    fn mapping_quality(&self) -> Option<io::Result<MappingQuality>> {
        self.mapping_quality.map(Ok)
    }

    fn cigar(&self) -> Box<dyn sam::alignment::record::Cigar + '_> {
        Box::new(Cigar::new(
            &self.features,
            self.bam_flags.is_unmapped(),
            self.read_length,
        ))
    }

    fn mate_reference_sequence_id<'r, 'h: 'r>(
        &'r self,
        _: &'h sam::Header,
    ) -> Option<io::Result<usize>> {
        self.mate_reference_sequence_id.map(Ok)
    }

    fn mate_alignment_start(&self) -> Option<io::Result<Position>> {
        self.mate_alignment_start.map(Ok)
    }

    fn template_length(&self) -> io::Result<i32> {
        Ok(self.template_length)
    }

    fn sequence(&self) -> Box<dyn sam::alignment::record::Sequence + '_> {
        match self.reconstructed_sequence() {
            Some(sequence) => Box::new(sequence),
            None => Box::new(Bases(&self.sequence[..])),
        }
    }

    fn quality_scores(&self) -> Box<dyn sam::alignment::record::QualityScores + '_> {
        if self.bam_flags.is_unmapped() || self.cram_flags.quality_scores_are_stored_as_array() {
            Box::new(Scores(&self.quality_scores[..]))
        } else {
            Box::new(QualityScores::new(&self.features, self.read_length))
        }
    }

    fn data(&self) -> Box<dyn sam::alignment::record::Data<'_> + '_> {
        if let Some(header) = self.header {
            Box::new(Data::new(header, &self.data, self.read_group_id))
        } else {
            Box::new(EmptyData::new())
        }
    }

    fn alignment_span(&self) -> Option<io::Result<usize>> {
        Some(Ok(self.alignment_span()))
    }

    fn alignment_end(&self) -> Option<io::Result<Position>> {
        self.alignment_end().map(Ok)
    }
}

struct Bases<'c>(&'c [u8]);

impl sam::alignment::record::Sequence for Bases<'_> {
    fn is_empty(&self) -> bool {
        self.0.is_empty()
    }

    fn len(&self) -> usize {
        self.0.len()
    }

    fn get(&self, i: usize) -> Option<u8> {
        self.0.get(i).copied()
    }

    fn iter(&self) -> Box<dyn Iterator<Item = u8> + '_> {
        Box::new(self.0.iter().copied())
    }
}

struct Scores<'c>(&'c [u8]);

impl sam::alignment::record::QualityScores for Scores<'_> {
    fn is_empty(&self) -> bool {
        self.0.is_empty()
    }

    fn len(&self) -> usize {
        self.0.len()
    }

    fn iter(&self) -> Box<dyn Iterator<Item = io::Result<u8>> + '_> {
        Box::new(self.0.iter().copied().map(Ok))
    }
}

struct EmptyData<'c>(PhantomData<&'c ()>);

impl EmptyData<'_> {
    fn new() -> Self {
        Self(PhantomData)
    }
}

impl<'c> sam::alignment::record::Data<'c> for EmptyData<'c> {
    fn is_empty(&self) -> bool {
        true
    }

    fn get(
        &self,
        _tag: &Tag,
    ) -> Option<io::Result<sam::alignment::record::data::field::Value<'c>>> {
        None
    }

    fn iter(
        &self,
    ) -> Box<
        dyn Iterator<Item = io::Result<(Tag, sam::alignment::record::data::field::Value<'c>)>> + 'c,
    > {
        Box::new(iter::empty())
    }
}

pub(crate) fn calculate_alignment_span(read_length: usize, features: &[Feature]) -> usize {
    features
        .iter()
        .fold(read_length, |alignment_span, feature| match feature {
            Feature::Insertion { bases, .. } => alignment_span - bases.len(),
            Feature::InsertBase { .. } => alignment_span - 1,
            Feature::Deletion { len, .. } => alignment_span + len,
            Feature::ReferenceSkip { len, .. } => alignment_span + len,
            Feature::SoftClip { bases, .. } => alignment_span - bases.len(),
            _ => alignment_span,
        })
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_try_into_alignment_record_buf() -> Result<(), Box<dyn std::error::Error>> {
        use std::sync::Arc;

        use noodles_fasta as fasta;

        let header = sam::Header::default();
        let reference = Arc::new(fasta::record::Sequence::from(
            b"ACGTACGTacgtNNACGTACGT".to_vec(),
        ));

        let features = vec![
            Feature::SoftClip {
                position: Position::try_from(1)?,
                bases: Cow::from(b"TT"),
            },
            Feature::Substitution {
                position: Position::try_from(4)?,
                code: 0b01,
            },
            Feature::Scores {
                position: Position::try_from(5)?,
                quality_scores: Cow::from(&[30, 31]),
            },
            Feature::Insertion {
                position: Position::try_from(7)?,
                bases: Cow::from(b"GG"),
            },
            Feature::Deletion {
                position: Position::try_from(10)?,
                len: 3,
            },
            Feature::QualityScore {
                position: Position::try_from(11)?,
                quality_score: 40,
            },
        ];

        let mapped = Record {
            name: Some(Cow::from(&b"r0"[..])),
            bam_flags: sam::alignment::record::Flags::empty(),
            reference_sequence: Some(ReferenceSequence::External {
                sequence: reference.clone(),
            }),
            alignment_start: Some(Position::try_from(3)?),
            reference_sequence_id: Some(0),
            mapping_quality: MappingQuality::new(60),
            read_length: 14,
            features: features.clone(),
            ..Default::default()
        };

        let with_stored_scores = Record {
            cram_flags: Flags::QUALITY_SCORES_ARE_STORED_AS_ARRAY,
            quality_scores: Cow::from((0..14).collect::<Vec<u8>>()),
            ..mapped.clone()
        };

        let embedded = Record {
            reference_sequence: Some(ReferenceSequence::Embedded {
                reference_start: Position::try_from(2)?,
                sequence: &b"CGTACGTacgtNNACGTACGT"[..],
            }),
            ..mapped.clone()
        };

        let unmapped = Record {
            bam_flags: sam::alignment::record::Flags::UNMAPPED,
            sequence: Cow::from(&b"ACGTN"[..]),
            quality_scores: Cow::from(&[1, 2, 3, 4, 5][..]),
            read_length: 5,
            features: Vec::new(),
            ..Default::default()
        };

        for record in [mapped, with_stored_scores, embedded, unmapped] {
            let expected = sam::alignment::RecordBuf::try_from_alignment_record(&header, &record)?;
            assert_eq!(record.try_into_alignment_record_buf(&header)?, expected);
        }

        Ok(())
    }

    #[test]
    fn test_calculate_alignment_span() -> Result<(), noodles_core::position::TryFromIntError> {
        let features = [];
        assert_eq!(calculate_alignment_span(4, &features), 4);

        let features = [Feature::HardClip {
            position: Position::try_from(1)?,
            len: 4,
        }];
        assert_eq!(calculate_alignment_span(4, &features), 4);

        let features = [
            Feature::Insertion {
                position: Position::try_from(1)?,
                bases: Cow::from(b"AC"),
            },
            Feature::InsertBase {
                position: Position::try_from(4)?,
                base: b'G',
            },
            Feature::Deletion {
                position: Position::try_from(6)?,
                len: 3,
            },
            Feature::ReferenceSkip {
                position: Position::try_from(10)?,
                len: 5,
            },
            Feature::SoftClip {
                position: Position::try_from(16)?,
                bases: Cow::from(b"ACGT"),
            },
        ];
        assert_eq!(calculate_alignment_span(20, &features), 21);

        Ok(())
    }
}
