//! A materialized sequence with borrowed record metadata and auxiliary values.
use super::*;

pub(super) struct SequenceRecord<'r, 'c> {
    pub(super) record: &'r Record<'c>,
    pub(super) sequence: &'r [u8],
}

impl sam::alignment::Record for SequenceRecord<'_, '_> {
    fn name(&self) -> Option<&BStr> {
        self.record.name.as_deref().map(|name| name.as_bstr())
    }

    fn flags(&self) -> io::Result<sam::alignment::record::Flags> {
        Ok(self.record.bam_flags)
    }

    fn reference_sequence_id<'r, 'h: 'r>(
        &'r self,
        _: &'h sam::Header,
    ) -> Option<io::Result<usize>> {
        self.record.reference_sequence_id.map(Ok)
    }

    fn alignment_start(&self) -> Option<io::Result<Position>> {
        self.record.alignment_start.map(Ok)
    }

    fn mapping_quality(&self) -> Option<io::Result<MappingQuality>> {
        self.record.mapping_quality.map(Ok)
    }

    fn cigar(&self) -> Box<dyn sam::alignment::record::Cigar + '_> {
        Box::new(Cigar::new(
            &self.record.features,
            self.record.bam_flags.is_unmapped(),
            self.record.read_length,
        ))
    }

    fn mate_reference_sequence_id<'r, 'h: 'r>(
        &'r self,
        _: &'h sam::Header,
    ) -> Option<io::Result<usize>> {
        self.record.mate_reference_sequence_id.map(Ok)
    }

    fn mate_alignment_start(&self) -> Option<io::Result<Position>> {
        self.record.mate_alignment_start.map(Ok)
    }

    fn template_length(&self) -> io::Result<i32> {
        Ok(self.record.template_length)
    }

    fn sequence(&self) -> Box<dyn sam::alignment::record::Sequence + '_> {
        Box::new(Bases(self.sequence))
    }
    fn sequence_ref(&self) -> sam::alignment::record::SequenceRef<'_> {
        sam::alignment::record::SequenceRef::Raw(self.sequence)
    }
    fn quality_scores_ref(&self) -> sam::alignment::record::QualityScoresRef<'_> {
        if self.record.bam_flags.is_unmapped()
            || self.record.cram_flags.quality_scores_are_stored_as_array()
        {
            sam::alignment::record::QualityScoresRef::Raw(&self.record.quality_scores)
        } else {
            sam::alignment::record::QualityScoresRef::QualityScores(self.quality_scores())
        }
    }

    fn quality_scores(&self) -> Box<dyn sam::alignment::record::QualityScores + '_> {
        if self.record.bam_flags.is_unmapped()
            || self.record.cram_flags.quality_scores_are_stored_as_array()
        {
            Box::new(Scores(&self.record.quality_scores[..]))
        } else {
            Box::new(QualityScores::new(
                &self.record.features,
                self.record.read_length,
            ))
        }
    }

    fn data(&self) -> Box<dyn sam::alignment::record::Data<'_> + '_> {
        if let Some(header) = self.record.header {
            Box::new(Data::new(
                header,
                &self.record.data,
                self.record.read_group_id,
            ))
        } else {
            Box::new(EmptyData::new())
        }
    }

    fn alignment_span(&self) -> Option<io::Result<usize>> {
        Some(Ok(self.record.alignment_span()))
    }

    fn alignment_end(&self) -> Option<io::Result<Position>> {
        self.record.alignment_end().map(Ok)
    }
}
