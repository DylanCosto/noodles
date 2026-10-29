//! CRAM record fields to decode.

use noodles_sam::alignment::record::data::field::Tag;

use crate::container::{
    CompressionHeader, block,
    compression_header::{
        Encoding,
        encoding::codec::{Byte, ByteArray, Integer},
        preservation_map::tag_sets,
    },
};

/// The record fields a CRAM reader decodes.
///
/// By default, every field is decoded. Fields that are left out are absent from the records the
/// reader returns. When the external blocks that store a left-out field hold no other field, they
/// are not decompressed at all, which is where the time is saved: quality scores, for example, are
/// often most of a CRAM file.
///
/// # Examples
///
/// ```
/// use noodles_cram::io::reader::Fields;
/// use noodles_sam::alignment::record::data::field::Tag;
///
/// let fields = Fields::default()
///     .set_quality_scores(false)
///     .set_tags([Tag::ALIGNMENT_HIT_COUNT, Tag::EDIT_DISTANCE]);
/// ```
#[derive(Clone, Debug, Default, Eq, PartialEq)]
pub struct Fields {
    skip_quality_scores: bool,
    tags: Option<Vec<Tag>>,
}

impl Fields {
    /// Sets whether quality scores are decoded.
    ///
    /// Records read without quality scores have none, or, when the file stores them as read
    /// features rather than as an array, all zero.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_cram::io::reader::Fields;
    /// let fields = Fields::default().set_quality_scores(false);
    /// assert!(!fields.quality_scores());
    /// ```
    pub fn set_quality_scores(mut self, decode: bool) -> Self {
        self.skip_quality_scores = !decode;
        self
    }

    /// Sets the data fields to decode; others are left out.
    ///
    /// This includes the read group (`RG`), which CRAM stores as a data series.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_cram::io::reader::Fields;
    /// use noodles_sam::alignment::record::data::field::Tag;
    ///
    /// let fields = Fields::default().set_tags([Tag::EDIT_DISTANCE]);
    /// assert!(fields.decodes_tag(Tag::EDIT_DISTANCE));
    /// assert!(!fields.decodes_tag(Tag::READ_GROUP));
    /// ```
    pub fn set_tags<I>(mut self, tags: I) -> Self
    where
        I: IntoIterator<Item = Tag>,
    {
        self.tags = Some(tags.into_iter().collect());
        self
    }

    /// Returns whether quality scores are decoded.
    pub fn quality_scores(&self) -> bool {
        !self.skip_quality_scores
    }

    /// Returns whether the data field with the given tag is decoded.
    pub fn decodes_tag(&self, tag: Tag) -> bool {
        self.tags.as_ref().is_none_or(|tags| tags.contains(&tag))
    }
}

/// What a record reader does with a field.
#[derive(Clone, Copy, Debug, Default, Eq, PartialEq)]
pub(crate) enum Action {
    /// Decodes the field.
    #[default]
    Decode,
    /// Decodes the field to advance its readers but leaves it out of the record.
    Discard,
    /// Leaves the field out without reading it; its blocks are not decompressed.
    Skip,
}

/// How to decode the slices of one container for a set of fields. The default decodes every
/// field.
#[derive(Debug, Default)]
pub(crate) struct Plan {
    quality_scores: Action,
    skip_read_group_ids: bool,
    tags: Vec<(tag_sets::Key, Action)>,
    skipped_block_content_ids: Vec<block::ContentId>,
}

impl Plan {
    pub(crate) fn new(fields: &Fields, compression_header: &CompressionHeader) -> Self {
        if *fields == Fields::default() {
            return Self::default();
        }

        let data_series_encodings = compression_header.data_series_encodings();

        // External blocks read by fields that are decoded.
        let mut needed = Vec::new();

        for encoding in [
            &data_series_encodings.bam_flags,
            &data_series_encodings.cram_flags,
            &data_series_encodings.reference_sequence_ids,
            &data_series_encodings.read_lengths,
            &data_series_encodings.alignment_starts,
            &data_series_encodings.read_group_ids,
            &data_series_encodings.mate_flags,
            &data_series_encodings.mate_reference_sequence_ids,
            &data_series_encodings.mate_alignment_starts,
            &data_series_encodings.template_lengths,
            &data_series_encodings.mate_distances,
            &data_series_encodings.tag_set_ids,
            &data_series_encodings.feature_counts,
            &data_series_encodings.feature_position_deltas,
            &data_series_encodings.deletion_lengths,
            &data_series_encodings.reference_skip_lengths,
            &data_series_encodings.padding_lengths,
            &data_series_encodings.hard_clip_lengths,
            &data_series_encodings.mapping_qualities,
        ]
        .into_iter()
        .flatten()
        {
            encoding.get().sources(&mut needed);
        }

        for encoding in [
            &data_series_encodings.feature_codes,
            &data_series_encodings.base_substitution_codes,
            &data_series_encodings.bases,
        ]
        .into_iter()
        .flatten()
        {
            encoding.get().sources(&mut needed);
        }

        for encoding in [
            &data_series_encodings.names,
            &data_series_encodings.stretches_of_bases,
            &data_series_encodings.insertion_bases,
            &data_series_encodings.soft_clip_bases,
        ]
        .into_iter()
        .flatten()
        {
            encoding.get().sources(&mut needed);
        }

        // Fields that are left out, with the blocks they read and whether they read the core data
        // block. The first is the quality scores, if they are left out.
        let mut left_out: Vec<(Option<tag_sets::Key>, Vec<block::ContentId>, bool, Action)> =
            Vec::new();

        let quality_score_encodings = (
            data_series_encodings.quality_scores.as_ref(),
            data_series_encodings.stretches_of_quality_scores.as_ref(),
        );

        if fields.quality_scores() {
            sources(quality_score_encodings.0, &mut needed);
            sources(quality_score_encodings.1, &mut needed);
        } else {
            let mut ids = Vec::new();
            let reads_core = sources(quality_score_encodings.0, &mut ids)
                | sources(quality_score_encodings.1, &mut ids);
            left_out.push((None, ids, reads_core, Action::Skip));
        }

        let mut tags = Vec::new();

        for &key in compression_header
            .preservation_map()
            .tag_sets()
            .iter()
            .flatten()
        {
            if tags.iter().any(|(k, _)| *k == key) {
                continue;
            }

            let Some(encoding) = compression_header
                .tag_encodings()
                .get(&block::ContentId::from(key))
            else {
                continue;
            };

            if fields.decodes_tag(key.tag()) {
                encoding.get().sources(&mut needed);
                tags.push((key, Action::Decode));
            } else {
                let mut ids = Vec::new();
                let reads_core = encoding.get().sources(&mut ids);
                left_out.push((Some(key), ids, reads_core, Action::Skip));
            }
        }

        // A left-out field must still be read if it reads the core data block or a block that is
        // decompressed anyway; its blocks are then decompressed too.
        loop {
            let mut changed = false;

            for (_, ids, reads_core, action) in &mut left_out {
                if *action == Action::Skip
                    && (*reads_core || ids.iter().any(|id| needed.contains(id)))
                {
                    *action = Action::Discard;
                    needed.extend_from_slice(ids);
                    changed = true;
                }
            }

            if !changed {
                break;
            }
        }

        let mut plan = Self {
            skip_read_group_ids: !fields.decodes_tag(Tag::READ_GROUP),
            tags,
            ..Default::default()
        };

        for (key, ids, _, action) in left_out {
            if action == Action::Skip {
                plan.skipped_block_content_ids.extend(ids);
            }

            match key {
                Some(key) => plan.tags.push((key, action)),
                None => plan.quality_scores = action,
            }
        }

        plan
    }

    pub(crate) fn quality_scores(&self) -> Action {
        self.quality_scores
    }

    pub(crate) fn read_group_ids(&self) -> bool {
        !self.skip_read_group_ids
    }

    pub(crate) fn tag(&self, key: tag_sets::Key) -> Action {
        self.tags
            .iter()
            .find(|(k, _)| *k == key)
            .map(|(_, action)| *action)
            .unwrap_or_default()
    }

    /// Returns the content IDs of the external blocks that no decoded field reads.
    pub(crate) fn skipped_block_content_ids(&self) -> &[block::ContentId] {
        &self.skipped_block_content_ids
    }
}

fn sources<C>(encoding: Option<&Encoding<C>>, ids: &mut Vec<block::ContentId>) -> bool
where
    C: Sources,
{
    encoding.is_some_and(|encoding| encoding.get().sources(ids))
}

/// The blocks an encoding reads values from.
trait Sources {
    /// Adds the content IDs of the external blocks the encoding reads and returns whether it also
    /// reads the core data block.
    fn sources(&self, ids: &mut Vec<block::ContentId>) -> bool;
}

impl Sources for Integer {
    fn sources(&self, ids: &mut Vec<block::ContentId>) -> bool {
        match self {
            Self::External { block_content_id } => {
                ids.push(*block_content_id);
                false
            }
            _ => true,
        }
    }
}

impl Sources for Byte {
    fn sources(&self, ids: &mut Vec<block::ContentId>) -> bool {
        match self {
            Self::External { block_content_id } => {
                ids.push(*block_content_id);
                false
            }
            Self::Huffman { .. } => true,
        }
    }
}

impl Sources for ByteArray {
    fn sources(&self, ids: &mut Vec<block::ContentId>) -> bool {
        match self {
            Self::ByteArrayLength {
                len_encoding,
                value_encoding,
            } => len_encoding.get().sources(ids) | value_encoding.get().sources(ids),
            Self::ByteArrayStop {
                block_content_id, ..
            } => {
                ids.push(*block_content_id);
                false
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use std::{io, num::NonZero};

    use noodles_core::Position;
    use noodles_fasta as fasta;
    use noodles_sam::{
        self as sam,
        alignment::{
            RecordBuf,
            io::Write as _,
            record::{
                Flags,
                cigar::{Op, op::Kind},
                data::field::Type,
            },
            record_buf::{Sequence, data::field::Value},
        },
        header::record::value::{
            Map,
            map::{ReadGroup, ReferenceSequence},
        },
    };

    use super::*;
    use crate::{
        container::compression_header::{
            PreservationMap, TagEncodings,
            data_series_encodings::{DataSeries, DataSeriesEncodings},
        },
        io::reader::Container,
    };

    const MODIFIED_BASES: Tag = Tag::new(b'M', b'M');

    fn build_cram() -> io::Result<(sam::Header, fasta::Repository, Vec<u8>)> {
        let repository = fasta::Repository::new(vec![fasta::Record::new(
            fasta::record::Definition::new("sq0", None),
            fasta::record::Sequence::from(b"ACGTACGTAC".to_vec()),
        )]);

        let header = sam::Header::builder()
            .set_header(Default::default())
            .add_reference_sequence(
                "sq0",
                Map::<ReferenceSequence>::new(NonZero::new(10).unwrap()),
            )
            .add_read_group("rg0", Map::<ReadGroup>::default())
            .build();

        let read_group = || (Tag::READ_GROUP, Value::from("rg0"));

        let records = [
            RecordBuf::builder()
                .set_name("r0")
                .set_flags(Flags::empty())
                .set_reference_sequence_id(0)
                .set_alignment_start(Position::MIN)
                .set_cigar([Op::new(Kind::Match, 4)].into_iter().collect())
                .set_sequence(Sequence::from(b"ACGA"))
                .set_quality_scores([10, 20, 30, 40].into_iter().collect())
                .set_data(
                    [
                        (Tag::EDIT_DISTANCE, Value::from(1)),
                        (MODIFIED_BASES, Value::from("C+m?;")),
                        read_group(),
                    ]
                    .into_iter()
                    .collect(),
                )
                .build(),
            RecordBuf::builder()
                .set_name("r1")
                .set_flags(Flags::empty())
                .set_reference_sequence_id(0)
                .set_alignment_start(Position::try_from(3).unwrap())
                .set_cigar(
                    [
                        Op::new(Kind::Match, 2),
                        Op::new(Kind::Insertion, 1),
                        Op::new(Kind::Match, 1),
                    ]
                    .into_iter()
                    .collect(),
                )
                .set_sequence(Sequence::from(b"GTTA"))
                .set_quality_scores([5, 6, 7, 8].into_iter().collect())
                .set_data(
                    [(Tag::EDIT_DISTANCE, Value::from(1)), read_group()]
                        .into_iter()
                        .collect(),
                )
                .build(),
            RecordBuf::builder()
                .set_name("r2")
                .set_flags(Flags::UNMAPPED)
                .set_sequence(Sequence::from(b"ACGT"))
                .set_quality_scores([1, 2, 3, 4].into_iter().collect())
                .set_data([read_group()].into_iter().collect())
                .build(),
        ];

        let mut writer = crate::io::writer::Builder::default()
            .set_reference_sequence_repository(repository.clone())
            .build_from_writer(Vec::new());

        writer.write_header(&header)?;

        for record in &records {
            writer.write_alignment_record(&header, record)?;
        }

        writer.try_finish(&header)?;

        Ok((header, repository, writer.into_inner()))
    }

    fn read_records(
        src: &[u8],
        repository: &fasta::Repository,
        fields: Fields,
    ) -> io::Result<Vec<RecordBuf>> {
        let mut reader = crate::io::reader::Builder::default()
            .set_reference_sequence_repository(repository.clone())
            .set_fields(fields)
            .build_from_reader(src);

        let header = reader.read_header()?;
        reader.records(&header).collect()
    }

    /// Reads records through the container API with the plan for `fields`, optionally turning
    /// every skipped field into a discarded one, and returns them with the IDs of the external
    /// blocks that were decompressed.
    fn read_records_with_plan(
        src: &[u8],
        repository: &fasta::Repository,
        fields: &Fields,
        discard_instead_of_skip: bool,
    ) -> io::Result<(Vec<RecordBuf>, Vec<block::ContentId>)> {
        let mut reader = crate::io::Reader::new(src);
        let header = reader.read_header()?;

        let mut container = Container::default();
        let mut records = Vec::new();
        let mut decoded_block_content_ids = Vec::new();

        while reader.read_container(&mut container)? != 0 {
            let compression_header = container.compression_header()?;
            let mut plan = Plan::new(fields, &compression_header);

            if discard_instead_of_skip {
                let discard = |action: &mut Action| {
                    if *action == Action::Skip {
                        *action = Action::Discard;
                    }
                };

                discard(&mut plan.quality_scores);
                plan.tags.iter_mut().for_each(|(_, action)| discard(action));
                plan.skipped_block_content_ids.clear();
            }

            for result in container.slices() {
                let slice = result?;
                let (core_data_src, external_data_srcs) = slice.decode_blocks_with_plan(&plan)?;

                decoded_block_content_ids.extend(external_data_srcs.iter().map(|(id, _)| *id));

                for record in slice.records_with_plan(
                    repository.clone(),
                    &header,
                    &compression_header,
                    &core_data_src,
                    &external_data_srcs,
                    &plan,
                )? {
                    records.push(RecordBuf::try_from_alignment_record(&header, &record)?);
                }
            }
        }

        Ok((records, decoded_block_content_ids))
    }

    /// The records as read with every field, keeping only the given tags and, if asked, dropping
    /// quality scores.
    fn without_fields(records: &[RecordBuf], quality_scores: bool, tags: &[Tag]) -> Vec<RecordBuf> {
        records
            .iter()
            .cloned()
            .map(|mut record| {
                if !quality_scores {
                    record.quality_scores_mut().as_mut().clear();
                }

                let data = record
                    .data()
                    .iter()
                    .filter(|(tag, _)| tags.contains(tag))
                    .map(|(tag, value)| (tag, value.clone()))
                    .collect();
                *record.data_mut() = data;

                record
            })
            .collect()
    }

    #[test]
    fn test_read_records_without_quality_scores_or_tags() -> io::Result<()> {
        let (_, repository, src) = build_cram()?;

        let all = read_records(&src, &repository, Fields::default())?;
        assert_eq!(all.len(), 3);
        assert_eq!(all[0].quality_scores().as_ref(), [10, 20, 30, 40]);
        assert!(all[0].data().get(&MODIFIED_BASES).is_some());

        let fields = Fields::default()
            .set_quality_scores(false)
            .set_tags([Tag::EDIT_DISTANCE]);
        let expected = without_fields(&all, false, &[Tag::EDIT_DISTANCE]);

        assert_eq!(read_records(&src, &repository, fields.clone())?, expected);

        // Skipped fields are not decompressed; discarded ones are, with the same records.
        let (records, decoded) = read_records_with_plan(&src, &repository, &fields, false)?;
        assert_eq!(records, expected);
        assert!(!decoded.contains(&block::ContentId::from(DataSeries::QualityScores)));
        assert!(
            !decoded.contains(&block::ContentId::from(tag_sets::Key::new(
                MODIFIED_BASES,
                Type::String
            )))
        );

        let (records, decoded) = read_records_with_plan(&src, &repository, &fields, true)?;
        assert_eq!(records, expected);
        assert!(decoded.contains(&block::ContentId::from(DataSeries::QualityScores)));

        // The read group is a data series, not a tag block.
        let fields = Fields::default().set_tags([Tag::READ_GROUP]);
        assert_eq!(
            read_records(&src, &repository, fields)?,
            without_fields(&all, true, &[Tag::READ_GROUP])
        );

        let fields = Fields::default().set_tags([]);
        assert_eq!(
            read_records(&src, &repository, fields)?,
            without_fields(&all, true, &[])
        );

        Ok(())
    }

    #[test]
    fn test_plan_discards_fields_sharing_a_decoded_block() {
        let bases = block::ContentId::from(DataSeries::Bases);
        let edit_distance = tag_sets::Key::new(Tag::EDIT_DISTANCE, Type::Int32);
        let modified_bases = tag_sets::Key::new(MODIFIED_BASES, Type::String);

        let mut data_series_encodings = DataSeriesEncodings::init();
        data_series_encodings.quality_scores = Some(Encoding::new(Byte::External {
            block_content_id: bases,
        }));

        let preservation_map = PreservationMap {
            tag_sets: vec![vec![edit_distance, modified_bases]],
            ..Default::default()
        };

        // MM's lengths are stored in NM's block.
        let mut tag_encodings = TagEncodings::new();
        tag_encodings.insert(
            block::ContentId::from(edit_distance),
            Encoding::new(ByteArray::ByteArrayStop {
                stop_byte: 0,
                block_content_id: block::ContentId::from(edit_distance),
            }),
        );
        tag_encodings.insert(
            block::ContentId::from(modified_bases),
            Encoding::new(ByteArray::ByteArrayLength {
                len_encoding: Encoding::new(Integer::External {
                    block_content_id: block::ContentId::from(edit_distance),
                }),
                value_encoding: Encoding::new(Byte::External {
                    block_content_id: block::ContentId::from(modified_bases),
                }),
            }),
        );

        let compression_header =
            CompressionHeader::new(preservation_map, data_series_encodings, tag_encodings);

        let fields = Fields::default()
            .set_quality_scores(false)
            .set_tags([Tag::EDIT_DISTANCE]);
        let plan = Plan::new(&fields, &compression_header);

        assert_eq!(plan.quality_scores(), Action::Discard);
        assert_eq!(plan.tag(edit_distance), Action::Decode);
        assert_eq!(plan.tag(modified_bases), Action::Discard);
        assert!(plan.skipped_block_content_ids().is_empty());

        // Without NM, nothing else reads MM's blocks.
        let plan = Plan::new(&Fields::default().set_tags([]), &compression_header);
        assert_eq!(plan.quality_scores(), Action::Decode);
        assert_eq!(plan.tag(edit_distance), Action::Skip);
        assert_eq!(plan.tag(modified_bases), Action::Skip);
        assert!(!plan.read_group_ids());

        let mut skipped = plan.skipped_block_content_ids().to_vec();
        skipped.sort_unstable();
        skipped.dedup();
        let mut expected = vec![
            block::ContentId::from(edit_distance),
            block::ContentId::from(modified_bases),
        ];
        expected.sort_unstable();
        assert_eq!(skipped, expected);
    }

    #[test]
    fn test_default_plan_decodes_everything() {
        let compression_header = CompressionHeader::new(
            PreservationMap::default(),
            DataSeriesEncodings::init(),
            TagEncodings::new(),
        );

        for plan in [
            Plan::default(),
            Plan::new(&Fields::default(), &compression_header),
        ] {
            assert_eq!(plan.quality_scores(), Action::Decode);
            assert!(plan.read_group_ids());
            assert_eq!(
                plan.tag(tag_sets::Key::new(MODIFIED_BASES, Type::String)),
                Action::Decode
            );
            assert!(plan.skipped_block_content_ids().is_empty());
        }
    }

    #[test]
    fn test_slice_records_decode_every_field() -> io::Result<()> {
        let (_, repository, src) = build_cram()?;

        let (records, _) = read_records_with_plan(&src, &repository, &Fields::default(), false)?;
        assert_eq!(records, read_records(&src, &repository, Fields::default())?);
        assert!(records[0].data().get(&Tag::READ_GROUP).is_some());

        Ok(())
    }
}
