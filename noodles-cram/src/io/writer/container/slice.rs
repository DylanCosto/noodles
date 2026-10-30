mod header;
pub mod records;

use std::{collections::HashMap, io};

use noodles_core::Position;
use noodles_fasta as fasta;
use noodles_sam as sam;

pub use self::header::write_header;
use self::records::ExternalDataWriters;
use crate::{
    calculate_normalized_sequence_digest,
    codecs::Encoder,
    container::{
        BlockContentEncoderMap, CompressionHeader, ReferenceSequenceContext,
        block::{self, CompressionMethod, ContentType},
        slice::Header,
    },
    io::{
        BitWriter,
        writer::{Context, Record, container::block::Block},
    },
    record::Flags,
};

pub struct Slice {
    pub header: Header,
    pub core_data_block: Block,
    pub external_data_blocks: Vec<Block>,
}

pub(super) fn build_slice(
    ctx: &Context,
    header: &sam::Header,
    record_counter: u64,
    compression_header: &CompressionHeader,
    records: &mut [Record],
) -> io::Result<Slice> {
    let reference_sequence_context = get_reference_sequence_context(records);

    set_mates(records);

    let (core_data_buf, external_data_bufs) =
        write_records(compression_header, reference_sequence_context, records)?;

    let (core_data_block, external_data_blocks) = build_blocks(
        &ctx.block_content_encoder_map,
        records,
        core_data_buf,
        external_data_bufs,
    )?;

    let mut block_content_ids = vec![core_data_block.content_id];
    block_content_ids.extend(external_data_blocks.iter().map(|block| block.content_id));

    // Without a reference (records whose sequence the reference lacks), there is nothing to
    // checksum.
    let (reference_sequence_context, reference_md5) = if compression_header
        .preservation_map()
        .external_reference_sequence_is_required()
    {
        calculate_reference_sequence_md5(
            &ctx.reference_sequence_repository,
            header,
            reference_sequence_context,
        )?
    } else {
        (reference_sequence_context, None)
    };

    let header = Header {
        reference_sequence_context,
        record_count: records.len(),
        record_counter,
        block_count: block_content_ids.len(),
        block_content_ids,
        embedded_reference_bases_block_content_id: None,
        reference_md5,
        optional_tags: Vec::new(),
    };

    Ok(Slice {
        header,
        core_data_block,
        external_data_blocks,
    })
}

fn get_reference_sequence_context(records: &[Record]) -> ReferenceSequenceContext {
    assert!(!records.is_empty());

    let record = &records[0];
    let mut reference_sequence_context = match (
        record.reference_sequence_id,
        record.alignment_start,
        record.alignment_end(),
    ) {
        (Some(id), Some(start), Some(end)) => ReferenceSequenceContext::some(id, start, end),
        _ => ReferenceSequenceContext::None,
    };

    for record in records.iter().skip(1) {
        reference_sequence_context.update(
            record.reference_sequence_id,
            record.alignment_start,
            record.alignment_end(),
        );
    }

    reference_sequence_context
}

fn set_mates(records: &mut [Record]) {
    assert!(!records.is_empty());

    let mut indices = HashMap::new();
    let mut i = records.len() - 1;

    loop {
        let record = &mut records[i];
        let flags = record.bam_flags;

        if flags.is_segmented() && !flags.is_secondary() && !flags.is_supplementary() {
            let name = record.name.as_ref().map(|name| name.to_owned());

            if let Some(j) = indices.insert(name, i) {
                let mid = i + 1;
                let (left, right) = records.split_at_mut(mid);

                let record = &mut left[i];
                let mate = &mut right[j - mid];

                if decodes_as_attached(record, mate) {
                    set_downstream_mate(i, record, j, mate);
                } else {
                    set_detached(record);
                }
            } else {
                set_detached(record);
            }
        } else {
            set_detached(record);
        }

        if i == 0 {
            break;
        }

        i -= 1;
    }
}

/// Returns whether `record` and its downstream `mate` decode to their own mate fields when they
/// are attached (not detached), as HTSlib decodes them (`cram_decode_slice_xref`).
///
/// Decoders take each record's mate reference sequence, position, and mate flags from the other
/// record, and compute template lengths, so a pair whose stored values differ from those must be
/// detached. Like HTSlib's encoder (`process_one_read`), pairs are also kept detached when they
/// are on different reference sequences, are both first or both last segments, or have a
/// template length of 0 (which includes pairs with an unmapped segment). Only pairs are attached:
/// a record that already has a downstream mate is not.
fn decodes_as_attached(record: &Record, mate: &Record) -> bool {
    let (flags, mate_flags) = (record.bam_flags, mate.bam_flags);

    let points_to = |a: &Record, b: &Record| {
        a.mate_reference_sequence_id == b.reference_sequence_id
            && a.mate_alignment_start == b.alignment_start
            && a.bam_flags.is_mate_unmapped() == b.bam_flags.is_unmapped()
            && a.bam_flags.is_mate_reverse_complemented() == b.bam_flags.is_reverse_complemented()
    };

    !mate.cram_flags.contains(Flags::MATE_IS_DOWNSTREAM)
        && record.reference_sequence_id == mate.reference_sequence_id
        && record.template_length != 0
        && !(flags.is_first_segment() && mate_flags.is_first_segment())
        && !(flags.is_last_segment() && mate_flags.is_last_segment())
        && points_to(record, mate)
        && points_to(mate, record)
        && decoded_template_lengths(record, mate) == (record.template_length, mate.template_length)
}

/// Returns the template lengths that HTSlib decodes for an attached pair.
fn decoded_template_lengths(record: &Record, mate: &Record) -> (i32, i32) {
    fn extent(record: &Record) -> (usize, usize) {
        let start = record.alignment_start.map(usize::from).unwrap_or_default();

        if record.bam_flags.is_unmapped() {
            (start, start)
        } else {
            let end = (start + record.alignment_span()).saturating_sub(1);
            (start, end.max(start))
        }
    }

    if record.bam_flags.is_unmapped()
        || mate.bam_flags.is_unmapped()
        || record.reference_sequence_id.is_none()
        || record.reference_sequence_id != mate.reference_sequence_id
    {
        return (0, 0);
    }

    let (start, end) = extent(record);
    let (mate_start, mate_end) = extent(mate);

    let (left, right) = (start.min(mate_start), end.max(mate_end));
    let left_count = usize::from(start == left) + usize::from(mate_start == left);
    let right_count = usize::from(end == right) + usize::from(mate_end == right);

    let Ok(len) = i32::try_from(right - left + 1) else {
        return (0, 0);
    };

    let is_leftmost = start == left && (end < right || left_count <= 1);
    let is_tied = start == left && end == right && left_count > 1 && right_count > 1;

    if is_leftmost || (is_tied && record.bam_flags.is_first_segment()) {
        (len, -len)
    } else {
        (-len, len)
    }
}

fn set_downstream_mate(i: usize, record: &mut Record, j: usize, mate: &mut Record) {
    record.mate_distance = Some(j - i - 1);
    record.cram_flags.insert(Flags::MATE_IS_DOWNSTREAM);
    mate.cram_flags.remove(Flags::IS_DETACHED);
}

fn set_detached(record: &mut Record) {
    record.cram_flags.insert(Flags::IS_DETACHED);
}

#[allow(clippy::type_complexity)]
fn write_records(
    compression_header: &CompressionHeader,
    reference_sequence_context: ReferenceSequenceContext,
    records: &[Record],
) -> io::Result<(Vec<u8>, Vec<(block::ContentId, Vec<u8>)>)> {
    use crate::container::compression_header::data_series_encodings::data_series::STANDARD_DATA_SERIES;

    let mut core_data_writer = BitWriter::default();
    let mut external_data_writers = ExternalDataWriters::default();

    for data_series in STANDARD_DATA_SERIES {
        let block_content_id = block::ContentId::from(*data_series);
        external_data_writers.insert(block_content_id, Vec::new());
    }

    for block_content_id in compression_header.tag_encodings.keys() {
        external_data_writers.insert(*block_content_id, Vec::new());
    }

    let mut writer = records::Writer::new(
        compression_header,
        &mut core_data_writer,
        &mut external_data_writers,
        reference_sequence_context,
    )?;

    for record in records {
        writer.write_record(record)?;
    }

    Ok((
        core_data_writer.finish()?,
        external_data_writers.into_iter().collect(),
    ))
}

fn build_blocks(
    block_content_encoder_map: &BlockContentEncoderMap,
    records: &[Record],
    core_data_buf: Vec<u8>,
    external_data_bufs: Vec<(block::ContentId, Vec<u8>)>,
) -> io::Result<(Block, Vec<Block>)> {
    const CORE_DATA_BLOCK_CONTENT_ID: block::ContentId = 0;

    let core_data_block = Block::encode(
        ContentType::CoreData,
        CORE_DATA_BLOCK_CONTENT_ID,
        block_content_encoder_map.core_data_encoder(),
        &core_data_buf,
    )?;

    let mut all_quality_scores_stored_as_arrays = true;

    for record in records.iter() {
        all_quality_scores_stored_as_arrays = all_quality_scores_stored_as_arrays
            && record.cram_flags.quality_scores_are_stored_as_array();
    }

    let external_data_blocks = external_data_bufs
        .into_iter()
        .filter(|(_, buf)| !buf.is_empty())
        .map(|(block_content_id, buf)| {
            encode_block(
                block_content_encoder_map,
                records,
                all_quality_scores_stored_as_arrays,
                block_content_id,
                &buf,
            )
        })
        .collect::<io::Result<_>>()?;

    Ok((core_data_block, external_data_blocks))
}

fn encode_block(
    block_content_encoder_map: &BlockContentEncoderMap,
    records: &[Record],
    all_quality_scores_stored_as_arrays: bool,
    block_content_id: block::ContentId,
    src: &[u8],
) -> io::Result<Block> {
    use crate::codecs::fqzcomp;

    let content_type = ContentType::ExternalData;

    if let Some(encoder) = block_content_encoder_map.get_data_series_encoder(block_content_id) {
        match encoder {
            Some(Encoder::Fqzcomp) => {
                if all_quality_scores_stored_as_arrays {
                    let lens: Vec<_> = records.iter().map(|r| r.read_length).collect();
                    let data = fqzcomp::encode(&lens, src)?;

                    Ok(Block {
                        compression_method: CompressionMethod::Fqzcomp,
                        content_type,
                        content_id: block_content_id,
                        uncompressed_size: data.len(),
                        src: data,
                    })
                } else {
                    Block::encode(
                        ContentType::ExternalData,
                        block_content_id,
                        block_content_encoder_map.default_encoder(),
                        src,
                    )
                }
            }
            _ => Block::encode(ContentType::ExternalData, block_content_id, encoder, src),
        }
    } else if let Some(encoder) =
        block_content_encoder_map.get_tag_values_encoders(block_content_id)
    {
        Block::encode(ContentType::ExternalData, block_content_id, encoder, src)
    } else {
        Block::encode(
            ContentType::ExternalData,
            block_content_id,
            block_content_encoder_map.default_encoder(),
            src,
        )
    }
}

/// Returns the slice's reference sequence context, clipped to the end of the reference sequence
/// as in HTSlib, and the MD5 checksum of the reference bases it covers.
///
/// Records can overhang the end of the reference sequence. A slice that starts past the end has
/// no bases to checksum, so its checksum is left unset.
fn calculate_reference_sequence_md5(
    reference_sequence_repository: &fasta::Repository,
    header: &sam::Header,
    reference_sequence_context: ReferenceSequenceContext,
) -> io::Result<(ReferenceSequenceContext, Option<[u8; 16]>)> {
    let ReferenceSequenceContext::Some(context) = reference_sequence_context else {
        return Ok((reference_sequence_context, None));
    };

    let reference_sequence_name = header
        .reference_sequences()
        .get_index(context.reference_sequence_id())
        .map(|(name, _)| name)
        .expect("invalid reference sequence ID");

    let reference_sequence = reference_sequence_repository
        .get(reference_sequence_name)
        .ok_or_else(|| {
            io::Error::new(
                io::ErrorKind::InvalidInput,
                format!("missing reference sequence: {reference_sequence_name}"),
            )
        })??;

    let start = context.alignment_start();

    let Some(reference_end) = Position::new(reference_sequence.len()).filter(|end| *end >= start)
    else {
        return Ok((reference_sequence_context, None));
    };

    let end = context.alignment_end().min(reference_end);
    let sequence = &reference_sequence[start..=end];

    Ok((
        ReferenceSequenceContext::some(context.reference_sequence_id(), start, end),
        Some(calculate_normalized_sequence_digest(sequence)),
    ))
}

#[cfg(test)]
mod tests {
    use noodles_sam::alignment::record::Flags as BamFlags;

    use super::*;

    fn segment(flags: BamFlags, start: usize, mate_start: usize, template_length: i32) -> Record {
        Record {
            bam_flags: BamFlags::SEGMENTED | flags,
            reference_sequence_id: Some(0),
            read_length: 10,
            alignment_start: Position::new(start),
            name: Some("r".into()),
            mate_reference_sequence_id: Some(0),
            mate_alignment_start: Position::new(mate_start),
            template_length,
            ..Default::default()
        }
    }

    fn is_attached(mut records: [Record; 2]) -> bool {
        set_mates(&mut records);
        records[0].cram_flags.contains(Flags::MATE_IS_DOWNSTREAM)
            && !records[1].cram_flags.contains(Flags::IS_DETACHED)
    }

    #[test]
    fn test_set_mates() {
        let r1 = BamFlags::FIRST_SEGMENT;
        let r2 = BamFlags::LAST_SEGMENT;

        // A consistent pair is attached.
        assert!(is_attached([
            segment(r1, 100, 150, 60),
            segment(r2, 150, 100, -60),
        ]));

        // Template lengths that decoders would not reproduce.
        assert!(!is_attached([
            segment(r1, 100, 150, 999),
            segment(r2, 150, 100, 5),
        ]));
        assert!(!is_attached([
            segment(r1, 100, 150, 0),
            segment(r2, 150, 100, 0)
        ]));

        // Identical extents: decoders give the first segment the positive length.
        assert!(is_attached([
            segment(r2, 100, 100, -10),
            segment(r1, 100, 100, 10),
        ]));
        assert!(!is_attached([
            segment(r2, 100, 100, 10),
            segment(r1, 100, 100, -10),
        ]));

        // A mate position that is not the mate's.
        assert!(!is_attached([
            segment(r1, 100, 105, 60),
            segment(r2, 150, 100, -60),
        ]));

        // Mates on different reference sequences.
        let mut mate = segment(r2, 150, 100, 0);
        mate.reference_sequence_id = Some(1);
        let mut record = segment(r1, 100, 150, 0);
        record.mate_reference_sequence_id = Some(1);
        assert!(!is_attached([record, mate]));

        // An unmapped segment (template length 0), as HTSlib keeps it.
        let mut record = segment(r1 | BamFlags::MATE_UNMAPPED, 100, 100, 0);
        record.mate_alignment_start = Position::new(100);
        let mut mate = segment(r2 | BamFlags::UNMAPPED, 100, 100, 0);
        mate.read_length = 0;
        assert!(!is_attached([record, mate]));

        // Supplementary alignments.
        assert!(!is_attached([
            segment(r1 | BamFlags::SUPPLEMENTARY, 100, 150, 60),
            segment(r2, 150, 100, -60),
        ]));
    }
}
