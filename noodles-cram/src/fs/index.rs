use std::{
    cmp,
    collections::HashMap,
    fs::File,
    io::{self, Read, Seek},
    path::Path,
};

use noodles_core::Position;

use crate::{
    container::{CompressionHeader, slice},
    crai,
    io::{
        Reader,
        reader::{Container, container::Slice},
    },
};

/// Indexes a CRAM file.
///
/// # Examples
///
/// ```no_run
/// use noodles_cram as cram;
/// let index = cram::fs::index("sample.cram")?;
/// # Ok::<(), std::io::Error>(())
/// ```
pub fn index<P>(src: P) -> io::Result<crai::Index>
where
    P: AsRef<Path>,
{
    let mut reader = File::open(src).map(Reader::new)?;
    index_reader(&mut reader)
}

fn index_reader<R>(reader: &mut Reader<R>) -> io::Result<crai::Index>
where
    R: Read + Seek,
{
    reader.read_header()?;

    let mut index = Vec::new();

    let mut container = Container::default();
    let mut container_position = reader.position()?;

    loop {
        let container_len = match reader.read_container(&mut container)? {
            0 => break,
            n => n,
        };

        let compression_header = container.compression_header()?;

        let landmarks = container.header().landmarks();
        let slice_count = landmarks.len();

        for (i, result) in container.slices().enumerate() {
            let slice = result?;
            let landmark = landmarks[i];

            let slice_length = if i < slice_count - 1 {
                landmarks[i + 1] - landmark
            } else {
                container_len - landmark
            };

            push_index_records(
                &mut index,
                &compression_header,
                &slice,
                container_position,
                landmark as u64,
                slice_length as u64,
            )?;
        }

        container_position = reader.position()?;
    }

    Ok(index)
}

fn push_index_records(
    index: &mut crai::Index,
    compression_header: &CompressionHeader,
    slice: &Slice,
    container_position: u64,
    landmark: u64,
    slice_length: u64,
) -> io::Result<()> {
    if slice.header().reference_sequence_context().is_many() {
        push_index_records_for_multi_reference_slice(
            index,
            compression_header,
            slice,
            container_position,
            landmark,
            slice_length,
        )
    } else {
        push_index_record_for_single_reference_slice(
            index,
            slice.header(),
            container_position,
            landmark,
            slice_length,
        )
    }
}

#[derive(Debug)]
struct SliceReferenceSequenceAlignmentRangeInclusive {
    start: Option<Position>,
    end: Option<Position>,
}

impl Default for SliceReferenceSequenceAlignmentRangeInclusive {
    fn default() -> Self {
        Self {
            start: Position::new(usize::MAX),
            end: None,
        }
    }
}

fn push_index_records_for_multi_reference_slice(
    index: &mut crai::Index,
    compression_header: &CompressionHeader,
    slice: &Slice,
    container_position: u64,
    landmark: u64,
    slice_length: u64,
) -> io::Result<()> {
    let mut reference_sequence_ids: HashMap<
        Option<usize>,
        SliceReferenceSequenceAlignmentRangeInclusive,
    > = HashMap::new();

    let (core_data_src, external_data_srcs) = slice.decode_blocks()?;

    // Alignment spans need positions and features, not reference sequences.
    for record in slice.records_without_reference_sequences(
        compression_header,
        &core_data_src,
        &external_data_srcs,
    )? {
        let range = reference_sequence_ids
            .entry(record.reference_sequence_id)
            .or_default();

        range.start = cmp::min(range.start, record.alignment_start);

        let alignment_end = record.alignment_end();
        range.end = cmp::max(range.end, alignment_end);
    }

    let mut sorted_reference_sequence_ids: Vec<_> =
        reference_sequence_ids.keys().copied().collect();
    // Unplaced records come last, as in a coordinate-sorted slice.
    sorted_reference_sequence_ids.sort_unstable_by_key(|id| (id.is_none(), *id));

    for reference_sequence_id in sorted_reference_sequence_ids {
        let (alignment_start, alignment_span) = if reference_sequence_id.is_some() {
            let range = &reference_sequence_ids[&reference_sequence_id];

            if let (Some(start), Some(end)) = (range.start, range.end) {
                let span = usize::from(end) - usize::from(start) + 1;
                (Some(start), span)
            } else {
                return Err(io::Error::new(
                    io::ErrorKind::InvalidData,
                    format!("invalid alignment range in slice: {range:?}"),
                ));
            }
        } else {
            // As in HTSlib (`cram_index_build_multiref`), which spans unplaced records from
            // position 0 to 0, inclusive. HTSlib does not find unplaced records through entries
            // with a span of 0 in multi-reference slices.
            (None, 1)
        };

        let record = crai::Record::new(
            reference_sequence_id,
            alignment_start,
            alignment_span,
            container_position,
            landmark,
            slice_length,
        );

        index.push(record);
    }

    Ok(())
}

fn push_index_record_for_single_reference_slice(
    index: &mut crai::Index,
    slice_header: &slice::Header,
    container_position: u64,
    landmark: u64,
    slice_length: u64,
) -> io::Result<()> {
    use crate::container::ReferenceSequenceContext;

    let (reference_sequence_id, alignment_start, alignment_span) =
        match slice_header.reference_sequence_context() {
            ReferenceSequenceContext::Some(context) => {
                let reference_sequence_id = Some(context.reference_sequence_id());
                let alignment_start = Some(context.alignment_start());
                let alignment_span = context.alignment_span();
                (reference_sequence_id, alignment_start, alignment_span)
            }
            ReferenceSequenceContext::None => (None, None, 0),
            ReferenceSequenceContext::Many => unreachable!(),
        };

    let record = crai::Record::new(
        reference_sequence_id,
        alignment_start,
        alignment_span,
        container_position,
        landmark,
        slice_length,
    );

    index.push(record);

    Ok(())
}

#[cfg(test)]
mod tests {
    use std::io::Cursor;

    use noodles_fasta as fasta;
    use noodles_sam::{
        self as sam,
        alignment::{
            RecordBuf,
            io::Write as _,
            record::{
                Flags,
                cigar::{Op, op::Kind},
            },
        },
    };

    use super::*;

    #[test]
    fn test_index_reader_with_multi_reference_slice() -> Result<(), Box<dyn std::error::Error>> {
        let header: sam::Header = "@SQ\tSN:sq0\tLN:8\n@SQ\tSN:sq1\tLN:8\n".parse()?;
        let repository = fasta::Repository::new(
            [b"sq0", b"sq1"]
                .map(|name| {
                    let definition = fasta::record::Definition::new(name, None);
                    fasta::Record::new(definition, b"ACGTACGT".to_vec().into())
                })
                .to_vec(),
        );

        // One slice with records on two reference sequences and an unplaced record.
        let mut writer = crate::io::writer::Builder::default()
            .set_reference_sequence_repository(repository)
            .build_from_writer(Vec::new());
        writer.write_header(&header)?;

        for (reference_sequence_id, start, len) in [(0, 2, 4), (0, 3, 5), (1, 1, 3)] {
            let record = RecordBuf::builder()
                .set_flags(Flags::empty())
                .set_reference_sequence_id(reference_sequence_id)
                .set_alignment_start(Position::try_from(start)?)
                .set_cigar([Op::new(Kind::Match, len)].into_iter().collect())
                .set_sequence(b"ACGTACGT"[..len].to_vec().into())
                .set_quality_scores(vec![30; len].into())
                .build();
            writer.write_alignment_record(&header, &record)?;
        }

        let unplaced = RecordBuf::builder()
            .set_flags(Flags::UNMAPPED)
            .set_sequence(b"ACGT".to_vec().into())
            .set_quality_scores(vec![30; 4].into())
            .build();
        writer.write_alignment_record(&header, &unplaced)?;

        writer.try_finish(&header)?;
        let src = writer.into_inner();

        // Indexing needs no reference sequences.
        let mut reader = Reader::new(Cursor::new(src));
        let index = index_reader(&mut reader)?;

        let actual: Vec<_> = index
            .iter()
            .map(|record| {
                (
                    record.reference_sequence_id(),
                    record.alignment_start(),
                    record.alignment_span(),
                )
            })
            .collect();

        let expected = [
            (Some(0), Position::new(2), 6),
            (Some(1), Position::new(1), 3),
            (None, None, 1),
        ];

        assert_eq!(actual, expected);

        // All entries are for the same slice.
        assert!(index.iter().all(|record| {
            record.offset() == index[0].offset() && record.landmark() == index[0].landmark()
        }));

        Ok(())
    }
}
