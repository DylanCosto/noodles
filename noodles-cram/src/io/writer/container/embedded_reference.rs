//! Reference sequences built from the records and embedded in slices, as HTSlib's `embed_ref=2`
//! does for sequences missing from the external reference.

use std::io;

use noodles_sam::alignment::{
    record::data::field::Tag, record_buf::data::field::Value as ValueBuf,
};

use crate::{
    container::{ReferenceSequenceContext, block},
    io::writer::{
        Record,
        record::{Feature, push_mismatches},
    },
};

/// The block content ID of an embedded reference sequence. It is past the data series' IDs, and
/// tag blocks' IDs are larger.
pub(crate) const BLOCK_CONTENT_ID: block::ContentId = 31;

/// HTSlib's CRAM flags tag: bit 1 means the record has no MD tag, bit 2 no NM tag, so readers
/// do not regenerate them from the embedded (consensus) reference.
const CRAM_FLAGS: Tag = Tag::new(b'c', b'F');

/// Builds a consensus reference sequence for a slice whose records store their aligned bases
/// verbatim and are all on one reference sequence, and re-encodes the records' bases as
/// differences from it.
///
/// This returns `None`, leaving the records unchanged, for slices on many (or no) reference
/// sequences, which cannot embed a reference.
pub(crate) fn embed(
    reference_sequence_context: ReferenceSequenceContext,
    records: &mut [Record],
) -> io::Result<Option<Vec<u8>>> {
    let ReferenceSequenceContext::Some(context) = reference_sequence_context else {
        return Ok(None);
    };

    let start = usize::from(context.alignment_start());
    let end = usize::from(context.alignment_end());
    let reference_sequence = consensus(records, start, end - start + 1);

    for record in records
        .iter_mut()
        .filter(|record| has_aligned_bases(record))
    {
        let offset = record_offset(record, start);
        record.features = reencode(record, &reference_sequence, offset)?;
        add_cram_flags(record);
    }

    Ok(Some(reference_sequence))
}

/// Whether a record has bases aligned to the reference, i.e., it is mapped and has a sequence.
fn has_aligned_bases(record: &Record) -> bool {
    !record.bam_flags.is_unmapped() && !record.sequence.is_empty()
}

/// The 0-based offset of a record's alignment start in the slice's reference sequence.
fn record_offset(record: &Record, slice_start: usize) -> usize {
    record
        .alignment_start
        .map(|start| usize::from(start) - slice_start)
        .unwrap_or_default()
}

/// Calls `f` with each stretch of aligned bases in verbatim features: its 0-based reference
/// offset from the alignment start, its 1-based read position, and its bases.
fn for_each_stretch<'f, F>(features: &'f [Feature], mut f: F) -> io::Result<()>
where
    F: FnMut(usize, usize, &'f [u8]) -> io::Result<()>,
{
    let mut read_position = 1;
    let mut reference_offset = 0;

    for feature in features {
        let position = usize::from(feature.position());

        // Aligned bases without features (none in verbatim records) still use the reference.
        reference_offset += position.saturating_sub(read_position);
        read_position = read_position.max(position);

        match feature {
            Feature::Bases { bases, .. } => {
                f(reference_offset, position, bases)?;
                reference_offset += bases.len();
                read_position += bases.len();
            }
            Feature::Insertion { bases, .. } | Feature::SoftClip { bases, .. } => {
                read_position += bases.len();
            }
            Feature::InsertBase { .. } => read_position += 1,
            Feature::Deletion { len, .. } | Feature::ReferenceSkip { len, .. } => {
                reference_offset += len;
            }
            _ => {}
        }
    }

    Ok(())
}

/// Builds the consensus of `len` reference positions from `start`, as HTSlib's
/// `cram_generate_reference`: the most frequent of `A`, `C`, `G`, and `T` (the first on ties), or
/// `N` where no base is aligned.
fn consensus(records: &[Record], start: usize, len: usize) -> Vec<u8> {
    let mut counts = vec![[0u32; 4]; len];

    for record in records.iter().filter(|record| has_aligned_bases(record)) {
        let offset = record_offset(record, start);

        // Verbatim features have no errors to report.
        let _ = for_each_stretch(&record.features, |reference_offset, _, bases| {
            let first = offset + reference_offset;

            for (count, base) in counts.iter_mut().skip(first).zip(bases) {
                let i = match base.to_ascii_uppercase() {
                    b'A' => 0,
                    b'C' => 1,
                    b'G' => 2,
                    b'T' => 3,
                    _ => continue,
                };

                count[i] += 1;
            }

            Ok(())
        });
    }

    counts
        .iter()
        .map(|count| {
            let mut best = (0, b'N');

            for (&n, base) in count.iter().zip(*b"ACGT") {
                if n > best.0 {
                    best = (n, base);
                }
            }

            best.1
        })
        .collect()
}

/// Returns the record's features with its stretches of bases replaced by differences from
/// `reference_sequence`, where the record starts at `offset`.
fn reencode(record: &Record, reference_sequence: &[u8], offset: usize) -> io::Result<Vec<Feature>> {
    let mut stretches = Vec::new();

    for_each_stretch(
        &record.features,
        |reference_offset, read_position, bases| {
            stretches.push((offset + reference_offset, read_position, bases.len()));
            Ok(())
        },
    )?;

    let mut features = Vec::with_capacity(record.features.len());
    let mut stretches = stretches.into_iter();

    for feature in &record.features {
        if let Feature::Bases { .. } = feature {
            let (reference_offset, read_position, len) =
                stretches.next().expect("one stretch per bases feature");
            let read_offset = read_position - 1;

            push_mismatches(
                &mut features,
                reference_sequence,
                reference_offset,
                &record.sequence[read_offset..read_offset + len],
                &record.quality_scores,
                read_offset,
            )?;
        } else {
            features.push(feature.clone());
        }
    }

    Ok(features)
}

/// Adds a `cF` tag for the MD and NM tags the record lacks.
fn add_cram_flags(record: &mut Record) {
    let has = |tag| record.data.iter().any(|(t, _)| *t == tag);

    let flags =
        u8::from(!has(Tag::MISMATCHED_POSITIONS)) | (u8::from(!has(Tag::EDIT_DISTANCE)) << 1);

    if flags != 0 {
        record.data.push((CRAM_FLAGS, ValueBuf::UInt8(flags)));
    }
}
