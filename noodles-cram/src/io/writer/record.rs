mod convert;
mod feature;

use bstr::BString;
use noodles_core::Position;
use noodles_sam::{
    self as sam,
    alignment::{
        record::{MappingQuality, data::field::Tag},
        record_buf::data::field::Value,
    },
};

pub(crate) use self::convert::push_mismatches;
pub use self::feature::Feature;
use crate::record::{Flags, MateFlags};

/// How a record's aligned bases are stored.
#[derive(Clone, Copy, Debug, Default, Eq, PartialEq)]
pub(crate) enum ReferenceBases {
    /// The record is not placed, so no reference is used.
    #[default]
    Unused,
    /// Differences from the external reference sequence.
    External,
    /// Verbatim, because the external reference lacks the record's sequence (or another record's
    /// in the same container).
    Verbatim,
}

#[derive(Clone, Debug, Default, PartialEq)]
pub struct Record {
    pub(crate) bam_flags: sam::alignment::record::Flags,
    pub(crate) cram_flags: Flags,
    pub(crate) reference_sequence_id: Option<usize>,
    pub(crate) read_length: usize,
    pub(crate) alignment_start: Option<Position>,
    pub(crate) read_group_id: Option<usize>,
    pub(crate) name: Option<BString>,
    pub(crate) mate_flags: MateFlags,
    pub(crate) mate_reference_sequence_id: Option<usize>,
    pub(crate) mate_alignment_start: Option<Position>,
    pub(crate) template_length: i32,
    pub(crate) mate_distance: Option<usize>,
    pub(crate) data: Vec<(Tag, Value)>,
    pub(crate) features: Vec<Feature>,
    pub(crate) reference_bases: ReferenceBases,
    pub(crate) mapping_quality: Option<MappingQuality>,
    pub(crate) sequence: Vec<u8>,
    pub(crate) quality_scores: Vec<u8>,
}

impl Record {
    pub fn alignment_span(&self) -> usize {
        calculate_alignment_span(self.read_length, &self.features)
    }

    pub fn alignment_end(&self) -> Option<Position> {
        self.alignment_start.and_then(|start| {
            let end = usize::from(start) + self.alignment_span() - 1;
            Position::new(end)
        })
    }

    /// Stores the aligned bases verbatim instead of as differences from the external reference.
    pub(crate) fn store_bases_verbatim(&mut self) {
        if self.reference_bases != ReferenceBases::External {
            return;
        }

        self.reference_bases = ReferenceBases::Verbatim;

        // A missing sequence (`*`) has no aligned bases to store.
        if !self.sequence.is_empty() {
            self.features = verbatim_features(&self.features, &self.sequence);
        }
    }
}

/// Returns `features` with substitutions and read bases replaced by stretches of `bases` for the
/// aligned read positions, i.e., those not in insertions or soft clips.
fn verbatim_features(features: &[Feature], bases: &[u8]) -> Vec<Feature> {
    // Pushes the aligned bases from read position `next` (1-based) up to `end`, exclusive.
    fn push_bases(dst: &mut Vec<Feature>, bases: &[u8], next: &mut usize, end: usize) {
        if end > *next
            && let Some(position) = Position::new(*next)
        {
            dst.push(Feature::Bases {
                position,
                bases: bases[*next - 1..end - 1].to_vec(),
            });

            *next = end;
        }
    }

    let mut verbatim = Vec::with_capacity(features.len());
    let mut next = 1;

    for feature in features {
        let position = usize::from(feature.position());

        match feature {
            Feature::Bases { .. } | Feature::ReadBase { .. } | Feature::Substitution { .. } => {
                continue;
            }
            Feature::Insertion {
                bases: inserted, ..
            }
            | Feature::SoftClip {
                bases: inserted, ..
            } => {
                push_bases(&mut verbatim, bases, &mut next, position);
                next = position + inserted.len();
            }
            Feature::InsertBase { .. } => {
                push_bases(&mut verbatim, bases, &mut next, position);
                next = position + 1;
            }
            _ => push_bases(&mut verbatim, bases, &mut next, position),
        }

        verbatim.push(feature.clone());
    }

    push_bases(&mut verbatim, bases, &mut next, bases.len() + 1);

    verbatim
}

fn calculate_alignment_span(read_length: usize, features: &[Feature]) -> usize {
    features
        .iter()
        .fold(read_length, |span, feature| match feature {
            Feature::Insertion { bases, .. } => span - bases.len(),
            Feature::InsertBase { .. } => span - 1,
            Feature::Deletion { len, .. } => span + len,
            Feature::ReferenceSkip { len, .. } => span + len,
            Feature::SoftClip { bases, .. } => span - bases.len(),
            _ => span,
        })
}
