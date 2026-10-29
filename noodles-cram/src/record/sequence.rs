mod iter;

use noodles_core::Position;
use noodles_sam as sam;

use self::iter::{Iter, with_positions::WithPositions};
use super::Feature;
use crate::container::compression_header::preservation_map::SubstitutionMatrix;

pub(super) struct Sequence<'r, 'c: 'r> {
    reference_sequence: Option<&'c [u8]>,
    substitution_matrix: SubstitutionMatrix,
    features: &'r [Feature<'c>],
    alignment_start: Position,
    read_length: usize,
}

impl<'r, 'c: 'r> Sequence<'r, 'c> {
    pub(super) fn new(
        reference_sequence: Option<&'c [u8]>,
        substitution_matrix: SubstitutionMatrix,
        features: &'r [Feature<'c>],
        alignment_start: Position,
        read_length: usize,
    ) -> Self {
        Self {
            reference_sequence,
            substitution_matrix,
            features,
            alignment_start,
            read_length,
        }
    }

    pub(super) fn to_vec(&self) -> Vec<u8> {
        let mut dst = Vec::with_capacity(self.read_length);
        self.extend_to(&mut dst);
        dst
    }

    pub(super) fn extend_to(&self, dst: &mut Vec<u8>) {
        use crate::container::compression_header::preservation_map::substitution_matrix::Base;

        dst.reserve(self.read_length);
        let mut features = WithPositions::new(self.features, self.alignment_start);
        let (mut previous_reference, mut previous_read) = features.positions();

        while let Some(((reference_position, read_position), feature)) = features.next() {
            if let Some(reference) = self.reference_sequence {
                let start = usize::from(previous_reference) - 1;
                let end = usize::from(reference_position) - 1;
                extend_reference(dst, reference, start, end - start);
            } else if read_position != previous_read {
                panic!("next: missing reference sequence");
            }

            match feature {
                Feature::Bases { bases, .. }
                | Feature::Insertion { bases, .. }
                | Feature::SoftClip { bases, .. } => dst.extend_from_slice(bases),
                Feature::ReadBase { base, .. } | Feature::InsertBase { base, .. } => {
                    dst.push(*base)
                }
                Feature::Substitution { code, .. } => {
                    let reference = self
                        .reference_sequence
                        .expect("missing reference sequence (substitution)");
                    // Past the end of the reference, HTSlib substitutes against N.
                    let raw_reference_base = reference
                        .get(usize::from(reference_position) - 1)
                        .copied()
                        .unwrap_or(b'N');
                    let reference_base = Base::try_from(raw_reference_base).unwrap_or(Base::N);
                    let mut base = u8::from(self.substitution_matrix.get(reference_base, *code));
                    if raw_reference_base.is_ascii_lowercase() {
                        base = base.to_ascii_lowercase();
                    }
                    dst.push(base);
                }
                _ => {}
            }

            (previous_reference, previous_read) = features.positions();
        }

        let last_read_position = usize::from(previous_read);
        if let Some(reference) = self.reference_sequence {
            if last_read_position <= self.read_length {
                let len = self.read_length - last_read_position + 1;
                let start = usize::from(previous_reference) - 1;
                extend_reference(dst, reference, start, len);
            }
        } else if last_read_position != self.read_length + 1 {
            panic!(
                "next: missing reference sequence: {} != {}",
                last_read_position, self.read_length
            );
        }
    }
}

/// Appends `len` reference bases from `start`, with `N`s for any past the end of the reference, as
/// HTSlib does for alignments that overhang it.
fn extend_reference(dst: &mut Vec<u8>, reference: &[u8], start: usize, len: usize) {
    let end = (start + len).min(reference.len());
    let available = reference.get(start..end).unwrap_or_default();
    dst.extend_from_slice(available);
    dst.resize(dst.len() + len - available.len(), b'N');
}

impl sam::alignment::record::Sequence for Sequence<'_, '_> {
    fn is_empty(&self) -> bool {
        self.len() == 0
    }

    fn len(&self) -> usize {
        self.read_length
    }

    fn get(&self, i: usize) -> Option<u8> {
        self.iter().nth(i)
    }

    fn iter(&self) -> Box<dyn Iterator<Item = u8> + '_> {
        Box::new(Iter::new(
            self.reference_sequence,
            self.substitution_matrix.clone(),
            self.features,
            self.alignment_start,
            self.read_length,
        ))
    }
}

#[cfg(test)]
mod tests {
    use std::borrow::Cow;

    use sam::alignment::record::Sequence as _;

    use super::*;

    #[test]
    fn test_to_vec_matches_iter() {
        let p = |n| Position::new(n).unwrap();
        let bases = |s: &'static [u8]| Cow::Borrowed(s);
        #[rustfmt::skip]
        let mut features = [
            Feature::SoftClip { position: p(1), bases: bases(b"Cg") },
            Feature::Bases { position: p(4), bases: bases(b"tG") },
            Feature::Substitution { position: p(7), code: 0 },
            Feature::InsertBase { position: p(8), base: b'Y' },
            Feature::Insertion { position: p(9), bases: bases(b"AN") },
            Feature::Deletion { position: p(11), len: 2 },
            Feature::ReadBase { position: p(12), base: b'N', quality_score: 3 },
            Feature::ReferenceSkip { position: p(13), len: 3 },
        ];
        let matrix = SubstitutionMatrix::default();
        for &base in b"ACGTNacgtnRY" {
            for code in 0..4 {
                features[2] = Feature::Substitution {
                    position: p(7),
                    code,
                };
                let reference = [base; 64];
                let reference = Some(&reference[..]);
                let sequence = Sequence::new(reference, matrix.clone(), &features, p(3), 20);
                assert_eq!(sequence.to_vec(), sequence.iter().collect::<Vec<_>>());
            }
        }
        // An empty read.
        let sequence = Sequence::new(None, matrix, &[], p(1), 0);
        assert!(sequence.to_vec().is_empty());
    }
}
