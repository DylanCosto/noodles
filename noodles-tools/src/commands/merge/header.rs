//! Merging input headers as `samtools merge` does (`trans_tbl_init`).
//!
//! The merged header takes `@HD` from the first header that has one, every new `@SQ` line in
//! input order, the `@RG` and `@PG` lines of each input (renaming IDs that collide with earlier
//! ones unless they are combined), and every `@CO` line. Each input gets a [`Translation`] from
//! its reference sequence, read group, and program IDs to the merged ones.

use std::collections::HashMap;

use bstr::{BStr, BString, ByteSlice};
use indexmap::IndexMap;
use noodles_sam::{
    self as sam,
    header::record::value::{
        Map,
        map::{
            header::tag::{GROUP_ORDER, SORT_ORDER, SUBSORT_ORDER},
            program::tag::PREVIOUS_PROGRAM_ID,
            read_group::tag::PROGRAM,
        },
    },
};

/// What to do with an `@RG` or `@PG` ID that an earlier header already has.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub enum Collisions {
    /// Give the later line a new ID with a random suffix (the default).
    Rename,
    /// Keep only the earlier line, and map the later ID to it (`-c`, `-p`).
    Combine,
}

/// The translation of one input's IDs to the merged header's.
#[derive(Debug, Default)]
pub struct Translation {
    /// Merged reference sequence IDs, by input ID.
    pub reference_sequence_ids: Vec<i32>,
    /// Whether the merged header orders this input's reference sequences differently.
    pub lost_coordinate_order: bool,
    /// Merged read group IDs, by input ID.
    pub read_groups: HashMap<BString, BString>,
    /// Merged program IDs, by input ID.
    pub programs: HashMap<BString, BString>,
}

/// Builds a merged header from input headers.
pub struct HeaderMerger {
    header: sam::Header,
    read_groups: Collisions,
    programs: Collisions,
    rng: Lrand48,
}

impl HeaderMerger {
    pub fn new(read_groups: Collisions, programs: Collisions, seed: i64) -> Self {
        Self {
            header: sam::Header::default(),
            read_groups,
            programs,
            rng: Lrand48::new(seed),
        }
    }

    /// Adds an input header. With `read_group`, every read group of the input gets that ID
    /// (`-r`), and an input without read groups gets an `@RG` line for it.
    pub fn add(&mut self, input: &sam::Header, read_group: Option<&BStr>) -> Translation {
        if self.header.header().is_none() {
            *self.header.header_mut() = input.header().cloned();
        }

        let (reference_sequence_ids, lost_coordinate_order) = self.add_reference_sequences(input);

        let (mut read_groups, new_read_groups) = translate_ids(
            self.header.read_groups_mut(),
            input.read_groups(),
            self.read_groups,
            read_group,
            &mut self.rng,
        );

        if let Some(id) = read_group
            && new_read_groups.is_empty()
        {
            let known = self.header.read_groups();
            let id = unique_id(
                id,
                |id| known.contains_key(id),
                Suffix::IfTaken,
                &mut self.rng,
            );
            self.header
                .read_groups_mut()
                .insert(id.clone(), Map::default());
            read_groups.insert(BString::default(), id);
        }

        let (programs, new_programs) = translate_ids(
            self.header.programs_mut().as_mut(),
            input.programs().as_ref(),
            self.programs,
            None,
            &mut self.rng,
        );

        for id in &new_read_groups {
            if let Some(read_group) = self.header.read_groups_mut().get_mut(id) {
                relink(
                    read_group.other_fields_mut().get_mut(&PROGRAM),
                    &programs,
                    "PG",
                );
            }
        }

        for id in &new_programs {
            if let Some(program) = self.header.programs_mut().as_mut().get_mut(id) {
                relink(
                    program.other_fields_mut().get_mut(&PREVIOUS_PROGRAM_ID),
                    &programs,
                    "PP",
                );
            }
        }

        self.header
            .comments_mut()
            .extend(input.comments().iter().cloned());

        Translation {
            reference_sequence_ids,
            lost_coordinate_order,
            read_groups,
            programs,
        }
    }

    /// Adds the input's new reference sequences, returning the merged ID of each input
    /// reference sequence and whether they are out of order.
    fn add_reference_sequences(&mut self, input: &sam::Header) -> (Vec<i32>, bool) {
        let merged = self.header.reference_sequences_mut();
        let mut ids = Vec::with_capacity(input.reference_sequences().len());
        let mut lost_order = false;

        for (name, reference_sequence) in input.reference_sequences() {
            let id = match merged.get_full(name) {
                Some((id, _, existing)) => {
                    // samtools matches reference sequences by name alone.
                    if existing.length() != reference_sequence.length() {
                        eprintln!(
                            "[W::merge] reference sequence {name} has length {} in an earlier \
                             header but {} here",
                            existing.length(),
                            reference_sequence.length()
                        );
                    }

                    id
                }
                None => {
                    merged.insert(name.clone(), reference_sequence.clone());
                    merged.len() - 1
                }
            };

            let id = i32::try_from(id).unwrap_or(i32::MAX);

            if ids.last().is_some_and(|&last| id <= last) {
                lost_order = true;
            }

            ids.push(id);
        }

        (ids, lost_order)
    }

    /// Returns the merged header.
    pub fn finish(self) -> sam::Header {
        self.header
    }
}

/// Marks `@HD` unsorted, without `SS` or `GO`, as samtools does when an input's reference
/// sequences are out of order.
pub fn mark_unsorted(header: &mut sam::Header) {
    if let Some(hd) = header.header_mut() {
        let fields = hd.other_fields_mut();
        fields.shift_remove(&SUBSORT_ORDER);
        fields.shift_remove(&GROUP_ORDER);
        fields.insert(SORT_ORDER, "unsorted".into());
    }
}

/// Adds an input's `@RG` or `@PG` lines (`trans_rg_pg`), returning the map from input IDs to
/// merged IDs and the merged IDs of the lines added.
fn translate_ids<V: Clone>(
    merged: &mut IndexMap<BString, V>,
    input: &IndexMap<BString, V>,
    collisions: Collisions,
    id_override: Option<&BStr>,
    rng: &mut Lrand48,
) -> (HashMap<BString, BString>, Vec<BString>) {
    let mut map = HashMap::with_capacity(input.len());
    let mut added = Vec::new();

    for (id, line) in input {
        let (new_id, is_new) = match id_override {
            Some(new_id) => (BString::from(new_id), !merged.contains_key(new_id)),
            None if collisions == Collisions::Combine || !merged.contains_key(id) => {
                (id.clone(), !merged.contains_key(id))
            }
            None => {
                let new_id = unique_id(
                    id.as_ref(),
                    |id| merged.contains_key(id),
                    Suffix::Always,
                    rng,
                );
                (new_id, true)
            }
        };

        if is_new {
            merged.insert(new_id.clone(), line.clone());
            added.push(new_id.clone());
        }

        map.insert(id.clone(), new_id);
    }

    (map, added)
}

/// Replaces a `PG` or `PP` value with its merged ID, warning when it names no program of its
/// input (`finish_rg_pg`).
fn relink(value: Option<&mut BString>, programs: &HashMap<BString, BString>, tag: &str) {
    let Some(value) = value else {
        return;
    };

    match programs.get(value) {
        Some(id) => value.clone_from(id),
        None => eprintln!("[W::finish_rg_pg] Tag {tag}:{value} not found in @PG records"),
    }
}

/// When [`unique_id`] adds a random suffix.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
enum Suffix {
    Always,
    IfTaken,
}

/// samtools' `gen_unique_id`: `prefix`, or `prefix-XXXXXXXX` with a random hexadecimal suffix
/// that `is_taken` rejects.
fn unique_id<F>(prefix: &BStr, is_taken: F, suffix: Suffix, rng: &mut Lrand48) -> BString
where
    F: Fn(&BStr) -> bool,
{
    if suffix == Suffix::IfTaken && !is_taken(prefix) {
        return prefix.into();
    }

    loop {
        let id = BString::from(format!("{prefix}-{:08X}", rng.next()));

        if !is_taken(id.as_bstr()) {
            return id;
        }
    }
}

/// glibc's `lrand48` (the 48-bit linear congruential generator samtools seeds with `-s`).
struct Lrand48(u64);

impl Lrand48 {
    const MASK: u64 = (1 << 48) - 1;

    /// `srand48`: the seed's low 32 bits become the state's high bits.
    fn new(seed: i64) -> Self {
        Self(((seed.cast_unsigned() & 0xffff_ffff) << 16) | 0x330e)
    }

    fn next(&mut self) -> u32 {
        self.0 = self.0.wrapping_mul(0x0005_deec_e66d).wrapping_add(0xb) & Self::MASK;
        u32::try_from(self.0 >> 17).unwrap_or(u32::MAX)
    }
}

#[cfg(test)]
mod tests {
    use std::collections::HashSet;

    use super::*;

    #[test]
    fn test_lrand48() {
        // samtools merge -s 42 names a clashing ID "g1-5F4C985F", then "p1-2BDDA45D".
        let mut rng = Lrand48::new(42);
        assert_eq!(rng.next(), 0x5F4C_985F);
        assert_eq!(rng.next(), 0x2BDD_A45D);
    }

    #[test]
    fn test_unique_id() {
        let mut rng = Lrand48::new(42);
        let taken: HashSet<BString> = [BString::from("g1")].into_iter().collect();
        let is_taken = |id: &BStr| taken.contains(id);

        assert_eq!(
            unique_id(b"g2".as_bstr(), is_taken, Suffix::IfTaken, &mut rng),
            "g2"
        );
        assert_eq!(
            unique_id(b"g1".as_bstr(), is_taken, Suffix::IfTaken, &mut rng),
            "g1-5F4C985F"
        );
        assert_eq!(
            unique_id(b"g2".as_bstr(), is_taken, Suffix::Always, &mut rng),
            "g2-2BDDA45D"
        );
    }
}
