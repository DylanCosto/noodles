//! Merged genomic intervals from BED files and region strings, for `-L` and `-M`.

use std::{fs, io, path::Path};

use anyhow::{Context, bail};
use noodles_core::{Position, Region};
use noodles_sam as sam;

/// Sorted, merged, 1-based closed intervals for each reference sequence (by header index).
#[derive(Clone, Debug, Default)]
pub struct RegionSet {
    by_reference: Vec<Vec<(usize, usize)>>,
}

impl RegionSet {
    pub fn new(header: &sam::Header) -> Self {
        Self {
            by_reference: vec![Vec::new(); header.reference_sequences().len()],
        }
    }

    /// Adds the intervals of a BED file (0-based, half-open). Names not in the header are
    /// ignored, as in samtools; comment, `track`, and `browser` lines are skipped.
    pub fn add_bed(&mut self, header: &sam::Header, path: &Path) -> anyhow::Result<()> {
        let text =
            fs::read_to_string(path).with_context(|| format!("reading {}", path.display()))?;

        for (i, line) in text.lines().enumerate() {
            let line = line.trim_end();

            if line.is_empty()
                || line.starts_with('#')
                || line.starts_with("track")
                || line.starts_with("browser")
            {
                continue;
            }

            let mut fields = line.split('\t');
            let (Some(name), Some(start), Some(end)) =
                (fields.next(), fields.next(), fields.next())
            else {
                bail!("{}:{}: expected at least 3 fields", path.display(), i + 1);
            };

            let start: usize = start
                .parse()
                .with_context(|| format!("{}:{}: invalid start", path.display(), i + 1))?;
            let end: usize = end
                .parse()
                .with_context(|| format!("{}:{}: invalid end", path.display(), i + 1))?;

            if let Some(id) = header.reference_sequences().get_index_of(name.as_bytes())
                && end > start
            {
                self.by_reference[id].push((start + 1, end));
            }
        }

        Ok(())
    }

    /// Adds a samtools region string (`chr1`, `chr1:100`, `chr1:100-200`).
    pub fn add_region(&mut self, header: &sam::Header, raw: &str) -> anyhow::Result<()> {
        let region: Region = raw
            .parse()
            .with_context(|| format!("invalid region: {raw}"))?;

        let Some((id, _, reference_sequence)) =
            header.reference_sequences().get_full(region.name())
        else {
            bail!("region {raw}: unknown reference sequence");
        };

        let interval = region.interval();
        let start = interval.start().map_or(1, usize::from);
        let end = interval
            .end()
            .map_or(reference_sequence.length().get(), usize::from);

        if end >= start {
            self.by_reference[id].push((start, end));
        }

        Ok(())
    }

    /// Sorts and merges overlapping or adjacent intervals.
    pub fn merge(&mut self) {
        for intervals in &mut self.by_reference {
            intervals.sort_unstable();

            let mut merged: Vec<(usize, usize)> = Vec::with_capacity(intervals.len());

            for &(start, end) in intervals.iter() {
                match merged.last_mut() {
                    Some(last) if start <= last.1.saturating_add(1) => last.1 = last.1.max(end),
                    _ => merged.push((start, end)),
                }
            }

            *intervals = merged;
        }
    }

    /// Merged intervals of one reference sequence.
    pub fn intervals(&self, reference_sequence_id: usize) -> &[(usize, usize)] {
        self.by_reference
            .get(reference_sequence_id)
            .map_or(&[], Vec::as_slice)
    }

    /// Returns the index of the first merged interval that `[start, end]` overlaps.
    pub fn first_overlap(
        &self,
        reference_sequence_id: usize,
        start: usize,
        end: usize,
    ) -> Option<usize> {
        let intervals = self.intervals(reference_sequence_id);
        let i = intervals.partition_point(|&(_, e)| e < start);
        intervals.get(i).filter(|&&(s, _)| s <= end).map(|_| i)
    }
}

/// Returns a record's reference sequence ID and the 1-based closed interval samtools uses for
/// overlap tests: the reference span of a mapped record (at least one base), or one base at the
/// position of an unmapped one. Returns `None` for records with no reference or position.
pub fn record_interval(
    record: &dyn sam::alignment::Record,
    header: &sam::Header,
) -> io::Result<Option<(usize, usize, usize)>> {
    let (Some(id), Some(start)) = (
        record.reference_sequence_id(header).transpose()?,
        record.alignment_start().transpose()?,
    ) else {
        return Ok(None);
    };

    let span = if record.flags()?.is_unmapped() {
        0
    } else {
        record.cigar().alignment_span()?
    };

    let start = usize::from(start);
    let end = start + span.max(1) - 1;

    Ok(Some((id, start, end)))
}

/// Formats a 1-based closed interval as a region query string.
pub fn region_string(name: &[u8], start: usize, end: usize) -> anyhow::Result<Region> {
    let name = std::str::from_utf8(name).context("reference sequence name is not UTF-8")?;
    let start = Position::try_from(start)?;
    let end = Position::try_from(end)?;
    Ok(Region::new(name, start..=end))
}

#[cfg(test)]
mod tests {
    use super::*;

    fn header() -> sam::Header {
        "@SQ\tSN:sq0\tLN:1000\n@SQ\tSN:sq1\tLN:500\n"
            .parse()
            .expect("valid header")
    }

    #[test]
    fn test_merge_and_first_overlap() -> anyhow::Result<()> {
        let header = header();
        let mut set = RegionSet::new(&header);

        set.add_region(&header, "sq0:100-200")?;
        set.add_region(&header, "sq0:150-250")?;
        set.add_region(&header, "sq0:251-260")?;
        set.add_region(&header, "sq0:400-500")?;
        set.add_region(&header, "sq1")?;
        set.merge();

        assert_eq!(set.intervals(0), [(100, 260), (400, 500)]);
        assert_eq!(set.intervals(1), [(1, 500)]);

        assert_eq!(set.first_overlap(0, 1, 99), None);
        assert_eq!(set.first_overlap(0, 90, 100), Some(0));
        assert_eq!(set.first_overlap(0, 250, 450), Some(0));
        assert_eq!(set.first_overlap(0, 261, 399), None);
        assert_eq!(set.first_overlap(0, 300, 400), Some(1));
        assert_eq!(set.first_overlap(1, 500, 600), Some(0));

        Ok(())
    }
}
