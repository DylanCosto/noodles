//! Positions and regions to pile up: BED and position files (`-l`) and `-r` regions.

use std::{fs, path::Path};

use anyhow::{Context, bail};
use noodles_core::Region;

/// `-l` intervals by sequence name: 0-based, half-open, sorted by start (`bed_read`), with the
/// running maximum of their ends.
#[derive(Debug, Default)]
pub(super) struct Bed {
    pub(super) by_name: std::collections::HashMap<Vec<u8>, Intervals>,
}

/// Intervals sorted by start, and the running maximum of their ends.
pub(super) type Intervals = (Vec<(i64, i64)>, Vec<i64>);

impl Bed {
    pub(super) fn read(path: &Path) -> anyhow::Result<Self> {
        let text = fs::read(path)
            .with_context(|| format!("Could not read file \"{}\"", path.display()))?;
        let mut bed = Self::default();

        for (n, line) in text.split(|&b| b == b'\n').enumerate() {
            let line = line.trim_ascii_start();

            if line.is_empty() || line[0] == b'#' {
                continue;
            }

            let name_end = line
                .iter()
                .position(u8::is_ascii_whitespace)
                .unwrap_or(line.len());
            let name = &line[..name_end];
            let rest = String::from_utf8_lossy(&line[name_end..]);
            let mut fields = rest.split_ascii_whitespace().map(str::parse::<u64>);

            let (num, mut beg, mut end) = match (fields.next(), fields.next()) {
                (Some(Ok(beg)), Some(Ok(end))) => (2, beg, end),
                (Some(Ok(beg)), _) => (1, beg, 0),
                _ => (0, 0, 0),
            };

            if num == 1 {
                // A position, counting from 1.
                end = beg;
                beg = beg.wrapping_sub(1);
            }

            if num < 1 || end < beg {
                if name == b"browser" || name == b"track" {
                    continue;
                }

                bail!(
                    "[bed_read] Parse error reading \"{}\" at line {}",
                    path.display(),
                    n + 1
                );
            }

            bed.by_name
                .entry(name.to_vec())
                .or_default()
                .0
                .push((beg as i64, end as i64));
        }

        bed.index();
        Ok(bed)
    }

    /// Sorts the intervals and computes the running maximum of their ends.
    pub(super) fn index(&mut self) {
        for (intervals, max_ends) in self.by_name.values_mut() {
            intervals.sort_unstable();

            *max_ends = intervals
                .iter()
                .scan(i64::MIN, |max, &(_, end)| {
                    *max = (*max).max(end);
                    Some(*max)
                })
                .collect();
        }
    }

    /// `bed_overlap`: whether `[beg, end)` overlaps an interval on `name`.
    pub(super) fn overlaps(&self, name: &[u8], beg: i64, end: i64) -> bool {
        let Some((intervals, max_ends)) = self.by_name.get(name) else {
            return false;
        };

        // Of the intervals starting before `end`, whether any ends after `beg`.
        let n = intervals.partition_point(|&(b, _)| b < end);
        n > 0 && max_ends[n - 1] > beg
    }
}

/// Parses `-r` as `sam_itr_querys` does: 0-based start and exclusive end.
pub(super) fn parse_region(raw: &str, names: &[Vec<u8>]) -> anyhow::Result<(i32, i64, i64)> {
    let parsed: Region = raw
        .parse()
        .with_context(|| format!("[E::mpileup] fail to parse region '{raw}'"))?;
    let name: &[u8] = parsed.name().as_ref();
    let tid = names
        .iter()
        .position(|n| n.as_slice() == name)
        .with_context(|| format!("[E::mpileup] fail to parse region '{raw}'"))?;
    let interval = parsed.interval();
    let beg = interval.start().map_or(0, |p| usize::from(p) as i64 - 1);
    let end = interval.end().map_or(i64::MAX, |p| usize::from(p) as i64);
    Ok((tid as i32, beg, end))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_bed_overlaps() {
        let mut bed = Bed::default();
        bed.by_name.insert(
            b"chr1".to_vec(),
            (vec![(200, 300), (10, 20), (0, 100_000)], Vec::new()),
        );
        bed.index();

        assert!(bed.overlaps(b"chr1", 15, 16));
        assert!(bed.overlaps(b"chr1", 99_999, 100_000));
        assert!(!bed.overlaps(b"chr1", 100_000, 100_001));
        assert!(!bed.overlaps(b"chr2", 15, 16));
    }
}
