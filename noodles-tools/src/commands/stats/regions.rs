//! Target regions (`-t`) and the regions given on the command line.

use std::{fs, path::Path};

use anyhow::{Context, bail};

use super::{Info, Regions, Stats, record::Bam};
use crate::regions::RegionSet;

impl Stats {
    /// `init_regions`: target regions from a file of `chr from to` lines (1-based, closed).
    pub(super) fn init_regions(&mut self, path: &Path, info: &Info) -> anyhow::Result<()> {
        let text = fs::read(path).with_context(|| format!("{}", path.display()))?;
        let mut regions: Vec<Regions> = Vec::new();
        let mut warned = false;
        let mut prev_tid: Option<usize> = None;
        let mut prev_pos = -1;
        let mut nchunks = 0;

        // Lines as `kfgetline` reads them: no empty line after a final newline, and no `\r`.
        let text = text.strip_suffix(b"\n").unwrap_or(&text);

        for line in text.split(|&b| b == b'\n') {
            let line = line.strip_suffix(b"\r").unwrap_or(line);

            if line.first() == Some(&b'#') {
                continue;
            }

            let name_end = line
                .iter()
                .position(|&b| b.is_ascii_whitespace())
                .with_context(|| {
                    format!(
                        "Could not parse the file: {} [{}]",
                        path.display(),
                        String::from_utf8_lossy(line)
                    )
                })?;

            let name = &line[..name_end];

            let Some(tid) = info.reference_names.iter().position(|n| n == name) else {
                if !warned {
                    eprintln!(
                        "Warning: Some sequences not present in the BAM, e.g. \"{}\". This message is printed only once.",
                        String::from_utf8_lossy(name)
                    );
                }

                warned = true;
                continue;
            };

            let rest = String::from_utf8_lossy(&line[name_end + 1..]);
            let mut fields = rest.split_ascii_whitespace();
            let (Some(Ok(beg)), Some(Ok(end))) = (
                fields.next().map(str::parse::<i64>),
                fields.next().map(str::parse::<i64>),
            ) else {
                bail!("Could not parse the region [{rest}]");
            };

            if regions.len() <= tid {
                regions.resize_with(tid + 100, Regions::default);
            }

            if prev_tid != Some(tid) {
                prev_tid = Some(tid);
                prev_pos = beg;
            }

            if prev_pos > beg {
                bail!(
                    "The positions are not in chromosomal order ({}:{beg} comes after {prev_pos})",
                    String::from_utf8_lossy(name)
                );
            }

            regions[tid].intervals.push((beg, end));
            nchunks = nchunks.max(regions[tid].intervals.len());
        }

        if regions.is_empty() {
            bail!("Unable to map the -t sequences to the BAM sequences.");
        }

        // Sort the intervals and merge overlapping ones (adjacent ones stay apart).
        for (tid, reg) in regions.iter_mut().enumerate() {
            if reg.intervals.len() > 1 {
                reg.intervals
                    .sort_by(|a, b| a.0.cmp(&b.0).then(a.1.cmp(&b.1)));

                let mut merged: Vec<(i64, i64)> = Vec::with_capacity(reg.intervals.len());

                for &(beg, end) in &reg.intervals {
                    match merged.last_mut() {
                        Some(last) if last.1 >= beg => last.1 = last.1.max(end),
                        _ => merged.push((beg, end)),
                    }
                }

                reg.intervals = merged;
            }

            for &(beg, end) in &reg.intervals {
                if end < i64::MAX {
                    self.target_count = self.target_count.wrapping_add_signed(end - beg + 1);
                } else if let Some(&len) = info.reference_lengths.get(tid)
                    && len > 0
                {
                    self.target_count = self
                        .target_count
                        .wrapping_add(len.wrapping_sub(beg as u64).wrapping_add(1));
                }
            }
        }

        self.regions = Some(regions);
        self.chunks = vec![(0, 0); nchunks];
        self.nchunks = 0;

        Ok(())
    }

    /// `replicate_regions`: target regions from the command-line regions.
    pub(super) fn replicate_regions(&mut self, set: &RegionSet, info: &Info, region_count: usize) {
        let mut regions = Vec::new();

        for tid in 0..info.reference_names.len() {
            let intervals = set.intervals(tid);

            if intervals.is_empty() {
                continue;
            }

            if regions.len() <= tid {
                regions.resize_with(tid + 10, Regions::default);
            }

            // htslib ends a region without an end (or ending at the sequence end here) at
            // `HTS_POS_MAX`.
            let length = info.reference_lengths[tid];

            regions[tid].intervals = intervals
                .iter()
                .map(|&(beg, end)| {
                    let end = if end as u64 == length {
                        i64::MAX
                    } else {
                        end as i64
                    };
                    (beg as i64, end)
                })
                .collect();

            for &(beg, end) in &regions[tid].intervals {
                let end = if end == i64::MAX { length as i64 } else { end };
                self.target_count = self.target_count.wrapping_add_signed(end - beg + 1);
            }
        }

        self.regions = Some(regions);
        self.chunks = vec![(0, 0); region_count.max(1)];
    }

    /// `is_in_regions`
    pub(super) fn is_in_regions(&mut self, bam: &Bam<'_>) -> anyhow::Result<bool> {
        let Some(regions) = &mut self.regions else {
            return Ok(true);
        };

        let Ok(tid) = usize::try_from(bam.tid) else {
            return Ok(false);
        };

        if tid >= regions.len() {
            return Ok(false);
        }

        if !self.is_sorted {
            bail!("The BAM must be sorted in order for -t to work.");
        }

        let reg = &mut regions[tid];

        if reg.cpos == reg.intervals.len() {
            return Ok(false);
        }

        let mut i = reg.cpos;

        while i < reg.intervals.len() && reg.intervals[i].1 <= bam.pos {
            i += 1;
        }

        if i >= reg.intervals.len() {
            reg.cpos = reg.intervals.len();
            return Ok(false);
        }

        let end_pos = bam.end_pos();

        if end_pos < reg.intervals[i].0 {
            return Ok(false);
        }

        reg.cpos = i;
        self.reg_from = reg.intervals[i].0;
        self.reg_to = reg.intervals[i].1;

        self.nchunks = 0;

        for &(beg, end) in &reg.intervals[i..] {
            if bam.pos < end && end_pos >= beg {
                if self.nchunks >= self.chunks.len() {
                    self.chunks.push((0, 0));
                }

                self.chunks[self.nchunks] = ((bam.pos + 1).max(beg), end_pos.min(end));
                self.nchunks += 1;
            }
        }

        Ok(true)
    }
}
