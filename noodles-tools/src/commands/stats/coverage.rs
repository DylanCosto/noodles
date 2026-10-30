//! Coverage and GC-depth: the round buffer of per-position depths and the reference sequence it is compared with.

use anyhow::bail;

use super::{GcDepth, Info, ReferenceStats, Regions, Stats};

/// `lidx2ridx`
pub(super) fn lidx2ridx(offset: i32, size: i32, refpos: i64, pos: i64) -> i32 {
    let size64 = i64::from(size);
    i32::try_from((i64::from(offset) + (pos - refpos) % size64) % size64).unwrap_or(0)
}

impl Stats {
    /// `coverage_idx`
    pub(super) fn coverage_index(&self, info: &Info, depth: i32) -> usize {
        if depth < info.cov_min {
            0
        } else if depth > info.cov_max {
            self.ncov - 1
        } else {
            usize::try_from(1 + (depth - info.cov_min) / info.cov_step).expect("nonnegative")
        }
    }

    /// `round_buffer_flush`
    pub(super) fn round_buffer_flush(&mut self, pos: i64, info: &Info) -> anyhow::Result<()> {
        if pos == self.cov_rbuf.pos {
            return Ok(());
        }

        let new_pos = pos;
        let size = i64::from(self.cov_rbuf.size);
        let mut pos = pos;

        if pos == -1 || pos - self.cov_rbuf.pos >= size {
            pos = self.cov_rbuf.pos + size - 1;
        }

        if pos < self.cov_rbuf.pos {
            bail!(
                "Expected coordinates in ascending order, got {pos} after {}",
                self.cov_rbuf.pos
            );
        }

        let mut ifrom = self.cov_rbuf.start;
        let ito = lidx2ridx(
            self.cov_rbuf.start,
            self.cov_rbuf.size,
            self.cov_rbuf.pos,
            pos - 1,
        );

        if ifrom > ito {
            for ibuf in ifrom..self.cov_rbuf.size {
                self.flush_depth(ibuf as usize, info);
            }

            ifrom = 0;
        }

        for ibuf in ifrom..=ito {
            self.flush_depth(ibuf as usize, info);
        }

        self.cov_rbuf.start = if new_pos == -1 {
            0
        } else {
            lidx2ridx(
                self.cov_rbuf.start,
                self.cov_rbuf.size,
                self.cov_rbuf.pos,
                pos,
            )
        };
        self.cov_rbuf.pos = new_pos;

        Ok(())
    }

    pub(super) fn flush_depth(&mut self, ibuf: usize, info: &Info) {
        let depth = self.cov_rbuf.buffer[ibuf];

        if depth != 0 {
            let i = self.coverage_index(info, depth);
            self.cov[i] += 1;
            self.cov_rbuf.buffer[ibuf] = 0;
        }
    }

    /// `round_buffer_insert_read`: adds the 0-based half-open `[from, to)`.
    pub(super) fn insert_read(&mut self, from: i64, to: i64) -> anyhow::Result<()> {
        let rbuf = &mut self.cov_rbuf;

        if to - from > i64::from(rbuf.size) {
            bail!(
                "The read length too big ({}), please increase the buffer length (currently {})",
                to - from,
                rbuf.size
            );
        }

        if from < rbuf.pos {
            bail!(
                "The reads are not sorted ({from} comes after {}).",
                rbuf.pos
            );
        }

        let mut ifrom = lidx2ridx(rbuf.start, rbuf.size, rbuf.pos, from);
        let ito = lidx2ridx(rbuf.start, rbuf.size, rbuf.pos, to);

        if ifrom > ito {
            for d in &mut rbuf.buffer[ifrom as usize..rbuf.size as usize] {
                *d += 1;
            }

            ifrom = 0;
        }

        for d in &mut rbuf.buffer[ifrom as usize..ito as usize] {
            *d += 1;
        }

        Ok(())
    }

    /// `read_ref_seq`
    pub(super) fn read_ref_seq(
        &mut self,
        info: &Info,
        tid: i32,
        pos: i64,
        end: i64,
    ) -> anyhow::Result<()> {
        let reference = info.reference.as_ref().expect("reference is loaded");

        if end >= pos + self.mrseq_buf - 1 && self.mrseq_buf < end - pos {
            self.mrseq_buf = end - pos;
            self.rseq_buf.resize(usize::try_from(self.mrseq_buf)?, 0);
        }

        let name = &info.reference_names[usize::try_from(tid)?];
        let mut sequence = Vec::new();

        if !reference.fetch(name, pos, pos + self.mrseq_buf - 1, &mut sequence)? {
            bail!(
                "Failed to fetch the sequence \"{}\"",
                String::from_utf8_lossy(name)
            );
        }

        for (dst, &b) in self.rseq_buf.iter_mut().zip(&sequence) {
            *dst = match b {
                b'A' | b'a' => 1,
                b'C' | b'c' => 2,
                b'G' | b'g' => 4,
                b'T' | b't' => 8,
                _ => 0,
            };
        }

        self.rseq_buf[sequence.len()..].fill(0);
        self.nrseq_buf = i64::try_from(sequence.len())?;
        self.rseq_pos = pos;
        self.tid = tid;

        Ok(())
    }

    /// `fai_gc_content`
    pub(super) fn gc_content(&self, pos: i64, len: i32) -> f32 {
        let i = usize::try_from(pos - self.rseq_pos).expect("in the buffer");
        let end = (pos - self.rseq_pos + i64::from(len)).min(self.nrseq_buf);
        let end = usize::try_from(end).unwrap_or(0);

        let (mut gc, mut count) = (0u32, 0u32);

        for &c in self.rseq_buf.get(i..end).unwrap_or(&[]) {
            match c {
                2 | 4 => {
                    gc += 1;
                    count += 1;
                }
                1 | 8 => count += 1,
                _ => {}
            }
        }

        if count > 0 {
            gc as f32 / count as f32
        } else {
            0.0
        }
    }

    /// `realloc_gcd_buffer`: makes `gcd[igcd]` valid, growing as `hts_expand0` does.
    pub(super) fn ensure_gcd(&mut self, info: &Info) {
        if self.igcd >= self.gcd.len() {
            let n = (self.igcd + 1).next_power_of_two();
            self.gcd.resize(n, GcDepth::default());
            self.realloc_rseq_buffer(info);
        }
    }

    /// `collect_refstats`
    pub(super) fn collect_reference_stats(&mut self, info: &Info) -> anyhow::Result<()> {
        let total_count = info.reference_names.len();
        let all_targets = self.regions.as_ref().is_none_or(Vec::is_empty);

        let regions: Vec<Regions> = if all_targets {
            (0..total_count)
                .map(|_| Regions {
                    intervals: vec![(1, i64::MAX)],
                    cpos: 0,
                })
                .collect()
        } else {
            self.regions.clone().unwrap_or_default()
        };

        let mut stats = ReferenceStats {
            total_count,
            count: 0,
            combined_length: 0,
            min_length: 0,
            max_length: 0,
            average_length: 0.0,
            average_gc: 0.0,
            rows: Vec::new(),
        };

        let mut gc_sum: f32 = 0.0;
        let mut warned = false;
        let mut buf = Vec::new();

        for (i, reg) in regions.iter().enumerate() {
            if reg.intervals.is_empty() {
                continue;
            }

            let Some(name) = info.reference_names.get(i) else {
                continue;
            };

            let header_length = i64::try_from(info.reference_lengths[i])?;

            for &(beg, original_end) in &reg.intervals {
                let start = beg;
                let end = if original_end == i64::MAX {
                    header_length
                } else {
                    original_end
                };

                if end < start {
                    stats.rows.push(None);
                    continue;
                }

                stats.count += 1;

                let (row_name, length) = if all_targets {
                    (String::from_utf8_lossy(name).into_owned(), header_length)
                } else {
                    let name = String::from_utf8_lossy(name);

                    // A whole sequence is named alone, a range as `name:start-end`.
                    let row_name = if beg == 1 && original_end == i64::MAX {
                        name.into_owned()
                    } else {
                        format!("{name}:{start}-{end}")
                    };

                    (row_name, (end - start + 1).min(header_length))
                };

                stats.combined_length += length;

                if stats.min_length == 0 || stats.min_length > length {
                    stats.min_length = length;
                }

                stats.max_length = stats.max_length.max(length);

                let Some(reference) = &info.reference else {
                    gc_sum = -1.0;
                    stats.rows.push(Some((row_name, length, -1.0, -1)));
                    continue;
                };

                let (mut gc, mut at, mut n) = (0i64, 0i64, 0i64);
                let mut start0 = start - 1;
                let mut rem = end - start0;

                while rem > 0 {
                    let chunk = rem.min(info.ref_chunk_size);

                    if !reference.fetch(name, start0, start0 + chunk - 1, &mut buf)? {
                        if !warned {
                            eprintln!(
                                "Warning: Some sequences not present in the reference, e.g. \"{}\".                             This message is printed only once.",
                                String::from_utf8_lossy(name)
                            );
                            warned = true;
                        }

                        break;
                    }

                    for &b in &buf {
                        match b {
                            b'G' | b'g' | b'C' | b'c' => gc += 1,
                            b'A' | b'a' | b'T' | b't' => at += 1,
                            b'N' | b'n' => n += 1,
                            _ => {}
                        }
                    }

                    let got = i64::try_from(buf.len())?;

                    if got < chunk {
                        break;
                    }

                    start0 += got;
                    rem -= got;
                }

                let gc_content = if gc + at > 0 {
                    gc as f32 / (gc + at) as f32
                } else {
                    0.0
                };

                gc_sum += gc_content;
                stats.rows.push(Some((row_name, length, gc_content, n)));
            }
        }

        if stats.count > 0 {
            stats.average_length = stats.combined_length as f32 / stats.count as f32;
            stats.average_gc = if gc_sum < 0.0 {
                -1.0
            } else {
                gc_sum / stats.count as f32
            };
        } else {
            stats.average_length = -1.0;
            stats.average_gc = -1.0;
        }

        self.reference_stats = Some(stats);

        Ok(())
    }
}
