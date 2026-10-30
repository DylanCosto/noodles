//! Per-cycle base and quality counts, gathered in batches and counted on worker threads.

use crate::sequence::Strand;

use super::{COMPLEMENT, CodeCounts, Stats};

/// `acgtno_count_t`
#[derive(Clone, Copy, Debug, Default)]
pub(super) struct BaseCounts {
    pub(super) a: u64,
    pub(super) c: u64,
    pub(super) g: u64,
    pub(super) t: u64,
    pub(super) n: u64,
    pub(super) other: u64,
}

impl From<&CodeCounts> for BaseCounts {
    fn from(counts: &CodeCounts) -> Self {
        let (a, c, g, t, n) = (counts[1], counts[2], counts[4], counts[8], counts[15]);
        let other = counts.iter().sum::<u64>() - a - c - g - t - n;
        Self {
            a,
            c,
            g,
            t,
            n,
            other,
        }
    }
}

/// A read whose per-cycle base and quality counts are pending, with its sequence and qualities in
/// [`PendingCycles::data`].
pub(super) struct PendingRead {
    pub(super) start: usize,
    pub(super) len: usize,
    pub(super) strand: Strand,
    pub(super) segment: Segment,
}

/// Which per-cycle tables a read is counted in.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub(super) enum Segment {
    /// The first segment of a template.
    First,
    /// Any other read.
    Other,
}

/// Reads whose per-cycle counts are deferred so that worker threads can count them together,
/// each over its own range of cycles. The counts are sums, so this gives the same tables as
/// counting each read as it is read.
#[derive(Default)]
pub(super) struct PendingCycles {
    pub(super) data: Vec<u8>,
    pub(super) reads: Vec<PendingRead>,
    pub(super) max_len: usize,
}

/// The per-cycle tables for cycles `lo..lo + revcomp.len()`.
pub(super) struct CycleTables<'a> {
    pub(super) lo: usize,
    pub(super) first: &'a mut [CodeCounts],
    pub(super) last: &'a mut [CodeCounts],
    pub(super) revcomp: &'a mut [CodeCounts],
    pub(super) quals_first: &'a mut [u64],
    pub(super) quals_last: &'a mut [u64],
    pub(super) rows: usize,
}

/// Counts the bases and qualities of one first or last read at the cycles `tables` covers
/// (`collect_orig_read_stats`). Cycle `c` holds base `c` of the read as sequenced: from the end of
/// the stored sequence for a reverse-strand read.
pub(super) fn count_cycles(
    seq: &[u8],
    qual: &[u8],
    strand: Strand,
    segment: Segment,
    tables: &mut CycleTables<'_>,
) {
    let len = qual.len();
    let lo = tables.lo;
    let hi = (lo + tables.revcomp.len()).min(len);

    if lo >= hi {
        return;
    }

    let (cycles, quals) = match segment {
        Segment::First => (&mut *tables.first, &mut *tables.quals_first),
        Segment::Other => (&mut *tables.last, &mut *tables.quals_last),
    };

    let reverse = strand == Strand::Reverse;
    let complement: [u8; 16] = if reverse {
        COMPLEMENT
    } else {
        std::array::from_fn(|i| i as u8)
    };

    let rows = tables.rows;

    for cycle in lo..hi {
        let i = if reverse { len - cycle - 1 } else { cycle };
        let code = usize::from((seq[i / 2] >> ((!i & 1) << 2)) & 0xf);
        let j = cycle - lo;

        cycles[j][code] += 1;
        tables.revcomp[j][usize::from(complement[code])] += 1;
        quals[j * rows + usize::from(qual[i])] += 1;
    }
}

/// Ends of up to `workers` cycle ranges covering `0..max_len` with about equal numbers of bases:
/// cycle `c` has a base from every read longer than `c`.
pub(super) fn balanced_bounds(reads: &[PendingRead], max_len: usize, workers: usize) -> Vec<usize> {
    let mut ending = vec![0u64; max_len + 1];

    for read in reads {
        ending[read.len] += 1;
    }

    let total: u64 = reads.iter().map(|read| read.len as u64).sum();
    let target = total.div_ceil(workers as u64).max(1);

    let mut bounds = Vec::with_capacity(workers);
    let mut active = reads.len() as u64;
    let mut work = 0;

    for (cycle, &n) in ending.iter().enumerate().take(max_len) {
        active -= n;
        work += active;

        if work >= target && bounds.len() + 1 < workers {
            bounds.push(cycle + 1);
            work = 0;
        }
    }

    bounds.push(max_len);
    bounds.dedup();
    bounds
}

impl Stats {
    /// The per-cycle tables for every cycle.
    pub(super) fn cycle_tables(&mut self) -> CycleTables<'_> {
        CycleTables {
            lo: 0,
            first: &mut self.acgtno_cycles_1st,
            last: &mut self.acgtno_cycles_2nd,
            revcomp: &mut self.acgtno_revcomp,
            quals_first: &mut self.quals_1st,
            quals_last: &mut self.quals_2nd,
            rows: self.quals_rows,
        }
    }

    /// Counts the pending reads on `workers` threads, splitting the cycles among them.
    pub(super) fn flush_cycles(&mut self) {
        let pending = std::mem::take(&mut self.pending);

        if pending.reads.is_empty() {
            self.pending = pending;
            return;
        }

        let workers = self.workers.max(1);
        let rows = self.quals_rows;
        let bounds = balanced_bounds(&pending.reads, pending.max_len, workers);

        let mut first = &mut self.acgtno_cycles_1st[..];
        let mut last = &mut self.acgtno_cycles_2nd[..];
        let mut revcomp = &mut self.acgtno_revcomp[..];
        let mut quals_first = &mut self.quals_1st[..];
        let mut quals_last = &mut self.quals_2nd[..];

        let mut parts = Vec::with_capacity(workers);
        let mut lo = 0;

        for &hi in &bounds {
            let n = hi - lo;

            let (a, rest) = first.split_at_mut(n);
            first = rest;
            let (b, rest) = last.split_at_mut(n);
            last = rest;
            let (c, rest) = revcomp.split_at_mut(n);
            revcomp = rest;
            let (d, rest) = quals_first.split_at_mut(n * rows);
            quals_first = rest;
            let (e, rest) = quals_last.split_at_mut(n * rows);
            quals_last = rest;

            parts.push(CycleTables {
                lo,
                first: a,
                last: b,
                revcomp: c,
                quals_first: d,
                quals_last: e,
                rows,
            });

            lo += n;
        }

        let reads = &pending.reads;
        let data = &pending.data;

        std::thread::scope(|scope| {
            for mut tables in parts {
                scope.spawn(move || {
                    for read in reads {
                        let seq_end = read.start + read.len.div_ceil(2);
                        let seq = &data[read.start..seq_end];
                        let qual = &data[seq_end..seq_end + read.len];
                        count_cycles(seq, qual, read.strand, read.segment, &mut tables);
                    }
                });
            }
        });

        let mut pending = pending;
        pending.data.clear();
        pending.reads.clear();
        pending.max_len = 0;
        self.pending = pending;
    }
}
