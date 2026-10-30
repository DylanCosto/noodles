//! Collecting statistics from each record.

use std::collections::HashMap;

use anyhow::bail;
use noodles_sam::alignment::record::Flags;

use super::{
    BARCODE_ERROR_LIMIT, BWA_MIN_READ_LENGTH, CIGAR_DELETION, CIGAR_HARD_CLIP, CIGAR_INSERTION,
    CIGAR_MATCH, CIGAR_PAD, CIGAR_SEQUENCE_MATCH, CIGAR_SEQUENCE_MISMATCH, CIGAR_SKIP,
    CIGAR_SOFT_CLIP, GC_IN_BYTE, INWARD, Info, OTHER, OUTWARD, PAIR_MAX, PENDING_BASES, Pair,
    QUALITY_BINS, Stats,
    cycles::{BaseCounts, PendingRead, count_cycles},
    record::{
        Bam, READ_ORDER_FIRST, READ_ORDER_LAST, consumes_reference, integer_field, string_value,
    },
};
use crate::bam_raw;

/// `cleanup_overlaps`: drops pairs whose chunks all end before `max`, returning how many.
pub(super) fn cleanup_overlaps(pairs: &mut HashMap<Vec<u8>, Pair>, max: i64) -> u32 {
    let before = pairs.len();
    pairs.retain(|_, pair| pair.chunks.last().is_some_and(|&(_, end)| end >= max));
    u32::try_from(before - pairs.len()).unwrap_or(u32::MAX)
}

/// `bwa_trim_read`: the number of bases BWA would trim from the 3' end.
pub(super) fn bwa_trim_read(trim_qual: i32, quals: &[u8], len: usize, reverse: bool) -> usize {
    if len < BWA_MIN_READ_LENGTH {
        return 0;
    }

    let max_trimmed = len - BWA_MIN_READ_LENGTH + 1;
    let (mut sum, mut max_sum, mut max_l) = (0i32, 0i32, 0usize);

    for l in 0..max_trimmed {
        sum += trim_qual - i32::from(quals[if reverse { l } else { len - 1 - l }]);

        if sum < 0 {
            break;
        }

        if sum > max_sum {
            max_sum = sum;
            // BWA clips one base less than the maximum.
            max_l = l;
        }
    }

    max_l
}

impl Stats {
    /// `count_indels`
    pub(super) fn count_indels(&mut self, bam: &Bam<'_>) -> anyhow::Result<()> {
        let is_fwd = !bam.has(Flags::REVERSE_COMPLEMENTED);
        let order = bam.order();
        let read_len = i64::try_from(bam.l_seq)?;
        let mut icycle: i64 = 0;

        for &op in &bam.cigar {
            let (cig, ncig) = (op & 0xf, i64::from(op >> 4));

            if ncig == 0 {
                continue;
            }

            if cig == CIGAR_INSERTION {
                let idx = if is_fwd {
                    icycle
                } else {
                    read_len - icycle - ncig
                };

                if idx < 0 {
                    bail!("FIXME: read_len={read_len} vs icycle={icycle}");
                }

                let idx = usize::try_from(idx)?;

                if idx >= self.nbases {
                    bail!("FIXME: {idx} vs {}", self.nbases);
                }

                if order == READ_ORDER_FIRST {
                    self.ins_cycles_1st[idx] += 1;
                }

                if order == READ_ORDER_LAST {
                    self.ins_cycles_2nd[idx] += 1;
                }

                icycle += ncig;

                if ncig as usize <= self.nindels {
                    self.insertions[ncig as usize - 1] += 1;
                }

                continue;
            }

            if cig == CIGAR_DELETION {
                let idx = if is_fwd {
                    icycle - 1
                } else {
                    read_len - icycle - 1
                };

                if idx < 0 {
                    continue;
                }

                let idx = usize::try_from(idx)?;

                if idx >= self.nbases {
                    bail!("FIXME: {idx} vs {}", self.nbases);
                }

                if order == READ_ORDER_FIRST {
                    self.del_cycles_1st[idx] += 1;
                }

                if order == READ_ORDER_LAST {
                    self.del_cycles_2nd[idx] += 1;
                }

                if ncig as usize <= self.nindels {
                    self.deletions[ncig as usize - 1] += 1;
                }

                continue;
            }

            if cig != CIGAR_SKIP && cig != CIGAR_HARD_CLIP && cig != CIGAR_PAD {
                icycle += ncig;
            }
        }

        Ok(())
    }

    /// `count_mismatches_per_cycle`
    pub(super) fn count_mismatches_per_cycle(
        &mut self,
        bam: &Bam<'_>,
        read_len: usize,
    ) -> anyhow::Result<()> {
        let is_fwd = !bam.has(Flags::REVERSE_COMPLEMENTED);
        let mut iread = 0usize;
        let mut icycle = 0usize;
        let mut iref = bam.pos - self.rseq_pos;
        let max_len = self.max_len;
        let nbases = self.nbases;
        let mpc = self.mpc_buf.as_mut().expect("reference is loaded");

        for &op in &bam.cigar {
            let (cig, ncig) = (op & 0xf, (op >> 4) as usize);

            match cig {
                CIGAR_INSERTION => {
                    iread += ncig;
                    icycle += ncig;
                    continue;
                }
                CIGAR_DELETION => {
                    iref += ncig as i64;
                    continue;
                }
                CIGAR_SOFT_CLIP => {
                    icycle += ncig;
                    iread += ncig;
                    continue;
                }
                CIGAR_HARD_CLIP => {
                    icycle += ncig;
                    continue;
                }
                CIGAR_SKIP | CIGAR_PAD => continue,
                CIGAR_MATCH | CIGAR_SEQUENCE_MATCH | CIGAR_SEQUENCE_MISMATCH => {}
                _ => bail!("TODO: cigar {cig}"),
            }

            if ncig as i64 + iref > self.nrseq_buf {
                bail!(
                    "FIXME: {ncig}+{iref} > {}, {}",
                    self.nrseq_buf,
                    String::from_utf8_lossy(bam.name)
                );
            }

            for _ in 0..ncig {
                let cread = bam.base(iread);
                let cref = self.rseq_buf[usize::try_from(iref)?];

                let cycle_index = |icycle: usize| {
                    if is_fwd {
                        icycle
                    } else {
                        read_len - icycle - 1
                    }
                };

                if cread == 15 {
                    let idx = cycle_index(icycle);

                    if idx > max_len {
                        bail!("mpc: {idx}>{max_len}");
                    }

                    let idx = idx * QUALITY_BINS;

                    if idx >= QUALITY_BINS * nbases {
                        bail!("FIXME: mpc_buf overflow");
                    }

                    mpc[idx] += 1;
                } else if cref != 0 && cread != 0 && cref != cread {
                    let qual = usize::from(bam.qual[iread]) + 1;
                    let idx = cycle_index(icycle);

                    if idx > max_len {
                        bail!("mpc: {idx}>{max_len}");
                    }

                    let idx = idx * QUALITY_BINS + qual;

                    if idx >= QUALITY_BINS * nbases {
                        bail!("FIXME: mpc_buf overflow");
                    }

                    mpc[idx] += 1;
                }

                iref += 1;
                iread += 1;
                icycle += 1;
            }
        }

        Ok(())
    }

    /// `update_checksum`
    pub(super) fn update_checksum(&mut self, bam: &Bam<'_>) {
        let crc = |bytes: &[u8]| {
            let mut crc = flate2::Crc::new();
            crc.update(bytes);
            crc.sum()
        };

        self.checksum_names = self.checksum_names.wrapping_add(crc(bam.name));

        if bam.l_seq == 0 {
            return;
        }

        self.checksum_reads = self.checksum_reads.wrapping_add(crc(bam.seq));
        self.checksum_quals = self.checksum_quals.wrapping_add(crc(bam.qual));
    }

    /// `collect_orig_read_stats`; returns the read's GC count.
    pub(super) fn collect_orig_read_stats(&mut self, bam: &Bam<'_>, info: &Info) -> usize {
        let seq_len = bam.l_seq;
        self.total_len += seq_len as u64;

        if bam.has(Flags::QC_FAIL) {
            self.nreads_qc_failed += 1;
        }

        if bam.has(Flags::SEGMENTED) {
            self.nreads_paired_tech += 1;
        }

        let order = bam.order();
        let reverse = bam.has(Flags::REVERSE_COMPLEMENTED);
        let mut gc_count = 0;

        if order == READ_ORDER_FIRST || order == READ_ORDER_LAST {
            gc_count = bam
                .seq
                .iter()
                .map(|&b| usize::from(GC_IN_BYTE[usize::from(b)]))
                .sum();
        }

        let gc_idx_min = gc_count * (self.ngc - 1) / seq_len;
        let gc_idx_max = ((gc_count + 1) * (self.ngc - 1) / seq_len).min(self.ngc - 1);

        let quals = match order {
            READ_ORDER_FIRST => {
                self.nreads_1st += 1;
                self.total_len_1st += seq_len as u64;

                for n in &mut self.gc_1st[gc_idx_min..gc_idx_max.max(gc_idx_min)] {
                    *n += 1;
                }

                Some(READ_ORDER_FIRST)
            }
            READ_ORDER_LAST => {
                self.nreads_2nd += 1;
                self.total_len_2nd += seq_len as u64;

                for n in &mut self.gc_2nd[gc_idx_min..gc_idx_max.max(gc_idx_min)] {
                    *n += 1;
                }

                Some(READ_ORDER_LAST)
            }
            _ => {
                self.nreads_other += 1;
                None
            }
        };

        if info.trim_qual > 0 {
            self.nbases_trimmed += bwa_trim_read(info.trim_qual, bam.qual, seq_len, reverse) as u64;
        }

        if let Some(order) = quals {
            let max_qual = usize::from(bam.qual.iter().copied().max().unwrap_or(0));

            if max_qual > self.max_qual {
                self.flush_cycles();
                self.realloc_qual_buffers(self.nbases, self.nbases, max_qual);
                self.max_qual = max_qual;
            }

            self.sum_qual += bam.qual.iter().map(|&q| u64::from(q)).sum::<u64>();

            let first = order == READ_ORDER_FIRST;

            if self.workers == 0 {
                let mut tables = self.cycle_tables();
                count_cycles(bam.seq, bam.qual, reverse, first, &mut tables);
            } else {
                let pending = &mut self.pending;
                let start = pending.data.len();
                pending.data.extend_from_slice(bam.seq);
                pending.data.extend_from_slice(bam.qual);
                pending.reads.push(PendingRead {
                    start,
                    len: seq_len,
                    reverse,
                    first,
                });
                pending.max_len = pending.max_len.max(seq_len);

                if pending.data.len() >= PENDING_BASES {
                    self.flush_cycles();
                }
            }
        }

        if order == READ_ORDER_FIRST {
            self.collect_barcode_stats(bam);
        }

        if bam.has(Flags::UNMAPPED) {
            self.nreads_unmapped += 1;
        } else {
            self.nbases_mapped += seq_len as u64;

            if bam.mapq == 0 {
                self.nreads_mq0 += 1;
            }

            if !bam.is_paired_and_mapped() {
                self.nreads_single_mapped += 1;
            } else {
                self.nreads_paired_and_mapped += 1;

                if bam.has(Flags::PROPERLY_SEGMENTED) {
                    self.nreads_properly_paired += 1;
                }

                if bam.tid != bam.mtid {
                    self.nreads_anomalous += 1;
                }
            }
        }

        gc_count
    }

    /// `collect_barcode_stats`
    pub(super) fn collect_barcode_stats(&mut self, bam: &Bam<'_>) {
        // `bam_aux2Z`: the value of a Z or H field.
        let text = |tag: [u8; 2]| -> Option<&[u8]> {
            let range = bam_raw::find_field(bam.raw, tag).ok()??;
            matches!(bam.raw[range.start + 2], b'Z' | b'H')
                .then(|| &bam.raw[range.start + 3..range.end - 1])
        };

        let name = String::from_utf8_lossy(bam.name);

        for barcode in &mut self.barcodes {
            let Some(sequence) = text(barcode.tag) else {
                continue;
            };

            // An empty barcode counts as none.
            if sequence.is_empty() {
                continue;
            }

            if barcode.len == 0 {
                barcode.len = sequence.len();
                barcode.counts = vec![BaseCounts::default(); sequence.len()];
                barcode.quals = vec![0; sequence.len() * QUALITY_BINS];
            }

            let tag = String::from_utf8_lossy(&barcode.tag);

            if sequence.len() > barcode.len {
                eprintln!("Barcodes with tag {tag} differ in length at sequence '{name}'");
                continue;
            }

            let mut error = false;

            for (i, &b) in sequence.iter().enumerate() {
                let counts = &mut barcode.counts[i];

                match b {
                    b'A' => counts.a += 1,
                    b'C' => counts.c += 1,
                    b'G' => counts.g += 1,
                    b'T' => counts.t += 1,
                    b'N' => counts.n += 1,
                    _ => match barcode.separator {
                        Some(sep) if sep != i => {
                            if self.barcode_errors < BARCODE_ERROR_LIMIT {
                                eprintln!(
                                    "Barcode separator for tag {tag} is in a different position or wrong barcode content('{}') at sequence '{name}'",
                                    String::from_utf8_lossy(sequence)
                                );
                                self.barcode_errors += 1;
                            }

                            error = true;
                        }
                        Some(_) => {}
                        None => barcode.separator = Some(i),
                    },
                }

                if error {
                    break;
                }
            }

            if error {
                continue;
            }

            let Some(qualities) = text(barcode.quality_tag) else {
                continue;
            };

            if qualities.len() == sequence.len() {
                for (i, &c) in qualities.iter().enumerate() {
                    if let Some(qual) = c.checked_sub(b'!').map(usize::from) {
                        barcode.quals[i * QUALITY_BINS + qual] += 1;
                        barcode.max_qual = barcode.max_qual.max(Some(qual));
                    }
                }
            } else {
                if self.barcode_errors < BARCODE_ERROR_LIMIT {
                    eprintln!(
                        "{tag} length and {} length don't match for sequence '{name}'",
                        String::from_utf8_lossy(&barcode.quality_tag)
                    );
                }

                self.barcode_errors += 1;
            }
        }
    }

    /// `remove_overlaps`: adds `[pmin, pmax)` to the coverage, less what the mate already covered.
    pub(super) fn remove_overlaps(
        &mut self,
        bam: &Bam<'_>,
        pairs: &mut HashMap<Vec<u8>, Pair>,
        pmin: i64,
        pmax: i64,
    ) -> anyhow::Result<()> {
        let order =
            u32::from(bam.has(Flags::FIRST_SEGMENT)) + 2 * u32::from(bam.has(Flags::LAST_SEGMENT));

        if !bam.has(Flags::SEGMENTED)
            || bam.has(Flags::MATE_UNMAPPED)
            || i64::from(bam.isize).abs() >= 2 * bam.l_seq as i64
            || (order != READ_ORDER_FIRST && order != READ_ORDER_LAST)
        {
            if pmin >= 0 {
                self.insert_read(pmin, pmax)?;
            }

            return Ok(());
        }

        let mut pmin = pmin;

        match pairs.get_mut(bam.name) {
            None => {
                if pmin == -1 {
                    return Ok(());
                }

                pairs.insert(
                    bam.name.to_vec(),
                    Pair {
                        first: order,
                        chunks: vec![(pmin, pmax)],
                    },
                );
                self.pair_count += 1;
            }
            Some(pair) if order == pair.first => {
                if pmin == -1 {
                    return Ok(());
                }

                pair.chunks.push((pmin, pmax));
            }
            Some(pair) => {
                if pmin == -1 {
                    pairs.remove(bam.name);
                    self.pair_count = self.pair_count.wrapping_sub(1);
                    return Ok(());
                }

                for &(beg, end) in &pair.chunks {
                    if pmin >= end {
                        continue;
                    }

                    if pmax <= beg {
                        break;
                    }

                    if pmin < beg {
                        self.insert_read(pmin, beg)?;
                        pmin = beg;
                    }

                    if pmax <= end {
                        self.nbases_mapped_cigar =
                            self.nbases_mapped_cigar.wrapping_add_signed(-(pmax - pmin));
                        return Ok(());
                    }

                    self.nbases_mapped_cigar =
                        self.nbases_mapped_cigar.wrapping_add_signed(-(end - pmin));
                    pmin = end;
                }
            }
        }

        self.insert_read(pmin, pmax)
    }

    /// `collect_stats`
    pub(super) fn collect(
        &mut self,
        bam: &Bam<'_>,
        info: &Info,
        pairs: &mut HashMap<Vec<u8>, Pair>,
    ) -> anyhow::Result<()> {
        if !self.is_in_regions(bam)? {
            return Ok(());
        }

        if let Some(read_groups) = &self.read_groups {
            let Some(rg) = string_value(bam.raw, *b"RG")? else {
                return Ok(());
            };

            if !read_groups.iter().any(|id| id == rg) {
                return Ok(());
            }
        }

        if info.flag_require != 0 && bam.flag & info.flag_require != info.flag_require {
            self.nreads_filtered += 1;
            return Ok(());
        }

        if info.flag_filter != 0 && bam.flag & info.flag_filter != 0 {
            self.nreads_filtered += 1;
            return Ok(());
        }

        if info.filter_readlen != -1 && bam.l_seq as i64 != info.filter_readlen {
            return Ok(());
        }

        self.update_checksum(bam);

        if bam.has(Flags::SECONDARY) {
            self.nreads_secondary += 1;
            return Ok(());
        }

        if bam.has(Flags::SUPPLEMENTARY) {
            self.nreads_supplementary += 1;
        }

        let seq_len = bam.l_seq;

        if seq_len == 0 {
            return Ok(());
        }

        if bam.has(Flags::DUPLICATE) {
            self.total_len_dup += seq_len as u64;
            self.nreads_dup += 1;
        }

        let order = bam.order();

        let read_len = seq_len
            + bam
                .cigar
                .iter()
                .filter(|&&op| op & 0xf == CIGAR_HARD_CLIP)
                .map(|&op| (op >> 4) as usize)
                .sum::<usize>();

        if read_len >= self.nbases {
            self.flush_cycles();
            self.realloc_buffers(read_len, info);
        }

        self.max_len = self.max_len.max(read_len);

        if order == READ_ORDER_FIRST {
            self.max_len_1st = self.max_len_1st.max(read_len);
        }

        if order == READ_ORDER_LAST {
            self.max_len_2nd = self.max_len_2nd.max(read_len);
        }

        let excluded = Flags::UNMAPPED
            | Flags::SECONDARY
            | Flags::SUPPLEMENTARY
            | Flags::QC_FAIL
            | Flags::DUPLICATE;

        if bam.flag & excluded.bits() == 0 {
            self.mapping_qualities[usize::from(bam.mapq)] += 1;
        }

        let mut gc_count = 0;

        if bam.is_original() {
            self.read_lengths[read_len] += 1;

            if order == READ_ORDER_FIRST {
                self.read_lengths_1st[read_len] += 1;
            }

            if order == READ_ORDER_LAST {
                self.read_lengths_2nd[read_len] += 1;
            }

            gc_count = self.collect_orig_read_stats(bam, info);
        }

        if bam.has(Flags::UNMAPPED) {
            return Ok(());
        }

        self.count_indels(bam)?;

        if bam.is_paired_and_mapped() && bam.is_original() {
            let mut isize = bam.isize.wrapping_abs();

            if info.nisize > 0 && isize > info.nisize {
                isize = info.nisize;
            }

            if isize > 0 || bam.tid == bam.mtid {
                let pos_fst = bam.mpos - bam.pos;
                let is_fst: i64 = if bam.has(Flags::FIRST_SEGMENT) { 1 } else { -1 };
                let is_fwd: i64 = if bam.has(Flags::REVERSE_COMPLEMENTED) {
                    -1
                } else {
                    1
                };
                let is_mfwd: i64 = if bam.has(Flags::MATE_REVERSE_COMPLEMENTED) {
                    -1
                } else {
                    1
                };

                let orientation = if is_fwd * is_mfwd > 0 {
                    OTHER
                } else if is_fst * pos_fst > 0 {
                    if is_fst * is_fwd > 0 { INWARD } else { OUTWARD }
                } else if is_fst * pos_fst < 0 {
                    if is_fst * is_fwd > 0 { OUTWARD } else { INWARD }
                } else {
                    // Exactly overlapping reads are taken as inward.
                    INWARD
                };

                self.isize.increment(isize, orientation);
            }
        }

        if let Some(nm) = integer_field(bam.raw, *b"NM")? {
            self.nmismatches = self.nmismatches.wrapping_add_signed(nm);
        }

        if bam.cigar.is_empty() {
            bail!("FIXME: mapped read with no cigar?");
        }

        let mut readlen = seq_len as i64;

        if self.regions.is_some() {
            // Count only on-target bases
            let mut iref = bam.pos + 1;

            for &op in &bam.cigar {
                let (cig, ncig_init) = (op & 0xf, i64::from(op >> 4));
                let mut ncig = ncig_init;

                if ncig == 0 {
                    continue;
                }

                match cig {
                    CIGAR_DELETION => readlen += ncig,
                    CIGAR_MATCH | CIGAR_SEQUENCE_MATCH | CIGAR_SEQUENCE_MISMATCH => {
                        if iref < self.reg_from {
                            ncig -= self.reg_from - iref;
                        }

                        if iref + ncig_init - 1 > self.reg_to {
                            ncig -= iref + ncig_init - 1 - self.reg_to;
                        }

                        self.nbases_mapped_cigar =
                            self.nbases_mapped_cigar.wrapping_add_signed(ncig.max(0));
                        iref += ncig_init;
                    }
                    CIGAR_INSERTION => {
                        iref += ncig;

                        if iref >= self.reg_from && iref <= self.reg_to {
                            self.nbases_mapped_cigar += ncig as u64;
                        }
                    }
                    _ => {}
                }
            }
        } else {
            for &op in &bam.cigar {
                let (cig, len) = (op & 0xf, u64::from(op >> 4));

                if matches!(
                    cig,
                    CIGAR_MATCH | CIGAR_INSERTION | CIGAR_SEQUENCE_MATCH | CIGAR_SEQUENCE_MISMATCH
                ) {
                    self.nbases_mapped_cigar += len;
                }

                if cig == CIGAR_DELETION {
                    readlen += len as i64;
                }
            }
        }

        if self.tid == bam.tid && bam.pos < self.pos {
            self.is_sorted = false;
        }

        self.pos = bam.pos;

        if !self.is_sorted {
            return Ok(());
        }

        if self.tid == -1 || self.tid != bam.tid {
            self.round_buffer_flush(-1, info)?;
        }

        // Drop pairs that can no longer overlap, to bound memory.
        self.last_read_flush += 1;

        if self.pair_count > PAIR_MAX && self.last_read_flush > PAIR_MAX {
            let removed = cleanup_overlaps(pairs, bam.pos);
            self.pair_count = self.pair_count.wrapping_sub(removed);
            self.last_read_flush = 0;
        }

        if self.last_pair_tid != bam.tid {
            let removed = cleanup_overlaps(pairs, i64::MAX - 1);
            self.pair_count = self.pair_count.wrapping_sub(removed);
            self.last_pair_tid = bam.tid;
            self.last_read_flush = 0;
        }

        let gcd_bin_size = i64::from(info.gcd_bin_size);

        if info.reference.is_some() {
            let mut inc_ref = 0;
            let mut inc_gcd = false;

            if self.rseq_pos == -1
                || self.tid != bam.tid
                || self.rseq_pos + self.nrseq_buf < bam.pos + readlen
            {
                inc_ref = bam.pos + readlen;
                inc_gcd = true;
            } else if self.gcd_pos + gcd_bin_size < bam.pos + readlen {
                inc_gcd = true;

                if self.rseq_pos + self.nrseq_buf < bam.pos + gcd_bin_size {
                    inc_ref = bam.pos + gcd_bin_size;
                }
            }

            if inc_gcd {
                self.igcd += 1;
                self.ensure_gcd(info);

                if inc_ref != 0 {
                    self.read_ref_seq(info, bam.tid, bam.pos, inc_ref)?;
                }

                self.gcd_pos = bam.pos;
                self.gcd[self.igcd].gc = self.gc_content(self.gcd_pos, info.gcd_bin_size);
            }

            self.count_mismatches_per_cycle(bam, read_len)?;
        } else if self.gcd_pos == -1 || self.tid != bam.tid || bam.pos - self.gcd_pos > gcd_bin_size
        {
            self.tid = bam.tid;
            self.gcd_pos = bam.pos;
            self.igcd += 1;
            self.ensure_gcd(info);
        }

        self.ensure_gcd(info);
        self.gcd[self.igcd].depth = self.gcd[self.igcd].depth.wrapping_add(1);

        if info.reference.is_none() {
            self.gcd[self.igcd].gc += gc_count as f32 / seq_len as f32;
        }

        // Coverage distribution
        self.round_buffer_flush(bam.pos, info)?;

        if self.regions.is_some() {
            let mut p = bam.pos;
            let (mut j, mut i) = (0, 0);

            while j < bam.cigar.len() && i < self.nchunks {
                let (op, oplen) = (bam.cigar[j] & 0xf, i64::from(bam.cigar[j] >> 4));

                if matches!(
                    op,
                    CIGAR_MATCH | CIGAR_SEQUENCE_MATCH | CIGAR_SEQUENCE_MISMATCH
                ) {
                    let pmin = p.max(self.chunks[i].0 - 1);
                    let pmax = (p + oplen).min(self.chunks[i].1);

                    if pmax > pmin {
                        if info.remove_overlaps {
                            self.remove_overlaps(bam, pairs, pmin, pmax)?;
                        } else {
                            self.insert_read(pmin, pmax)?;
                        }
                    }
                }

                let pnew = p + if consumes_reference(op) { oplen } else { 0 };

                if pnew >= self.chunks[i].1 {
                    i += 1;
                } else {
                    j += 1;
                    p = pnew;
                }
            }
        } else {
            let mut p = bam.pos;

            for &op in &bam.cigar {
                let (op, oplen) = (op & 0xf, i64::from(op >> 4));

                if matches!(
                    op,
                    CIGAR_MATCH | CIGAR_SEQUENCE_MATCH | CIGAR_SEQUENCE_MISMATCH
                ) {
                    if info.remove_overlaps {
                        self.remove_overlaps(bam, pairs, p, p + oplen)?;
                    } else {
                        self.insert_read(p, p + oplen)?;
                    }
                }

                if consumes_reference(op) {
                    p += oplen;
                }
            }
        }

        if info.remove_overlaps {
            self.remove_overlaps(bam, pairs, -1, -1)?;
        }

        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_bwa_trim_read() {
        let quals = [30; 50];
        assert_eq!(bwa_trim_read(20, &quals, 50, false), 0);

        let mut quals = [30; 50];
        quals[45..].fill(2);
        assert_eq!(bwa_trim_read(20, &quals, 50, false), 4);
        assert_eq!(bwa_trim_read(20, &quals, 30, false), 0);
    }
}
