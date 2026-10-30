//! Writing the report.

use std::io::{self, Write};

use noodles_sam as sam;

use super::{GcDepth, Info, InsertSizeRows, QUALITY_BINS, Stats, cycles::BaseCounts};

/// `gcd_percentile`, in single precision.
pub(super) fn gcd_percentile(gcd: &[GcDepth], n: usize, p: u32) -> f32 {
    let x = p as f32 * (n as u32 + 1) as f32 / 100.0;
    #[expect(
        clippy::cast_possible_truncation,
        reason = "samtools truncates the float to an int"
    )]
    let k = x as i32;

    if k <= 0 {
        return gcd[0].depth as f32;
    }

    let k = k as usize;

    if k >= n {
        return gcd[n - 1].depth as f32;
    }

    let d = x - k as f32;
    gcd[k - 1].depth as f32 + d * gcd[k].depth.wrapping_sub(gcd[k - 1].depth) as f32
}

pub(super) fn write_quality_rows<W: Write>(
    out: &mut W,
    label: &str,
    quals: &[u64],
    rows: usize,
    max_len: usize,
    max_qual: usize,
) -> io::Result<()> {
    let mut line = Vec::new();

    for (ibase, counts) in quals.chunks_exact(rows).take(max_len).enumerate() {
        line.clear();
        line.extend_from_slice(label.as_bytes());
        line.push(b'\t');
        push_u64(&mut line, ibase as u64 + 1);

        for &n in &counts[..=max_qual] {
            line.push(b'\t');
            push_u64(&mut line, n);
        }

        line.push(b'\n');
        out.write_all(&line)?;
    }

    Ok(())
}

/// Writes `label`, the cycle, and each count as a percentage of `sum` (`%.2f`).
pub(super) fn write_percent_row<W: Write>(
    out: &mut W,
    label: &str,
    cycle: usize,
    counts: &[u64],
    sum: u64,
) -> io::Result<()> {
    let mut line = Vec::with_capacity(64);
    line.extend_from_slice(label.as_bytes());
    line.push(b'\t');
    push_u64(&mut line, cycle as u64);

    for &n in counts {
        line.push(b'\t');
        push_fixed2(&mut line, 100.0 * n as f64 / sum as f64);
    }

    line.push(b'\n');
    out.write_all(&line)
}

pub(super) fn push_u64(dst: &mut Vec<u8>, n: u64) {
    let mut buf = [0; lexical_core::BUFFER_SIZE];
    dst.extend_from_slice(lexical_core::write(n, &mut buf));
}

/// Appends `printf("%.2f", value)` for a nonnegative value.
pub(super) fn push_fixed2(dst: &mut Vec<u8>, value: f64) {
    // `value * 100` is within an ulp of the exact product, so unless it lies near a rounding
    // boundary, rounding it gives the digits printf would; otherwise defer to exact formatting.
    let scaled = value * 100.0;
    let rounded = scaled.round();

    if value.is_finite()
        && (0.0..1e15).contains(&value)
        && ((scaled - rounded).abs() - 0.5).abs() > 1e-6
    {
        let n = rounded as u64;
        push_u64(dst, n / 100);
        dst.push(b'.');
        dst.push(b'0' + (n % 100 / 10) as u8);
        dst.push(b'0' + (n % 10) as u8);
    } else {
        dst.extend_from_slice(format!("{value:.2}").as_bytes());
    }
}

pub(super) fn write_gc_rows<W: Write>(
    out: &mut W,
    label: &str,
    counts: &[u64],
    ngc: usize,
) -> io::Result<()> {
    let mut prev = 0;

    for ibase in 0..ngc {
        if counts[ibase] == counts[prev] {
            continue;
        }

        writeln!(
            out,
            "{label}\t{:.2}\t{}",
            (ibase + prev) as f64 * 0.5 * 100.0 / (ngc - 1) as f64,
            counts[prev]
        )?;
        prev = ibase;
    }

    Ok(())
}

/// The read group IDs to include for `-I`: those with this ID or sample name.
pub(super) fn read_groups_for(header: &sam::Header, id: &str) -> Vec<Vec<u8>> {
    use sam::header::record::value::map::read_group::tag::SAMPLE;

    header
        .read_groups()
        .iter()
        .filter(|(key, map)| {
            key.as_slice() == id.as_bytes()
                || map
                    .other_fields()
                    .get(&SAMPLE)
                    .is_some_and(|sm| sm.as_slice() == id.as_bytes())
        })
        .map(|(key, _)| key.to_vec())
        .collect()
}

/// Formats like C's `printf("%e", value)`.
pub(super) fn format_e(value: f64) -> String {
    let s = format!("{value:.6e}");
    let (mantissa, exponent) = s.split_once('e').expect("scientific notation");
    let exponent: i32 = exponent.parse().expect("valid exponent");
    let sign = if exponent < 0 { '-' } else { '+' };
    format!("{mantissa}e{sign}{:02}", exponent.abs())
}

/// Insert size statistics for the summary, from [`Stats::insert_size_summary`].
#[derive(Clone, Copy, Debug)]
struct InsertSizeSummary {
    /// One past the last insert size in the main bulk.
    ibulk: i32,
    average: f64,
    sd: f64,
    inward: u64,
    outward: u64,
    other: u64,
}

impl Stats {
    /// `output_stats`
    pub(super) fn output<W: Write>(
        &mut self,
        out: &mut W,
        info: &Info,
        command_line: &str,
        split_tag: Option<&str>,
    ) -> io::Result<()> {
        self.flush_cycles();

        let insert_sizes = self.insert_size_summary(info);

        self.write_preamble(out, command_line, split_tag)?;
        self.write_summary(out, info, &insert_sizes)?;
        self.round_up_maxima();
        self.write_quality_tables(out)?;
        self.write_base_content(out)?;
        self.write_barcodes(out)?;
        self.write_distributions(out, info, insert_sizes.ibulk)?;

        if self.is_sorted {
            self.write_coverage(out, info)?;
            self.write_gc_depth(out, info, self.average_read_length())?;
        }

        self.write_reference_stats(out)
    }

    /// The insert size average and standard deviation, from the main bulk of insert sizes only,
    /// and the pair counts by orientation.
    fn insert_size_summary(&mut self, info: &Info) -> InsertSizeSummary {
        self.isize.halve();

        let nitems = self.isize.len();
        let (mut nisize, mut nisize_inward, mut nisize_outward, mut nisize_other) =
            (0u64, 0u64, 0u64, 0u64);

        for isize in 0..nitems {
            let [i, o, x] = self.isize.get(isize);
            nisize_inward += i;
            nisize_outward += o;
            nisize_other += x;
            nisize += i + o + x;
        }

        let mut ibulk = 0;
        let mut bulk = 0.0f64;
        let mut avg_isize = 0.0f64;

        for isize in 0..nitems {
            let [i, o, x] = self.isize.get(isize);
            let num = i + o + x;

            if num > 0 {
                ibulk = isize + 1;
            }

            bulk += num as f64;
            avg_isize += (isize as u64).wrapping_mul(num) as f64;

            if bulk / nisize as f64 > f64::from(info.isize_main_bulk) {
                ibulk = isize + 1;
                nisize = bulk as u64;
                break;
            }
        }

        let divisor = if nisize > 0 { nisize } else { 1 };
        avg_isize /= divisor as f64;

        let mut sd_isize = 0.0f64;

        for isize in 1..ibulk {
            let [i, o, x] = self.isize.get(isize);
            let d = f64::from(isize) - avg_isize;
            sd_isize += (i + o + x) as f64 * d * d / divisor as f64;
        }

        InsertSizeSummary {
            ibulk,
            average: avg_isize,
            sd: sd_isize.sqrt(),
            inward: nisize_inward,
            outward: nisize_outward,
            other: nisize_other,
        }
    }

    /// The average read length.
    fn average_read_length(&self) -> f32 {
        let nreads = self.nreads_1st + self.nreads_2nd + self.nreads_other;

        if nreads > 0 {
            self.total_len as f32 / nreads as f32
        } else {
            0.0
        }
    }

    /// The file header and the checksums (CHK).
    fn write_preamble<W: Write>(
        &self,
        out: &mut W,
        command_line: &str,
        split_tag: Option<&str>,
    ) -> io::Result<()> {
        writeln!(
            out,
            "# This file was produced by ntools stats ({}) and can be plotted using plot-bamstats",
            env!("CARGO_PKG_VERSION")
        )?;

        match (&self.split_name, split_tag) {
            (Some(name), Some(tag)) => writeln!(
                out,
                "# This file contains statistics only for reads with tag: {tag}={}",
                String::from_utf8_lossy(name)
            )?,
            _ => writeln!(out, "# This file contains statistics for all reads.")?,
        }

        writeln!(out, "# The command line was:  {command_line}")?;
        writeln!(
            out,
            "# CHK, Checksum\t[2]Read Names\t[3]Sequences\t[4]Qualities"
        )?;
        writeln!(
            out,
            "# CHK, CRC32 of reads which passed filtering followed by addition (32bit overflow)"
        )?;
        writeln!(
            out,
            "CHK\t{:08x}\t{:08x}\t{:08x}",
            self.checksum_names, self.checksum_reads, self.checksum_quals
        )?;

        Ok(())
    }

    /// The summary numbers (SN).
    fn write_summary<W: Write>(
        &self,
        out: &mut W,
        info: &Info,
        insert_sizes: &InsertSizeSummary,
    ) -> io::Result<()> {
        let InsertSizeSummary {
            average: avg_isize,
            sd: sd_isize,
            inward: nisize_inward,
            outward: nisize_outward,
            other: nisize_other,
            ..
        } = *insert_sizes;

        let nreads = self.nreads_1st + self.nreads_2nd + self.nreads_other;

        writeln!(
            out,
            "# Summary Numbers. Use `grep ^SN | cut -f 2-` to extract this part."
        )?;
        writeln!(
            out,
            "SN\traw total sequences:\t{}\t# excluding supplementary and secondary reads",
            self.nreads_filtered + nreads
        )?;
        writeln!(out, "SN\tfiltered sequences:\t{}", self.nreads_filtered)?;
        writeln!(out, "SN\tsequences:\t{nreads}")?;
        writeln!(
            out,
            "SN\tis sorted:\t{}\t# {} by coordinate",
            u8::from(self.is_sorted),
            if self.is_sorted {
                "sorted"
            } else {
                "not sorted"
            }
        )?;
        writeln!(out, "SN\t1st fragments:\t{}", self.nreads_1st)?;
        writeln!(out, "SN\tlast fragments:\t{}", self.nreads_2nd)?;
        writeln!(
            out,
            "SN\treads mapped:\t{}",
            self.nreads_paired_and_mapped + self.nreads_single_mapped
        )?;
        writeln!(
            out,
            "SN\treads mapped and paired:\t{}\t# paired-end technology bit set + both mates mapped",
            self.nreads_paired_and_mapped
        )?;
        writeln!(out, "SN\treads unmapped:\t{}", self.nreads_unmapped)?;
        writeln!(
            out,
            "SN\treads properly paired:\t{}\t# proper-pair bit set",
            self.nreads_properly_paired
        )?;
        writeln!(
            out,
            "SN\treads paired:\t{}\t# paired-end technology bit set",
            self.nreads_paired_tech
        )?;
        writeln!(
            out,
            "SN\treads duplicated:\t{}\t# PCR or optical duplicate bit set",
            self.nreads_dup
        )?;
        writeln!(
            out,
            "SN\treads MQ0:\t{}\t# mapped and MQ=0",
            self.nreads_mq0
        )?;
        writeln!(out, "SN\treads QC failed:\t{}", self.nreads_qc_failed)?;
        writeln!(
            out,
            "SN\tnon-primary alignments:\t{}",
            self.nreads_secondary
        )?;
        writeln!(
            out,
            "SN\tsupplementary alignments:\t{}",
            self.nreads_supplementary
        )?;
        writeln!(
            out,
            "SN\ttotal length:\t{}\t# ignores clipping",
            self.total_len
        )?;
        writeln!(
            out,
            "SN\ttotal first fragment length:\t{}\t# ignores clipping",
            self.total_len_1st
        )?;
        writeln!(
            out,
            "SN\ttotal last fragment length:\t{}\t# ignores clipping",
            self.total_len_2nd
        )?;
        writeln!(
            out,
            "SN\tbases mapped:\t{}\t# ignores clipping",
            self.nbases_mapped
        )?;
        writeln!(
            out,
            "SN\tbases mapped (cigar):\t{}\t# more accurate",
            self.nbases_mapped_cigar.cast_signed()
        )?;
        writeln!(out, "SN\tbases trimmed:\t{}", self.nbases_trimmed)?;
        writeln!(out, "SN\tbases duplicated:\t{}", self.total_len_dup)?;
        writeln!(
            out,
            "SN\tmismatches:\t{}\t# from NM fields",
            self.nmismatches.cast_signed()
        )?;

        let error_rate = if self.nbases_mapped_cigar > 0 {
            f64::from(self.nmismatches as f32 / self.nbases_mapped_cigar as f32)
        } else {
            0.0
        };

        writeln!(
            out,
            "SN\terror rate:\t{}\t# mismatches / bases mapped (cigar)",
            format_e(error_rate)
        )?;

        let avg_read_length = self.average_read_length();

        let ratio = |a: u64, b: u64| -> f64 {
            if b > 0 {
                f64::from(a as f32 / b as f32)
            } else {
                0.0
            }
        };

        writeln!(
            out,
            "SN\taverage length:\t{:.0}",
            f64::from(avg_read_length)
        )?;
        writeln!(
            out,
            "SN\taverage first fragment length:\t{:.0}",
            ratio(self.total_len_1st, self.nreads_1st)
        )?;
        writeln!(
            out,
            "SN\taverage last fragment length:\t{:.0}",
            ratio(self.total_len_2nd, self.nreads_2nd)
        )?;
        writeln!(out, "SN\tmaximum length:\t{}", self.max_len)?;
        writeln!(
            out,
            "SN\tmaximum first fragment length:\t{}",
            self.max_len_1st
        )?;
        writeln!(
            out,
            "SN\tmaximum last fragment length:\t{}",
            self.max_len_2nd
        )?;
        writeln!(
            out,
            "SN\taverage quality:\t{:.1}",
            if self.total_len > 0 {
                self.sum_qual as f64 / self.total_len as f64
            } else {
                0.0
            }
        )?;
        writeln!(out, "SN\tinsert size average:\t{avg_isize:.1}")?;
        writeln!(out, "SN\tinsert size standard deviation:\t{sd_isize:.1}")?;
        writeln!(out, "SN\tinward oriented pairs:\t{nisize_inward}")?;
        writeln!(out, "SN\toutward oriented pairs:\t{nisize_outward}")?;
        writeln!(out, "SN\tpairs with other orientation:\t{nisize_other}")?;
        writeln!(
            out,
            "SN\tpairs on different chromosomes:\t{}",
            self.nreads_anomalous / 2
        )?;
        writeln!(
            out,
            "SN\tpercentage of properly paired reads (%):\t{:.1}",
            ratio(100 * self.nreads_properly_paired, nreads)
        )?;

        if self.target_count > 0 {
            writeln!(out, "SN\tbases inside the target:\t{}", self.target_count)?;

            let threshold = usize::try_from(info.cov_threshold).unwrap_or(0);
            let cov_sum: u64 = self.cov.iter().skip(threshold + 1).sum();

            writeln!(
                out,
                "SN\tpercentage of target genome with coverage > {} (%):\t{:.2}",
                info.cov_threshold,
                f64::from((100 * cov_sum) as f32 / self.target_count as f32)
            )?;
        }

        Ok(())
    }

    /// Includes the longest read and the highest quality in the per-cycle tables, as samtools
    /// does after the summary.
    fn round_up_maxima(&mut self) {
        if self.max_len < self.nbases {
            self.max_len += 1;
        }

        if self.max_qual + 1 < QUALITY_BINS && self.max_qual < 255 {
            self.max_qual += 1;
        }
    }

    /// Base qualities per cycle (FFQ, LFQ) and mismatches per cycle (MPC).
    fn write_quality_tables<W: Write>(&self, out: &mut W) -> io::Result<()> {
        writeln!(
            out,
            "# First Fragment Qualities. Use `grep ^FFQ | cut -f 2-` to extract this part."
        )?;
        writeln!(
            out,
            "# Columns correspond to qualities and rows to cycles. First column is the cycle number."
        )?;
        write_quality_rows(
            out,
            "FFQ",
            &self.quals_1st,
            self.quals_rows,
            self.max_len_1st,
            self.max_qual,
        )?;

        writeln!(
            out,
            "# Last Fragment Qualities. Use `grep ^LFQ | cut -f 2-` to extract this part."
        )?;
        writeln!(
            out,
            "# Columns correspond to qualities and rows to cycles. First column is the cycle number."
        )?;
        write_quality_rows(
            out,
            "LFQ",
            &self.quals_2nd,
            self.quals_rows,
            self.max_len_2nd,
            self.max_qual,
        )?;

        if let Some(mpc) = &self.mpc_buf {
            writeln!(
                out,
                "# Mismatches per cycle and quality. Use `grep ^MPC | cut -f 2-` to extract this part."
            )?;
            writeln!(
                out,
                "# Columns correspond to qualities, rows to cycles. First column is the cycle number, second"
            )?;
            writeln!(
                out,
                "# is the number of N's and the rest is the number of mismatches"
            )?;

            write_quality_rows(out, "MPC", mpc, QUALITY_BINS, self.max_len, self.max_qual)?;
        }

        Ok(())
    }

    /// GC content (GCF, GCL) and base content per cycle (GCC, GCT, and by fragment).
    fn write_base_content<W: Write>(&self, out: &mut W) -> io::Result<()> {
        writeln!(
            out,
            "# GC Content of first fragments. Use `grep ^GCF | cut -f 2-` to extract this part."
        )?;
        write_gc_rows(out, "GCF", &self.gc_1st, self.ngc)?;
        writeln!(
            out,
            "# GC Content of last fragments. Use `grep ^GCL | cut -f 2-` to extract this part."
        )?;
        write_gc_rows(out, "GCL", &self.gc_2nd, self.ngc)?;

        writeln!(
            out,
            "# ACGT content per cycle. Use `grep ^GCC | cut -f 2-` to extract this part. The columns are: cycle; A,C,G,T base counts as a percentage of all A/C/G/T bases [%]; and N and O counts as a percentage of all A/C/G/T bases [%]"
        )?;

        for ibase in 0..self.max_len {
            let (c1, c2) = (
                BaseCounts::from(&self.acgtno_cycles_1st[ibase]),
                BaseCounts::from(&self.acgtno_cycles_2nd[ibase]),
            );
            let sum = c1.a + c1.c + c1.g + c1.t + c2.a + c2.c + c2.g + c2.t;

            if sum == 0 {
                continue;
            }

            write_percent_row(
                out,
                "GCC",
                ibase + 1,
                &[
                    c1.a + c2.a,
                    c1.c + c2.c,
                    c1.g + c2.g,
                    c1.t + c2.t,
                    c1.n + c2.n,
                    c1.other + c2.other,
                ],
                sum,
            )?;
        }

        writeln!(
            out,
            "# ACGT content per cycle, read oriented. Use `grep ^GCT | cut -f 2-` to extract this part. The columns are: cycle; A,C,G,T base counts as a percentage of all A/C/G/T bases [%]"
        )?;

        for ibase in 0..self.max_len {
            let c = BaseCounts::from(&self.acgtno_revcomp[ibase]);
            let sum = c.a + c.c + c.g + c.t;

            if sum == 0 {
                continue;
            }

            write_percent_row(out, "GCT", ibase + 1, &[c.a, c.c, c.g, c.t], sum)?;
        }

        for (label, total_label, ordinal, cycles) in [
            ("FBC", "FTC", "first", &self.acgtno_cycles_1st),
            ("LBC", "LTC", "last", &self.acgtno_cycles_2nd),
        ] {
            writeln!(
                out,
                "# ACGT content per cycle for {ordinal} fragments. Use `grep ^{label} | cut -f 2-` to extract this part. The columns are: cycle; A,C,G,T base counts as a percentage of all A/C/G/T bases [%]; and N and O counts as a percentage of all A/C/G/T bases [%]"
            )?;

            let mut totals = BaseCounts::default();

            for (ibase, c) in cycles.iter().take(self.max_len).enumerate() {
                let c = BaseCounts::from(c);
                let sum = c.a + c.c + c.g + c.t;
                totals.a += c.a;
                totals.c += c.c;
                totals.g += c.g;
                totals.t += c.t;
                totals.n += c.n;

                if sum > 0 {
                    write_percent_row(
                        out,
                        label,
                        ibase + 1,
                        &[c.a, c.c, c.g, c.t, c.n, c.other],
                        sum,
                    )?;
                }
            }

            writeln!(
                out,
                "# ACGT raw counters for {ordinal} fragments. Use `grep ^{total_label} | cut -f 2-` to extract this part. The columns are: A,C,G,T,N base counters"
            )?;
            writeln!(
                out,
                "{total_label}\t{}\t{}\t{}\t{}\t{}",
                totals.a, totals.c, totals.g, totals.t, totals.n
            )?;
        }

        Ok(())
    }

    /// Barcode base content and qualities per cycle.
    fn write_barcodes<W: Write>(&self, out: &mut W) -> io::Result<()> {
        for barcode in &self.barcodes {
            if barcode.len == 0 {
                continue;
            }

            let (tag, quality_tag) = (
                String::from_utf8_lossy(&barcode.tag),
                String::from_utf8_lossy(&barcode.quality_tag),
            );

            writeln!(
                out,
                "# ACGT content per cycle for barcodes. Use `grep ^{tag}C | cut -f 2-` to extract this part. The columns are: cycle; A,C,G,T base counts as a percentage of all A/C/G/T bases [%]; and N counts as a percentage of all A/C/G/T bases [%]"
            )?;

            for (i, c) in barcode.counts.iter().enumerate() {
                if barcode.separator == Some(i) {
                    continue;
                }

                let sum = c.a + c.c + c.g + c.t;

                if sum > 0 {
                    let (segment, cycle) = barcode.segment(i);
                    let label = format!("{tag}C{segment}");
                    write_percent_row(out, &label, cycle, &[c.a, c.c, c.g, c.t, c.n], sum)?;
                }
            }

            writeln!(
                out,
                "# Barcode Qualities. Use `grep ^{quality_tag}Q | cut -f 2-` to extract this part."
            )?;
            writeln!(
                out,
                "# Columns correspond to qualities and rows to barcode cycles. First column is the cycle number."
            )?;

            for i in 0..barcode.len {
                if barcode.separator == Some(i) {
                    continue;
                }

                let (segment, cycle) = barcode.segment(i);
                write!(out, "{quality_tag}Q{segment}\t{cycle}")?;

                if let Some(max_qual) = barcode.max_qual {
                    for &n in &barcode.quals[i * QUALITY_BINS..=i * QUALITY_BINS + max_qual] {
                        write!(out, "\t{n}")?;
                    }
                }

                writeln!(out)?;
            }
        }

        Ok(())
    }

    /// Insert sizes (IS) below `ibulk`, read lengths (RL, FRL, LRL), mapping qualities (MAPQ), and indels (ID, IC).
    fn write_distributions<W: Write>(
        &self,
        out: &mut W,
        info: &Info,
        ibulk: i32,
    ) -> io::Result<()> {
        writeln!(
            out,
            "# Insert sizes. Use `grep ^IS | cut -f 2-` to extract this part. The columns are: insert size, pairs total, inward oriented pairs, outward oriented pairs, other pairs"
        )?;

        for isize in 0..ibulk {
            let [i, o, x] = self.isize.get(isize);

            if info.insert_size_rows == InsertSizeRows::All || i + o + x > 0 {
                writeln!(out, "IS\t{isize}\t{}\t{i}\t{o}\t{x}", i + o + x)?;
            }
        }

        for (label, description, lengths, max) in [
            ("RL", "", &self.read_lengths, self.max_len),
            (
                "FRL",
                " - first fragments",
                &self.read_lengths_1st,
                self.max_len_1st,
            ),
            (
                "LRL",
                " - last fragments",
                &self.read_lengths_2nd,
                self.max_len_2nd,
            ),
        ] {
            writeln!(
                out,
                "# Read lengths{description}. Use `grep ^{label} | cut -f 2-` to extract this part. The columns are: read length, count"
            )?;

            for ilen in 0..max {
                if lengths[ilen + 1] > 0 {
                    writeln!(out, "{label}\t{}\t{}", ilen + 1, lengths[ilen + 1])?;
                }
            }
        }

        writeln!(
            out,
            "# Mapping qualities for reads !(UNMAP|SECOND|SUPPL|QCFAIL|DUP). Use `grep ^MAPQ | cut -f 2-` to extract this part. The columns are: mapq, count"
        )?;

        for (mapq, &n) in self.mapping_qualities.iter().enumerate() {
            if n > 0 {
                writeln!(out, "MAPQ\t{mapq}\t{n}")?;
            }
        }

        writeln!(
            out,
            "# Indel distribution. Use `grep ^ID | cut -f 2-` to extract this part. The columns are: length, number of insertions, number of deletions"
        )?;

        for ilen in 0..self.nindels {
            if self.insertions[ilen] > 0 || self.deletions[ilen] > 0 {
                writeln!(
                    out,
                    "ID\t{}\t{}\t{}",
                    ilen + 1,
                    self.insertions[ilen],
                    self.deletions[ilen]
                )?;
            }
        }

        writeln!(
            out,
            "# Indels per cycle. Use `grep ^IC | cut -f 2-` to extract this part. The columns are: cycle, number of insertions (fwd), .. (rev) , number of deletions (fwd), .. (rev)"
        )?;

        for ilen in 0..=self.nbases {
            let counts = [
                self.ins_cycles_1st[ilen],
                self.ins_cycles_2nd[ilen],
                self.del_cycles_1st[ilen],
                self.del_cycles_2nd[ilen],
            ];

            if counts.iter().any(|&n| n > 0) {
                writeln!(
                    out,
                    "IC\t{}\t{}\t{}\t{}\t{}",
                    ilen + 1,
                    counts[0],
                    counts[1],
                    counts[2],
                    counts[3]
                )?;
            }
        }

        Ok(())
    }

    /// Reference statistics (RFS).
    fn write_reference_stats<W: Write>(&self, out: &mut W) -> io::Result<()> {
        if let Some(rs) = &self.reference_stats {
            writeln!(
                out,
                "# Reference statistics. Use `grep ^RFS | cut -f 2-` to extract this part."
            )?;
            writeln!(
                out,
                "# Total count, Output count, Average GC, Min length, Max length, Average length, Total length in first row."
            )?;
            writeln!(
                out,
                "# Sequence name, Length, GC content, Unknown count in following rows."
            )?;
            writeln!(
                out,
                "RFS\t{}\t{}\t{:.2}\t{}\t{}\t{:.2}\t{}",
                rs.total_count,
                rs.count,
                f64::from(rs.average_gc),
                rs.min_length,
                rs.max_length,
                f64::from(rs.average_length),
                rs.combined_length
            )?;

            for (name, length, gc, n) in rs.rows.iter().flatten() {
                writeln!(out, "RFS\t{name}\t{length}\t{:.2}\t{n}", f64::from(*gc))?;
            }
        }

        Ok(())
    }

    pub(super) fn write_coverage<W: Write>(&self, out: &mut W, info: &Info) -> io::Result<()> {
        writeln!(
            out,
            "# Coverage distribution. Use `grep ^COV | cut -f 2-` to extract this part."
        )?;

        let (min, step) = (info.cov_min, info.cov_step);

        if self.cov[0] > 0 {
            writeln!(out, "COV\t[<{min}]\t{}\t{}", min - 1, self.cov[0])?;
        }

        for icov in 1..self.ncov - 1 {
            if self.cov[icov] > 0 {
                let icov_i32 = i32::try_from(icov).unwrap_or(i32::MAX);
                let lo = min + (icov_i32 - 1) * step;
                let hi = min + icov_i32 * step - 1;
                writeln!(out, "COV\t[{lo}-{hi}]\t{hi}\t{}", self.cov[icov])?;
            }
        }

        if self.cov[self.ncov - 1] > 0 {
            let last = min + (i32::try_from(self.ncov).unwrap_or(i32::MAX) - 2) * step - 1;
            writeln!(out, "COV\t[{last}<]\t{last}\t{}", self.cov[self.ncov - 1])?;
        }

        Ok(())
    }

    pub(super) fn write_gc_depth<W: Write>(
        &mut self,
        out: &mut W,
        info: &Info,
        avg_read_length: f32,
    ) -> io::Result<()> {
        writeln!(
            out,
            "# GC-depth. Use `grep ^GCD | cut -f 2-` to extract this part. The columns are: GC%, unique sequence percentiles, 10th, 25th, 50th, 75th and 90th depth percentile"
        )?;

        // As in samtools, the last bin keeps its unrounded GC, and is sorted but not reported.
        let reported = self.igcd.min(self.gcd.len());

        for bin in &mut self.gcd[..reported] {
            if info.reference.is_some() {
                bin.gc = (100.0 * f64::from(bin.gc)).round_ties_even() as f32;
            } else if bin.depth > 0 {
                bin.gc =
                    (100.0 * f64::from(bin.gc) / f64::from(bin.depth)).round_ties_even() as f32;
            }
        }

        if !self.gcd.is_empty() {
            self.gcd[..=self.igcd].sort_by(|a, b| {
                a.gc.partial_cmp(&b.gc)
                    .unwrap_or(std::cmp::Ordering::Equal)
                    .then(a.depth.cmp(&b.depth))
            });
        }

        let bin_size = info.gcd_bin_size as f32;
        let mut igcd = 0;

        while igcd < self.igcd {
            let gc = self.gcd[igcd].gc;
            let mut nbins = 0;

            while igcd + nbins < self.igcd && f64::from(self.gcd[igcd + nbins].gc - gc).abs() < 0.1
            {
                nbins += 1;
            }

            let bins = &self.gcd[igcd..];
            let p = |percentile: u32| {
                f64::from(gcd_percentile(bins, nbins, percentile) * avg_read_length / bin_size)
            };

            writeln!(
                out,
                "GCD\t{:.1}\t{:.3}\t{:.3}\t{:.3}\t{:.3}\t{:.3}\t{:.3}",
                f64::from(gc),
                ((igcd + nbins + 1) as u32) as f64 * 100.0 / (self.igcd + 1) as f64,
                p(10),
                p(25),
                p(50),
                p(75),
                p(90)
            )?;

            igcd += nbins;
        }

        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_format_e() {
        assert_eq!(format_e(0.0), "0.000000e+00");
        assert_eq!(format_e(0.001_234_567_8), "1.234568e-03");
        assert_eq!(format_e(12345.0), "1.234500e+04");
    }
}
