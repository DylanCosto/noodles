//! Per-position read depth over coordinate-sorted alignments, as samtools depth counts it.

use std::{
    cmp::Reverse,
    collections::{BinaryHeap, VecDeque},
    io,
};

use noodles_sam::{
    self as sam,
    alignment::record::{Flags, cigar::op::Kind},
};

use crate::{alignment::Rec, bam_raw::Record};

// BAM CIGAR operation codes.
const CIGAR_MATCH: u32 = 0;
const CIGAR_INSERTION: u32 = 1;
const CIGAR_DELETION: u32 = 2;
const CIGAR_SKIP: u32 = 3;
const CIGAR_SOFT_CLIP: u32 = 4;
const CIGAR_SEQUENCE_MATCH: u32 = 7;
const CIGAR_SEQUENCE_MISMATCH: u32 = 8;

/// samtools' default filter for depth and coverage.
pub const DEFAULT_EXCLUDE_FLAGS: Flags = Flags::UNMAPPED
    .union(Flags::SECONDARY)
    .union(Flags::QC_FAIL)
    .union(Flags::DUPLICATE);

/// Which bases count toward depth.
#[derive(Clone, Copy, Debug, Default)]
pub struct CountOptions {
    /// Bases below this quality are not counted (when positive).
    pub min_base_quality: u8,
    /// Count deletions (`-J`).
    pub include_deletions: bool,
    /// Sum base qualities into [`Cell::quality_sum`] (for coverage).
    pub track_quality: bool,
    /// Count reference skips (`N`), as samtools bedcov does by default.
    pub include_skips: bool,
}

/// Depth at one position, and the summed base qualities of the bases counted there.
#[derive(Clone, Copy, Debug, Default, Eq, PartialEq)]
pub struct Cell {
    pub depth: u32,
    /// Deletions and reference skips at this position, whether or not `depth` counts them.
    pub gaps: u32,
    pub quality_sum: u64,
    /// Counted bases without a quality (no sequence, or `0xff`), so not in `quality_sum`.
    pub missing_quality: u32,
    /// Whether any read aligns here, counted or not (e.g., a filtered base or a deletion).
    /// samtools reports touched positions even at depth 0.
    pub touched: bool,
}

/// A sliding window of depths for one reference sequence at a time.
///
/// Records must arrive sorted by reference sequence and start, so positions before a record's
/// start are final and can be emitted.
#[derive(Debug, Default)]
pub struct DepthWindow {
    reference_sequence_id: Option<usize>,
    /// 1-based position of `cells[0]`.
    start: usize,
    cells: VecDeque<Cell>,
    options: CountOptions,
    /// Reused per record: the CIGAR, BAM-encoded (length << 4 | operation).
    ops: Vec<u32>,
    quality_scores: Vec<u8>,
    cap: Option<PileupCap>,
}

impl DepthWindow {
    pub fn new(options: CountOptions) -> Self {
        Self {
            options,
            ..Self::default()
        }
    }

    /// A window that leaves out reads as HTSlib's pileup does with a maximum count (`maxcnt`,
    /// samtools coverage's `-d`).
    pub fn with_max_count(options: CountOptions, max_count: usize) -> Self {
        Self {
            options,
            cap: Some(PileupCap::new(max_count)),
            ..Self::default()
        }
    }

    /// Adds a record's aligned bases, first emitting every position before its start.
    ///
    /// `emit` receives `(reference sequence ID, position, cell)` for every touched position.
    pub fn add<F>(&mut self, record: &Rec<'_>, header: &sam::Header, emit: &mut F) -> io::Result<()>
    where
        F: FnMut(usize, usize, Cell) -> io::Result<()>,
    {
        self.add_clipped(record, header, 0, emit)
    }

    /// [`Self::add`], leaving out the aligned bases (and counted deletions) before the 1-based
    /// position `clip`, where they overlap the read's mate (samtools depth `-s`).
    pub fn add_clipped<F>(
        &mut self,
        record: &Rec<'_>,
        header: &sam::Header,
        clip: usize,
        emit: &mut F,
    ) -> io::Result<()>
    where
        F: FnMut(usize, usize, Cell) -> io::Result<()>,
    {
        // BAM records are read raw; others through the record trait.
        let raw = match record {
            Rec::Bam(b) => Some(Record::new(b.as_ref())?),
            _ => None,
        };

        let Some((id, start)) = record_position(record, raw.as_ref(), header)? else {
            return Ok(());
        };

        if self.reference_sequence_id != Some(id) {
            self.finish(emit)?;
            self.reference_sequence_id = Some(id);
            self.start = start;
        } else if start < self.start {
            return Err(io::Error::new(
                io::ErrorKind::InvalidData,
                "input is not coordinate-sorted",
            ));
        }

        self.flush_before(start, emit)?;

        let needs_quality = self.load_record(record, raw.as_ref())?;

        let span: usize = self
            .ops
            .iter()
            .filter(|&&op| consumes_reference(op & 0xf))
            .map(|&op| (op >> 4) as usize)
            .sum();

        if let Some(cap) = &mut self.cap
            && !cap.admit(id, start - 1, start - 1 + span)
        {
            return Ok(());
        }

        self.count_bases(start, span, needs_quality, clip);

        Ok(())
    }

    /// Reads a record's CIGAR, BAM-encoded, and its qualities if counting needs them; returns
    /// whether it does.
    fn load_record(&mut self, record: &Rec<'_>, raw: Option<&Record<'_>>) -> io::Result<bool> {
        let r = record.as_dyn();

        match raw {
            Some(raw) => raw.read_cigar(&mut self.ops)?,
            None => {
                self.ops.clear();

                for result in r.cigar().iter() {
                    let op = result?;
                    let len = u32::try_from(op.len())
                        .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))?;
                    self.ops.push(len << 4 | kind_code(op.kind()));
                }
            }
        }

        let needs_quality = self.options.min_base_quality > 0 || self.options.track_quality;

        self.quality_scores.clear();

        if needs_quality {
            match (&raw, record) {
                (Some(raw), _) => self.quality_scores.extend_from_slice(raw.quality_scores()),
                (None, Rec::Buf(r)) => self
                    .quality_scores
                    .extend_from_slice(r.quality_scores().as_ref()),
                (None, _) => {
                    for result in r.quality_scores().iter() {
                        self.quality_scores.push(result?);
                    }
                }
            }
        }

        Ok(needs_quality)
    }

    /// Counts the bases of the record in [`Self::ops`] and [`Self::quality_scores`], which
    /// starts at the 1-based `start` and spans `span` reference bases.
    fn count_bases(&mut self, start: usize, span: usize, needs_quality: bool, clip: usize) {
        let CountOptions {
            min_base_quality,
            include_deletions,
            include_skips,
            ..
        } = self.options;

        // Size the window once for the record's whole reference span.
        let needed = start + span - self.start;
        if needed > self.cells.len() {
            self.cells.resize(needed, Cell::default());
        }

        // The reusable buffers are moved out while cells are updated, then put back.
        let ops = std::mem::take(&mut self.ops);
        let quality_scores = std::mem::take(&mut self.quality_scores);

        // One slice for the whole record (rotating the ring only when it has wrapped), so the
        // loops below run over plain slices.
        let window = self.cells.make_contiguous();

        let mut reference_position = start;
        let mut read_position = 0;

        for &op in &ops {
            let len = (op >> 4) as usize;
            let offset = reference_position - self.start;
            // Bases before `clip` overlap the mate and are not counted.
            let clipped = clip.saturating_sub(reference_position).min(len);

            match op & 0xf {
                CIGAR_MATCH | CIGAR_SEQUENCE_MATCH | CIGAR_SEQUENCE_MISMATCH => {
                    // Clipped positions are still reported, at the depth of the other reads.
                    for cell in &mut window[offset..offset + clipped] {
                        cell.touched = true;
                    }

                    let cells = &mut window[offset + clipped..offset + len];
                    let first_base = read_position + clipped;

                    if needs_quality {
                        for (i, cell) in cells.iter_mut().enumerate() {
                            cell.touched = true;

                            // Missing qualities (0xff) count, as in samtools.
                            let q = match quality_scores.get(first_base + i) {
                                Some(&q) if q != 0xff => q,
                                _ => {
                                    cell.depth += 1;
                                    cell.missing_quality += 1;
                                    continue;
                                }
                            };

                            if q < min_base_quality {
                                continue;
                            }

                            cell.depth += 1;
                            cell.quality_sum += u64::from(q);
                        }
                    } else {
                        for cell in cells {
                            cell.touched = true;
                            cell.depth += 1;
                        }
                    }

                    reference_position += len;
                    read_position += len;
                }
                CIGAR_DELETION if include_deletions => {
                    for cell in &mut window[offset..offset + clipped] {
                        cell.touched = true;
                    }

                    // As in samtools, a deletion takes the quality of the base after it.
                    let passes = !needs_quality
                        || quality_scores
                            .get(read_position)
                            .is_none_or(|&q| q == 0xff || q >= min_base_quality);

                    for cell in &mut window[offset + clipped..offset + len] {
                        cell.touched = true;
                        cell.gaps += 1;

                        if passes {
                            cell.depth += 1;
                        }
                    }

                    reference_position += len;
                }
                CIGAR_DELETION => {
                    for cell in &mut window[offset..offset + len] {
                        cell.touched = true;
                        cell.gaps += 1;
                    }

                    reference_position += len;
                }
                CIGAR_SKIP => {
                    // Reference skips are within the read's span, so they are reported too.
                    for cell in &mut window[offset..offset + len] {
                        cell.touched = true;
                        cell.gaps += 1;

                        if include_skips {
                            cell.depth += 1;
                        }
                    }

                    reference_position += len;
                }
                CIGAR_INSERTION | CIGAR_SOFT_CLIP => read_position += len,
                _ => {}
            }
        }

        self.ops = ops;
        self.quality_scores = quality_scores;
    }

    fn flush_before<F>(&mut self, position: usize, emit: &mut F) -> io::Result<()>
    where
        F: FnMut(usize, usize, Cell) -> io::Result<()>,
    {
        let Some(id) = self.reference_sequence_id else {
            return Ok(());
        };

        while self.start < position {
            match self.cells.pop_front() {
                Some(cell) => {
                    if cell.touched {
                        emit(id, self.start, cell)?;
                    }
                    self.start += 1;
                }
                None => {
                    self.start = position;
                    break;
                }
            }
        }

        Ok(())
    }

    /// Emits the remaining positions of the current reference sequence.
    pub fn finish<F>(&mut self, emit: &mut F) -> io::Result<()>
    where
        F: FnMut(usize, usize, Cell) -> io::Result<()>,
    {
        self.flush_before(usize::MAX, emit)?;
        self.reference_sequence_id = None;
        Ok(())
    }
}

/// HTSlib's pileup depth limit (`maxcnt`) for one input: a read that starts where the previous
/// one did is left out once `maxcnt` reads are buffered. The buffer holds the reads that have
/// not ended before the position just before the new read's start.
#[derive(Clone, Debug)]
struct PileupCap {
    max_count: usize,
    /// The pileup's reference sequence and 0-based position: the last added read's start.
    id: usize,
    position: usize,
    /// The 0-based exclusive ends of the buffered reads.
    ends: BinaryHeap<Reverse<usize>>,
}

impl PileupCap {
    fn new(max_count: usize) -> Self {
        Self {
            max_count,
            id: 0,
            position: 0,
            ends: BinaryHeap::new(),
        }
    }

    /// Returns whether a read at the 0-based `start` with the exclusive `end` enters the pileup
    /// (`bam_plp_push`).
    fn admit(&mut self, id: usize, start: usize, end: usize) -> bool {
        let new_reference = id != self.id;

        if new_reference {
            self.ends.clear();
        } else {
            // The pileup has dropped the reads that end before `start - 1`.
            while self.ends.peek().is_some_and(|&Reverse(e)| e < start) {
                self.ends.pop();
            }

            if start == self.position && self.ends.len() >= self.max_count {
                return false;
            }
        }

        // An empty read starting at the pileup's position is not buffered.
        if new_reference || end > self.position {
            self.ends.push(Reverse(end));
        }

        self.id = id;
        self.position = start;

        true
    }
}

/// A record's reference sequence ID and 1-based start, or `None` if it is unplaced; `raw` is
/// the record's raw view if it is BAM.
fn record_position(
    record: &Rec<'_>,
    raw: Option<&Record<'_>>,
    header: &sam::Header,
) -> io::Result<Option<(usize, usize)>> {
    let r = record.as_dyn();

    let position = match raw {
        Some(raw) => {
            let (Ok(id), Ok(position)) = (
                usize::try_from(raw.reference_sequence_id()),
                usize::try_from(raw.position()),
            ) else {
                return Ok(None);
            };

            (id, position + 1)
        }
        None => {
            let (Some(id), Some(start)) = (
                r.reference_sequence_id(header).transpose()?,
                r.alignment_start().transpose()?,
            ) else {
                return Ok(None);
            };

            (id, usize::from(start))
        }
    };

    Ok(Some(position))
}

/// samtools depth's `qlen_used`: the sequence length without leading and trailing soft clips,
/// or, without a sequence, the query length of the CIGAR.
pub fn used_query_length(record: &Rec<'_>) -> io::Result<usize> {
    let r = record.as_dyn();
    let ops = r.cigar().iter().collect::<io::Result<Vec<_>>>()?;
    let sequence_len = r.sequence().len();

    if sequence_len == 0 {
        return Ok(query_length(&ops));
    }

    let leading = ops
        .iter()
        .take_while(|op| op.kind() == Kind::SoftClip)
        .count();
    let soft_clipped: usize = ops[..leading]
        .iter()
        .chain(
            ops[leading..]
                .iter()
                .rev()
                .take_while(|op| op.kind() == Kind::SoftClip),
        )
        .map(|op| op.len())
        .sum();

    Ok(sequence_len.saturating_sub(soft_clipped))
}

/// HTSlib's `bam_cigar2qlen`: the query length of a record's CIGAR, soft clips included.
pub fn cigar_query_length(record: &Rec<'_>) -> io::Result<usize> {
    let ops = record
        .as_dyn()
        .cigar()
        .iter()
        .collect::<io::Result<Vec<_>>>()?;
    Ok(query_length(&ops))
}

fn query_length(ops: &[noodles_sam::alignment::record::cigar::Op]) -> usize {
    ops.iter()
        .filter(|op| op.kind().consumes_read())
        .map(|op| op.len())
        .sum()
}

/// The BAM code of a CIGAR operation.
fn kind_code(kind: Kind) -> u32 {
    match kind {
        Kind::Match => CIGAR_MATCH,
        Kind::Insertion => CIGAR_INSERTION,
        Kind::Deletion => CIGAR_DELETION,
        Kind::Skip => CIGAR_SKIP,
        Kind::SoftClip => CIGAR_SOFT_CLIP,
        Kind::HardClip => 5,
        Kind::Pad => 6,
        Kind::SequenceMatch => CIGAR_SEQUENCE_MATCH,
        Kind::SequenceMismatch => CIGAR_SEQUENCE_MISMATCH,
    }
}

/// Whether a BAM CIGAR operation code consumes the reference.
fn consumes_reference(code: u32) -> bool {
    matches!(
        code,
        CIGAR_MATCH | CIGAR_DELETION | CIGAR_SKIP | CIGAR_SEQUENCE_MATCH | CIGAR_SEQUENCE_MISMATCH
    )
}

/// Formats like C's `printf("%.*g", precision, value)`.
pub fn format_g(value: f64, precision: usize) -> String {
    let precision = precision.max(1);

    if value == 0.0 {
        return "0".to_owned();
    }

    if !value.is_finite() {
        return if value.is_nan() {
            "nan".to_owned()
        } else if value > 0.0 {
            "inf".to_owned()
        } else {
            "-inf".to_owned()
        };
    }

    // The exponent after rounding to `precision` significant digits decides the style.
    let scientific = format!("{:.*e}", precision - 1, value);
    let (mantissa, exponent) = scientific
        .split_once('e')
        .expect("scientific notation has an exponent");
    let exponent: i32 = exponent.parse().expect("valid exponent");

    let precision_i32 = i32::try_from(precision).unwrap_or(i32::MAX);

    if exponent < -4 || exponent >= precision_i32 {
        let mantissa = trim_fraction(mantissa);
        let sign = if exponent < 0 { '-' } else { '+' };
        format!("{mantissa}e{sign}{:02}", exponent.abs())
    } else {
        let decimals = usize::try_from(precision_i32 - 1 - exponent).unwrap_or(0);
        trim_fraction(&format!("{value:.decimals$}")).to_owned()
    }
}

fn trim_fraction(s: &str) -> &str {
    if s.contains('.') {
        s.trim_end_matches('0').trim_end_matches('.')
    } else {
        s
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_pileup_cap() {
        let mut cap = PileupCap::new(2);

        assert!(cap.admit(0, 10, 20));
        assert!(cap.admit(0, 10, 20));
        // Two reads are buffered at the pileup's position.
        assert!(!cap.admit(0, 10, 20));
        // The first read at a new position is always added.
        assert!(cap.admit(0, 11, 12));
        assert!(!cap.admit(0, 11, 30));
        // The reads before have ended; an empty read is buffered only when it comes first.
        assert!(cap.admit(0, 40, 40));
        assert!(cap.admit(0, 40, 40));
        assert!(cap.admit(0, 40, 50));
        assert!(!cap.admit(0, 40, 50));
        assert!(cap.admit(1, 40, 50));
    }

    #[test]
    fn test_format_g() {
        let cases = [
            (50.001_775, 6, "50.0018"),
            (17.718_2, 6, "17.7182"),
            (25.0, 3, "25"),
            (35.2, 3, "35.2"),
            (0.652_32, 6, "0.65232"),
            (100.0, 6, "100"),
            (0.000_012_34, 6, "1.234e-05"),
            (1_234_567.0, 6, "1.23457e+06"),
            (99.999_99, 3, "100"),
            (0.0, 6, "0"),
        ];

        for (value, precision, expected) in cases {
            assert_eq!(
                format_g(value, precision),
                expected,
                "{value} with {precision}"
            );
        }
    }
}
