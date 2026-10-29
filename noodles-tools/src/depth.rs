//! Per-position read depth over coordinate-sorted alignments, as samtools depth counts it.

use std::{collections::VecDeque, io};

use noodles_sam::{
    self as sam,
    alignment::record::{
        Flags,
        cigar::{Op, op::Kind},
    },
};

use crate::alignment::Rec;

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
    pub quality_sum: u64,
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
    /// Reused per record.
    ops: Vec<Op>,
    quality_scores: Vec<u8>,
}

impl DepthWindow {
    pub fn new(options: CountOptions) -> Self {
        Self {
            options,
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
        let r = record.as_dyn();

        let (Some(id), Some(start)) = (
            r.reference_sequence_id(header).transpose()?,
            r.alignment_start().transpose()?,
        ) else {
            return Ok(());
        };

        let start = usize::from(start);

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

        self.ops.clear();
        match record {
            Rec::Bam(r) => {
                for result in r.cigar().iter() {
                    self.ops.push(result?);
                }
            }
            _ => {
                for result in r.cigar().iter() {
                    self.ops.push(result?);
                }
            }
        }

        let CountOptions {
            min_base_quality,
            include_deletions,
            track_quality,
            include_skips,
        } = self.options;

        let needs_quality = min_base_quality > 0 || track_quality;

        self.quality_scores.clear();

        if needs_quality {
            match record {
                Rec::Bam(r) => self
                    .quality_scores
                    .extend_from_slice(r.quality_scores().as_ref()),
                Rec::Buf(r) => self
                    .quality_scores
                    .extend_from_slice(r.quality_scores().as_ref()),
                _ => {
                    for result in r.quality_scores().iter() {
                        self.quality_scores.push(result?);
                    }
                }
            }
        }

        // Size the window once for the record's whole reference span.
        let span: usize = self
            .ops
            .iter()
            .filter(|op| op.kind().consumes_reference())
            .map(|op| op.len())
            .sum();
        let needed = start + span - self.start;
        if needed > self.cells.len() {
            self.cells.resize(needed, Cell::default());
        }

        // The reusable buffers are moved out while cells are updated, then put back.
        let ops = std::mem::take(&mut self.ops);
        let quality_scores = std::mem::take(&mut self.quality_scores);

        let mut reference_position = start;
        let mut read_position = 0;

        for op in &ops {
            let len = op.len();
            let offset = reference_position - self.start;

            match op.kind() {
                Kind::Match | Kind::SequenceMatch | Kind::SequenceMismatch => {
                    let cells = self.cells.range_mut(offset..offset + len);

                    if needs_quality {
                        for (i, cell) in cells.enumerate() {
                            cell.touched = true;

                            // Missing qualities (0xff) count, as in samtools.
                            let q = match quality_scores.get(read_position + i) {
                                Some(&q) if q != 0xff => q,
                                _ => {
                                    cell.depth += 1;
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
                Kind::Deletion => {
                    for cell in self.cells.range_mut(offset..offset + len) {
                        cell.touched = true;

                        if include_deletions {
                            cell.depth += 1;
                        }
                    }

                    reference_position += len;
                }
                Kind::Skip => {
                    if include_skips {
                        for cell in self.cells.range_mut(offset..offset + len) {
                            cell.depth += 1;
                        }
                    }

                    reference_position += len;
                }
                Kind::Insertion | Kind::SoftClip => read_position += len,
                Kind::HardClip | Kind::Pad => {}
            }
        }

        self.ops = ops;
        self.quality_scores = quality_scores;

        Ok(())
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
