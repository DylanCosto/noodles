//! samtools coverage's text histograms (`-m`, `-D`), drawn as samtools draws them.

use std::io::{self, Write};

use super::Tally;
use crate::depth::format_g;

/// Rows of the histogram.
const ROWS: u8 = 10;

/// What each histogram bin shows.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub enum Plot {
    /// The percentage of the bin's bases that are covered (`-m`).
    Breadth,
    /// The mean depth in the bin (`-D`).
    Depth,
}

/// The characters bars are drawn with.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub enum Blocks {
    /// Eighth blocks and box-drawing lines.
    Utf8,
    /// `.` and `:` (`-A`).
    Ascii,
}

impl Blocks {
    fn bars(self) -> &'static [&'static str] {
        match self {
            Self::Utf8 => &["▁", "▂", "▃", "▄", "▅", "▆", "▇", "█"],
            Self::Ascii => &[".", ":"],
        }
    }

    fn vertical_line(self) -> &'static str {
        match self {
            Self::Utf8 => "│",
            Self::Ascii => "|",
        }
    }
}

/// Writes one reference sequence's histogram (samtools' `print_hist`).
pub(super) fn write<W: Write>(
    out: &mut W,
    name: &[u8],
    length: usize,
    tally: &Tally,
    plot: Plot,
    blocks: Blocks,
) -> io::Result<()> {
    let scale = match plot {
        Plot::Breadth => 100.0,
        Plot::Depth => 1.0,
    };

    let bin_width = tally.bin_width as f64;
    let data: Vec<f64> = tally
        .bins
        .iter()
        .map(|&n| scale * n as f64 / bin_width)
        .collect();
    let max = data.iter().copied().fold(0.0, f64::max);

    out.write_all(name)?;
    writeln!(out, " ({}bp)", readable_bps(length as f64))?;

    let bars = blocks.bars();
    let bar_levels = bars.len() as f64;
    let row_height = max / f64::from(ROWS);
    let len = tally.len();
    let counts = &tally.counts;

    for row in (0..ROWS).rev() {
        let floor = row_height * f64::from(row);

        match plot {
            Plot::Breadth => write!(out, ">{floor:7.2}% ")?,
            Plot::Depth => write!(out, ">{floor:8.1} ")?,
        }

        out.write_all(blocks.vertical_line().as_bytes())?;

        for &value in &data {
            // The bar's height above this row, in bar levels.
            let level = c_int((bar_levels * (value - floor) / row_height).round() - 1.0);

            match usize::try_from(level) {
                Ok(level) => out.write_all(bars[level.min(bars.len() - 1)].as_bytes())?,
                Err(_) => out.write_all(b" ")?,
            }
        }

        out.write_all(blocks.vertical_line().as_bytes())?;
        out.write_all(b" ")?;

        match row {
            9 => write!(out, "Number of reads: {}", counts.selected)?,
            8 if counts.reads > counts.selected => {
                write!(out, "    ({} filtered)", counts.reads - counts.selected)?;
            }
            7 => write!(
                out,
                "Covered bases:   {}bp",
                readable_bps(tally.covered_bases as f64)
            )?,
            6 => write!(
                out,
                "Percent covered: {}%",
                format_g(100.0 * tally.covered_bases as f64 / len, 4)
            )?,
            5 => write!(
                out,
                "Mean coverage:   {}x",
                format_g(tally.depth_sum as f64 / len, 3)
            )?,
            4 => write!(
                out,
                "Mean baseQ:      {}",
                format_g(tally.mean_base_quality(), 3)
            )?,
            3 => {
                // samtools divides by the read count unguarded: 0/0 prints as -nan.
                let mean = if counts.selected > 0 {
                    format_g(
                        counts.mapping_quality_sum as f64 / counts.selected as f64,
                        3,
                    )
                } else {
                    "-nan".to_owned()
                };
                write!(out, "Mean mapQ:       {mean}")?;
            }
            1 => write!(out, "Histo bin width: {}bp", readable_bps(bin_width))?,
            0 => match plot {
                Plot::Breadth => write!(out, "Histo max bin:   {}%", format_g(max, 5))?,
                Plot::Depth => write!(out, "Histo max cov:   {}", format_g(max, 5))?,
            },
            _ => {}
        }

        writeln!(out)?;
    }

    // The x axis: a label every 10 bins.
    let n_bins = i64::try_from(tally.bins.len()).unwrap_or(i64::MAX);

    write!(
        out,
        "     {}",
        center(&readable_bps((tally.start + 1) as f64))
    )?;

    for bin in (10..10 * (n_bins / 10)).step_by(10) {
        let position = tally.start + tally.bin_width * bin;
        write!(out, "{}", center(&readable_bps(position as f64)))?;
    }

    let padding = usize::try_from(n_bins % 10).unwrap_or(0);
    writeln!(
        out,
        "{:>padding$}{}",
        " ",
        center(&readable_bps(tally.end as f64))
    )?;

    Ok(())
}

/// Converts like C's `(int)` cast on x86-64: NaN and out-of-range values become `INT_MIN`.
fn c_int(value: f64) -> i64 {
    if value.is_nan() || value < f64::from(i32::MIN) || value >= f64::from(i32::MAX) + 1.0 {
        i64::from(i32::MIN)
    } else {
        value as i64
    }
}

/// Formats a base count with a unit (`1.2K`, `3.45M`), with as many decimals as the unit's
/// rank (samtools' `readable_bps`).
fn readable_bps(mut n: f64) -> String {
    const UNITS: [&str; 5] = ["", "K", "M", "G", "T"];

    let mut i = 0;

    while n >= 1000.0 && i < UNITS.len() - 1 {
        n /= 1000.0;
        i += 1;
    }

    format!("{n:.i$}{}", UNITS[i])
}

/// Centers an axis label in 10 columns (samtools' `center_text`).
fn center(text: &str) -> String {
    const WIDTH: usize = 10;

    let len = text.len();
    let padding = WIDTH.saturating_sub(len) / 2;
    let odd = WIDTH.saturating_sub(len) % 2;

    if padding >= 1 {
        format!(
            " {text:>w1$}{:>w2$}",
            " ",
            w1 = len + padding,
            w2 = padding - 1 + odd
        )
    } else {
        text.to_owned()
    }
}

/// The width of the terminal: `COLUMNS`, or the terminal on stderr (0 if neither).
pub(super) fn terminal_columns() -> i64 {
    match std::env::var("COLUMNS") {
        Ok(columns) => atoi(&columns),
        Err(_) => stderr_columns(),
    }
}

#[cfg(unix)]
fn stderr_columns() -> i64 {
    let mut size = libc::winsize {
        ws_row: 0,
        ws_col: 0,
        ws_xpixel: 0,
        ws_ypixel: 0,
    };

    // SAFETY: `TIOCGWINSZ` writes a `winsize` through the pointer, which is valid for writes.
    let ret = unsafe { libc::ioctl(libc::STDERR_FILENO, libc::TIOCGWINSZ, &raw mut size) };

    if ret == 0 { i64::from(size.ws_col) } else { 0 }
}

#[cfg(not(unix))]
fn stderr_columns() -> i64 {
    0
}

/// C's `atoi`: an optionally signed decimal prefix, or 0.
fn atoi(s: &str) -> i64 {
    let s = s.trim_start();
    let (sign, digits) = match s.as_bytes().first() {
        Some(b'-') => (-1, &s[1..]),
        Some(b'+') => (1, &s[1..]),
        _ => (1, s),
    };

    let end = digits
        .find(|c: char| !c.is_ascii_digit())
        .unwrap_or(digits.len());

    sign * digits[..end].parse::<i64>().unwrap_or(0)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_readable_bps() {
        assert_eq!(readable_bps(1.0), "1");
        assert_eq!(readable_bps(999.0), "999");
        assert_eq!(readable_bps(20_000.0), "20.0K");
        assert_eq!(readable_bps(1_250.0), "1.2K");
        assert_eq!(readable_bps(8_000_000.0), "8.00M");
    }

    #[test]
    fn test_center() {
        assert_eq!(center("1"), "     1    ");
        assert_eq!(center("250"), "    250   ");
        assert_eq!(center("1.0K"), "    1.0K  ");
        assert_eq!(center("123.456G"), "  123.456G ");
        assert_eq!(center("1000.000T"), "1000.000T");
    }

    #[test]
    fn test_atoi() {
        assert_eq!(atoi("120"), 120);
        assert_eq!(atoi(" -5x"), -5);
        assert_eq!(atoi("abc"), 0);
        assert_eq!(atoi(""), 0);
    }
}
