//! Read names: tile and coordinate parsing for optical duplicates (`--read-coords`), barcodes (`--barcode-tag/-name/-rgx`), and UMIs moved to a tag (`--move-umi-to-tag`).

use std::ops::Range;

use anyhow::{Context, bail};
use regex::bytes::{Regex, RegexBuilder};

use super::{MAX_WARNINGS, UMI, record::name, warn};
use crate::bam_raw;

/// Flowcell coordinates parsed from a read name.
pub(super) struct Coordinates {
    /// The part of the name (such as the lane and tile) that must match.
    pub(super) prefix: Range<usize>,
    pub(super) x: i64,
    pub(super) y: i64,
}

/// How flowcell coordinates are read from read names.
pub(super) enum CoordinateParser {
    /// Illumina's colon-separated fields (`get_coordinates_colons`).
    Colons,
    /// `--read-coords`, with the groups of the x and y coordinates and the part that must match.
    Regex {
        regex: Regex,
        x: usize,
        y: usize,
        t: Option<usize>,
    },
}

/// Why coordinates could not be read from a read name.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub(super) enum CoordinateError {
    /// The `--read-coords` regex did not match (samtools fails if this is a kept read's name).
    Unmatched,
    /// The coordinates were not numbers, or the name had the wrong fields; a warning is counted.
    Unreadable,
}

/// Where barcodes that duplicates must match come from.
pub(super) enum Barcodes {
    None,
    Tag([u8; 2]),
    /// The first group of a regex on the read name.
    Name(Regex),
}

/// `--move-umi-to-tag`: the barcode regex and the separator removed with the UMI.
pub(super) struct UmiMover {
    pub(super) regex: Regex,
    pub(super) separator: u8,
}

/// `get_coordinates_colons`: the x and y coordinates from an Illumina-style read name
/// (`machine:run:flowcell:lane:tile:x:y`, with an optional `:UMI`, or older forms with 3 or 4
/// colons).
pub(super) fn colon_coordinates(name: &[u8], warnings: &mut u64) -> Option<Coordinates> {
    let mut separators = 0;
    let (mut x_start, mut y_start) = (0, 0);

    for (i, &b) in name.iter().enumerate() {
        if b == b':' {
            separators += 1;

            match separators {
                2 => x_start = i + 1,
                3 => y_start = i + 1,
                // HiSeq-style names
                4 => {
                    x_start = y_start;
                    y_start = i + 1;
                }
                // Newer Illumina names
                5 => x_start = i + 1,
                6 => y_start = i + 1,
                _ => {}
            }
        }
    }

    let display = String::from_utf8_lossy(name);

    if !matches!(separators, 3 | 4 | 6 | 7) {
        warn(
            warnings,
            format_args!("cannot decipher read name {display} for optical duplicate marking."),
        );
        return None;
    }

    let Some(x) = strtol(&name[x_start..]) else {
        warn(
            warnings,
            format_args!("cannot decipher x coordinate in {display} ."),
        );
        return None;
    };

    let Some(y) = strtol(&name[y_start..]) else {
        warn(
            warnings,
            format_args!("cannot decipher y coordinate in {display} ."),
        );
        return None;
    };

    Some(Coordinates {
        prefix: 0..x_start,
        x,
        y,
    })
}

/// `--read-coords` and `--coords-order`.
pub(super) fn coordinate_regex(pattern: &str, order: &str) -> anyhow::Result<CoordinateParser> {
    let (t, x, y) = if order.starts_with("txy") || order.starts_with("tyx") {
        (Some(1), 2, 3)
    } else if order.starts_with("xyt") || order.starts_with("yxt") {
        (Some(3), 1, 2)
    } else if order.starts_with("xty") || order.starts_with("ytx") {
        (Some(2), 1, 3)
    } else if order.starts_with("xy") || order.starts_with("yx") {
        (None, 1, 2)
    } else {
        bail!("[markdup] error,  could not recognise regex coordinate order \"{order}\".");
    };

    let regex = compile_regex(pattern).context("[markdup] error, regex fail")?;
    Ok(CoordinateParser::Regex { regex, x, y, t })
}

/// A regex over the bytes of read names.
pub(super) fn compile_regex(pattern: &str) -> Result<Regex, regex::Error> {
    RegexBuilder::new(pattern).unicode(false).build()
}

impl CoordinateParser {
    /// `get_coordinates`: the coordinates in a read name, counting a warning if they cannot be
    /// read. A name that a `--read-coords` regex does not match is not warned about.
    pub(super) fn parse(
        &self,
        name: &[u8],
        warnings: &mut u64,
    ) -> Result<Coordinates, CoordinateError> {
        let Self::Regex { regex, x, y, t } = self else {
            return colon_coordinates(name, warnings).ok_or(CoordinateError::Unreadable);
        };

        let groups = || {
            let captures = regex.captures(name)?;
            let prefix = match t {
                Some(t) => captures.get(*t)?.range(),
                None => 0..0,
            };
            Some((captures.get(*x)?, captures.get(*y)?, prefix))
        };

        let (x, y, prefix) = groups().ok_or(CoordinateError::Unmatched)?;

        let display = String::from_utf8_lossy(name);
        let mut coordinate = |axis: &str, m: regex::bytes::Match<'_>| {
            // samtools copies each coordinate into a 255-byte buffer.
            if m.len() > 254 {
                warn(
                    warnings,
                    format_args!(
                        "{axis} coordinate string longer than allowed qname length in {display} ({} long).",
                        m.len()
                    ),
                );
                return None;
            }

            let n = strtol(m.as_bytes());

            if n.is_none() {
                warn(
                    warnings,
                    format_args!(
                        "cannot decipher {axis} coordinate in {display} ({}).",
                        String::from_utf8_lossy(m.as_bytes())
                    ),
                );
            }

            n
        };

        let x = coordinate("x", x).ok_or(CoordinateError::Unreadable)?;
        let y = coordinate("y", y).ok_or(CoordinateError::Unreadable)?;

        Ok(Coordinates { prefix, x, y })
    }
}

impl Barcodes {
    /// The barcode's hash for a read's keys (0 if it has none), counting a warning if it cannot
    /// be read.
    pub(super) fn hash(&self, record: &[u8], warnings: &mut u64) -> i32 {
        let incoming_warnings = *warnings;

        let barcode = match self {
            Self::None => None,
            Self::Tag(tag) => match bam_raw::find_field(record, *tag) {
                Ok(Some(range)) if matches!(record[range.start + 2], b'Z' | b'H') => {
                    Some(&record[range.start + 3..range.end - 1])
                }
                Ok(None) => None,
                _ => {
                    warn(
                        warnings,
                        format_args!(
                            "{} tag wrong type. Aux tag needs to be a string type.",
                            String::from_utf8_lossy(tag).trim_end_matches('\0')
                        ),
                    );
                    None
                }
            },
            Self::Name(regex) => {
                let qname = name(record);

                match regex.captures(qname) {
                    Some(captures) => {
                        if captures.get(1).is_none() {
                            warn(
                                warnings,
                                format_args!(
                                    "barcode regex unable to match substring on {}.",
                                    String::from_utf8_lossy(qname)
                                ),
                            );
                        }

                        captures.get(1).map(|m| m.as_bytes())
                    }
                    None => {
                        warn(
                            warnings,
                            format_args!(
                                "barcode regex match error \"No match\" on {}.",
                                String::from_utf8_lossy(qname)
                            ),
                        );
                        None
                    }
                }
            }
        };

        if *warnings == MAX_WARNINGS && incoming_warnings != *warnings {
            eprintln!(
                "[markdup] warning, {} barcode read warnings.  New warnings will not be reported.",
                *warnings
            );
        }

        barcode.map_or(0, |b| do_hash(b) as i32)
    }
}

/// samtools markdup's `do_hash` (Jenkins' one-at-a-time hash).
pub(super) fn do_hash(key: &[u8]) -> u32 {
    let mut hash = 0u32;

    for &b in key {
        hash = hash.wrapping_add(u32::from(b));
        hash = hash.wrapping_add(hash << 10);
        hash ^= hash >> 6;
    }

    hash = hash.wrapping_add(hash << 3);
    hash ^= hash >> 11;
    hash.wrapping_add(hash << 15)
}

impl UmiMover {
    /// `move_umi_to_tag`: moves the UMI (the regex's first group) from the read name to an RX
    /// tag, removing one separator next to it.
    pub(super) fn move_to_tag(&self, record: &mut Vec<u8>) -> anyhow::Result<()> {
        let qname = name(record);

        let Some(umi) = self.regex.captures(qname).and_then(|c| c.get(1)) else {
            return Ok(());
        };

        let Range { mut start, mut end } = umi.range();
        let umi = umi.as_bytes().to_vec();

        if start > 0 && qname[start - 1] == self.separator {
            start -= 1;
        } else if end < qname.len() && qname[end] == self.separator {
            end += 1;
        }

        let mut new_name = Vec::with_capacity(qname.len() - (end - start));
        new_name.extend_from_slice(&qname[..start]);
        new_name.extend_from_slice(&qname[end..]);

        bam_raw::update_string_field(record, UMI, &umi)
            .context("[markdup] error, unable to add RX tag.")?;

        // l_read_name, then the NUL-terminated name at 32.
        let old_len = usize::from(record[8]);
        record.splice(32..32 + old_len - 1, new_name.iter().copied());
        record[8] = u8::try_from(new_name.len() + 1).expect("a shorter name fits");

        Ok(())
    }
}

/// C's `strtol(s, &end, 10)`, or `None` if it reads no digits. Out-of-range values saturate.
pub(super) fn strtol(s: &[u8]) -> Option<i64> {
    let s = s.trim_ascii_start();

    let (negative, digits) = match s.first() {
        Some(b'-') => (true, &s[1..]),
        Some(b'+') => (false, &s[1..]),
        _ => (false, s),
    };

    let len = digits.iter().take_while(|b| b.is_ascii_digit()).count();

    if len == 0 {
        return None;
    }

    let n = digits[..len].iter().fold(0i64, |n, &d| {
        n.saturating_mul(10).saturating_add(i64::from(d - b'0'))
    });

    Some(if negative { n.saturating_neg() } else { n })
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_coordinates() {
        let mut warnings = 0;

        let c = colon_coordinates(b"M1:7:FC1:1:1101:1234:5678", &mut warnings).expect("parsed");
        assert_eq!((c.prefix, c.x, c.y), (0..16, 1234, 5678));

        let c =
            colon_coordinates(b"M1:7:FC1:1:1101:1234:5678:ACGT", &mut warnings).expect("parsed");
        assert_eq!((c.prefix, c.x, c.y), (0..16, 1234, 5678));

        // HiSeq style: machine:lane:tile:x:y
        let c = colon_coordinates(b"HS:1:1101:12:34", &mut warnings).expect("parsed");
        assert_eq!((c.prefix, c.x, c.y), (0..10, 12, 34));

        assert_eq!(warnings, 0);
        assert!(colon_coordinates(b"read1", &mut warnings).is_none());
        assert!(colon_coordinates(b"a:b:c:x:y", &mut warnings).is_none());
        assert_eq!(warnings, 2);
    }

    #[test]
    fn test_coordinate_regex() -> anyhow::Result<()> {
        let mut warnings = 0;

        let parser = coordinate_regex(r"^(\w+)_(\d+)_(\d+)$", "txy")?;
        let c = parser
            .parse(b"tile7_120_340", &mut warnings)
            .expect("parsed");
        assert_eq!((c.prefix, c.x, c.y), (0..5, 120, 340));

        let parser = coordinate_regex(r"_(\d+)_(\d+)$", "yx")?;
        let c = parser.parse(b"a_1_2", &mut warnings).expect("parsed");
        assert_eq!((c.prefix, c.x, c.y), (0..0, 1, 2));

        // A name the regex does not match is skipped silently.
        assert!(matches!(
            parser.parse(b"a", &mut warnings),
            Err(CoordinateError::Unmatched)
        ));
        assert_eq!(warnings, 0);

        assert!(coordinate_regex("(x)", "abc").is_err());
        Ok(())
    }

    #[test]
    fn test_do_hash() {
        assert_eq!(do_hash(b""), 0);
        assert_eq!(do_hash(b"a"), 0xca2e_9442);
    }

    #[test]
    fn test_strtol() {
        assert_eq!(strtol(b"123:4"), Some(123));
        assert_eq!(strtol(b" -5x"), Some(-5));
        assert_eq!(strtol(b"x1"), None);
        assert_eq!(strtol(b"99999999999999999999"), Some(i64::MAX));
    }
}
