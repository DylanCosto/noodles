//! Region strings (`hts_parse_region`) and BED region lists (`bed_read`, `bed_reglist`) as
//! samtools consensus reads them.

use std::{collections::HashMap, fs, io::Read, path::Path};

use anyhow::{Context, bail};
use flate2::read::MultiGzDecoder;

use super::pileup::POS_MAX;
use crate::khash::StringSet;

/// `hts_parse_decimal` with `HTS_PARSE_THOUSANDS_SEP`: returns the value and the number of
/// bytes used (0 if there were no digits).
fn parse_decimal(s: &[u8]) -> (i64, usize) {
    let mut n: u64 = 0;
    let mut digits = 0;
    let mut decimals = 0;
    let mut e: i64 = 0;

    let mut i = s.iter().take_while(|c| c.is_ascii_whitespace()).count();
    let mut sign = b'+';

    if let Some(&c @ (b'+' | b'-')) = s.get(i) {
        sign = c;
        i += 1;
    }

    let push = |n: u64, c: u8| n.wrapping_mul(10).wrapping_add(u64::from(c - b'0'));

    while let Some(&c) = s.get(i) {
        if c.is_ascii_digit() {
            digits += 1;
            n = push(n, c);
        } else if c != b',' {
            break;
        }

        i += 1;
    }

    if s.get(i) == Some(&b'.') {
        i += 1;

        while let Some(&c) = s.get(i).filter(|c| c.is_ascii_digit()) {
            decimals += 1;
            digits += 1;
            n = push(n, c);
            i += 1;
        }
    }

    match s.get(i) {
        Some(b'e' | b'E') => {
            i += 1;
            let mut esign = b'+';

            if let Some(&c @ (b'+' | b'-')) = s.get(i) {
                esign = c;
                i += 1;
            }

            while let Some(&c) = s.get(i).filter(|c| c.is_ascii_digit()) {
                e = e.wrapping_mul(10).wrapping_add(i64::from(c - b'0'));
                i += 1;
            }

            if esign == b'-' {
                e = -e;
            }
        }
        Some(b'k' | b'K') => {
            e += 3;
            i += 1;
        }
        Some(b'm' | b'M') => {
            e += 6;
            i += 1;
        }
        Some(b'g' | b'G') => {
            e += 9;
            i += 1;
        }
        _ => {}
    }

    e -= decimals;

    while e > 0 && n != 0 {
        n = n.wrapping_mul(10);
        e -= 1;
    }

    while e < 0 && n != 0 {
        n /= 10;
        e += 1;
    }

    let n = n as i64;
    let value = if sign == b'+' { n } else { n.wrapping_neg() };
    (value, if digits > 0 { i } else { 0 })
}

/// `hts_parse_region` (samtools conventions: `chr:100` runs to the end): returns the reference
/// sequence ID and the 0-based, half-open interval, or `None` if the region is invalid.
pub fn parse_region(s: &str, names: &[Vec<u8>]) -> Option<(usize, i64, i64)> {
    let getid = |name: &[u8]| names.iter().position(|n| n.as_slice() == name);
    let b = s.as_bytes();

    let (body, colon, quoted) = if let Some(rest) = b.strip_prefix(b"{") {
        let close = rest.iter().position(|&c| c == b'}')?;
        let colon = (rest.get(close + 1) == Some(&b':')).then_some(close + 1);
        (rest, colon, 1)
    } else {
        (b, b.iter().rposition(|&c| c == b':'), 0)
    };

    let Some(colon) = colon else {
        let tid = getid(&body[..body.len() - quoted])?;
        return Some((tid, 0, POS_MAX));
    };

    if quoted == 0
        && let Some(tid) = getid(body)
    {
        // The whole string is a name; it is ambiguous if the part before the colon is too.
        return match getid(&body[..colon]) {
            Some(_) => None,
            None => Some((tid, 0, POS_MAX)),
        };
    }

    let tid = getid(&body[..colon - quoted])?;
    let coords = &body[colon + 1..];

    let (beg, used) = parse_decimal(coords);
    let beg = beg - 1;
    let hyphen = &coords[used..];

    if beg < 0 {
        if beg != -1 && hyphen.first() == Some(&b'-') && !coords.is_empty() {
            return None;
        }

        match hyphen.first() {
            None | Some(b',') => {
                return Some((tid, 0, if beg == -1 { POS_MAX } else { -(beg + 1) }));
            }
            Some(c) if c.is_ascii_digit() => {
                return Some((tid, 0, if beg == -1 { POS_MAX } else { -(beg + 1) }));
            }
            _ if beg < -1 => return None,
            _ => {}
        }
    }

    let mut end = match hyphen.first() {
        None => POS_MAX,
        Some(b'-') => {
            let (end, used) = parse_decimal(&hyphen[1..]);

            match hyphen[1 + used..].first() {
                None | Some(b',') => end,
                _ => return None,
            }
        }
        _ => return None,
    };

    if end == 0 {
        end = POS_MAX;
    }

    if beg >= end {
        return None;
    }

    // `hts_itr_query` (and the threaded code) raise a negative start to 0.
    Some((tid, beg.max(0), end))
}

/// A `bed_read` region list: intervals by name, with names in khash order and each name's
/// intervals sorted (not merged).
#[derive(Debug, Default)]
pub struct Bed {
    order: StringSet,
    intervals: HashMap<Vec<u8>, Vec<(i64, i64)>>,
}

impl Bed {
    /// `bed_read`: reads a BED file (plain or gzip-compressed); 2-column lines are 1-based
    /// positions.
    pub fn read(path: &Path) -> anyhow::Result<Self> {
        let raw = fs::read(path)?;

        let text = if raw.starts_with(&[0x1f, 0x8b]) {
            let mut text = Vec::new();
            MultiGzDecoder::new(&raw[..]).read_to_end(&mut text)?;
            text
        } else {
            raw
        };

        let mut bed = Self::default();

        for (n, line) in text
            .split(|&c| c == b'\n')
            .filter(|line| !line.is_empty())
            .enumerate()
        {
            let line = line.trim_ascii_start();

            if line.is_empty() || line[0] == b'#' {
                continue;
            }

            let name_len = line
                .iter()
                .position(u8::is_ascii_whitespace)
                .unwrap_or(line.len());
            let name = &line[..name_len];

            let (mut beg, mut end, num) = if name_len < line.len() {
                scan_two_u64(&line[name_len + 1..])
            } else {
                (0, 0, 0)
            };

            if num == 1 {
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

            bed.order.insert(name);
            bed.intervals
                .entry(name.to_vec())
                .or_default()
                .push((beg as i64, end as i64));
        }

        for intervals in bed.intervals.values_mut() {
            intervals.sort_unstable();
        }

        Ok(bed)
    }

    /// `bed_reglist` and `fill_reglist_tid`: the intervals of the names in the header, in
    /// khash order, with their reference sequence IDs. Unknown names are skipped with a warning.
    pub fn reglist(&self, names: &[Vec<u8>]) -> anyhow::Result<Vec<(usize, i64, i64)>> {
        if self.order.is_empty() {
            bail!("the regions file has no regions");
        }

        let mut list = Vec::new();
        let mut known = 0;

        for name in self.order.iter() {
            match names.iter().position(|n| n.as_slice() == name) {
                Some(tid) => {
                    known += 1;
                    let intervals = self.intervals.get(name).map_or(&[][..], Vec::as_slice);
                    list.extend(intervals.iter().map(|&(beg, end)| (tid, beg, end)));
                }
                None => eprintln!(
                    "[W::fill_reglist_tid] Region '{}' specifies an unknown reference name",
                    String::from_utf8_lossy(name)
                ),
            }
        }

        if known == 0 {
            bail!("no regions file region names are in the header");
        }

        Ok(list)
    }
}

/// `sscanf(s, "%" SCNu64 " %" SCNu64, ...)`: returns the values and how many were read.
fn scan_two_u64(s: &[u8]) -> (u64, u64, u32) {
    fn scan(s: &[u8]) -> Option<(u64, &[u8])> {
        let s = s.trim_ascii_start();
        let (neg, s) = match s.first() {
            Some(b'-') => (true, &s[1..]),
            Some(b'+') => (false, &s[1..]),
            _ => (false, s),
        };

        let digits = s.iter().take_while(|c| c.is_ascii_digit()).count();

        if digits == 0 {
            return None;
        }

        let n = std::str::from_utf8(&s[..digits])
            .ok()
            .and_then(|d| d.parse::<u64>().ok())
            .unwrap_or(u64::MAX);

        Some((if neg { n.wrapping_neg() } else { n }, &s[digits..]))
    }

    match scan(s) {
        None => (0, 0, 0),
        Some((beg, rest)) => match scan(rest) {
            None => (beg, 0, 1),
            Some((end, _)) => (beg, end, 2),
        },
    }
}

/// Reads a BED file, with the error message samtools gives.
pub fn read_bed(path: &Path) -> anyhow::Result<Bed> {
    Bed::read(path).with_context(|| format!("Could not read file \"{}\"", path.display()))
}

#[cfg(test)]
mod tests {
    use super::*;

    fn names() -> Vec<Vec<u8>> {
        vec![b"chr1".to_vec(), b"chr2".to_vec(), b"HLA:1".to_vec()]
    }

    #[test]
    fn test_parse_region() {
        let names = names();

        assert_eq!(parse_region("chr1", &names), Some((0, 0, POS_MAX)));
        assert_eq!(parse_region("chr2:100", &names), Some((1, 99, POS_MAX)));
        assert_eq!(parse_region("chr2:1,000-2k", &names), Some((1, 999, 2000)));
        assert_eq!(parse_region("chr1:-100", &names), Some((0, 0, 100)));
        assert_eq!(parse_region("chr1:100-", &names), Some((0, 99, POS_MAX)));
        assert_eq!(parse_region("HLA:1", &names), Some((2, 0, POS_MAX)));
        assert_eq!(parse_region("{HLA:1}:5-6", &names), Some((2, 4, 6)));
        assert_eq!(parse_region("chr1:200-100", &names), None);
        assert_eq!(parse_region("chr1:1-2x", &names), None);
        assert_eq!(parse_region("chr3", &names), None);
    }

    #[test]
    fn test_scan_two_u64() {
        assert_eq!(scan_two_u64(b"\t10\t20\tname"), (10, 20, 2));
        assert_eq!(scan_two_u64(b" 7"), (7, 0, 1));
        assert_eq!(scan_two_u64(b"x"), (0, 0, 0));
    }
}
