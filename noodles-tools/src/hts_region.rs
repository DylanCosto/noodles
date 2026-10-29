//! HTSlib's region strings (`hts_parse_region`, `hts_parse_decimal`).

/// `HTS_POS_MAX`: the end of a region without one.
pub const POS_MAX: i64 = (1 << 62) - 1;

/// `hts_parse_decimal` with `HTS_PARSE_THOUSANDS_SEP`: returns the value and the number of
/// bytes used (0 if there were no digits).
pub fn parse_decimal(s: &[u8]) -> (i64, usize) {
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

/// `hts_parse_region` (samtools conventions: `chr:100` runs to the end), with `getid` looking up
/// names: returns the reference sequence ID and the 0-based, half-open interval, or `None` if the region is invalid.
pub fn parse_region<F>(s: &str, getid: F) -> Option<(usize, i64, i64)>
where
    F: Fn(&[u8]) -> Option<usize>,
{
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
