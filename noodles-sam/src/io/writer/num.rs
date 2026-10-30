use std::io::{self, Write};

use lexical_core::FormattedSize;

pub fn write_i8<W>(writer: &mut W, n: i8) -> io::Result<()>
where
    W: Write,
{
    let mut dst = [0; i8::FORMATTED_SIZE_DECIMAL];
    let buf = lexical_core::write(n, &mut dst);
    writer.write_all(buf)
}

pub fn write_u8<W>(writer: &mut W, n: u8) -> io::Result<()>
where
    W: Write,
{
    let mut dst = [0; u8::FORMATTED_SIZE_DECIMAL];
    let buf = lexical_core::write(n, &mut dst);
    writer.write_all(buf)
}

pub fn write_i16<W>(writer: &mut W, n: i16) -> io::Result<()>
where
    W: Write,
{
    let mut dst = [0; i16::FORMATTED_SIZE_DECIMAL];
    let buf = lexical_core::write(n, &mut dst);
    writer.write_all(buf)
}

pub fn write_u16<W>(writer: &mut W, n: u16) -> io::Result<()>
where
    W: Write,
{
    let mut dst = [0; u16::FORMATTED_SIZE_DECIMAL];
    let buf = lexical_core::write(n, &mut dst);
    writer.write_all(buf)
}

pub fn write_i32<W>(writer: &mut W, n: i32) -> io::Result<()>
where
    W: Write,
{
    let mut dst = [0; i32::FORMATTED_SIZE_DECIMAL];
    let buf = lexical_core::write(n, &mut dst);
    writer.write_all(buf)
}

pub fn write_u32<W>(writer: &mut W, n: u32) -> io::Result<()>
where
    W: Write,
{
    let mut dst = [0; u32::FORMATTED_SIZE_DECIMAL];
    let buf = lexical_core::write(n, &mut dst);
    writer.write_all(buf)
}

pub fn write_usize<W>(writer: &mut W, n: usize) -> io::Result<()>
where
    W: Write,
{
    let mut dst = [0; usize::FORMATTED_SIZE_DECIMAL];
    let buf = lexical_core::write(n, &mut dst);
    writer.write_all(buf)
}

/// Writes a float as HTSlib does (`kputd`): like `printf("%g")`, with 6 significant digits,
/// trailing zeros removed, and an exponent (`1e+30`) only for magnitudes below 1e-4 or of 1e6
/// and above. Non-finite values are written as `nan`, `inf`, and `-inf`.
pub fn write_f32<W>(writer: &mut W, n: f32) -> io::Result<()>
where
    W: Write,
{
    let mut buf = [0; 32];
    let len = format_g(&mut buf, f64::from(n));
    writer.write_all(&buf[..len])
}

/// Formats `n` as `printf("%g")` does into `buf` and returns the length written.
fn format_g(buf: &mut [u8; 32], n: f64) -> usize {
    use std::io::Write as _;

    const PRECISION: i32 = 6;

    let mut dst = &mut buf[..];
    let capacity = dst.len();

    if n.is_nan() {
        let s: &[u8] = if n.is_sign_negative() {
            b"-nan"
        } else {
            b"nan"
        };
        dst.write_all(s).ok();
        return s.len();
    } else if n.is_infinite() {
        let s: &[u8] = if n.is_sign_negative() {
            b"-inf"
        } else {
            b"inf"
        };
        dst.write_all(s).ok();
        return s.len();
    } else if n == 0.0 {
        let s: &[u8] = if n.is_sign_negative() { b"-0" } else { b"0" };
        dst.write_all(s).ok();
        return s.len();
    }

    // The exponent of `n` after rounding to the precision.
    let mut scientific = [0; 32];
    let scientific_len = {
        let mut dst = &mut scientific[..];
        write!(dst, "{:.*e}", (PRECISION - 1) as usize, n).ok();
        32 - dst.len()
    };
    let scientific = &scientific[..scientific_len];
    let e = scientific
        .iter()
        .position(|&b| b == b'e')
        .unwrap_or(scientific_len);
    let exponent: i32 = std::str::from_utf8(&scientific[e + 1..])
        .ok()
        .and_then(|s| s.parse().ok())
        .unwrap_or_default();

    if !(-4..PRECISION).contains(&exponent) {
        let mantissa = trim_fraction(&scientific[..e]);
        dst.write_all(mantissa).ok();
        let sign = if exponent < 0 { '-' } else { '+' };
        write!(dst, "e{sign}{:02}", exponent.unsigned_abs()).ok();
    } else {
        let mut fixed = [0; 32];
        let fixed_len = {
            let mut dst = &mut fixed[..];
            write!(dst, "{:.*}", (PRECISION - 1 - exponent) as usize, n).ok();
            32 - dst.len()
        };
        dst.write_all(trim_fraction(&fixed[..fixed_len])).ok();
    }

    capacity - dst.len()
}

/// Removes trailing zeros from a decimal fraction, and the decimal point if nothing follows it.
fn trim_fraction(s: &[u8]) -> &[u8] {
    if !s.contains(&b'.') {
        return s;
    }

    let s = s.trim_ascii_end();
    let end = s.iter().rposition(|&b| b != b'0').map_or(0, |i| i + 1);
    let s = &s[..end];
    s.strip_suffix(b".").unwrap_or(s)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_write_f32() -> io::Result<()> {
        fn t(buf: &mut Vec<u8>, n: f32, expected: &[u8]) -> io::Result<()> {
            buf.clear();
            write_f32(buf, n)?;
            assert_eq!(buf, expected);
            Ok(())
        }

        let mut buf = Vec::new();

        t(&mut buf, -1.0, b"-1")?;
        t(&mut buf, -0.5, b"-0.5")?;
        t(&mut buf, 0.0, b"0")?;
        t(&mut buf, 0.5, b"0.5")?;
        t(&mut buf, 1.0, b"1")?;
        t(&mut buf, -0.0, b"-0")?;
        t(&mut buf, 2.345_67, b"2.34567")?;
        t(&mut buf, 1.234_567_9, b"1.23457")?;
        t(&mut buf, 0.1, b"0.1")?;
        t(&mut buf, 0.000_1, b"0.0001")?;
        t(&mut buf, 0.000_012_5, b"1.25e-05")?;
        t(&mut buf, 999_999.0, b"999999")?;
        t(&mut buf, 1_000_000.0, b"1e+06")?;
        t(&mut buf, 1e30, b"1e+30")?;
        t(&mut buf, -1e-30, b"-1e-30")?;
        t(&mut buf, 1e10, b"1e+10")?;
        t(&mut buf, 123_456.7, b"123457")?;
        t(&mut buf, f32::NAN, b"nan")?;
        t(&mut buf, f32::INFINITY, b"inf")?;
        t(&mut buf, f32::NEG_INFINITY, b"-inf")?;
        t(&mut buf, f32::MAX, b"3.40282e+38")?;
        t(&mut buf, f32::MIN_POSITIVE, b"1.17549e-38")?;

        Ok(())
    }
}
