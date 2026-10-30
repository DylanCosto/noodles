use std::io::{self, Write};

use crate::io::writer::num;

pub(super) fn write_float<W>(writer: &mut W, n: f32) -> io::Result<()>
where
    W: Write,
{
    // § 1.5 "The alignment section: optional fields" (2023-05-24) only allows finite values,
    // but HTSlib writes `nan` and `inf`, so they are written rather than rejected.
    num::write_f32(writer, n)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_write_float() -> io::Result<()> {
        let mut buf = Vec::new();

        buf.clear();
        write_float(&mut buf, 0.0)?;
        assert_eq!(buf, b"0");

        buf.clear();
        write_float(&mut buf, f32::NAN)?;
        assert_eq!(buf, b"nan");

        Ok(())
    }
}
