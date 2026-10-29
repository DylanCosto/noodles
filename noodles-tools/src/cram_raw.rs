//! Parsing of raw CRAM 3 structures that noodles does not expose: ITF-8 and LTF-8 integers and
//! slice headers.

use std::io;

/// The parts of a slice header that the block-level commands use.
pub struct SliceHeader {
    /// The reference sequence ID: -1 for unmapped and -2 for multiple reference sequences.
    pub reference_sequence_id: i32,
    /// The 1-based alignment start.
    pub alignment_start: i32,
    /// The number of blocks following the slice header.
    pub block_count: usize,
    /// The content ID of the embedded reference block, if any.
    pub embedded_reference: Option<i32>,
}

/// Reads a CRAM 3 slice header.
pub fn read_slice_header(src: &mut &[u8]) -> io::Result<SliceHeader> {
    let reference_sequence_id = read_itf8(src)?;
    let alignment_start = read_itf8(src)?;
    read_itf8(src)?; // alignment span
    read_itf8(src)?; // record count
    read_ltf8(src)?; // record counter

    let block_count = read_len(src)?;

    let id_count = read_len(src)?;
    for _ in 0..id_count {
        read_itf8(src)?;
    }

    let id = read_itf8(src)?;
    let embedded_reference = (id >= 0).then_some(id);

    Ok(SliceHeader {
        reference_sequence_id,
        alignment_start,
        block_count,
        embedded_reference,
    })
}

pub fn invalid_data<E: Into<Box<dyn std::error::Error + Send + Sync>>>(e: E) -> io::Error {
    io::Error::new(io::ErrorKind::InvalidData, e)
}

/// Splits off the first `n` bytes.
pub fn take<'a>(src: &mut &'a [u8], n: usize) -> io::Result<&'a [u8]> {
    src.split_off(..n)
        .ok_or_else(|| io::Error::from(io::ErrorKind::UnexpectedEof))
}

fn read_u8(src: &mut &[u8]) -> io::Result<u8> {
    src.split_off_first()
        .copied()
        .ok_or_else(|| io::Error::from(io::ErrorKind::UnexpectedEof))
}

/// Reads an ITF-8 integer.
pub fn read_itf8(src: &mut &[u8]) -> io::Result<i32> {
    let b0 = read_u8(src)?;
    let extra = b0.leading_ones().min(4);

    let mut n = if extra == 4 {
        u32::from(b0 & 0x0f)
    } else {
        u32::from(b0 & (0xff >> (extra + 1)))
    };

    for i in 0..extra {
        let b = u32::from(read_u8(src)?);
        n = if i == 3 {
            (n << 4) | (b & 0x0f)
        } else {
            (n << 8) | b
        };
    }

    Ok(n as i32)
}

/// Reads a non-negative ITF-8 length or count.
pub fn read_len(src: &mut &[u8]) -> io::Result<usize> {
    usize::try_from(read_itf8(src)?).map_err(invalid_data)
}

/// Reads an LTF-8 integer.
pub fn read_ltf8(src: &mut &[u8]) -> io::Result<i64> {
    let b0 = read_u8(src)?;
    let extra = b0.leading_ones();
    let mut n = if extra >= 7 {
        0
    } else {
        u64::from(b0 & (0xff >> (extra + 1)))
    };

    for _ in 0..extra {
        n = (n << 8) | u64::from(read_u8(src)?);
    }

    Ok(n as i64)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_read_itf8() -> io::Result<()> {
        assert_eq!(read_itf8(&mut &[0x7f][..])?, 0x7f);
        assert_eq!(read_itf8(&mut &[0x80, 0xff][..])?, 0xff);
        assert_eq!(read_itf8(&mut &[0xc1, 0x00, 0x00][..])?, 0x1_0000);
        assert_eq!(read_itf8(&mut &[0xe1, 0x00, 0x00, 0x00][..])?, 0x100_0000);
        assert_eq!(read_itf8(&mut &[0xff, 0xff, 0xff, 0xff, 0x0f][..])?, -1);
        Ok(())
    }

    #[test]
    fn test_read_ltf8() -> io::Result<()> {
        assert_eq!(read_ltf8(&mut &[0x7f][..])?, 0x7f);
        assert_eq!(read_ltf8(&mut &[0x80, 0xff][..])?, 0xff);
        assert_eq!(read_ltf8(&mut &[0xff, 0, 0, 0, 0, 0, 0, 1, 0][..])?, 0x100);
        Ok(())
    }
}
