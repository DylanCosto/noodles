//! Rewriting BAM headers and concatenating BAM files by copying compressed BGZF blocks.
//!
//! Only the block holding the end of the header is decompressed and recompressed; every later
//! block is copied as is, as in `samtools reheader` and `samtools cat`.

use std::{
    fs::File,
    io::{self, Read, Seek, SeekFrom, Write},
    path::Path,
};

use anyhow::{Context, bail};
use flate2::read::DeflateDecoder;
use noodles_bam as bam;
use noodles_bgzf as bgzf;
use noodles_sam as sam;

/// The BGZF end-of-file marker block.
const BGZF_EOF: [u8; 28] = [
    0x1f, 0x8b, 0x08, 0x04, 0x00, 0x00, 0x00, 0x00, 0x00, 0xff, 0x06, 0x00, 0x42, 0x43, 0x02, 0x00,
    0x1b, 0x00, 0x03, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00,
];

/// Length of a BGZF block header, up to and including `BSIZE`.
const BLOCK_HEADER_LEN: usize = 18;

/// Length of a BGZF block footer (`CRC32` and `ISIZE`).
const BLOCK_FOOTER_LEN: usize = 8;

/// A BAM file split at the end of its header.
pub struct Body {
    file: File,
    /// Uncompressed record bytes that follow the header in its last block.
    leftover: Vec<u8>,
    /// Compressed offset of the first block after the header's last block.
    next_block: u64,
    /// Compressed offset where the copy ends: before the EOF marker, if there is one.
    end: u64,
}

/// Reads a BAM file's header and locates the records after it.
pub fn split(path: &Path) -> anyhow::Result<(sam::Header, Body)> {
    let mut reader = File::open(path)
        .map(bam::io::Reader::new)
        .with_context(|| format!("opening {}", path.display()))?;

    let header = reader
        .read_header()
        .with_context(|| format!("reading header of {}", path.display()))?;

    let (block_start, offset) = reader.get_ref().virtual_position().into();
    let offset = usize::from(offset);

    let mut file = reader.into_inner().into_inner();
    let len = file.metadata()?.len();

    let (data, block_len) = read_block(&mut file, block_start)
        .with_context(|| format!("reading BGZF block at {block_start} of {}", path.display()))?;

    let leftover = data
        .get(offset..)
        .with_context(|| format!("{}: header ends past its block", path.display()))?
        .to_vec();

    let has_eof = ends_with_eof(&mut file, len)?;
    let eof_len = if has_eof { BGZF_EOF.len() as u64 } else { 0 };

    Ok((
        header,
        Body {
            file,
            leftover,
            next_block: block_start + block_len,
            end: len - eof_len,
        },
    ))
}

impl Body {
    /// Writes the records: the leftover bytes recompressed, then the remaining blocks copied.
    /// The input's EOF marker is not copied.
    pub fn copy_to<W: Write>(mut self, dst: &mut W) -> anyhow::Result<()> {
        if !self.leftover.is_empty() {
            let mut writer = bgzf::io::Writer::new(&mut *dst);
            writer.write_all(&self.leftover)?;
            // Flushing writes the block; taking the inner writer skips the EOF marker that
            // finishing (or dropping) would write.
            writer.flush()?;
            writer.into_inner();
        }

        if self.end > self.next_block {
            self.file.seek(SeekFrom::Start(self.next_block))?;
            let mut rest = (&mut self.file).take(self.end - self.next_block);
            io::copy(&mut rest, dst)?;
        }

        Ok(())
    }
}

/// Writes a BAM header in its own BGZF blocks, without an EOF marker.
pub fn write_header<W: Write>(dst: &mut W, header: &sam::Header) -> io::Result<()> {
    let mut bgzf_writer = bgzf::io::Writer::new(&mut *dst);
    bam::io::Writer::from(&mut bgzf_writer).write_header(header)?;
    bgzf_writer.flush()?;
    bgzf_writer.into_inner();
    Ok(())
}

pub fn write_eof<W: Write>(dst: &mut W) -> io::Result<()> {
    dst.write_all(&BGZF_EOF)
}

/// Returns whether the file ends with a BGZF EOF marker.
pub fn ends_with_eof(file: &mut File, len: u64) -> io::Result<bool> {
    let eof_len = BGZF_EOF.len() as u64;

    if len < eof_len {
        return Ok(false);
    }

    let mut buf = [0; BGZF_EOF.len()];
    file.seek(SeekFrom::Start(len - eof_len))?;
    file.read_exact(&mut buf)?;

    Ok(buf == BGZF_EOF)
}

/// Reads and inflates the BGZF block at `offset`, returning its data and compressed length.
fn read_block(file: &mut File, offset: u64) -> anyhow::Result<(Vec<u8>, u64)> {
    file.seek(SeekFrom::Start(offset))?;

    let mut header = [0; BLOCK_HEADER_LEN];

    match file.read_exact(&mut header) {
        Ok(()) => {}
        // The header ends at the end of the file.
        Err(e) if e.kind() == io::ErrorKind::UnexpectedEof => return Ok((Vec::new(), 0)),
        Err(e) => return Err(e.into()),
    }

    if header[..4] != [0x1f, 0x8b, 0x08, 0x04] || header[12..14] != *b"BC" {
        bail!("not a BGZF block");
    }

    let block_len = usize::from(u16::from_le_bytes([header[16], header[17]])) + 1;

    let payload_len = block_len
        .checked_sub(BLOCK_HEADER_LEN + BLOCK_FOOTER_LEN)
        .context("invalid BGZF block size")?;

    let mut payload = vec![0; payload_len + BLOCK_FOOTER_LEN];
    file.read_exact(&mut payload)?;

    let mut data = Vec::new();
    DeflateDecoder::new(&payload[..payload_len]).read_to_end(&mut data)?;

    Ok((data, block_len as u64))
}
