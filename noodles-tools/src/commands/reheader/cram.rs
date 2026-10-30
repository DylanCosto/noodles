//! Replacing a CRAM file's header container, copying every other container as is.

use std::{
    fs::{File, OpenOptions},
    io::{self, Read, Seek, SeekFrom, Write},
    path::Path,
};

use anyhow::{Context, bail};
use flate2::{Compression, write::GzEncoder};
use noodles_sam as sam;

use crate::cram_raw::{read_itf8, read_ltf8, take};

/// The length of the file definition: magic number, version, and file ID.
const FILE_DEFINITION_LEN: u64 = 26;

/// The header container's start is read in one go; container headers are much smaller.
const CONTAINER_HEADER_MAX_LEN: usize = 1 << 12;

/// The block content type of the SAM header.
const FILE_HEADER_CONTENT_TYPE: u8 = 0;

const RAW: u8 = 0;
const GZIP: u8 = 1;

/// The parsed start of a CRAM file: its version and its header container.
struct FileStart {
    major_version: u8,
    /// The file definition.
    definition: [u8; FILE_DEFINITION_LEN as usize],
    container: ContainerHeader,
}

/// The header of the header container.
struct ContainerHeader {
    /// The bytes from the reference sequence ID through the base count, kept as they are.
    fields: Vec<u8>,
    landmark_count: usize,
    /// The encoded length of the container header.
    len: usize,
    /// The length of the blocks that follow it.
    data_len: usize,
}

impl FileStart {
    fn read(file: &mut File) -> anyhow::Result<Self> {
        let mut definition = [0; FILE_DEFINITION_LEN as usize];
        file.read_exact(&mut definition)
            .context("reading the CRAM file definition")?;

        if &definition[..4] != b"CRAM" {
            bail!("not a CRAM file");
        }

        let major_version = definition[4];

        if !(2..=3).contains(&major_version) {
            bail!("unsupported CRAM version {major_version}");
        }

        let mut buf = Vec::with_capacity(CONTAINER_HEADER_MAX_LEN);
        file.take(CONTAINER_HEADER_MAX_LEN as u64)
            .read_to_end(&mut buf)?;

        let container = ContainerHeader::parse(&buf, major_version)
            .context("reading the CRAM header container")?;

        Ok(Self {
            major_version,
            definition,
            container,
        })
    }

    fn has_checksums(&self) -> bool {
        self.major_version >= 3
    }

    /// The offset of the first container after the header container.
    fn data_start(&self) -> u64 {
        FILE_DEFINITION_LEN + (self.container.len + self.container.data_len) as u64
    }
}

impl ContainerHeader {
    fn parse(buf: &[u8], major_version: u8) -> io::Result<Self> {
        let mut src = buf;

        let data_len = take(&mut src, 4)?;
        let data_len = i32::from_le_bytes([data_len[0], data_len[1], data_len[2], data_len[3]]);
        let data_len = usize::try_from(data_len).map_err(crate::cram_raw::invalid_data)?;

        let fields_start = buf.len() - src.len();

        // Reference sequence ID, start, span, record count.
        for _ in 0..4 {
            read_itf8(&mut src)?;
        }

        // Record counter and base count.
        if major_version >= 3 {
            read_ltf8(&mut src)?;
        } else {
            read_itf8(&mut src)?;
        }

        read_ltf8(&mut src)?;

        let fields = buf[fields_start..buf.len() - src.len()].to_vec();

        // The block count.
        read_itf8(&mut src)?;

        let landmark_count = to_count(read_itf8(&mut src)?)?;

        for _ in 0..landmark_count {
            read_itf8(&mut src)?;
        }

        if major_version >= 3 {
            take(&mut src, 4)?;
        }

        Ok(Self {
            fields,
            landmark_count,
            len: buf.len() - src.len(),
            data_len,
        })
    }
}

fn to_count(n: i32) -> io::Result<usize> {
    usize::try_from(n).map_err(crate::cram_raw::invalid_data)
}

/// Writes the input with `header` in place of its header (`cram_reheader`). The version and the
/// containers are copied as they are.
pub fn reheader<W: Write>(input: &Path, header: &sam::Header, dst: &mut W) -> anyhow::Result<()> {
    let mut file = File::open(input).with_context(|| format!("opening {}", input.display()))?;
    let start = FileStart::read(&mut file)?;

    let block = encode_block(&header_data(header)?, GZIP, start.has_checksums())?;
    let container = encode_container_header(&[0; 6], 1, &[], block.len(), start.has_checksums())?;

    dst.write_all(&start.definition)?;
    dst.write_all(&container)?;
    dst.write_all(&block)?;

    file.seek(SeekFrom::Start(start.data_start()))?;
    io::copy(&mut file, dst)?;

    Ok(())
}

/// Replaces the header in the file itself, if it fits in the space of the old header container
/// (`cram_reheader_inplace`).
pub fn reheader_in_place(input: &Path, header: &sam::Header) -> anyhow::Result<()> {
    let mut file = OpenOptions::new()
        .read(true)
        .write(true)
        .open(input)
        .with_context(|| format!("opening {}", input.display()))?;

    let start = FileStart::read(&mut file)?;
    let data = header_data(header)?;

    if start.major_version == 2 {
        reheader_in_place_2(&mut file, &start, &data)
    } else {
        reheader_in_place_3(&mut file, &start, &data)
    }
}

/// CRAM 2.1 keeps the header in one uncompressed block, padded with NULs for growth; the new
/// header must fit in that block.
fn reheader_in_place_2(file: &mut File, start: &FileStart, data: &[u8]) -> anyhow::Result<()> {
    let block_start = FILE_DEFINITION_LEN + start.container.len as u64;

    let mut buf = vec![0; start.container.data_len.min(CONTAINER_HEADER_MAX_LEN)];
    file.seek(SeekFrom::Start(block_start))?;
    file.read_exact(&mut buf)?;

    let mut src = buf.as_slice();
    let method = take(&mut src, 2)?[0];
    read_itf8(&mut src)?;
    let compressed_len = to_count(read_itf8(&mut src)?)?;
    let len = to_count(read_itf8(&mut src)?)?;
    let data_start = block_start + (buf.len() - src.len()) as u64;

    if method != RAW || compressed_len != len {
        bail!("the CRAM header block is compressed; use the non-inplace version");
    }

    if len < data.len() {
        bail!(
            "New header will not fit. Use non-inplace version ({} > {len})",
            data.len()
        );
    }

    let mut block_data = data.to_vec();
    block_data.resize(len, 0);

    file.seek(SeekFrom::Start(data_start))?;
    file.write_all(&block_data)?;

    Ok(())
}

/// CRAM 3 header containers hold a (possibly compressed) header block and spare space. The new
/// header is written as one gzip block, the container keeps its size, and the rest is zeroed.
fn reheader_in_place_3(file: &mut File, start: &FileStart, data: &[u8]) -> anyhow::Result<()> {
    let old = &start.container;
    let total_len = old.len + old.data_len;

    // One block, and one landmark (at 0) if there were any.
    let landmarks: &[i32] = if old.landmark_count > 0 { &[0] } else { &[] };
    let header_len = encode_container_header(&old.fields, 1, landmarks, 0, true)?.len();

    let Some(data_len) = total_len.checked_sub(header_len) else {
        bail!("Quirk of fate makes this troublesome! Please use non-inplace version.");
    };

    let container = encode_container_header(&old.fields, 1, landmarks, data_len, true)?;
    let block = encode_block(data, GZIP, true)?;

    if block.len() > data_len {
        bail!(
            "New header will not fit. Use non-inplace version ({} > {data_len})",
            block.len()
        );
    }

    let mut buf = Vec::with_capacity(total_len);
    buf.extend_from_slice(&container);
    buf.extend_from_slice(&block);
    buf.resize(total_len, 0);

    file.seek(SeekFrom::Start(FILE_DEFINITION_LEN))?;
    file.write_all(&buf)?;

    Ok(())
}

/// The header block's data: the length of the SAM header text, then the text.
fn header_data(header: &sam::Header) -> io::Result<Vec<u8>> {
    let mut text = Vec::new();
    sam::io::Writer::new(&mut text).write_header(header)?;

    let len = i32::try_from(text.len()).map_err(crate::cram_raw::invalid_data)?;

    let mut data = Vec::with_capacity(4 + text.len());
    data.extend_from_slice(&len.to_le_bytes());
    data.extend_from_slice(&text);

    Ok(data)
}

/// Encodes a file header block, with a CRC32 for CRAM 3.
fn encode_block(data: &[u8], method: u8, with_checksum: bool) -> io::Result<Vec<u8>> {
    let compressed = match method {
        GZIP => {
            let mut encoder = GzEncoder::new(Vec::new(), Compression::best());
            encoder.write_all(data)?;
            encoder.finish()?
        }
        _ => data.to_vec(),
    };

    let mut block = vec![method, FILE_HEADER_CONTENT_TYPE];
    write_itf8(&mut block, 0);
    write_itf8(&mut block, to_i32(compressed.len())?);
    write_itf8(&mut block, to_i32(data.len())?);
    block.extend_from_slice(&compressed);

    if with_checksum {
        let crc = crc32fast::hash(&block);
        block.extend_from_slice(&crc.to_le_bytes());
    }

    Ok(block)
}

/// Encodes a container header from its raw fields (reference sequence ID through base count).
fn encode_container_header(
    fields: &[u8],
    block_count: i32,
    landmarks: &[i32],
    data_len: usize,
    with_checksum: bool,
) -> io::Result<Vec<u8>> {
    let mut buf = to_i32(data_len)?.to_le_bytes().to_vec();
    buf.extend_from_slice(fields);
    write_itf8(&mut buf, block_count);
    write_itf8(&mut buf, to_i32(landmarks.len())?);

    for &landmark in landmarks {
        write_itf8(&mut buf, landmark);
    }

    if with_checksum {
        let crc = crc32fast::hash(&buf);
        buf.extend_from_slice(&crc.to_le_bytes());
    }

    Ok(buf)
}

fn to_i32(n: usize) -> io::Result<i32> {
    i32::try_from(n).map_err(crate::cram_raw::invalid_data)
}

/// Writes an ITF-8 integer.
fn write_itf8(dst: &mut Vec<u8>, n: i32) {
    let n = n.cast_unsigned();
    let [b0, b1, b2, b3] = n.to_be_bytes();

    if n < 0x80 {
        dst.push(b3);
    } else if n < 0x4000 {
        dst.extend_from_slice(&[0x80 | b2, b3]);
    } else if n < 0x20_0000 {
        dst.extend_from_slice(&[0xc0 | b1, b2, b3]);
    } else if n < 0x1000_0000 {
        dst.extend_from_slice(&[0xe0 | b0, b1, b2, b3]);
    } else {
        dst.extend_from_slice(&[
            0xf0 | (b0 >> 4),
            (b0 << 4) | (b1 >> 4),
            (b1 << 4) | (b2 >> 4),
            (b2 << 4) | (b3 >> 4),
            b3 & 0x0f,
        ]);
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_write_itf8() -> io::Result<()> {
        for n in [
            0,
            1,
            0x7f,
            0x80,
            0x3fff,
            0x4000,
            0x1f_ffff,
            0x20_0000,
            0x0fff_ffff,
            0x1000_0000,
            -1,
            i32::MAX,
        ] {
            let mut buf = Vec::new();
            write_itf8(&mut buf, n);
            let mut src = buf.as_slice();
            assert_eq!(read_itf8(&mut src)?, n);
            assert!(src.is_empty());
        }

        Ok(())
    }
}
