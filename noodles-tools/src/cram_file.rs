//! CRAM file structure at the container level: the file definition, container headers, and
//! header containers, for copying containers without decoding them (`reheader`, `cat`).
//!
//! CRAM 2.1 and 3.x are supported.

use std::io::{self, Read, Write};

use anyhow::{Context, bail};
use flate2::{Compression, write::GzEncoder};
use noodles_sam as sam;

use crate::cram_raw::{invalid_data, read_itf8, read_ltf8};

/// The length of the file definition: magic number, version, and file ID.
pub const FILE_DEFINITION_LEN: usize = 26;

/// The block content type of the SAM header.
const FILE_HEADER_CONTENT_TYPE: u8 = 0;

/// Block compression methods.
pub const RAW: u8 = 0;
pub const GZIP: u8 = 1;

/// The reference sequence ID of multi-reference containers.
pub const MULTI_REFERENCE: i32 = -2;

/// The start position of EOF containers ("EOF" in ASCII).
const EOF_START: i32 = 0x0045_4f46;

/// The CRAM 3 EOF container.
pub const EOF_3: [u8; 38] = [
    0x0f, 0x00, 0x00, 0x00, 0xff, 0xff, 0xff, 0xff, 0x0f, 0xe0, 0x45, 0x4f, 0x46, 0x00, 0x00, 0x00,
    0x00, 0x01, 0x00, 0x05, 0xbd, 0xd9, 0x4f, 0x00, 0x01, 0x00, 0x06, 0x06, 0x01, 0x00, 0x01, 0x00,
    0x01, 0x00, 0xee, 0x63, 0x01, 0x4b,
];

/// The CRAM 2.1 EOF container.
pub const EOF_2: [u8; 30] = [
    0x0b, 0x00, 0x00, 0x00, 0xff, 0xff, 0xff, 0xff, 0xff, 0xe0, 0x45, 0x4f, 0x46, 0x00, 0x00, 0x00,
    0x00, 0x01, 0x00, 0x00, 0x01, 0x00, 0x06, 0x06, 0x01, 0x00, 0x01, 0x00, 0x01, 0x00,
];

/// A CRAM file definition.
#[derive(Clone, Copy)]
pub struct FileDefinition(pub [u8; FILE_DEFINITION_LEN]);

impl FileDefinition {
    pub fn read<R: Read>(reader: &mut R) -> anyhow::Result<Self> {
        let mut buf = [0; FILE_DEFINITION_LEN];
        reader
            .read_exact(&mut buf)
            .context("reading the CRAM file definition")?;

        if &buf[..4] != b"CRAM" {
            bail!("not a CRAM file");
        }

        if !(2..=3).contains(&buf[4]) {
            bail!("unsupported CRAM version {}", buf[4]);
        }

        Ok(Self(buf))
    }

    pub fn major_version(&self) -> u8 {
        self.0[4]
    }

    /// Whether containers and blocks carry CRC32 checksums (CRAM 3).
    pub fn checksums(&self) -> Checksums {
        if self.major_version() >= 3 {
            Checksums::Crc32
        } else {
            Checksums::None
        }
    }
}

/// Whether CRAM containers and blocks end with a checksum.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub enum Checksums {
    /// CRAM 3: a CRC32 of the preceding bytes.
    Crc32,
    /// CRAM 2.
    None,
}

/// A container header.
pub struct ContainerHeader {
    pub reference_sequence_id: i32,
    pub start: i32,
    pub span: i32,
    pub record_count: i32,
    /// The bytes from the reference sequence ID through the base count, as read.
    pub fields: Vec<u8>,
    pub landmark_count: usize,
    /// The encoded header.
    pub raw: Vec<u8>,
    /// The length of the blocks that follow the header.
    pub data_len: usize,
}

impl ContainerHeader {
    /// Reads a container header, or returns `None` at the end of the stream.
    pub fn read<R: Read>(reader: &mut R, definition: &FileDefinition) -> io::Result<Option<Self>> {
        let mut raw = Vec::new();

        let mut len = [0; 4];

        match reader.read(&mut len[..1])? {
            0 => return Ok(None),
            _ => reader.read_exact(&mut len[1..])?,
        }

        raw.extend_from_slice(&len);

        // Reference sequence ID, start, span, record count.
        for _ in 0..4 {
            read_varint(reader, &mut raw, Varint::Itf8)?;
        }

        // Record counter and base count.
        let counter = if definition.major_version() >= 3 {
            Varint::Ltf8
        } else {
            Varint::Itf8
        };

        read_varint(reader, &mut raw, counter)?;
        read_varint(reader, &mut raw, Varint::Ltf8)?;

        let fields_end = raw.len();

        // Block count and landmark count.
        read_varint(reader, &mut raw, Varint::Itf8)?;
        read_varint(reader, &mut raw, Varint::Itf8)?;

        let mut src = &raw[fields_end..];
        read_itf8(&mut src)?;
        let landmark_count = to_count(read_itf8(&mut src)?)?;

        for _ in 0..landmark_count {
            read_varint(reader, &mut raw, Varint::Itf8)?;
        }

        if definition.checksums() == Checksums::Crc32 {
            let mut crc = [0; 4];
            reader.read_exact(&mut crc)?;
            raw.extend_from_slice(&crc);
        }

        let data_len = to_count(i32::from_le_bytes(len))?;

        let mut src = &raw[4..fields_end];
        let reference_sequence_id = read_itf8(&mut src)?;
        let start = read_itf8(&mut src)?;
        let span = read_itf8(&mut src)?;
        let record_count = read_itf8(&mut src)?;

        Ok(Some(Self {
            reference_sequence_id,
            start,
            span,
            record_count,
            fields: raw[4..fields_end].to_vec(),
            landmark_count,
            raw,
            data_len,
        }))
    }

    /// The encoded length of the header and its blocks.
    pub fn total_len(&self) -> usize {
        self.raw.len() + self.data_len
    }

    /// Whether this is an EOF container (HTSlib's `empty_container`).
    pub fn is_eof(&self) -> bool {
        self.record_count == 0 && self.reference_sequence_id == -1 && self.start == EOF_START
    }
}

#[derive(Clone, Copy)]
enum Varint {
    Itf8,
    Ltf8,
}

/// Reads one ITF-8 or LTF-8 integer's bytes into `raw`.
fn read_varint<R: Read>(reader: &mut R, raw: &mut Vec<u8>, kind: Varint) -> io::Result<()> {
    let mut b0 = [0; 1];
    reader.read_exact(&mut b0)?;
    raw.push(b0[0]);

    let ones = usize::try_from(b0[0].leading_ones()).unwrap_or(0);

    let extra = match kind {
        Varint::Itf8 => ones.min(4),
        Varint::Ltf8 => ones,
    };

    let start = raw.len();
    raw.resize(start + extra, 0);
    reader.read_exact(&mut raw[start..])?;

    // LTF-8 values are only checked to be well-formed.
    if let Varint::Ltf8 = kind {
        read_ltf8(&mut &raw[start - 1..])?;
    }

    Ok(())
}

fn to_count(n: i32) -> io::Result<usize> {
    usize::try_from(n).map_err(invalid_data)
}

fn to_i32(n: usize) -> io::Result<i32> {
    i32::try_from(n).map_err(invalid_data)
}

/// Encodes a header container (one gzip block) holding `header`.
pub fn encode_header_container(
    header: &sam::Header,
    definition: &FileDefinition,
) -> io::Result<Vec<u8>> {
    // Reference sequence ID, start, span, record count, record counter, and base count: all 0.
    const FIELDS: [u8; 6] = [0; 6];

    let block = encode_block(&header_data(header)?, GZIP, definition.checksums())?;
    let mut container =
        encode_container_header(&FIELDS, 1, &[], block.len(), definition.checksums())?;
    container.extend_from_slice(&block);

    Ok(container)
}

/// The header block's data: the length of the SAM header text, then the text.
pub fn header_data(header: &sam::Header) -> io::Result<Vec<u8>> {
    let mut text = Vec::new();
    sam::io::Writer::new(&mut text).write_header(header)?;

    let mut data = Vec::with_capacity(4 + text.len());
    data.extend_from_slice(&to_i32(text.len())?.to_le_bytes());
    data.extend_from_slice(&text);

    Ok(data)
}

/// Encodes a file header block, with a CRC32 for CRAM 3.
pub fn encode_block(data: &[u8], method: u8, checksums: Checksums) -> io::Result<Vec<u8>> {
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

    if checksums == Checksums::Crc32 {
        let crc = crc32fast::hash(&block);
        block.extend_from_slice(&crc.to_le_bytes());
    }

    Ok(block)
}

/// Encodes a container header from its raw fields (reference sequence ID through base count).
pub fn encode_container_header(
    fields: &[u8],
    block_count: i32,
    landmarks: &[i32],
    data_len: usize,
    checksums: Checksums,
) -> io::Result<Vec<u8>> {
    let mut buf = to_i32(data_len)?.to_le_bytes().to_vec();
    buf.extend_from_slice(fields);
    write_itf8(&mut buf, block_count);
    write_itf8(&mut buf, to_i32(landmarks.len())?);

    for &landmark in landmarks {
        write_itf8(&mut buf, landmark);
    }

    if checksums == Checksums::Crc32 {
        let crc = crc32fast::hash(&buf);
        buf.extend_from_slice(&crc.to_le_bytes());
    }

    Ok(buf)
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

    #[test]
    fn test_read_container_header() -> io::Result<()> {
        let definition = FileDefinition(*b"CRAM\x03\x00xxxxxxxxxxxxxxxxxxxx");
        let raw = encode_container_header(
            &[0x81, 0x00, 5, 7, 3, 0, 0x80, 0x80],
            2,
            &[0, 300],
            9,
            Checksums::Crc32,
        )?;

        let mut src = raw.as_slice();
        let header = ContainerHeader::read(&mut src, &definition)?.expect("header");

        assert!(src.is_empty());
        assert_eq!(header.reference_sequence_id, 0x100);
        assert_eq!(header.start, 5);
        assert_eq!(header.span, 7);
        assert_eq!(header.record_count, 3);
        assert_eq!(header.landmark_count, 2);
        assert_eq!(header.data_len, 9);
        assert_eq!(header.raw, raw);

        Ok(())
    }
}
