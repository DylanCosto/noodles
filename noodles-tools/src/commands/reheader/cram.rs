//! Replacing a CRAM file's header container, copying every other container as is.

use std::{
    fs::{File, OpenOptions},
    io::{self, Read, Seek, SeekFrom, Write},
    path::Path,
};

use anyhow::{Context, bail};
use noodles_sam as sam;

use crate::{
    cram_file::{
        ContainerHeader, FILE_DEFINITION_LEN, FileDefinition, GZIP, RAW, encode_block,
        encode_container_header, encode_header_container, header_data,
    },
    cram_raw::{read_itf8, read_len, take},
};

/// The start of a CRAM file: its file definition and header container header.
struct FileStart {
    definition: FileDefinition,
    container: ContainerHeader,
}

impl FileStart {
    fn read(file: &mut File) -> anyhow::Result<Self> {
        let definition = FileDefinition::read(file)?;
        let container = ContainerHeader::read(file, &definition)
            .context("reading the CRAM header container")?
            .context("missing CRAM header container")?;

        Ok(Self {
            definition,
            container,
        })
    }

    /// The offset of the first container after the header container.
    fn data_start(&self) -> u64 {
        (FILE_DEFINITION_LEN + self.container.total_len()) as u64
    }
}

/// Writes the input with `header` in place of its header (`cram_reheader`). The version and the
/// containers are copied as they are.
pub fn reheader<W: Write>(input: &Path, header: &sam::Header, dst: &mut W) -> anyhow::Result<()> {
    let mut file = File::open(input).with_context(|| format!("opening {}", input.display()))?;
    let start = FileStart::read(&mut file)?;

    dst.write_all(&start.definition.0)?;
    dst.write_all(&encode_header_container(header, &start.definition)?)?;

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

    if start.definition.major_version() == 2 {
        reheader_in_place_2(&mut file, &start, &data)
    } else {
        reheader_in_place_3(&mut file, &start, &data)
    }
}

/// CRAM 2.1 keeps the header in one uncompressed block, padded with NULs for growth; the new
/// header must fit in that block.
fn reheader_in_place_2(file: &mut File, start: &FileStart, data: &[u8]) -> anyhow::Result<()> {
    /// Enough for the block's method, content type, content ID, and sizes.
    const BLOCK_HEADER_MAX_LEN: usize = 17;

    let block_start = (FILE_DEFINITION_LEN + start.container.raw.len()) as u64;

    let mut buf = vec![0; start.container.data_len.min(BLOCK_HEADER_MAX_LEN)];
    file.read_exact(&mut buf)?;

    let mut src = buf.as_slice();
    let method = take(&mut src, 2)?[0];
    read_itf8(&mut src)?;
    let compressed_len = read_len(&mut src)?;
    let len = read_len(&mut src)?;
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
    let total_len = old.total_len();

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

    file.seek(SeekFrom::Start(FILE_DEFINITION_LEN as u64))?;
    file.write_all(&buf)?;

    Ok(())
}
