//! Opening alignment inputs from paths or stdin and detecting their format.

use std::{
    fs::File,
    io::{self, BufRead, BufReader, Read, Seek, SeekFrom},
    path::PathBuf,
};

use anyhow::Context;
use flate2::bufread::MultiGzDecoder;

/// Large enough to hold a whole BGZF block, so format detection can always see the first one.
const BUFFER_CAPACITY: usize = 1 << 17;

/// The path that names stdin, as in samtools.
pub const STDIN: &str = "-";

/// An alignment file format.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub enum Format {
    Sam,
    Bam,
    Cram,
}

/// Whether a stream is BGZF-compressed.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub enum Compression {
    None,
    Bgzf,
}

/// A raw input stream: a file, or stdin.
pub enum Stream {
    File(File),
    Stdin(io::Stdin),
}

impl Read for Stream {
    fn read(&mut self, buf: &mut [u8]) -> io::Result<usize> {
        match self {
            Self::File(f) => f.read(buf),
            Self::Stdin(s) => s.read(buf),
        }
    }
}

impl Seek for Stream {
    fn seek(&mut self, pos: SeekFrom) -> io::Result<u64> {
        match self {
            Self::File(f) => f.seek(pos),
            Self::Stdin(_) => Err(io::Error::new(
                io::ErrorKind::Unsupported,
                "stdin is not seekable",
            )),
        }
    }
}

/// An opened input with its detected format.
pub struct Input {
    pub path: Option<PathBuf>,
    pub reader: BufReader<Stream>,
    pub format: Format,
    pub compression: Compression,
}

impl Input {
    /// Opens `src` (`-` for stdin) and detects its format from its leading bytes.
    pub fn open(src: &str) -> anyhow::Result<Self> {
        let (path, stream) = if src == STDIN {
            (None, Stream::Stdin(io::stdin()))
        } else {
            let file = File::open(src).with_context(|| format!("opening {src}"))?;
            (Some(PathBuf::from(src)), Stream::File(file))
        };

        let mut reader = BufReader::with_capacity(BUFFER_CAPACITY, stream);
        let (format, compression) =
            detect(&mut reader).with_context(|| format!("detecting the format of {src}"))?;

        Ok(Self {
            path,
            reader,
            format,
            compression,
        })
    }
}

/// Detects the format and compression of a buffered stream without consuming it.
pub fn detect<R>(reader: &mut R) -> io::Result<(Format, Compression)>
where
    R: BufRead,
{
    const GZIP_MAGIC_NUMBER: [u8; 2] = [0x1f, 0x8b];
    const BAM_MAGIC_NUMBER: [u8; 4] = *b"BAM\x01";
    const CRAM_MAGIC_NUMBER: [u8; 4] = *b"CRAM";

    let src = reader.fill_buf()?;

    if src.starts_with(&GZIP_MAGIC_NUMBER) {
        let mut decoder = MultiGzDecoder::new(src);
        let mut buf = [0; BAM_MAGIC_NUMBER.len()];

        let format = match decoder.read_exact(&mut buf) {
            Ok(()) if buf == BAM_MAGIC_NUMBER => Format::Bam,
            _ => Format::Sam,
        };

        Ok((format, Compression::Bgzf))
    } else if src.starts_with(&CRAM_MAGIC_NUMBER) {
        Ok((Format::Cram, Compression::None))
    } else if src.starts_with(&BAM_MAGIC_NUMBER) {
        bail_io("uncompressed (raw) BAM streams are not supported")
    } else {
        Ok((Format::Sam, Compression::None))
    }
}

fn bail_io<T>(msg: &str) -> io::Result<T> {
    Err(io::Error::new(io::ErrorKind::InvalidData, msg.to_owned()))
}

#[cfg(test)]
mod tests {
    use std::io::Write;

    use noodles_bgzf as bgzf;

    use super::*;

    fn bgzf_of(data: &[u8]) -> io::Result<Vec<u8>> {
        let mut writer = bgzf::io::Writer::new(Vec::new());
        writer.write_all(data)?;
        writer.finish()
    }

    #[test]
    fn test_detect() -> io::Result<()> {
        let cases: [(Vec<u8>, Format, Compression); 5] = [
            (b"@HD\tVN:1.6\n".to_vec(), Format::Sam, Compression::None),
            (Vec::new(), Format::Sam, Compression::None),
            (b"CRAM\x03\x01".to_vec(), Format::Cram, Compression::None),
            (bgzf_of(b"BAM\x01\0\0\0\0")?, Format::Bam, Compression::Bgzf),
            (bgzf_of(b"@HD\tVN:1.6\n")?, Format::Sam, Compression::Bgzf),
        ];

        for (src, format, compression) in cases {
            let got = detect(&mut &src[..])?;
            assert_eq!(got, (format, compression), "input {src:?}");
        }

        Ok(())
    }
}
