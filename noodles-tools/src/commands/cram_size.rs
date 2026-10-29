//! `ntools cram-size`: CRAM block sizes per content ID and data series (samtools cram-size).
//!
//! This follows samtools 1.24's `cram_size.c`: the compressed and uncompressed sizes of every
//! slice data block are summed by content ID and by compression method, expanded with the
//! parameters stored at the start of the compressed data (gzip level, rANS order, ...). Content
//! IDs are labelled with the data series and tags whose encodings use them.

use std::{
    collections::{BTreeMap, HashMap},
    fs::File,
    io::{self, BufReader, BufWriter, Read, Write},
    path::PathBuf,
};

use anyhow::Context;
use noodles_cram::{
    self as cram,
    container::block::{CompressionMethod, ContentType},
};

use crate::cram_raw::{invalid_data, read_itf8, read_len, read_slice_header, take};

#[derive(Debug, clap::Args)]
pub struct Args {
    /// Input CRAM file [default: stdin].
    input: Option<PathBuf>,

    /// Write the report to this file.
    #[arg(short = 'o', long = "output")]
    output: Option<PathBuf>,

    /// Report each compression method on its own line.
    #[arg(short = 'v', long = "verbose")]
    verbose: bool,

    /// Describe each container's encodings.
    #[arg(short = 'e', long = "encodings")]
    encodings: bool,
}

/// The number of expanded compression methods (samtools' `COMP_MAX`).
const METHOD_COUNT: usize = 53;

/// One-character method names, indexed by expanded method (samtools' `comp_method2char`).
const METHOD_CHARS: &[u8; METHOD_COUNT] = b".gblr0afn_GbbbbbbbbBrR010101014545454582aAaAaAaAaaanN";

/// Method names, indexed by expanded method (samtools' `comp_method2str`).
const METHOD_NAMES: [&str; METHOD_COUNT] = [
    "raw",
    "gzip",
    "bzip2",
    "lzma",
    "r4x8",
    "rNx16",
    "arith",
    "fqzcomp",
    "tok3",
    "gzip-min",
    "gzip-max",
    "bzip2-1",
    "bzip2-2",
    "bzip2-3",
    "bzip2-4",
    "bzip2-5",
    "bzip2-6",
    "bzip2-7",
    "bzip2-8",
    "bzip2-9",
    "r4x8-o0",
    "r4x8-o1",
    "r4x16-o0",
    "r4x16-o1",
    "r4x16-o0R",
    "r4x16-o1R",
    "r4x16-o0P",
    "r4x16-o1P",
    "r4x16-o0PR",
    "r4x16-o1PR",
    "r32x16-o0",
    "r32x16-o1",
    "r32x16-o0R",
    "r32x16-o1R",
    "r32x16-o0P",
    "r32x16-o1P",
    "r32x16-o0PR",
    "r32x16-o1PR",
    "rNx16-xo0",
    "rNx16-cat",
    "arith-o0",
    "arith-o1",
    "arith-o0R",
    "arith-o1R",
    "arith-o0P",
    "arith-o1P",
    "arith-o0PR",
    "arith-o1PR",
    "arith-stripe",
    "arith-cat",
    "arith-ext",
    "tok3-rans",
    "tok3-arith",
];

// Expanded method indices (samtools' `enum comp_expanded`).
const GZIP_1: usize = 9;
const GZIP_9: usize = 10;
const BZIP2_1: usize = 11;
const RANS4X8_O0: usize = 20;
const RANS4X16_O0: usize = 22;
const RANSNX16_STRIPE: usize = 38;
const RANSNX16_CAT: usize = 39;
const ARITH_O0: usize = 40;
const ARITH_STRIPE: usize = 48;
const ARITH_CAT: usize = 49;
const ARITH_EXT: usize = 50;
const TOK3_RANS: usize = 51;
const TOK3_ARITH: usize = 52;

// rANS Nx16 and arithmetic coder format flags (htscodecs' `RANS_ORDER_*`).
const ORDER_X32: u8 = 0x04;
const ORDER_STRIPE: u8 = 0x08;
const ORDER_CAT: u8 = 0x20;
const ORDER_RLE: u8 = 0x40;
const ORDER_PACK: u8 = 0x80;
/// The arithmetic coder's external entropy coder flag (`X_EXT`).
const ARITH_EXTERNAL: u8 = 0x04;

/// The content ID reported for core data blocks.
const CORE: i32 = -1;

pub fn run(args: Args) -> anyhow::Result<()> {
    let input: Box<dyn Read> = match &args.input {
        Some(path) if path.as_os_str() != "-" => {
            Box::new(File::open(path).with_context(|| format!("opening {}", path.display()))?)
        }
        _ => Box::new(io::stdin().lock()),
    };

    let mut reader = cram::io::Reader::new(BufReader::new(CountingReader::new(input)));
    reader.read_header().context("reading CRAM header")?;

    let mut out: Box<dyn Write> = match &args.output {
        Some(path) => Box::new(BufWriter::new(
            File::create(path).with_context(|| format!("creating {}", path.display()))?,
        )),
        None => Box::new(BufWriter::new(io::stdout().lock())),
    };

    let mut stats = Stats::default();
    let mut container = cram::io::reader::Container::default();

    loop {
        match reader.read_container(&mut container) {
            Ok(0) => break,
            Ok(_) => {}
            Err(e) if e.kind() == io::ErrorKind::UnexpectedEof => {
                eprintln!("[W::cram-size] EOF marker is absent; the input may be truncated");
                break;
            }
            Err(e) => return Err(e).context("reading container"),
        }

        let encodings = stats.add_container(&container)?;

        if args.encodings {
            writeln!(out, "Container encodings")?;
            for (key, codec) in &encodings {
                writeln!(out, "\t{}\t{codec}", DataSeries(*key))?;
            }
            writeln!(out)?;
        }
    }

    // Count the rest of the input (the EOF container) to report the file size.
    io::copy(reader.get_mut(), &mut io::sink())?;
    let file_size = reader.get_ref().get_ref().count;

    let compressed_size = stats.write_report(&mut out, args.verbose)?;

    writeln!(out)?;
    writeln!(out, "Number of containers  {:18}", stats.containers)?;
    writeln!(out, "Number of slices      {:18}", stats.slices)?;
    writeln!(out, "Number of sequences   {:18}", stats.records)?;
    writeln!(out, "Number of bases       {:18}", stats.bases)?;
    writeln!(out, "Total file size       {file_size:18}")?;
    writeln!(
        out,
        "Format overhead size  {:18}",
        i128::from(file_size) - i128::from(compressed_size)
    )?;

    out.flush()?;

    Ok(())
}

/// A reader that counts the bytes read through it.
struct CountingReader<R> {
    inner: R,
    count: u64,
}

impl<R> CountingReader<R> {
    fn new(inner: R) -> Self {
        Self { inner, count: 0 }
    }
}

impl<R: Read> Read for CountingReader<R> {
    fn read(&mut self, buf: &mut [u8]) -> io::Result<usize> {
        let n = self.inner.read(buf)?;
        self.count += n as u64;
        Ok(n)
    }
}

/// Summed compressed and uncompressed sizes, by expanded method.
struct Sizes {
    compressed: [u64; METHOD_COUNT],
    uncompressed: [u64; METHOD_COUNT],
}

impl Default for Sizes {
    fn default() -> Self {
        Self {
            compressed: [0; METHOD_COUNT],
            uncompressed: [0; METHOD_COUNT],
        }
    }
}

impl Sizes {
    /// Returns the method indices ordered by compressed size, largest first, then by index.
    fn methods_by_size(&self) -> [usize; METHOD_COUNT] {
        let mut order: [usize; METHOD_COUNT] = std::array::from_fn(|i| i);
        order.sort_by(|&a, &b| self.compressed[b].cmp(&self.compressed[a]).then(a.cmp(&b)));
        order
    }
}

#[derive(Default)]
struct Stats {
    containers: u64,
    slices: u64,
    records: u64,
    bases: u64,
    /// Block sizes by content ID (`CORE` for core data blocks).
    sizes: BTreeMap<i32, Sizes>,
    /// Data series keys by content ID, in first-seen order.
    data_series: HashMap<i32, Vec<u32>>,
    /// The content ID of the embedded reference, if any.
    embedded_reference: Option<i32>,
}

impl Stats {
    /// Adds a container's blocks and returns its encodings in samtools' iteration order.
    fn add_container(
        &mut self,
        container: &cram::io::reader::Container,
    ) -> anyhow::Result<Vec<(u32, Codec)>> {
        let header = container.header();
        self.containers += 1;
        self.records += header.record_count() as u64;
        self.bases += header.base_count();

        let mut src = container.src();

        let block = cram::io::reader::read_block(&mut src).context("reading compression header")?;
        let data = block.decode().context("decoding compression header")?;
        let encodings =
            read_encodings(&mut &data[..]).context("reading compression header encodings")?;

        for (key, codec) in &encodings {
            for id in codec.content_ids().into_iter().flatten() {
                let keys = self.data_series.entry(id).or_default();
                if !keys.contains(key) {
                    keys.push(*key);
                }
            }
        }

        // The container block count can be wrong (htslib bug), so the slice headers' block
        // counts are used instead.
        for _ in header.landmarks() {
            self.slices += 1;

            let block = cram::io::reader::read_block(&mut src).context("reading slice header")?;
            let data = block.decode().context("decoding slice header")?;
            let slice_header = read_slice_header(&mut &data[..]).context("reading slice header")?;

            if let Some(id) = slice_header.embedded_reference {
                match self.embedded_reference {
                    None => self.embedded_reference = Some(id),
                    Some(used) if used != id => eprintln!(
                        "Embedded reference is not consistently using the same Content-Id.\n\
                         Reported figures for reference will be invalid."
                    ),
                    Some(_) => {}
                }
            }

            for _ in 0..slice_header.block_count {
                let block = cram::io::reader::read_block(&mut src).context("reading block")?;

                let id = if block.content_type == ContentType::CoreData {
                    CORE
                } else {
                    block.content_id
                };

                let method = expand_method(block.compression_method, block.src);
                let sizes = self.sizes.entry(id).or_default();
                sizes.compressed[method] += block.src.len() as u64;
                sizes.uncompressed[method] += block.uncompressed_size as u64;
            }
        }

        Ok(encodings)
    }

    /// Writes the per-content-ID table and returns the total compressed size.
    fn write_report<W: Write>(&self, out: &mut W, verbose: bool) -> io::Result<u64> {
        writeln!(
            out,
            "#   Content_ID  Uncomp.size    Comp.size   Ratio Method{}  Data_series",
            if verbose { "    " } else { "" }
        )?;

        let mut total = 0;

        for (&id, sizes) in &self.sizes {
            let order = sizes.methods_by_size();

            if verbose {
                for (i, &method) in order.iter().enumerate() {
                    if i > 0 && sizes.compressed[method] == 0 {
                        break;
                    }

                    if i > 0 {
                        writeln!(out)?;
                    }

                    write_sizes(
                        out,
                        id,
                        sizes.uncompressed[method],
                        sizes.compressed[method],
                    )?;
                    write!(out, " {:<11}", METHOD_NAMES[method])?;
                    self.write_data_series(out, id)?;
                }
            } else {
                let methods: String = order
                    .iter()
                    .take_while(|&&method| sizes.compressed[method] != 0)
                    .map(|&method| char::from(METHOD_CHARS[method]))
                    .collect();

                write_sizes(
                    out,
                    id,
                    sizes.uncompressed.iter().sum(),
                    sizes.compressed.iter().sum(),
                )?;
                let methods = if methods.is_empty() { "." } else { &methods };
                write!(out, " {methods:<7}")?;
                self.write_data_series(out, id)?;
            }

            if id >= 0 && self.embedded_reference == Some(id) {
                write!(out, " embedded_ref")?;
            }

            writeln!(out)?;

            total += sizes.compressed.iter().sum::<u64>();
        }

        Ok(total)
    }

    /// Writes the data series using a content ID, most recently seen first.
    fn write_data_series<W: Write>(&self, out: &mut W, id: i32) -> io::Result<()> {
        if let Some(keys) = self.data_series.get(&id) {
            for &key in keys.iter().rev() {
                write!(out, " {}", DataSeries(key))?;
            }
        }

        Ok(())
    }
}

/// Writes the content ID, sizes, and compression ratio columns.
fn write_sizes<W: Write>(
    out: &mut W,
    id: i32,
    uncompressed: u64,
    compressed: u64,
) -> io::Result<()> {
    if id < 0 {
        write!(out, "BLOCK {:>8}", "CORE")?;
    } else {
        write!(out, "BLOCK {id:8}")?;
    }

    write!(out, " {uncompressed:12} {compressed:12}")?;

    let ratio = 100.0 * (compressed as f64 + 0.0001) / (uncompressed as f64 + 0.0001);

    if ratio > 999.0 {
        write!(out, "   >999%")
    } else {
        write!(out, " {ratio:6.2}%")
    }
}

/// Expands a block's compression method with the parameters at the start of its data
/// (htslib's `cram_expand_method` and samtools' `comp_method2expanded`).
fn expand_method(method: CompressionMethod, data: &[u8]) -> usize {
    match method {
        CompressionMethod::Gzip => match data.get(8) {
            Some(4) => GZIP_1,
            Some(2) => GZIP_9,
            _ => 1,
        },
        CompressionMethod::Bzip2 => match data.get(3) {
            Some(&level @ b'1'..=b'9') => BZIP2_1 + usize::from(level - b'1'),
            _ => 2,
        },
        CompressionMethod::Rans4x8 => RANS4X8_O0 + usize::from(data.first() == Some(&1)),
        CompressionMethod::RansNx16 => {
            let flags = data.first().copied().unwrap_or(0);

            if flags & ORDER_STRIPE != 0 {
                RANSNX16_STRIPE
            } else if flags & ORDER_CAT != 0 {
                RANSNX16_CAT
            } else {
                RANS4X16_O0
                    + usize::from(flags & 1)
                    + transform_offset(flags)
                    + if flags & ORDER_X32 != 0 { 8 } else { 0 }
            }
        }
        CompressionMethod::AdaptiveArithmeticCoding => {
            let flags = data.first().copied().unwrap_or(0);

            if flags & ORDER_STRIPE != 0 {
                ARITH_STRIPE
            } else if flags & ORDER_CAT != 0 {
                ARITH_CAT
            } else if flags & ARITH_EXTERNAL != 0 {
                ARITH_EXT
            } else {
                // As in samtools, both order bits are added, so order 2 and 3 alias other methods.
                ARITH_O0 + usize::from(flags & 3) + transform_offset(flags)
            }
        }
        CompressionMethod::NameTokenizer => match data.get(8) {
            Some(1) => TOK3_ARITH,
            _ => TOK3_RANS,
        },
        CompressionMethod::Lzma => 3,
        CompressionMethod::Fqzcomp => 7,
        _ => 0,
    }
}

/// Returns the method index offset of the RLE and PACK transforms.
fn transform_offset(flags: u8) -> usize {
    (if flags & ORDER_RLE != 0 { 2 } else { 0 }) + if flags & ORDER_PACK != 0 { 4 } else { 0 }
}

/// A data series or tag key, displayed as its two or three characters.
#[derive(Clone, Copy)]
struct DataSeries(u32);

impl std::fmt::Display for DataSeries {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        let [_, a, b, c] = self.0.to_be_bytes();

        if a != 0 {
            write!(f, "{}", char::from(a))?;
        }

        write!(f, "{}{}", char::from(b), char::from(c))
    }
}

/// A decoded encoding (htslib's `cram_codec`), as far as cram-size describes it.
#[derive(Debug)]
enum Codec {
    External(i32),
    Huffman(Vec<HuffmanCode>),
    ByteArrayLen { len: Box<Codec>, value: Box<Codec> },
    ByteArrayStop { stop: u8, id: i32 },
    Beta { offset: i32, bits: i32 },
    Subexp { offset: i32, k: i32 },
    Gamma { offset: i32 },
}

#[derive(Debug)]
struct HuffmanCode {
    symbol: i32,
    len: i32,
}

impl Codec {
    /// Returns the content IDs this codec reads (`CORE` for the core block; htslib's
    /// `cram_codec_to_id`).
    fn content_ids(&self) -> [Option<i32>; 2] {
        match self {
            Self::External(id) | Self::ByteArrayStop { id, .. } => [Some(*id), None],
            Self::Huffman(codes) if codes.len() == 1 => [None, None],
            Self::Huffman(_) | Self::Beta { .. } | Self::Subexp { .. } | Self::Gamma { .. } => {
                [Some(CORE), None]
            }
            Self::ByteArrayLen { len, value } => [len.content_ids()[0], value.content_ids()[0]],
        }
    }
}

impl std::fmt::Display for Codec {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Self::External(id) => write!(f, "EXTERNAL(id={id})"),
            // htslib leaves an empty HUFFMAN codec undescribed.
            Self::Huffman(codes) if codes.is_empty() => f.write_str("?"),
            Self::Huffman(codes) => {
                f.write_str("HUFFMAN(codes={")?;
                for (i, code) in codes.iter().enumerate() {
                    let sep = if i > 0 { "," } else { "" };
                    write!(f, "{sep}{}", code.symbol)?;
                }
                f.write_str("},lengths={")?;
                for (i, code) in codes.iter().enumerate() {
                    let sep = if i > 0 { "," } else { "" };
                    write!(f, "{sep}{}", code.len)?;
                }
                f.write_str("})")
            }
            // htslib omits the closing parenthesis.
            Self::ByteArrayLen { len, value } => {
                write!(
                    f,
                    "BYTE_ARRAY_LEN(len_codec={{{len}}},val_codec={{{value}}}"
                )
            }
            Self::ByteArrayStop { stop, id } => write!(f, "BYTE_ARRAY_STOP(stop={stop},id={id})"),
            Self::Beta { offset, bits } => write!(f, "BETA(offset={offset}, nbits={bits})"),
            Self::Subexp { offset, k } => write!(f, "SUBEXP(offset={offset},k={k})"),
            Self::Gamma { offset } => write!(f, "GAMMA(offset={offset})"),
        }
    }
}

/// Data series in htslib's `enum cram_DS_ID` order, which is the order samtools reports them.
const DATA_SERIES: [&[u8; 2]; 29] = [
    b"RN", b"QS", b"IN", b"SC", b"BF", b"CF", b"AP", b"RG", b"MQ", b"NS", b"MF", b"TS", b"NP",
    b"NF", b"RL", b"FN", b"FC", b"FP", b"DL", b"BA", b"BS", b"TL", b"RI", b"RS", b"PD", b"HC",
    b"BB", b"QQ", b"TN",
];

/// The tag key for CRAM 1.0's TC data series, reported after the others.
const LEGACY_DATA_SERIES: &[u8; 2] = b"TC";

/// The number of htslib's tag encoding hash buckets (`CRAM_MAP_HASH`).
const TAG_BUCKETS: usize = 32;

/// Reads the data series and tag encodings of a compression header, in htslib's iteration
/// order: data series in `enum cram_DS_ID` order (the last encoding of a repeated key wins),
/// then tags by hash bucket, each bucket in reverse file order.
fn read_encodings(src: &mut &[u8]) -> io::Result<Vec<(u32, Codec)>> {
    // Preservation map
    let len = read_len(src)?;
    take(src, len)?;

    let len = read_len(src)?;
    let mut map = take(src, len)?;
    let count = read_len(&mut map)?;

    let mut data_series: [Option<Codec>; DATA_SERIES.len() + 1] = Default::default();

    for _ in 0..count {
        let key = take(&mut map, 2)?;
        let Some(codec) = read_codec(&mut map)? else {
            continue;
        };

        let i = if key == LEGACY_DATA_SERIES {
            Some(DATA_SERIES.len())
        } else {
            DATA_SERIES.iter().position(|ds| ds[..] == *key)
        };

        if let Some(i) = i {
            data_series[i] = Some(codec);
        }
    }

    let len = read_len(src)?;
    let mut map = take(src, len)?;
    let count = read_len(&mut map)?;

    let mut tags: [Vec<(u32, Codec)>; TAG_BUCKETS] = Default::default();

    for _ in 0..count {
        let key = read_itf8(&mut map)? as u32;
        let codec = read_codec(&mut map)?
            .ok_or_else(|| invalid_data("tag encoding uses the NULL codec"))?;
        let [_, a, b, _] = key.to_be_bytes();
        let bucket = (usize::from(a) * 3 + usize::from(b)) % TAG_BUCKETS;
        tags[bucket].push((key, codec));
    }

    let keys = DATA_SERIES.iter().copied().chain([LEGACY_DATA_SERIES]);
    let mut encodings: Vec<(u32, Codec)> = keys
        .zip(data_series)
        .filter_map(|(key, codec)| Some((u32::from(u16::from_be_bytes(*key)), codec?)))
        .collect();

    for bucket in tags {
        encodings.extend(bucket.into_iter().rev());
    }

    Ok(encodings)
}

/// Reads an encoding; the NULL codec reads as `None`.
fn read_codec(src: &mut &[u8]) -> io::Result<Option<Codec>> {
    let id = read_itf8(src)?;
    let len = read_len(src)?;
    let mut params = take(src, len)?;
    let src = &mut params;

    let codec = match id {
        0 => return Ok(None),
        1 => Codec::External(read_itf8(src)?),
        3 => {
            let n = read_len(src)?;
            let symbols = (0..n)
                .map(|_| read_itf8(src))
                .collect::<io::Result<Vec<_>>>()?;

            if read_len(src)? != n {
                return Err(invalid_data("HUFFMAN symbol and length counts differ"));
            }

            let mut codes = symbols
                .into_iter()
                .map(|symbol| read_itf8(src).map(|len| HuffmanCode { symbol, len }))
                .collect::<io::Result<Vec<_>>>()?;
            codes.sort_by_key(|code| (code.len, code.symbol));
            Codec::Huffman(codes)
        }
        4 => {
            let len = read_codec(src)?.ok_or_else(|| invalid_data("NULL length codec"))?;
            let value = read_codec(src)?.ok_or_else(|| invalid_data("NULL value codec"))?;
            Codec::ByteArrayLen {
                len: Box::new(len),
                value: Box::new(value),
            }
        }
        5 => {
            let stop = take(src, 1)?[0];
            let id = read_itf8(src)?;
            Codec::ByteArrayStop { stop, id }
        }
        6 => {
            let offset = read_itf8(src)?;
            let bits = read_itf8(src)?;
            Codec::Beta { offset, bits }
        }
        7 => {
            let offset = read_itf8(src)?;
            let k = read_itf8(src)?;
            Codec::Subexp { offset, k }
        }
        9 => Codec::Gamma {
            offset: read_itf8(src)?,
        },
        _ => return Err(invalid_data(format!("unsupported codec {id}"))),
    };

    Ok(Some(codec))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_expand_method() {
        let name = |method, data: &[u8]| METHOD_NAMES[expand_method(method, data)];

        assert_eq!(name(CompressionMethod::None, &[]), "raw");
        assert_eq!(name(CompressionMethod::Gzip, &[0; 9]), "gzip");
        assert_eq!(
            name(CompressionMethod::Gzip, &[0, 0, 0, 0, 0, 0, 0, 0, 4]),
            "gzip-min"
        );
        assert_eq!(name(CompressionMethod::Bzip2, b"BZh9"), "bzip2-9");
        assert_eq!(name(CompressionMethod::Rans4x8, &[1]), "r4x8-o1");
        assert_eq!(name(CompressionMethod::RansNx16, &[0x45]), "r32x16-o1R");
        assert_eq!(name(CompressionMethod::RansNx16, &[0x08]), "rNx16-xo0");
        assert_eq!(
            name(CompressionMethod::AdaptiveArithmeticCoding, &[0x81]),
            "arith-o1P"
        );
        assert_eq!(name(CompressionMethod::NameTokenizer, &[0; 9]), "tok3-rans");
    }

    #[test]
    fn test_read_encodings_order() -> io::Result<()> {
        fn external(key: &[u8], id: u8) -> Vec<u8> {
            [key, &[1, 1, id]].concat()
        }

        let ds = [&[2u8][..], &external(b"RL", 1), &external(b"BF", 2)].concat();
        // Both XAZ encodings share a hash bucket, so they come out in reverse file order.
        let tags = [
            &[3u8][..],
            &[0xe0, 0x58, 0x41, 0x5a, 1, 1, 3],
            &[0xe0, 0x58, 0x42, 0x5a, 1, 1, 4],
            &[0xe0, 0x58, 0x41, 0x5a, 1, 1, 5],
        ]
        .concat();
        let src = [&[0, ds.len() as u8][..], &ds, &[tags.len() as u8], &tags].concat();

        let encodings = read_encodings(&mut &src[..])?;
        let keys: Vec<_> = encodings
            .iter()
            .map(|(key, codec)| format!("{}:{codec}", DataSeries(*key)))
            .collect();

        assert_eq!(
            keys,
            [
                "BF:EXTERNAL(id=2)",
                "RL:EXTERNAL(id=1)",
                "XAZ:EXTERNAL(id=5)",
                "XAZ:EXTERNAL(id=3)",
                "XBZ:EXTERNAL(id=4)",
            ]
        );

        Ok(())
    }
}
