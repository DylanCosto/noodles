//! `ntools fastq` and `ntools fasta`: convert alignments to reads (samtools fastq/fasta).

use std::{
    fs::File,
    io::{self, BufWriter, Write},
    path::{Path, PathBuf},
};

use anyhow::Context;
use noodles_bgzf::io::writer::CompressionLevel;
use noodles_sam::alignment::record::{
    Flags,
    data::field::{Tag, Value},
};

use crate::{
    alignment::{AlignmentReader, ReadOptions, Rec},
    bgzf_stream::BgzfWriter,
    filter::Filter,
    flags,
    sequence::reverse_complement,
};

const BUFFER_CAPACITY: usize = 1 << 17;

/// samtools' default exclusion: secondary and supplementary alignments.
const DEFAULT_EXCLUDE_FLAGS: Flags = Flags::SECONDARY.union(Flags::SUPPLEMENTARY);

/// Tags copied by `-t`.
const DEFAULT_COPIED_TAGS: [[u8; 2]; 3] = [*b"RG", *b"BC", *b"QT"];

/// Which read format to write.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub enum ReadFormat {
    Fastq,
    Fasta,
}

#[derive(Debug, clap::Args)]
pub struct Args {
    /// Input SAM, BAM, or CRAM file (`-` for stdin).
    input: String,

    /// Write READ1 or READ2 reads (and others not sent elsewhere) to this file.
    #[arg(short = 'o')]
    output: Option<PathBuf>,

    /// Write reads with neither or both of READ1 and READ2 to this file.
    #[arg(short = '0')]
    other: Option<PathBuf>,

    /// Write READ1 reads to this file.
    #[arg(short = '1')]
    read1: Option<PathBuf>,

    /// Write READ2 reads to this file.
    #[arg(short = '2')]
    read2: Option<PathBuf>,

    /// Write singletons (reads whose mate is not adjacent) to this file.
    #[arg(short = 's')]
    singletons: Option<PathBuf>,

    /// Never append /1 and /2 to read names.
    #[arg(short = 'n', conflicts_with = "suffix_always")]
    suffix_never: bool,

    /// Always append /1 and /2 to read names.
    #[arg(short = 'N')]
    suffix_always: bool,

    /// Copy RG, BC, and QT tags to the name line.
    #[arg(short = 't')]
    copy_default_tags: bool,

    /// Copy these tags to the name line (comma-separated).
    #[arg(short = 'T', value_delimiter = ',')]
    copy_tags: Vec<String>,

    /// Only include reads with all of these flags set.
    #[arg(short = 'f', value_parser = flags::parse)]
    require_flags: Option<Flags>,

    /// Exclude reads with any of these flags set [default: 0x900].
    #[arg(short = 'F', value_parser = flags::parse)]
    exclude_flags: Option<Flags>,

    /// Exclude reads with all of these flags set.
    #[arg(short = 'G', value_parser = flags::parse)]
    exclude_all_flags: Option<Flags>,

    /// Quality score for reads without qualities.
    #[arg(short = 'v', default_value_t = 1)]
    default_quality: u8,

    /// Compression level for `.gz` outputs.
    #[arg(short = 'c', default_value_t = 1)]
    level: u8,

    /// Reference FASTA (indexed) for CRAM input.
    #[arg(long = "reference")]
    reference: Option<PathBuf>,

    /// Additional threads for decompression and compression.
    #[arg(short = '@', long = "threads", default_value_t = 0)]
    threads: usize,
}

/// A read held until its mate can be identified.
#[derive(Debug, Default)]
struct Read {
    name: Vec<u8>,
    flags: Flags,
    sequence: Vec<u8>,
    quality_scores: Vec<u8>,
    tags: Vec<u8>,
}

impl Read {
    fn fill(
        &mut self,
        record: &Rec<'_>,
        copied_tags: &[[u8; 2]],
        default_quality: u8,
    ) -> io::Result<()> {
        let r = record.as_dyn();

        self.name.clear();
        self.name
            .extend_from_slice(r.name().map_or(b"*", |name| name.as_ref()));

        self.flags = r.flags()?;

        self.sequence.clear();
        self.quality_scores.clear();

        match record {
            Rec::Bam(r) => {
                self.sequence.extend(r.sequence().iter());
                self.quality_scores
                    .extend_from_slice(r.quality_scores().as_ref());
            }
            Rec::Buf(r) => {
                self.sequence.extend_from_slice(r.sequence().as_ref());
                self.quality_scores
                    .extend_from_slice(r.quality_scores().as_ref());
            }
            _ => {
                self.sequence.extend(r.sequence().iter());

                for result in r.quality_scores().iter() {
                    self.quality_scores.push(result?);
                }
            }
        }

        // BAM marks missing qualities with 0xff in every position.
        if self.quality_scores.len() != self.sequence.len()
            || self.quality_scores.first() == Some(&0xff)
        {
            self.quality_scores.clear();
            self.quality_scores
                .resize(self.sequence.len(), default_quality);
        }

        if self.flags.is_reverse_complemented() {
            reverse_complement(&mut self.sequence);
            self.quality_scores.reverse();
        }

        self.tags.clear();

        if !copied_tags.is_empty() {
            let data = r.data();

            for &tag in copied_tags {
                if let Some(value) = data.get(&Tag::from(tag)).transpose()? {
                    write_tag(&mut self.tags, tag, &value)?;
                }
            }
        }

        Ok(())
    }

    fn segment(&self) -> Segment {
        match (self.flags.is_first_segment(), self.flags.is_last_segment()) {
            (true, false) => Segment::First,
            (false, true) => Segment::Last,
            _ => Segment::Other,
        }
    }
}

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
enum Segment {
    First,
    Last,
    Other,
}

/// A read output: plain or BGZF-compressed, by file extension.
enum Sink {
    Plain(BufWriter<Box<dyn Write + Send>>),
    Bgzf(BgzfWriter<File>),
}

impl Sink {
    fn create(
        path: Option<&Path>,
        level: CompressionLevel,
        threads: usize,
    ) -> anyhow::Result<Self> {
        let Some(path) = path.filter(|p| *p != Path::new("-")) else {
            return Ok(Self::Plain(BufWriter::with_capacity(
                BUFFER_CAPACITY,
                crate::output::raw_stdout()?,
            )));
        };

        let file = File::create(path).with_context(|| format!("creating {}", path.display()))?;
        let is_compressed = path
            .extension()
            .is_some_and(|ext| ext == "gz" || ext == "bgz" || ext == "bgzf");

        Ok(if is_compressed {
            Self::Bgzf(BgzfWriter::new(file, level, threads))
        } else {
            Self::Plain(BufWriter::with_capacity(BUFFER_CAPACITY, Box::new(file)))
        })
    }

    fn writer(&mut self) -> &mut dyn Write {
        match self {
            Self::Plain(w) => w,
            Self::Bgzf(w) => w,
        }
    }

    fn finish(self) -> io::Result<()> {
        match self {
            Self::Plain(mut w) => w.flush(),
            Self::Bgzf(w) => w.finish().map(|_| ()),
        }
    }
}

/// Indices into the list of sinks for each kind of read.
struct Routes {
    first: usize,
    last: usize,
    other: usize,
    singleton: Option<usize>,
}

pub fn run(args: Args, format: ReadFormat) -> anyhow::Result<()> {
    let options = ReadOptions {
        threads: args.threads,
        reference: args.reference.clone(),
    };

    let (mut reader, _) = AlignmentReader::open(&args.input, &options)?;
    let header = reader.read_header()?;

    let filter = Filter {
        require: args.require_flags.unwrap_or_else(Flags::empty),
        exclude: args.exclude_flags.unwrap_or(DEFAULT_EXCLUDE_FLAGS),
        exclude_all: args.exclude_all_flags,
        ..Filter::default()
    };

    let copied_tags = copied_tags(&args)?;

    if (args.read1.is_some() || args.read2.is_some()) && is_coordinate_sorted(&header) {
        eprintln!("samtools fastq: Coordinate sorted file.  Read pairs may be out of order");
    }

    // As in samtools, /1 and /2 are appended by default only when writing a single stream.
    let has_output_files = args.output.is_some()
        || args.other.is_some()
        || args.read1.is_some()
        || args.read2.is_some()
        || args.singletons.is_some();
    let append_suffix = args.suffix_always || (!args.suffix_never && !has_output_files);

    let level = crate::output::parse_level(args.level)?;
    let mut sinks = vec![Sink::create(args.output.as_deref(), level, args.threads)?];
    let mut add = |path: &Option<PathBuf>| -> anyhow::Result<Option<usize>> {
        match path {
            Some(path) => {
                sinks.push(Sink::create(Some(path), level, args.threads)?);
                Ok(Some(sinks.len() - 1))
            }
            None => Ok(None),
        }
    };

    let routes = Routes {
        other: add(&args.other)?.unwrap_or(0),
        first: add(&args.read1)?.unwrap_or(0),
        last: add(&args.read2)?.unwrap_or(0),
        singleton: add(&args.singletons)?,
    };

    let mut writer = ReadWriter {
        format,
        append_suffix,
        sinks,
        routes,
        processed: 0,
        singletons: 0,
    };

    let mut pending = Read::default();
    let mut has_pending = false;
    let mut current = Read::default();

    reader.for_each(&header, |record| {
        if !filter.matches(&record, &header)? {
            return Ok(());
        }

        writer.processed += 1;
        current.fill(&record, &copied_tags, args.default_quality)?;

        if current.segment() == Segment::Other {
            writer.write(writer.routes.other, &current)?;
            return Ok(());
        }

        if has_pending && pending.name == current.name && pending.segment() != current.segment() {
            let (first, last) = if pending.segment() == Segment::First {
                (&pending, &current)
            } else {
                (&current, &pending)
            };

            writer.write(writer.routes.first, first)?;
            writer.write(writer.routes.last, last)?;
            has_pending = false;
        } else {
            if has_pending {
                writer.write_singleton(&pending)?;
            }

            std::mem::swap(&mut pending, &mut current);
            has_pending = true;
        }

        Ok(())
    })?;

    if has_pending {
        writer.write_singleton(&pending)?;
    }

    writer.finish()
}

struct ReadWriter {
    format: ReadFormat,
    append_suffix: bool,
    sinks: Vec<Sink>,
    routes: Routes,
    processed: u64,
    singletons: u64,
}

impl ReadWriter {
    fn write_singleton(&mut self, read: &Read) -> io::Result<()> {
        let sink = match self.routes.singleton {
            Some(i) => {
                self.singletons += 1;
                i
            }
            None if read.segment() == Segment::First => self.routes.first,
            None => self.routes.last,
        };

        self.write(sink, read)
    }

    fn write(&mut self, sink: usize, read: &Read) -> io::Result<()> {
        let w = self.sinks[sink].writer();

        let marker = match self.format {
            ReadFormat::Fastq => b'@',
            ReadFormat::Fasta => b'>',
        };

        w.write_all(&[marker])?;
        w.write_all(&read.name)?;

        if self.append_suffix {
            match read.segment() {
                Segment::First => w.write_all(b"/1")?,
                Segment::Last => w.write_all(b"/2")?,
                Segment::Other => {}
            }
        }

        w.write_all(&read.tags)?;
        w.write_all(b"\n")?;
        w.write_all(&read.sequence)?;
        w.write_all(b"\n")?;

        if self.format == ReadFormat::Fastq {
            w.write_all(b"+\n")?;

            // Encode in chunks to keep writes few.
            let mut buf = [0; 256];
            for chunk in read.quality_scores.chunks(buf.len()) {
                for (d, &q) in buf.iter_mut().zip(chunk) {
                    *d = q.saturating_add(b'!');
                }
                w.write_all(&buf[..chunk.len()])?;
            }

            w.write_all(b"\n")?;
        }

        Ok(())
    }

    fn finish(self) -> anyhow::Result<()> {
        for sink in self.sinks {
            sink.finish()?;
        }

        eprintln!(
            "[M::bam2fq_mainloop] discarded {} singletons",
            self.singletons
        );
        eprintln!("[M::bam2fq_mainloop] processed {} reads", self.processed);

        Ok(())
    }
}

fn is_coordinate_sorted(header: &noodles_sam::Header) -> bool {
    use noodles_sam::header::record::value::map::header::{
        sort_order::COORDINATE, tag::SORT_ORDER,
    };

    header
        .header()
        .and_then(|hd| hd.other_fields().get(&SORT_ORDER))
        .is_some_and(|so| so == COORDINATE)
}

fn copied_tags(args: &Args) -> anyhow::Result<Vec<[u8; 2]>> {
    let mut tags = Vec::new();

    if args.copy_default_tags {
        tags.extend(DEFAULT_COPIED_TAGS);
    }

    for raw in args.copy_tags.iter().filter(|s| !s.is_empty()) {
        let tag: [u8; 2] = raw
            .as_bytes()
            .try_into()
            .with_context(|| format!("invalid tag: {raw}"))?;

        if !tags.contains(&tag) {
            tags.push(tag);
        }
    }

    Ok(tags)
}

/// Appends `\tTAG:TYPE:VALUE`, as samtools prints copied tags.
fn write_tag(dst: &mut Vec<u8>, tag: [u8; 2], value: &Value<'_>) -> io::Result<()> {
    let (ty, text): (u8, Vec<u8>) = match value {
        Value::Character(c) => (b'A', vec![*c]),
        Value::String(s) => (b'Z', s.to_vec()),
        Value::Hex(s) => (b'H', s.to_vec()),
        Value::Int8(n) => (b'i', n.to_string().into_bytes()),
        Value::UInt8(n) => (b'i', n.to_string().into_bytes()),
        Value::Int16(n) => (b'i', n.to_string().into_bytes()),
        Value::UInt16(n) => (b'i', n.to_string().into_bytes()),
        Value::Int32(n) => (b'i', n.to_string().into_bytes()),
        Value::UInt32(n) => (b'i', n.to_string().into_bytes()),
        Value::Float(n) => (b'f', format!("{n}").into_bytes()),
        // Arrays are not copied by samtools either.
        Value::Array(_) => return Ok(()),
    };

    dst.push(b'\t');
    dst.extend_from_slice(&tag);
    dst.push(b':');
    dst.push(ty);
    dst.push(b':');
    dst.extend_from_slice(&text);

    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_write_tag() -> io::Result<()> {
        let mut dst = Vec::new();
        write_tag(&mut dst, *b"RG", &Value::String(b"rg0".into()))?;
        write_tag(&mut dst, *b"NH", &Value::UInt8(2))?;
        assert_eq!(dst, b"\tRG:Z:rg0\tNH:i:2");
        Ok(())
    }
}
