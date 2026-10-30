//! `ntools fastq` and `ntools fasta`: convert alignments to reads (samtools fastq/fasta).
//!
//! This follows `bam2fq_mainloop` in samtools 1.24 and HTSlib's `fastq_format1`. Records are
//! grouped by name; of each group, the first READ1, READ2, and other read is kept (a later copy
//! with base qualities replaces one without), and the group is written when the next name
//! starts: READ1 and READ2 as a pair or as singletons, then the other read.

use std::{
    collections::HashSet,
    fs::{self, File},
    io::{self, Write},
    path::{Path, PathBuf},
};

use anyhow::{Context, bail};
use noodles_bgzf::io::writer::CompressionLevel;
use noodles_pileup::baq::SEQ_NT16_TABLE;
use noodles_sam::alignment::record::Flags;

use crate::{
    alignment::{AlignmentReader, RawEncoder, ReadOptions},
    bam_raw::{self, BASES, Record, field, le_u16, le_u32},
    bgzf_stream::BgzfWriter,
    commands::mpileup::push_kputd,
    flags,
};

/// Bytes buffered per output before they are written.
const BUFFER_CAPACITY: usize = 1 << 17;

/// samtools' default exclusion: secondary and supplementary alignments.
const DEFAULT_EXCLUDE_FLAGS: Flags = Flags::SECONDARY.union(Flags::SUPPLEMENTARY);

/// Tags copied by `-t`.
const DEFAULT_COPIED_TAGS: &str = "RG,BC,QT";

/// The most UMI tags HTSlib searches (`UMI_TAGS`).
const MAX_UMI_TAGS: usize = 5;

/// The largest `-v` quality.
const MAX_DEFAULT_QUALITY: u8 = 93;

/// Missing base qualities, as BAM stores them.
const MISSING_QUALITY: u8 = 0xff;

/// The quality character HTSlib writes for bases without qualities.
const MISSING_QUALITY_CHAR: u8 = b'B';

const SOFT_CLIP: u32 = 4;
const REVERSE_COMPLEMENTED: u16 = 0x10;
const UNMAPPED: u16 = 0x04;

/// The complement of each 4-bit base code, as HTSlib writes reverse-strand reads (`=` becomes
/// `!`).
const COMPLEMENTS: &[u8; 16] = b"!TGKCYSBAWRDMHVN";

/// CIGAR operation characters by code (`bam_cigar_opchr`).
const CIGAR_OPS: &[u8; 16] = b"MIDNSHP=XB??????";

/// Decodes a byte of two 4-bit bases to the reverse complement of the pair.
static REVERSE_COMPLEMENT_PAIRS: [[u8; 2]; 256] = {
    let mut table = [[0; 2]; 256];
    let mut i = 0;
    while i < 256 {
        table[i] = [COMPLEMENTS[i & 0x0f], COMPLEMENTS[i >> 4]];
        i += 1;
    }
    table
};

/// Which read format to write.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub enum ReadFormat {
    Fastq,
    Fasta,
}

impl ReadFormat {
    fn command(self) -> &'static str {
        match self {
            Self::Fastq => "fastq",
            Self::Fasta => "fasta",
        }
    }
}

#[derive(Debug, clap::Args)]
pub struct Args {
    /// Input SAM, BAM, or CRAM file (`-` for stdin).
    #[arg(default_value = "-")]
    input: String,

    /// Write READ1 and READ2 reads (and others not sent elsewhere) to this file.
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

    /// Write singletons (READ1 or READ2 reads without their mate) to this file.
    #[arg(short = 's')]
    singletons: Option<PathBuf>,

    /// Write first index reads (from the barcode tag) to this file.
    #[arg(long = "i1", visible_alias = "I1")]
    index1: Option<PathBuf>,

    /// Write second index reads to this file.
    #[arg(long = "i2", visible_alias = "I2")]
    index2: Option<PathBuf>,

    /// How to split the barcode and quality tags into index reads (e.g. `i*i*`, `n8i14`).
    #[arg(long = "index-format", visible_aliases = ["if", "IF"])]
    index_format: Option<String>,

    /// The barcode tag.
    #[arg(long = "barcode-tag", default_value = "BC")]
    barcode_tag: String,

    /// The barcode quality tag.
    #[arg(long = "quality-tag", default_value = "QT")]
    quality_tag: String,

    /// Add an Illumina CASAVA 1.8 entry (e.g. `1:N:0:ATCACG`) to the name line.
    #[arg(short = 'i')]
    casava: bool,

    /// Add the UMI to the read name.
    #[arg(short = 'U', long = "UMI")]
    umi: bool,

    /// Tags to search for the UMI, in order.
    #[arg(long = "UMI-tag", default_value = "OX,RX")]
    umi_tags: String,

    /// Never append /1 and /2 to read names.
    #[arg(short = 'n')]
    suffix_never: bool,

    /// Always append /1 and /2 to read names.
    #[arg(short = 'N')]
    suffix_always: bool,

    /// Copy RG, BC, and QT tags to the name line.
    #[arg(short = 't')]
    copy_default_tags: bool,

    /// Copy these tags to the name line (comma-separated; `*` for all).
    #[arg(short = 'T')]
    copy_tags: Option<String>,

    /// Write qualities from the OQ tag when present.
    #[arg(short = 'O')]
    original_qualities: bool,

    /// Only include reads with this tag, and one of these values if given (`TAG[:VALUE]`).
    #[arg(short = 'd', long = "tag")]
    tag_values: Vec<String>,

    /// Only include reads with this tag and a value listed in this file (`TAG:FILE`).
    #[arg(short = 'D', long = "tag-file")]
    tag_files: Vec<String>,

    /// Only include reads with all of these flags set.
    #[arg(short = 'f', long = "require-flags", value_parser = flags::parse)]
    require_flags: Vec<Flags>,

    /// Exclude reads with any of these flags set [default: 0x900].
    #[arg(
        short = 'F',
        long = "excl-flags",
        visible_alias = "exclude-flags",
        value_parser = flags::parse
    )]
    exclude_flags: Option<Flags>,

    /// Only include reads with any of these flags set.
    #[arg(
        long = "rf",
        visible_aliases = ["incl-flags", "include-flags"],
        value_parser = flags::parse
    )]
    include_flags: Vec<Flags>,

    /// Exclude reads with all of these flags set.
    #[arg(short = 'G', value_parser = flags::parse)]
    exclude_all_flags: Vec<Flags>,

    /// Remove soft-clipped bases.
    #[arg(long = "no-sc")]
    no_soft_clips: bool,

    /// Do not keep removed soft clips in a tag.
    #[arg(long = "no-sc-bkp")]
    no_soft_clip_backup: bool,

    /// The tag keeping removed soft clips (CIGAR, bases, and qualities).
    #[arg(long = "sc-aux", default_value = "s0")]
    soft_clip_tag: String,

    /// Quality of index read bases without qualities.
    #[arg(short = 'v', default_value_t = 1, allow_negative_numbers = true)]
    default_quality: i32,

    /// Compression level for `.gz` outputs [0..9].
    #[arg(short = 'c', default_value_t = 1)]
    level: i32,

    /// Reference FASTA (indexed) for CRAM input.
    #[arg(long = "reference")]
    reference: Option<PathBuf>,

    /// Additional threads for decompression and compression.
    #[arg(short = '@', long = "threads", default_value_t = 0)]
    threads: usize,
}

pub fn run(args: Args, format: ReadFormat) -> anyhow::Result<()> {
    let config = Config::new(&args, format)?;

    let options = ReadOptions {
        threads: args.threads,
        reference: args.reference.clone(),
    };

    let (mut reader, _) = AlignmentReader::open(&args.input, &options)?;
    let header = reader.read_header()?;

    if is_coordinate_sorted(&header) {
        eprintln!(
            "samtools {}: Coordinate sorted file.  Read pairs may be out of order",
            format.command()
        );
    }

    let outputs = Outputs::open(&args)?;
    let mut converter = Converter {
        config,
        outputs,
        processed: 0,
        singletons: 0,
        index_record: Vec::new(),
    };

    let mut encoder = RawEncoder::new();
    let mut slots: [Vec<u8>; 4] = Default::default();
    let mut group = Group::default();
    let mut scratch = Vec::new();
    let mut cigar = Vec::new();

    reader.for_each(&header, |record| {
        let config = &converter.config;

        if !config.flags.keeps(record.flags()?.bits()) {
            return Ok(());
        }

        let slot = &mut slots[group.free];
        slot.clear();
        slot.extend_from_slice(record.encode_raw(&header, &mut encoder)?);

        if let Some(filter) = &config.tag_filter
            && !filter.keeps(slot)?
        {
            return Ok(());
        }

        if config.original_qualities {
            apply_original_qualities(slot)?;
        }

        if let Some(clips) = config.soft_clips {
            if remove_soft_clips(slot, clips, &mut cigar, &mut scratch)? {
                std::mem::swap(slot, &mut scratch);
            }

            if Record::new(slot)?.sequence_len() == 0 {
                return Ok(());
            }
        }

        converter.processed += 1;
        group.add(&mut slots, &mut converter)
    })?;

    if group.name_slot.is_some() {
        converter.flush(&mut slots, &group)?;
    }

    converter.finish()
}

/// The reads kept for the current name: indices into the record slots by part.
#[derive(Debug, Default)]
struct Group {
    /// The kept record of each part (other, READ1, READ2), if any.
    best: [Option<usize>; 3],
    /// Whether each kept record has base qualities.
    has_qualities: [bool; 3],
    /// A kept record, holding the group's name.
    name_slot: Option<usize>,
    /// The slot the next record is read into.
    free: usize,
}

impl Group {
    /// Adds the record in the free slot, first writing the previous group if its name differs.
    fn add(&mut self, slots: &mut [Vec<u8>; 4], converter: &mut Converter) -> anyhow::Result<()> {
        let n = self.free;
        let record = Record::new(&slots[n])?;

        let is_same_name = self
            .name_slot
            .is_some_and(|i| Record::new(&slots[i]).is_ok_and(|r| r.name() == record.name()));

        let part = Part::of(record.flags()) as usize;
        let has_qualities = record
            .quality_scores()
            .first()
            .is_none_or(|&q| q != MISSING_QUALITY);

        if !is_same_name {
            if self.name_slot.is_some() {
                converter.flush(slots, self)?;
            }

            self.best = [None; 3];
        }

        // Prefer a copy of the read that has base qualities.
        let replaces = match self.best[part] {
            None => true,
            Some(_) => has_qualities && !self.has_qualities[part],
        };

        if replaces {
            self.best[part] = Some(n);
            self.has_qualities[part] = has_qualities;
            self.name_slot = Some(n);
            self.free = (0..slots.len())
                .find(|i| !self.best.contains(&Some(*i)))
                .expect("4 slots hold at most 3 kept records");
        }

        Ok(())
    }
}

/// samtools' read designation.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
enum Part {
    Other = 0,
    First = 1,
    Last = 2,
}

impl Part {
    fn of(flags: u16) -> Self {
        const FIRST: u16 = 0x40;
        const LAST: u16 = 0x80;

        match (flags & FIRST != 0, flags & LAST != 0) {
            (true, false) => Self::First,
            (false, true) => Self::Last,
            _ => Self::Other,
        }
    }
}

/// Settings derived from the arguments.
struct Config {
    flags: FlagFilter,
    tag_filter: Option<TagFilter>,
    original_qualities: bool,
    soft_clips: Option<SoftClips>,
    formatter: Formatter,
    /// `--index-format`, if given.
    index_format: Option<Vec<u8>>,
    quality_tag: [u8; 2],
    default_quality: u8,
}

impl Config {
    fn new(args: &Args, format: ReadFormat) -> anyhow::Result<Self> {
        let index_count = args
            .index_format
            .as_deref()
            .map_or(0, |s| s.bytes().filter(|&b| b == b'i').count());

        if index_count > 2 {
            bail!("Invalid index format: more than 2 indexes");
        }

        if args.index2.is_some() && args.index1.is_none() {
            bail!("Index one specified, but index two not given");
        }

        if args.casava && index_count == 0 {
            bail!(
                "You must specify an index format (--index-format) with the Illumina Casava (-i) option"
            );
        }

        if index_count == 0 && args.index1.is_some() {
            bail!("index_format not specified, but index file given");
        }

        let default_quality = u8::try_from(args.default_quality)
            .ok()
            .filter(|&q| q <= MAX_DEFAULT_QUALITY)
            .with_context(|| {
                format!(
                    "Invalid -v default quality {}, allowed range 0 to {MAX_DEFAULT_QUALITY}",
                    args.default_quality
                )
            })?;

        let soft_clips = if args.no_soft_clips {
            let backup = if args.no_soft_clip_backup {
                None
            } else {
                Some(
                    <[u8; 2]>::try_from(args.soft_clip_tag.as_bytes())
                        .ok()
                        .context("Invalid backup tag")?,
                )
            };

            Some(SoftClips { backup })
        } else {
            None
        };

        // As in samtools, /1 and /2 are appended by default unless READ1 or READ2 files are given.
        let rnum = args.suffix_always
            || (!args.suffix_never
                && args.output.is_none()
                && args.read1.is_none()
                && args.read2.is_none());

        let formatter = Formatter {
            fastq: format == ReadFormat::Fastq,
            rnum,
            casava: args.casava,
            barcode_tag: short_tag(&args.barcode_tag),
            umi_tags: if args.umi {
                parse_umi_tags(&args.umi_tags)
            } else {
                Vec::new()
            },
            aux: AuxTags::new(
                args.copy_default_tags.then_some(DEFAULT_COPIED_TAGS),
                args.copy_tags.as_deref(),
            ),
        };

        let mut flags = FlagFilter {
            require: 0,
            exclude: args.exclude_flags.unwrap_or(DEFAULT_EXCLUDE_FLAGS).bits(),
            exclude_all: 0,
            include_any: 0,
        };

        for f in &args.require_flags {
            flags.require |= f.bits();
        }

        for f in &args.exclude_all_flags {
            flags.exclude_all |= f.bits();
        }

        for f in &args.include_flags {
            flags.include_any |= f.bits();
        }

        Ok(Self {
            flags,
            tag_filter: TagFilter::new(&args.tag_values, &args.tag_files)?,
            original_qualities: args.original_qualities,
            soft_clips,
            formatter,
            index_format: args.index_format.as_ref().map(|s| s.as_bytes().to_vec()),
            quality_tag: short_tag(&args.quality_tag),
            default_quality,
        })
    }
}

/// The first two bytes of a tag argument, as HTSlib looks tags up.
fn short_tag(s: &str) -> [u8; 2] {
    let b = s.as_bytes();
    [
        b.first().copied().unwrap_or(0),
        b.get(1).copied().unwrap_or(0),
    ]
}

/// `-f`, `-F`, `--rf`, and `-G`.
#[derive(Clone, Copy, Debug)]
struct FlagFilter {
    require: u16,
    exclude: u16,
    exclude_all: u16,
    include_any: u16,
}

impl FlagFilter {
    #[inline]
    fn keeps(&self, flags: u16) -> bool {
        flags & self.require == self.require
            && flags & self.exclude == 0
            && (self.include_any == 0 || flags & self.include_any != 0)
            && !(flags & self.exclude_all != 0 && flags & self.exclude_all == self.exclude_all)
    }
}

/// `-d TAG[:VALUE]` and `-D TAG:FILE`: keep records with the tag, and with one of the values if
/// any are given. This is not `view`'s [`crate::filter::TagFilter`]: samtools fastq reads whole
/// lines from `-D` files and matches floats as `%f` prints them, where samtools view reads
/// whitespace-separated words and never matches floats.
#[derive(Debug)]
struct TagFilter {
    tag: [u8; 2],
    values: Option<HashSet<Vec<u8>>>,
}

impl TagFilter {
    fn new(tag_values: &[String], tag_files: &[String]) -> anyhow::Result<Option<Self>> {
        let mut filter: Option<Self> = None;

        for arg in tag_values {
            let b = arg.as_bytes();

            if b.len() < 2 || (b.len() > 2 && b[2] != b':') {
                bail!("Invalid \"tag:value\" option: \"{arg}\"");
            }

            let filter = Self::for_tag(&mut filter, arg)?;

            if let Some(value) = b.get(3..) {
                filter
                    .values
                    .get_or_insert_with(HashSet::new)
                    .insert(value.to_vec());
            }
        }

        for arg in tag_files {
            let b = arg.as_bytes();

            // samtools also accepts `;`, for MinGW path translation.
            if b.len() < 4 || (b[2] != b':' && b[2] != b';') {
                bail!("Invalid \"tag:file\" option: \"{arg}\"");
            }

            let filter = Self::for_tag(&mut filter, arg)?;
            let path = &arg[3..];
            let text =
                fs::read(path).with_context(|| format!("failed to open \"{path}\" for reading"))?;
            let values = filter.values.get_or_insert_with(HashSet::new);

            for line in text.split_inclusive(|&c| c == b'\n') {
                let line = line.strip_suffix(b"\n").unwrap_or(line);
                let line = line.strip_suffix(b"\r").unwrap_or(line);
                values.insert(line.to_vec());
            }
        }

        Ok(filter)
    }

    /// The filter for the tag starting `arg`. As in samtools, every argument must name the same
    /// tag.
    fn for_tag<'a>(filter: &'a mut Option<Self>, arg: &str) -> anyhow::Result<&'a mut Self> {
        let tag = [arg.as_bytes()[0], arg.as_bytes()[1]];
        let filter = filter.get_or_insert(Self { tag, values: None });

        if filter.tag != tag {
            bail!("Different tag type specified to before");
        }

        Ok(filter)
    }

    fn keeps(&self, record: &[u8]) -> io::Result<bool> {
        let Some(range) = bam_raw::find_field(record, self.tag)? else {
            return Ok(false);
        };

        let Some(values) = &self.values else {
            return Ok(true);
        };

        let ty = record[range.start + 2];
        let raw = &record[range.start + 3..range.end];

        // Values compare as text: integers in decimal and floats as `%f` prints them.
        let text: Vec<u8> = match ty {
            b'c' | b'C' | b's' | b'S' | b'i' | b'I' => integer_value(ty, raw)
                .unwrap_or_default()
                .to_string()
                .into_bytes(),
            b'f' => format_f(f64::from(f32::from_le_bytes(
                raw[..4].try_into().expect("a float field has 4 bytes"),
            )))
            .into_bytes(),
            b'A' => raw.to_vec(),
            b'Z' | b'H' => raw[..raw.len() - 1].to_vec(),
            // Anything else fails the match.
            _ => return Ok(false),
        };

        Ok(values.contains(&text))
    }
}

/// C's `%f`.
fn format_f(n: f64) -> String {
    if n.is_nan() {
        if n.is_sign_negative() { "-nan" } else { "nan" }.to_owned()
    } else {
        format!("{n:.6}")
    }
}

fn integer_value(ty: u8, raw: &[u8]) -> Option<i64> {
    Some(match ty {
        b'c' => i64::from(raw.first()?.cast_signed()),
        b'C' => i64::from(*raw.first()?),
        b's' => i64::from(i16::from_le_bytes(raw.get(..2)?.try_into().ok()?)),
        b'S' => i64::from(u16::from_le_bytes(raw.get(..2)?.try_into().ok()?)),
        b'i' => i64::from(i32::from_le_bytes(raw.get(..4)?.try_into().ok()?)),
        b'I' => i64::from(u32::from_le_bytes(raw.get(..4)?.try_into().ok()?)),
        _ => return None,
    })
}

/// The byte offset of the base qualities in an encoded record.
fn quality_offset(record: &[u8]) -> io::Result<usize> {
    let r = Record::new(record)?;
    Ok(record.len() - r.data().len() - r.sequence_len())
}

/// `-O`: replaces base qualities with those of a string `OQ` tag, as far as both go.
fn apply_original_qualities(record: &mut [u8]) -> io::Result<()> {
    let Some(range) = bam_raw::find_field(record, *b"OQ")? else {
        return Ok(());
    };

    if record[range.start + 2] != b'Z' {
        return Ok(());
    }

    let start = quality_offset(record)?;
    let len = Record::new(record)?.sequence_len();
    let n = (range.len() - 4).min(len);

    for i in 0..n {
        record[start + i] = record[range.start + 3 + i].wrapping_sub(b'!');
    }

    Ok(())
}

/// `--no-sc`, and the tag keeping what is removed unless `--no-sc-bkp`.
#[derive(Clone, Copy, Debug)]
struct SoftClips {
    backup: Option<[u8; 2]>,
}

/// `removesoftclips`: writes `record` without its soft-clipped bases and qualities to `dst`,
/// returning whether there were any. The backup tag holds the CIGAR, the removed bases, and
/// their qualities (`10S90M:ACGT...:IIII...`), reverse complemented for reverse-strand reads as
/// they are written. The rewritten record has no CIGAR, which reads do not need.
fn remove_soft_clips(
    record: &[u8],
    clips: SoftClips,
    cigar: &mut Vec<u32>,
    dst: &mut Vec<u8>,
) -> io::Result<bool> {
    /// Query-consuming operations (bit 0 of `bam_cigar_type`): M, I, S, =, and X.
    const CONSUMES_QUERY: u16 = 0b1_1001_0011;

    let r = Record::new(record)?;
    r.read_cigar(cigar)?;

    if !cigar.iter().any(|&op| op & 0x0f == SOFT_CLIP) {
        return Ok(false);
    }

    let len = r.sequence_len();
    let qualities = r.quality_scores();
    let mut kept = Vec::with_capacity(len);
    let mut kept_qualities = Vec::with_capacity(len);
    let mut clipped = Vec::new();
    let mut clipped_qualities = Vec::new();
    let mut pos = 0;

    for &op in cigar.iter() {
        let kind = op & 0x0f;

        if CONSUMES_QUERY & (1 << kind) == 0 {
            continue;
        }

        let end = (pos + (op >> 4) as usize).min(len);

        let (bases, quals) = if kind == SOFT_CLIP {
            (&mut clipped, &mut clipped_qualities)
        } else {
            (&mut kept, &mut kept_qualities)
        };

        bases.extend((pos..end).map(|i| r.base(i)));
        quals.extend_from_slice(&qualities[pos..end]);

        pos = end;
    }

    let name_len = usize::from(record[field::NAME_LENGTH]);
    let new_len = u32::try_from(kept.len()).expect("fewer bases than before");

    dst.clear();
    dst.extend_from_slice(&record[..field::NAME + name_len]);
    dst[field::CIGAR_OP_COUNT].copy_from_slice(&0u16.to_le_bytes());
    dst[field::SEQUENCE_LENGTH].copy_from_slice(&new_len.to_le_bytes());

    for pair in kept.chunks(2) {
        dst.push(pair[0] << 4 | pair.get(1).copied().unwrap_or(0));
    }

    dst.extend_from_slice(&kept_qualities);

    let Some(tag) = clips.backup else {
        dst.extend_from_slice(r.data());
        return Ok(true);
    };

    // The existing backup tag is replaced.
    dst.extend_from_slice(r.data());

    if let Some(range) = bam_raw::find_field(dst, tag)? {
        dst.drain(range);
    }

    let reverse = r.flags() & REVERSE_COMPLEMENTED != 0;

    dst.extend_from_slice(&tag);
    dst.push(b'Z');

    let push_op = |dst: &mut Vec<u8>, op: u32| {
        push_int(dst, i64::from(op >> 4));
        dst.push(CIGAR_OPS[(op & 0x0f) as usize]);
    };

    if reverse {
        cigar.iter().rev().for_each(|&op| push_op(dst, op));
    } else {
        cigar.iter().for_each(|&op| push_op(dst, op));
    }

    dst.push(b':');

    if reverse {
        dst.extend(clipped.iter().rev().map(|&c| COMPLEMENTS[usize::from(c)]));
        dst.push(b':');
        dst.extend(clipped_qualities.iter().rev().map(|q| q.wrapping_add(b'!')));
    } else {
        dst.extend(clipped.iter().map(|&c| BASES[usize::from(c)]));
        dst.push(b':');
        dst.extend(clipped_qualities.iter().map(|q| q.wrapping_add(b'!')));
    }

    dst.push(0);

    Ok(true)
}

/// Which tags `-t`/`-T` copy to the name line.
#[derive(Debug)]
enum AuxTags {
    None,
    All,
    Only(Vec<[u8; 2]>),
}

impl AuxTags {
    /// As samtools sets HTSlib's `FASTQ_OPT_AUX`: `*` or an empty `-T` copies every tag;
    /// otherwise `-t`'s tags (`default_tags`) and then `-T`'s.
    fn new(default_tags: Option<&str>, copy_tags: Option<&str>) -> Self {
        if let Some(tags) = copy_tags
            && (tags.is_empty() || tags.starts_with('*'))
        {
            return Self::All;
        }

        let mut list = String::from(default_tags.unwrap_or_default());

        if let Some(tags) = copy_tags {
            if !list.is_empty() {
                list.push(',');
            }

            list.push_str(tags);
        }

        if list.is_empty() {
            Self::None
        } else if list == "1" {
            Self::All
        } else {
            Self::Only(parse_tag_list(&list))
        }
    }

    #[inline]
    fn includes(&self, tag: [u8; 2]) -> bool {
        match self {
            Self::None => false,
            Self::All => true,
            Self::Only(tags) => tags.contains(&tag),
        }
    }
}

/// HTSlib's parsing of a `FASTQ_OPT_AUX` list: two-character tags separated by commas, up to the
/// first malformed entry.
fn parse_tag_list(list: &str) -> Vec<[u8; 2]> {
    let b = list.as_bytes();
    let at = |i: usize| b.get(i).copied().unwrap_or(0);
    let mut tags = Vec::new();
    let mut i = 0;

    while i + 3 <= b.len() + 1 {
        if at(i) == b',' || at(i + 1) == b',' || !(at(i + 2) == b',' || at(i + 2) == 0) {
            let end = (i + 3).min(b.len());
            eprintln!(
                "[W::fastq_state_set] Bad tag format '{}'; skipping option",
                String::from_utf8_lossy(&b[i..end])
            );
            break;
        }

        tags.push([at(i), at(i + 1)]);
        i += 3;
    }

    tags
}

/// HTSlib's parsing of `FASTQ_OPT_UMI`: up to five tags separated by commas, up to the first
/// malformed entry.
fn parse_umi_tags(list: &str) -> Vec<[u8; 2]> {
    let list = if list == "1" { "RX" } else { list };
    let b = list.as_bytes();
    let at = |i: usize| b.get(i).copied().unwrap_or(0);
    let mut tags = Vec::new();
    let mut i = 0;

    while i < b.len() && tags.len() < MAX_UMI_TAGS {
        if !at(i).is_ascii_alphabetic()
            || !at(i + 1).is_ascii_alphanumeric()
            || (i + 2 < b.len() && b[i + 2] != b',')
        {
            eprintln!("[W::fastq_state_set] Bad UMI tag list '{list}'");
            break;
        }

        tags.push([b[i], b[i + 1]]);
        i += 3;
    }

    tags
}

/// HTSlib's `fastq_format1` with samtools' options.
#[derive(Debug)]
struct Formatter {
    fastq: bool,
    /// Append /1 and /2 to the names of paired reads.
    rnum: bool,
    casava: bool,
    barcode_tag: [u8; 2],
    umi_tags: Vec<[u8; 2]>,
    aux: AuxTags,
}

impl Formatter {
    fn write(&self, record: &[u8], dst: &mut Vec<u8>) -> io::Result<()> {
        const PAIRED: u16 = 0x01;
        const FIRST: u16 = 0x40;
        const LAST: u16 = 0x80;
        const QC_FAIL: u16 = 0x200;

        let r = Record::new(record)?;
        let flags = r.flags();

        dst.push(if self.fastq { b'@' } else { b'>' });

        if self.umi_tags.is_empty() {
            dst.extend_from_slice(r.name());
        } else {
            self.write_name_with_umi(&r, dst)?;
        }

        if self.rnum && flags & PAIRED != 0 {
            match flags & (FIRST | LAST) {
                FIRST => dst.extend_from_slice(b"/1"),
                LAST => dst.extend_from_slice(b"/2"),
                _ => {}
            }
        }

        if self.casava {
            let rnum = if flags & FIRST != 0 {
                b'1'
            } else if flags & LAST != 0 {
                b'2'
            } else {
                b'0'
            };

            let filtered = if flags & QC_FAIL != 0 { b'Y' } else { b'N' };
            dst.extend_from_slice(&[b' ', rnum, b':', filtered, b':', b'0', b':']);

            match string_field(&r, self.barcode_tag)? {
                // A barcode must start with a base; non-letters become `+` (`ACGT-TTGA` is
                // written `ACGT+TTGA`).
                Some((b'Z', bc)) if bc.first().is_some_and(u8::is_ascii_alphabetic) => {
                    dst.extend(bc.iter().map(|&c| {
                        if c.is_ascii_alphabetic() {
                            c.to_ascii_uppercase()
                        } else {
                            b'+'
                        }
                    }));
                }
                Some(_) => {
                    eprintln!("[W::fastq_format1] BC tag starts with non-sequence base; using '0'");
                    dst.push(b'0');
                }
                None => dst.push(b'0'),
            }
        }

        if !matches!(self.aux, AuxTags::None) {
            // HTSlib moves a long CIGAR's `CG` field into the CIGAR when it reads a record, so
            // that field is not a tag.
            let data = r.data();
            let data_start = r.as_bytes().len() - data.len();
            let long_cigar = r.long_cigar_field()?.map(|range| range.start - data_start);
            self.write_tags(data, long_cigar, dst)?;
        }

        dst.push(b'\n');

        let len = r.sequence_len();
        let packed = r.packed_sequence();
        let reverse = flags & REVERSE_COMPLEMENTED != 0;

        if reverse {
            // Writing into a presized slice avoids a capacity check per base pair.
            let start = dst.len();
            dst.resize(start + len, 0);

            // An odd last base is written first.
            let (head, pairs) = dst[start..].split_at_mut(len % 2);

            if let [first] = head {
                *first = COMPLEMENTS[usize::from(packed[len / 2] >> 4)];
            }

            for (pair, &byte) in pairs
                .as_chunks_mut::<2>()
                .0
                .iter_mut()
                .zip(packed[..len / 2].iter().rev())
            {
                *pair = REVERSE_COMPLEMENT_PAIRS[usize::from(byte)];
            }
        } else {
            r.decode_sequence(dst);
        }

        if self.fastq {
            dst.extend_from_slice(b"\n+\n");

            let qualities = r.quality_scores();

            if qualities.first() == Some(&MISSING_QUALITY) {
                dst.resize(dst.len() + len, MISSING_QUALITY_CHAR);
            } else if reverse {
                dst.extend(qualities.iter().rev().map(|q| q.wrapping_add(b'!')));
            } else {
                dst.extend(qualities.iter().map(|q| q.wrapping_add(b'!')));
            }
        }

        dst.push(b'\n');

        Ok(())
    }

    /// The name with `:UMI` inserted before any `#` suffix (the last `:` or `#` decides), the
    /// UMI uppercased and with non-letters as `+`.
    fn write_name_with_umi(&self, r: &Record<'_>, dst: &mut Vec<u8>) -> io::Result<()> {
        /// HTSlib's limit on a moved `#` suffix.
        const MAX_SUFFIX_LEN: usize = 255;

        let name = r.name();
        let suffix_start = name
            .iter()
            .rposition(|&c| c == b':' || c == b'#')
            .filter(|&i| name[i] == b'#' && name.len() - i < MAX_SUFFIX_LEN)
            .unwrap_or(name.len());

        dst.extend_from_slice(&name[..suffix_start]);

        let mut umi = None;

        for &tag in &self.umi_tags {
            if let Some(field) = string_field(r, tag)? {
                umi = Some(field);
                break;
            }
        }

        if let Some((b'Z', umi)) = umi {
            dst.push(b':');
            dst.extend(umi.iter().map(|&c| {
                if c.is_ascii_alphabetic() {
                    c.to_ascii_uppercase()
                } else {
                    b'+'
                }
            }));
        }

        dst.extend_from_slice(&name[suffix_start..]);

        Ok(())
    }

    /// Appends the selected tags, in record order, as `\tTAG:TYPE:VALUE`, leaving out the field
    /// at offset `skip` in `data`.
    fn write_tags(
        &self,
        mut data: &[u8],
        skip: Option<usize>,
        dst: &mut Vec<u8>,
    ) -> io::Result<()> {
        let mut offset = 0;

        // As HTSlib, stop when a field cannot fit.
        while data.len() >= 4 {
            let tag = [data[0], data[1]];

            let len = if self.aux.includes(tag) && skip != Some(offset) {
                dst.push(b'\t');
                write_field(data, dst)?
            } else {
                match field_len(data) {
                    Some(len) => len,
                    None => break,
                }
            };

            data = &data[len..];
            offset += len;
        }

        Ok(())
    }
}

/// The type and value of a record's first field with `tag`, the value of a string (`Z` or `H`)
/// without its NUL.
fn string_field<'a>(r: &Record<'a>, tag: [u8; 2]) -> io::Result<Option<(u8, &'a [u8])>> {
    let Some(range) = r.find_field(tag)? else {
        return Ok(None);
    };

    let buf = r.as_bytes();
    let ty = buf[range.start + 2];
    let value = &buf[range.start + 3..range.end];

    Ok(Some(match ty {
        b'Z' | b'H' => (ty, &value[..value.len() - 1]),
        _ => (ty, value),
    }))
}

/// The byte length of the field at the start of `data`, including its tag and type.
fn field_len(data: &[u8]) -> Option<usize> {
    let ty = *data.get(2)?;
    let rest = data.get(3..)?;

    let value_len = match ty {
        b'A' | b'c' | b'C' => 1,
        b's' | b'S' => 2,
        b'i' | b'I' | b'f' => 4,
        b'd' => 8,
        b'Z' | b'H' => memchr::memchr(0, rest)? + 1,
        b'B' => {
            let width = array_width(*rest.first()?)?;
            let count = usize::try_from(le_u32(rest.get(1..5)?)).ok()?;
            5 + count.checked_mul(width)?
        }
        _ => return None,
    };

    (3 + value_len <= data.len()).then_some(3 + value_len)
}

fn array_width(subtype: u8) -> Option<usize> {
    match subtype {
        b'A' | b'c' | b'C' => Some(1),
        b's' | b'S' => Some(2),
        b'i' | b'I' | b'f' => Some(4),
        _ => None,
    }
}

/// HTSlib's `sam_format_aux1`: writes the field at the start of `data` as `TAG:TYPE:VALUE`,
/// returning its byte length.
fn write_field(data: &[u8], dst: &mut Vec<u8>) -> io::Result<usize> {
    let invalid = || io::Error::new(io::ErrorKind::InvalidData, "invalid auxiliary field");

    let len = field_len(data).ok_or_else(invalid)?;
    let ty = data[2];
    let value = &data[3..len];

    dst.extend_from_slice(&data[..2]);
    dst.push(b':');

    match ty {
        b'c' | b'C' | b's' | b'S' | b'i' | b'I' => {
            dst.extend_from_slice(b"i:");
            push_int(dst, integer_value(ty, value).ok_or_else(invalid)?);
        }
        b'A' => dst.extend_from_slice(&[b'A', b':', value[0]]),
        b'f' => {
            let n = f32::from_le_bytes(value.try_into().map_err(|_| invalid())?);
            dst.extend_from_slice(b"f:");
            // HTSlib's `kputd`, which keeps the sign of -0.
            noodles_sam::io::writer::write_f32(dst, n)?;
        }
        b'd' => {
            let n = f64::from_le_bytes(value.try_into().map_err(|_| invalid())?);
            dst.extend_from_slice(b"d:");
            dst.extend_from_slice(crate::depth::format_g(n, 6).as_bytes());
        }
        b'Z' | b'H' => {
            dst.extend_from_slice(&[ty, b':']);
            dst.extend_from_slice(&value[..value.len() - 1]);
        }
        b'B' => {
            let subtype = value[0];
            let width = array_width(subtype).ok_or_else(invalid)?;

            dst.extend_from_slice(&[b'B', b':', subtype]);

            for item in value[5..].chunks_exact(width) {
                dst.push(b',');

                match subtype {
                    b'f' => push_kputd(
                        dst,
                        f64::from(f32::from_le_bytes(item.try_into().map_err(|_| invalid())?)),
                    ),
                    // HTSlib fails on `A` arrays.
                    b'A' => return Err(invalid()),
                    _ => push_int(dst, integer_value(subtype, item).ok_or_else(invalid)?),
                }
            }
        }
        _ => return Err(invalid()),
    }

    Ok(len)
}

fn push_int(dst: &mut Vec<u8>, n: i64) {
    let mut buf = [0; lexical_core::BUFFER_SIZE];
    dst.extend_from_slice(lexical_core::write(n, &mut buf));
}

/// A read output: plain or BGZF-compressed, by file extension.
struct Sink {
    buf: Vec<u8>,
    out: Out,
}

enum Out {
    Plain(crate::output::Sink),
    Bgzf(BgzfWriter<File>),
}

impl Sink {
    fn create(
        path: Option<&Path>,
        level: CompressionLevel,
        threads: usize,
    ) -> anyhow::Result<Self> {
        let out = match path.filter(|p| *p != Path::new("-")) {
            None => Out::Plain(crate::output::raw_stdout()?),
            Some(path) => {
                let file =
                    File::create(path).with_context(|| format!("creating {}", path.display()))?;

                let is_compressed = path
                    .extension()
                    .is_some_and(|ext| ext == "gz" || ext == "bgz" || ext == "bgzf");

                if is_compressed {
                    Out::Bgzf(BgzfWriter::new(file, level, threads))
                } else {
                    Out::Plain(Box::new(file))
                }
            }
        };

        Ok(Self {
            buf: Vec::with_capacity(BUFFER_CAPACITY + (1 << 12)),
            out,
        })
    }

    fn write_buffered(&mut self) -> io::Result<()> {
        match &mut self.out {
            Out::Plain(w) => w.write_all(&self.buf)?,
            Out::Bgzf(w) => w.write_all(&self.buf)?,
        }

        self.buf.clear();
        Ok(())
    }

    #[inline]
    fn flush_if_full(&mut self) -> io::Result<()> {
        if self.buf.len() >= BUFFER_CAPACITY {
            self.write_buffered()?;
        }

        Ok(())
    }

    fn finish(mut self) -> io::Result<()> {
        self.write_buffered()?;

        match self.out {
            Out::Plain(mut w) => w.flush(),
            Out::Bgzf(w) => w.finish().map(|_| ()),
        }
    }
}

/// The open outputs, and which one each kind of read goes to. Outputs named by the same path
/// are shared.
struct Outputs {
    sinks: Vec<Sink>,
    paths: Vec<Option<PathBuf>>,
    /// Other, READ1, and READ2 reads.
    reads: [usize; 3],
    singletons: Option<usize>,
    index: [Option<usize>; 2],
}

impl Outputs {
    fn open(args: &Args) -> anyhow::Result<Self> {
        let level = crate::output::parse_level(u8::try_from(args.level.clamp(0, 9))?)?;

        let mut outputs = Self {
            sinks: Vec::new(),
            paths: Vec::new(),
            reads: [0; 3],
            singletons: None,
            index: [None; 2],
        };

        if let Some(path) = &args.singletons {
            outputs.singletons = Some(outputs.open_one(Some(path), level, args.threads)?);
        }

        let read_paths = [
            args.other.as_deref(),
            args.read1.as_deref().or(args.output.as_deref()),
            args.read2.as_deref().or(args.output.as_deref()),
        ];

        for (i, path) in read_paths.into_iter().enumerate() {
            outputs.reads[i] = outputs.open_one(path, level, args.threads)?;
        }

        for (i, path) in [&args.index1, &args.index2].into_iter().enumerate() {
            if let Some(path) = path {
                outputs.index[i] = Some(outputs.open_one(Some(path), level, args.threads)?);
            }
        }

        Ok(outputs)
    }

    /// Opens an output, or returns the one already open for the path (`None` is stdout).
    fn open_one(
        &mut self,
        path: Option<&Path>,
        level: CompressionLevel,
        threads: usize,
    ) -> anyhow::Result<usize> {
        let path = path.filter(|p| *p != Path::new("-")).map(Path::to_path_buf);

        if let Some(i) = self.paths.iter().position(|p| *p == path) {
            return Ok(i);
        }

        self.sinks
            .push(Sink::create(path.as_deref(), level, threads)?);
        self.paths.push(path);

        Ok(self.sinks.len() - 1)
    }
}

/// Writes groups of reads.
struct Converter {
    config: Config,
    outputs: Outputs,
    processed: u64,
    singletons: u64,
    /// A reusable buffer for index reads.
    index_record: Vec<u8>,
}

impl Converter {
    fn write(&mut self, sink: usize, record: &[u8]) -> io::Result<()> {
        let sink = &mut self.outputs.sinks[sink];
        self.config.formatter.write(record, &mut sink.buf)?;
        sink.flush_if_full()
    }

    /// `flush_rec`: writes a group's reads.
    fn flush(&mut self, slots: &mut [Vec<u8>; 4], group: &Group) -> anyhow::Result<()> {
        let [other, first, last] = group.best;

        match (first, last) {
            (Some(i), Some(j)) => {
                // With CASAVA entries, READ2 gets READ1's barcode and barcode qualities.
                if self.config.formatter.casava {
                    let barcode_tag = self.config.formatter.barcode_tag;
                    let quality_tag = self.config.quality_tag;

                    for tag in [barcode_tag, quality_tag] {
                        let value = match string_field(&Record::new(&slots[i])?, tag)? {
                            Some((b'Z' | b'H', value)) => value.to_vec(),
                            _ => continue,
                        };

                        bam_raw::update_string_field(&mut slots[j], tag, &value)?;
                    }
                }

                self.write(self.outputs.reads[1], &slots[i])?;
                self.write(self.outputs.reads[2], &slots[j])?;
                self.write_index(Some(&slots[i]), Some(&slots[j]))?;
            }
            (Some(_), None) | (None, Some(_)) => {
                let (i, part) = match first {
                    Some(i) => (i, 1),
                    None => (last.expect("one read is present"), 2),
                };

                let sink = match self.outputs.singletons {
                    Some(sink) => {
                        self.singletons += 1;
                        sink
                    }
                    None => self.outputs.reads[part],
                };

                self.write(sink, &slots[i])?;

                let (r1, r2) = if part == 1 {
                    (Some(&slots[i][..]), None)
                } else {
                    (None, Some(&slots[i][..]))
                };

                self.write_index(r1, r2)?;
            }
            (None, None) => {}
        }

        if let Some(i) = other {
            self.write(self.outputs.reads[0], &slots[i])?;
            self.write_index(Some(&slots[i]), None)?;
        }

        Ok(())
    }

    /// `output_index`: splits the barcode and its qualities (from READ1, else READ2) into index
    /// reads, as `--index-format` describes.
    fn write_index(&mut self, r1: Option<&[u8]>, r2: Option<&[u8]>) -> anyhow::Result<()> {
        // Without --index-format there are no index outputs, and the default format is valid.
        if self.config.index_format.is_none() {
            return Ok(());
        }

        let r1 = r1.map(Record::new).transpose()?;
        let r2 = r2.map(Record::new).transpose()?;
        let records = [r1, r2];

        let Some(src) = r1.or(r2) else {
            return Ok(());
        };

        let find = |tag: [u8; 2]| -> io::Result<Option<&[u8]>> {
            for r in records.iter().flatten() {
                if let Some((ty, value)) = string_field(r, tag)? {
                    return Ok(matches!(ty, b'Z' | b'H').then_some(value));
                }
            }

            Ok(None)
        };

        let Some(barcode) = find(self.config.formatter.barcode_tag)? else {
            return Ok(());
        };

        let qualities = find(self.config.quality_tag)?.filter(|q| q.len() == barcode.len());

        let format = self.config.index_format.as_deref().unwrap_or(b"i*i*");
        let mut format = format.iter().copied().peekable();
        let mut index = 0;
        let mut pos = 0;

        while index < 2 {
            let Some(code) = format.next() else {
                break;
            };

            // A length, or `*` (any character) for "up to the next non-letter".
            let mut len = 0usize;
            let mut has_len = false;

            while let Some(d) = format.next_if(u8::is_ascii_digit) {
                len = len.saturating_mul(10).saturating_add(usize::from(d - b'0'));
                has_len = true;
            }

            // Any other character stands for `*`.
            if !has_len {
                format.next();
            }

            let end = if len > 0 {
                (pos + len).min(barcode.len())
            } else {
                pos + barcode[pos..]
                    .iter()
                    .position(|c| !c.is_ascii_alphabetic())
                    .unwrap_or(barcode.len() - pos)
            };

            match code {
                b'n' => {}
                b'i' => {
                    if let Some(sink) = self.outputs.index[index]
                        && end > pos
                    {
                        build_index_record(
                            &src,
                            &barcode[pos..end],
                            qualities.map(|q| &q[pos..end]),
                            self.config.default_quality,
                            &mut self.index_record,
                        );

                        let sink = &mut self.outputs.sinks[sink];
                        self.config
                            .formatter
                            .write(&self.index_record, &mut sink.buf)?;
                        sink.flush_if_full()?;
                    }

                    index += 1;
                }
                _ => bail!("Unknown index-format code"),
            }

            // `*` skips the separator. samtools reads past the end of the tag here.
            pos = if len == 0 {
                (end + 1).min(barcode.len())
            } else {
                end
            };
        }

        Ok(())
    }

    fn finish(self) -> anyhow::Result<()> {
        for sink in self.outputs.sinks {
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

/// `write_index_rec`: an unmapped forward read of the index bases, with `src`'s name, flags,
/// and tags, and with the barcode qualities or else `default_quality`.
fn build_index_record(
    src: &Record<'_>,
    bases: &[u8],
    qualities: Option<&[u8]>,
    default_quality: u8,
    dst: &mut Vec<u8>,
) {
    const BIN_UNMAPPED: u16 = 4680;

    let buf = src.as_bytes();
    let name_len = usize::from(buf[field::NAME_LENGTH]);
    let flags = (src.flags() | UNMAPPED) & !REVERSE_COMPLEMENTED;
    let len = u32::try_from(bases.len()).expect("an index is shorter than its tag");

    dst.clear();
    dst.extend_from_slice(&(-1i32).to_le_bytes());
    dst.extend_from_slice(&(-1i32).to_le_bytes());
    dst.push(buf[field::NAME_LENGTH]);
    dst.push(0);
    dst.extend_from_slice(&BIN_UNMAPPED.to_le_bytes());
    dst.extend_from_slice(&0u16.to_le_bytes());
    dst.extend_from_slice(&flags.to_le_bytes());
    dst.extend_from_slice(&len.to_le_bytes());
    dst.extend_from_slice(&(-1i32).to_le_bytes());
    dst.extend_from_slice(&(-1i32).to_le_bytes());
    dst.extend_from_slice(&0i32.to_le_bytes());
    dst.extend_from_slice(&buf[field::NAME..field::NAME + name_len]);

    for pair in bases.chunks(2) {
        let code = |c: u8| SEQ_NT16_TABLE[usize::from(c)];
        dst.push(code(pair[0]) << 4 | pair.get(1).map_or(0, |&c| code(c)));
    }

    match qualities {
        Some(q) => dst.extend(q.iter().map(|c| c.wrapping_sub(b'!'))),
        None => dst.resize(dst.len() + bases.len(), default_quality),
    }

    dst.extend_from_slice(src.data());

    debug_assert_eq!(le_u16(&dst[field::FLAGS]), flags);
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

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_aux_tags() {
        assert!(matches!(AuxTags::new(None, None), AuxTags::None));
        assert!(matches!(
            AuxTags::new(Some(DEFAULT_COPIED_TAGS), Some("*")),
            AuxTags::All
        ));
        assert!(matches!(AuxTags::new(None, Some("")), AuxTags::All));
        assert!(matches!(AuxTags::new(None, Some("1")), AuxTags::All));

        let AuxTags::Only(tags) = AuxTags::new(Some(DEFAULT_COPIED_TAGS), Some("XA,XB")) else {
            panic!("expected a tag list");
        };
        assert_eq!(tags, [*b"RG", *b"BC", *b"QT", *b"XA", *b"XB"]);

        // Parsing stops at the first malformed entry.
        assert_eq!(parse_tag_list("XA,XBC,XD"), [*b"XA"]);
    }

    #[test]
    fn test_parse_umi_tags() {
        assert_eq!(parse_umi_tags("OX,RX"), [*b"OX", *b"RX"]);
        assert_eq!(parse_umi_tags("1"), [*b"RX"]);
        assert_eq!(parse_umi_tags("RX,1A"), [*b"RX"]);
        assert_eq!(parse_umi_tags("RXY"), Vec::<[u8; 2]>::new());
    }

    #[test]
    fn test_flag_filter() {
        let filter = FlagFilter {
            require: 0x1,
            exclude: 0x900,
            exclude_all: 0xc0,
            include_any: 0x0,
        };

        assert!(filter.keeps(0x41));
        assert!(!filter.keeps(0x40));
        assert!(!filter.keeps(0x101));
        assert!(!filter.keeps(0xc1));
    }

    #[test]
    fn test_write_field() -> io::Result<()> {
        let mut dst = Vec::new();
        write_field(b"XFf\x00\x00\xc0\x3f", &mut dst)?;
        assert_eq!(dst, b"XF:f:1.5");

        dst.clear();
        write_field(b"XBBc\x02\x00\x00\x00\x01\xfe", &mut dst)?;
        assert_eq!(dst, b"XB:B:c,1,-2");

        Ok(())
    }
}
