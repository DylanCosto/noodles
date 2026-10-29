//! `ntools markdup`: mark duplicate alignments (samtools markdup).
//!
//! This follows `bam_mark_duplicates` in samtools 1.24: reads are held in a coordinate window,
//! keyed by the unclipped 5' positions of the read (and its mate, from `MC`), and the read with
//! the highest base quality score (plus the `ms` mate score from `fixmate -m`) is kept. Hash
//! entries are removed only by the read that inserted them, as in samtools.
//!
//! With `-t` or `-d`, the duplicates of each kept read form a chain that is checked when its
//! first read leaves the window: `do` tags are pointed at the kept read, and duplicates close to
//! each other on the flowcell (from Illumina read names) are retagged as optical. With `-S`,
//! output goes through a temporary file, and a second pass marks the supplementary, secondary,
//! and unmapped records of primary duplicates.
//!
//! Regexes (`--read-coords`, `--barcode-rgx`) use Rust's `regex` syntax rather than POSIX
//! extended regular expressions; for the usual patterns of literals, classes, and groups, they
//! match the same text.

use std::{
    collections::{HashMap, VecDeque},
    fs::File,
    io::{self, BufRead, Read as _, Write},
    ops::Range,
    path::{Path, PathBuf},
};

use anyhow::{Context, bail};
use noodles_bam as bam;
use noodles_sam::{
    self as sam,
    alignment::{
        io::Write as _,
        record::{
            Flags,
            cigar::op::Kind,
            data::field::{Tag, Value},
        },
    },
    header::record::value::map::header::{sort_order::QUERY_NAME, tag::SORT_ORDER},
};
use regex::bytes::{Regex, RegexBuilder};

use crate::{
    alignment::{AlignmentReader, ReadOptions, Rec},
    bam_raw::{self, field, le_u16, le_u32},
    bgzf_stream::{BgzfReader, BgzfWriter},
    output::{self, AlignmentWriter, OutputFormat, WriteOptions},
    program,
};

/// Base qualities below this are not scored (`MD_MIN_QUALITY`).
const MIN_SCORED_QUALITY: u8 = 15;

/// Read name coordinate (or barcode) warnings reported before going quiet (`BMD_WARNING_MAX`).
const MAX_WARNINGS: u64 = 10;

// Orientations and read sides, the values samtools hashes.
const O_FF: i8 = 2;
const O_RR: i8 = 3;
const O_FR: i8 = 5;
const O_RF: i8 = 7;
const R_LE: i8 = 11;
const R_RI: i8 = 13;

const MATE_CIGAR: Tag = Tag::new(b'M', b'C');
const MATE_SCORE: Tag = Tag::new(b'm', b's');
const DUPLICATE_COUNT: Tag = Tag::new(b'd', b'c');

const UMI: [u8; 2] = *b"RX";
const ORIGINAL_NAME: [u8; 2] = *b"do";
const DUPLICATE_TYPE: [u8; 2] = *b"dt";
const OPTICAL: &[u8] = b"SQ";
const LIBRARY: &[u8] = b"LB";

#[derive(Clone, Copy, Debug, Eq, PartialEq, clap::ValueEnum)]
enum Mode {
    /// Template start and end positions.
    T,
    /// Sequence start positions.
    S,
}

#[derive(Debug, clap::Args)]
pub struct Args {
    /// Coordinate-sorted input, after `fixmate -m`.
    input: String,

    /// Output file (`-` for stdout).
    output: PathBuf,

    /// Remove duplicate reads.
    #[arg(short = 'r')]
    remove: bool,

    /// Maximum read length.
    #[arg(short = 'l', default_value_t = 300)]
    max_length: i64,

    /// Report statistics on stderr.
    #[arg(short = 's')]
    stats: bool,

    /// Write statistics to this file (implies -s).
    #[arg(short = 'f')]
    stats_file: Option<PathBuf>,

    /// Write statistics as JSON (implies -s).
    #[arg(long = "json")]
    json: bool,

    /// Also mark the supplementary, secondary, and unmapped records of duplicates.
    #[arg(short = 'S')]
    supplementary: bool,

    /// Tag duplicates with the name of the read they duplicate (do), and with dt when using -d.
    #[arg(short = 't')]
    tag: bool,

    /// Mark duplicates within this many pixels of each other (from Illumina read names) as
    /// optical (dt:Z:SQ; others dt:Z:LB). 0 disables.
    #[arg(short = 'd', default_value_t = 0)]
    optical_distance: i64,

    /// Do not check duplicates of the same read against each other (with -t or -d).
    #[arg(long = "no-multi-dup")]
    no_multi_dup: bool,

    /// Prefix for the temporary file -S writes [default: next to the output].
    #[arg(short = 'T')]
    temp_prefix: Option<PathBuf>,

    /// Clear existing duplicate flags and dt/do tags.
    #[arg(short = 'c')]
    clear: bool,

    /// Duplicate decision method for pairs.
    #[arg(short = 'm', long = "mode", value_enum, default_value_t = Mode::T)]
    mode: Mode,

    /// Output uncompressed BAM.
    #[arg(short = 'u')]
    uncompressed: bool,

    /// Include QC-failed reads.
    #[arg(long = "include-fails")]
    include_fails: bool,

    /// Use read groups in duplicate matching.
    #[arg(long = "use-read-groups")]
    use_read_groups: bool,

    /// Record each kept read's duplicate count (including itself) in a dc tag.
    #[arg(long = "duplicate-count")]
    duplicate_count: bool,

    /// Output format [default: from the output extension, else BAM].
    #[arg(short = 'O', long = "output-fmt", ignore_case = true)]
    output_format: Option<OutputFormat>,

    /// Reference FASTA (indexed) for CRAM input or output.
    #[arg(long = "reference")]
    reference: Option<PathBuf>,

    /// Additional threads for compression and decompression.
    #[arg(short = '@', long = "threads", default_value_t = 0)]
    threads: usize,

    /// Do not add an @PG line to the header.
    #[arg(long = "no-PG")]
    no_pg: bool,

    /// Regex for the coordinates in read names, with groups for the x and y coordinates and
    /// (optionally) a part, such as the lane and tile, that must match [default: Illumina names].
    #[arg(long = "read-coords")]
    read_coords: Option<String>,

    /// Order of the --read-coords groups: t (the part that must match), x, and y; txy, tyx, xyt,
    /// yxt, xty, ytx, xy, or yx.
    #[arg(long = "coords-order", default_value = "txy")]
    coords_order: String,

    /// Tag with a barcode that duplicates must match.
    #[arg(long = "barcode-tag")]
    barcode_tag: Option<String>,

    /// Use the UMI/barcode in the read name (the eighth colon-delimited part).
    #[arg(long = "barcode-name")]
    barcode_name: bool,

    /// Regex for the barcode in the read name, as its first group (instead of --barcode-name).
    #[arg(long = "barcode-rgx")]
    barcode_rgx: Option<String>,

    /// Move the UMI from the read name to an RX tag (with --barcode-name or --barcode-rgx).
    #[arg(long = "move-umi-to-tag")]
    move_umi_to_tag: bool,

    /// The separator removed with the UMI from read names.
    #[arg(long = "umi-separator", default_value = ":")]
    umi_separator: String,
}

/// `key_data_t`. Single-read keys leave the mate fields zero, so derived equality matches
/// samtools' `key_equal`.
#[derive(Clone, Copy, Debug, Default, Eq, Hash, PartialEq)]
struct Key {
    this_coord: i64,
    other_coord: i64,
    this_ref: i32,
    other_ref: i32,
    /// The barcode's hash (`do_hash`), as samtools compares barcodes.
    barcode: i32,
    read_group: i32,
    single: bool,
    leftmost: i8,
    orientation: i8,
}

/// `read_queue_t`
struct Read {
    record: bam::Record,
    /// The window coordinate: the 5' unclipped position for examined reads, else `pos`.
    pos: i64,
    tid: i32,
    /// Keys this read inserted, which it removes when written.
    pair_key: Option<Key>,
    single_key: Option<Key>,
    duplicate_count: i64,
    read_group: usize,
    /// The next read in this read's duplicate chain, by ID (`-t`/`-d`).
    duplicate: Option<u64>,
    /// The read this one was found to duplicate, by ID.
    original: Option<u64>,
    /// Whether this read's chain has been checked.
    dup_checked: bool,
}

/// `dup_map_t`: a primary duplicate whose other records `-S` marks too.
struct DuplicateName {
    /// The kept read's name, for `do` tags (with `-t`).
    original: Option<Vec<u8>>,
    optical: bool,
    read_group: usize,
}

/// `check_t`: a duplicate in a chain being checked for optical duplicates.
struct Check {
    id: u64,
    name: Vec<u8>,
    /// The part of the name (lane and tile) that must match for an optical duplicate.
    prefix: Vec<u8>,
    x: i64,
    y: i64,
    /// Whether `x`, `y`, and `prefix` were read from the name.
    has_coordinates: bool,
    optical: bool,
    paired: bool,
    qc_fail: bool,
    score: i64,
    mate_score: i64,
}

/// Flowcell coordinates parsed from a read name.
struct Coordinates {
    /// The part of the name (such as the lane and tile) that must match.
    prefix: Range<usize>,
    x: i64,
    y: i64,
}

/// How flowcell coordinates are read from read names.
enum CoordinateParser {
    /// Illumina's colon-separated fields (`get_coordinates_colons`).
    Colons,
    /// `--read-coords`, with the groups of the x and y coordinates and the part that must match.
    Regex {
        regex: Regex,
        x: usize,
        y: usize,
        t: Option<usize>,
    },
}

/// Why coordinates could not be read from a read name.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
enum CoordinateError {
    /// The `--read-coords` regex did not match (samtools fails if this is a kept read's name).
    Unmatched,
    /// The coordinates were not numbers, or the name had the wrong fields; a warning is counted.
    Unreadable,
}

/// Where barcodes that duplicates must match come from.
enum Barcodes {
    None,
    Tag([u8; 2]),
    /// The first group of a regex on the read name.
    Name(Regex),
}

/// `--move-umi-to-tag`: the barcode regex and the separator removed with the UMI.
struct UmiMover {
    regex: Regex,
    separator: u8,
}

#[derive(Debug, Default, Clone, Copy)]
struct Stats {
    reading: i64,
    writing: i64,
    excluded: i64,
    examined: i64,
    pair: i64,
    single: i64,
    duplicate: i64,
    single_duplicate: i64,
    optical: i64,
    single_optical: i64,
    non_primary_duplicate: i64,
    non_primary_optical: i64,
}

impl Stats {
    fn add(mut self, other: &Self) -> Self {
        self.reading += other.reading;
        self.writing += other.writing;
        self.excluded += other.excluded;
        self.examined += other.examined;
        self.pair += other.pair;
        self.single += other.single;
        self.duplicate += other.duplicate;
        self.single_duplicate += other.single_duplicate;
        self.optical += other.optical;
        self.single_optical += other.single_optical;
        self.non_primary_duplicate += other.non_primary_duplicate;
        self.non_primary_optical += other.non_primary_optical;
        self
    }
}

/// Which optical duplicate count a marked duplicate adds to.
#[derive(Clone, Copy, Debug)]
enum DuplicateKind {
    Pair,
    Single,
}

/// Where written reads go.
enum Output {
    Writer(AlignmentWriter),
    /// With `-S`, `block_size`-prefixed BAM records in a temporary file, for the second pass.
    Temp {
        writer: AlignmentWriter,
        temp: BgzfWriter<File>,
    },
}

impl Output {
    fn write(&mut self, header: &sam::Header, record: &[u8]) -> io::Result<()> {
        match self {
            Self::Writer(writer) => writer.write_raw(header, record),
            Self::Temp { temp, .. } => {
                let block_size = u32::try_from(record.len())
                    .map_err(|e| io::Error::new(io::ErrorKind::InvalidInput, e))?;
                temp.write_all(&block_size.to_le_bytes())?;
                temp.write_all(record)
            }
        }
    }
}

struct Markdup {
    options: Options,
    reads: VecDeque<Read>,
    /// ID of `reads[0]`; IDs increase with input order.
    front_id: u64,
    pair_hash: HashMap<Key, u64>,
    single_hash: HashMap<Key, u64>,
    /// Primary duplicates that `-S` looks up by name.
    duplicates: HashMap<Vec<u8>, DuplicateName>,
    stats: Vec<Stats>,
    /// Failed attempts to read coordinates from read names.
    warnings: u64,
    coordinates: CoordinateParser,
    barcodes: Barcodes,
    /// Failed attempts to read barcodes.
    barcode_warnings: u64,
    umi: Option<UmiMover>,
}

#[derive(Clone, Copy, Debug)]
struct Options {
    remove: bool,
    max_length: i64,
    mode: Mode,
    duplicate_count: bool,
    tag: bool,
    optical_distance: i64,
    supplementary: bool,
    /// Whether duplicate chains are built and checked (`-t` or `-d`, without `--no-multi-dup`).
    check_chain: bool,
}

pub fn run(args: Args) -> anyhow::Result<()> {
    let barcode_name = args.barcode_name || args.barcode_rgx.is_some();

    if args.barcode_tag.is_some() && barcode_name {
        bail!(
            "[markdup] error, cannot specify --barcode-tag and --barcode-name (or --barcode-rgx) at same time."
        );
    }

    if args.move_umi_to_tag && !barcode_name {
        bail!("[markdup] error, --move-umi-to-tag requires --barcode-name or --barcode-rgx.");
    }

    let coordinates = match &args.read_coords {
        Some(pattern) => coordinate_regex(pattern, &args.coords_order)?,
        None => CoordinateParser::Colons,
    };

    let barcodes = if let Some(tag) = &args.barcode_tag {
        let tag = tag.as_bytes();
        Barcodes::Tag([
            tag.first().copied().unwrap_or(0),
            tag.get(1).copied().unwrap_or(0),
        ])
    } else if barcode_name {
        // From Illumina's UMI documentation: "The UMI sequence is located in the eighth
        // colon-delimited field of the read name (QNAME)".
        const ILLUMINA_UMI: &str = "[0-9A-Za-z]+:[0-9A-Za-z]+:[0-9A-Za-z]+:[0-9A-Za-z]+:[0-9A-Za-z]+:[0-9A-Za-z]+:[0-9A-Za-z]+:([!-?A-~]+)";
        let pattern = args.barcode_rgx.as_deref().unwrap_or(ILLUMINA_UMI);
        Barcodes::Name(
            compile_regex(pattern)
                .with_context(|| "[markdup] error, barcode regex fail".to_string())?,
        )
    } else {
        Barcodes::None
    };

    let umi = match &barcodes {
        Barcodes::Name(regex) if args.move_umi_to_tag => Some(UmiMover {
            regex: regex.clone(),
            separator: args.umi_separator.bytes().next().unwrap_or(0),
        }),
        _ => None,
    };

    let read_options = ReadOptions {
        threads: args.threads,
        reference: args.reference.clone(),
    };

    let (mut reader, _) = AlignmentReader::open(&args.input, &read_options)?;
    let mut header = reader.read_header()?;

    let is_queryname_sorted = header
        .header()
        .and_then(|hd| hd.other_fields().get(&SORT_ORDER))
        .is_some_and(|so| so == QUERY_NAME);

    if is_queryname_sorted {
        bail!("[markdup] error, queryname sorted, must be sorted by coordinate.");
    }

    let read_groups: HashMap<Vec<u8>, usize> = if args.use_read_groups {
        header
            .read_groups()
            .keys()
            .enumerate()
            .map(|(i, id)| (id.to_vec(), i + 1))
            .collect()
    } else {
        HashMap::new()
    };

    let group_count = read_groups.len();

    if !args.no_pg {
        program::add(&mut header)?;
    }

    let format = args
        .output_format
        .or_else(|| OutputFormat::from_path(&args.output))
        .unwrap_or(OutputFormat::Bam);

    let write_options = WriteOptions {
        format,
        level: output::parse_level(if args.uncompressed { 0 } else { 6 })?,
        threads: args.threads,
        reference: args.reference.clone(),
    };

    let mut writer = AlignmentWriter::create(Some(&args.output), &write_options, &header)?;
    writer.write_header(&header)?;

    let optical_distance = args.optical_distance.max(0);

    let mut out = if args.supplementary {
        let dir = temp_dir(args.temp_prefix.as_deref(), &args.output);
        let file = tempfile::tempfile_in(dir)
            .with_context(|| format!("creating a temporary file in {}", dir.display()))?;

        Output::Temp {
            writer,
            temp: BgzfWriter::new(file, output::parse_level(1)?, args.threads),
        }
    } else {
        Output::Writer(writer)
    };

    let mut state = Markdup {
        options: Options {
            remove: args.remove,
            max_length: args.max_length,
            mode: args.mode,
            duplicate_count: args.duplicate_count,
            tag: args.tag,
            optical_distance,
            supplementary: args.supplementary,
            check_chain: !args.no_multi_dup && (args.tag || optical_distance > 0),
        },
        reads: VecDeque::new(),
        front_id: 0,
        pair_hash: HashMap::new(),
        single_hash: HashMap::new(),
        duplicates: HashMap::new(),
        stats: vec![Stats::default(); group_count + 1],
        warnings: 0,
        coordinates,
        barcodes,
        barcode_warnings: 0,
        umi,
    };

    let exclude = if args.include_fails {
        Flags::SECONDARY | Flags::SUPPLEMENTARY | Flags::UNMAPPED
    } else {
        Flags::SECONDARY | Flags::SUPPLEMENTARY | Flags::UNMAPPED | Flags::QC_FAIL
    };

    let mut prev_tid: i32 = 0;
    let mut prev_coord: i64 = 0;
    let mut encoder = bam::io::Writer::from(Vec::new());
    let mut cleared = Vec::new();

    reader.for_each(&header.clone(), |rec| {
        let mut raw: Vec<u8> = match &rec {
            Rec::Bam(r) => r.as_ref().to_vec(),
            _ => {
                encoder.get_mut().clear();
                encoder.write_alignment_record(&header, rec.as_dyn())?;
                encoder.get_ref()[size_of::<u32>()..].to_vec()
            }
        };

        if args.clear && flags_of(&raw).is_duplicate() {
            let flags = flags_of(&raw) - Flags::DUPLICATE;
            set_flags(&mut raw, flags);
            bam_raw::strip_fields(&raw, &mut cleared, |tag| {
                tag == DUPLICATE_TYPE || tag == ORIGINAL_NAME
            })?;
            std::mem::swap(&mut raw, &mut cleared);
        }

        let record = bam::Record::try_from(raw)?;

        let tid = i32_field(record.as_ref(), field::REFERENCE_SEQUENCE_ID);
        let pos = i64::from(i32_field(record.as_ref(), field::POSITION));

        if tid >= 0 && (tid < prev_tid || (tid == prev_tid && pos < prev_coord)) {
            bail!("[markdup] error, input is not coordinate sorted.");
        }

        prev_coord = pos;
        prev_tid = tid;

        let read_group = if args.use_read_groups {
            read_group_of(&record, &read_groups)?
        } else {
            0
        };

        state.stats[read_group].reading += 1;

        let id = state.front_id + state.reads.len() as u64;
        state.reads.push_back(Read {
            record,
            pos,
            tid,
            pair_key: None,
            single_key: None,
            duplicate_count: 1,
            read_group,
            duplicate: None,
            original: None,
            dup_checked: false,
        });

        let flags = flags_of(state.get(id).record.as_ref());

        if flags.intersects(exclude) {
            state.stats[read_group].excluded += 1;
        } else {
            state.stats[read_group].examined += 1;
            state.examine(id, read_group)?;
        }

        state.flush(prev_tid, prev_coord, &header, &mut out)?;

        Ok(())
    })?;

    state.flush_all(&header, &mut out)?;

    let writer = match out {
        Output::Writer(writer) => writer,
        Output::Temp { mut writer, temp } => {
            let file = temp.finish()?;
            state.mark_non_primary(file, args.threads, &header, &mut writer)?;
            writer
        }
    };

    writer.finish(&header)?;

    if state.warnings > 0 {
        eprintln!(
            "[markdup] warning, number of failed attempts to get coordinates from read names = {}",
            state.warnings
        );
    }

    if state.barcode_warnings > 0 {
        eprintln!(
            "[markdup] warning, number of failed attempts to get barcodes = {}",
            state.barcode_warnings
        );
    }

    if args.stats || args.stats_file.is_some() || args.json {
        let mut out: Box<dyn Write> = match &args.stats_file {
            Some(path) => Box::new(
                File::create(path).with_context(|| format!("creating {}", path.display()))?,
            ),
            None => Box::new(io::stderr()),
        };

        let command_line = std::env::args().collect::<Vec<_>>().join(" ");
        let total = state.stats.iter().fold(Stats::default(), Stats::add);
        let group_names: Vec<String> = header
            .read_groups()
            .keys()
            .map(|id| id.to_string())
            .collect();

        if args.json {
            write_json_report(
                &mut out,
                &command_line,
                &total,
                args.use_read_groups,
                &state.stats,
                &group_names,
            )?;
        } else {
            write_stats(&mut out, "COMMAND: ", &command_line, &total)?;
            writeln!(out)?;

            if args.use_read_groups {
                if state.stats[0].reading > 0 {
                    write_stats(&mut out, "READ GROUP: ", "ungrouped", &state.stats[0])?;
                    writeln!(out)?;
                }

                for (id, stats) in group_names.iter().zip(&state.stats[1..]) {
                    write_stats(&mut out, "READ GROUP: ", id, stats)?;
                    writeln!(out)?;
                }
            }
        }
    }

    Ok(())
}

/// The directory for `-S`'s temporary file: the `-T` prefix's, else the output's, else the current
/// directory (for stdout).
fn temp_dir<'a>(prefix: Option<&'a Path>, output: &'a Path) -> &'a Path {
    let is_stdout = output == Path::new("-") || output.starts_with("/dev");

    let dir = match prefix {
        Some(prefix) if prefix.is_dir() => Some(prefix),
        Some(prefix) => prefix.parent(),
        None if is_stdout => None,
        None => output.parent(),
    };

    dir.filter(|p| !p.as_os_str().is_empty())
        .unwrap_or(Path::new("."))
}

impl Markdup {
    fn read(&self, id: u64) -> Option<&Read> {
        id.checked_sub(self.front_id)
            .and_then(|i| self.reads.get(usize::try_from(i).ok()?))
    }

    fn read_mut(&mut self, id: u64) -> Option<&mut Read> {
        id.checked_sub(self.front_id)
            .and_then(|i| self.reads.get_mut(usize::try_from(i).ok()?))
    }

    /// A read that hashes or duplicate chains refer to, which stays buffered until written.
    fn get(&self, id: u64) -> &Read {
        self.read(id)
            .expect("hashed and chained reads are buffered")
    }

    fn get_mut(&mut self, id: u64) -> &mut Read {
        self.read_mut(id)
            .expect("hashed and chained reads are buffered")
    }

    /// Edits a buffered read's encoded record.
    fn edit<F>(&mut self, id: u64, f: F) -> io::Result<()>
    where
        F: FnOnce(&mut Vec<u8>) -> io::Result<()>,
    {
        let read = self.get_mut(id);
        let mut raw = Vec::from(std::mem::take(&mut read.record));
        let result = f(&mut raw);
        read.record = bam::Record::try_from(raw)?;
        result
    }

    /// The barcode hash of a buffered read, counting a warning if it cannot be read.
    fn barcode(&mut self, id: u64) -> i32 {
        let Self {
            reads,
            front_id,
            barcodes,
            barcode_warnings,
            ..
        } = self;

        let i = usize::try_from(id - *front_id).expect("buffered read");
        barcodes.hash(reads[i].record.as_ref(), barcode_warnings)
    }

    /// The last read in the duplicate chain starting at `id`.
    fn chain_end(&self, mut id: u64) -> u64 {
        while let Some(next) = self.get(id).duplicate {
            id = next;
        }

        id
    }

    /// Keys the read and resolves duplicates against earlier reads with the same key.
    fn examine(&mut self, id: u64, read_group: usize) -> anyhow::Result<()> {
        let read_group_key = i32::try_from(read_group).unwrap_or(i32::MAX);
        let chain = self.options.check_chain;

        if has_mate(&self.get(id).record) {
            // samtools reads the barcode for each key, warning twice.
            let pair_barcode = self.barcode(id);
            let single_barcode = self.barcode(id);
            let record = &self.get(id).record;
            let pair_key = make_pair_key(record, self.options.mode, read_group_key, pair_barcode)?;
            let single_key = make_single_key(record, read_group_key, single_barcode);

            self.stats[read_group].pair += 1;
            self.get_mut(id).pos = single_key.this_coord;

            // A single read at a paired read's position is a duplicate, whatever its score.
            match self.existing(&self.single_hash, &single_key) {
                None => {
                    self.single_hash.insert(single_key, id);
                    self.get_mut(id).single_key = Some(single_key);
                }
                Some(other) if !has_mate(&self.get(other).record) => {
                    if chain {
                        self.get_mut(id).duplicate = Some(other);
                        self.get_mut(other).original = Some(id);
                    }

                    self.single_hash.insert(single_key, id);
                    self.get_mut(id).duplicate_count += 1;
                    self.mark_duplicate(id, other, read_group, DuplicateKind::Single)?;
                    self.stats[read_group].single_duplicate += 1;
                }
                Some(_) => {}
            }

            match self.existing(&self.pair_hash, &pair_key) {
                None => {
                    self.pair_hash.insert(pair_key, id);
                    self.get_mut(id).pair_key = Some(pair_key);
                }
                Some(other) => {
                    let old = &self.get(other).record;
                    let new = &self.get(id).record;

                    let (old_score, new_score) = if flags_of(old.as_ref()).is_qc_fail()
                        != flags_of(new.as_ref()).is_qc_fail()
                    {
                        if flags_of(old.as_ref()).is_qc_fail() {
                            (0, 1)
                        } else {
                            (1, 0)
                        }
                    } else {
                        (
                            read_score(old.as_ref()) + mate_score(old)?,
                            read_score(new.as_ref()) + mate_score(new)?,
                        )
                    };

                    let tie_add = if new_score == old_score {
                        if name(new.as_ref()) < name(old.as_ref()) {
                            1
                        } else {
                            -1
                        }
                    } else {
                        0
                    };

                    let (original, duplicate) = if new_score + tie_add > old_score {
                        let other_count = self.get(other).duplicate_count;
                        self.get_mut(id).duplicate_count += other_count;

                        if chain {
                            match self.get(id).duplicate {
                                Some(first) => {
                                    let last = self.chain_end(first);
                                    self.get_mut(last).duplicate = Some(other);
                                }
                                None => self.get_mut(id).duplicate = Some(other),
                            }

                            self.get_mut(other).original = Some(id);
                        }

                        self.pair_hash.insert(pair_key, id);
                        (id, other)
                    } else {
                        if chain {
                            // The new read goes first in the kept read's chain, followed by the
                            // kept read's duplicates and then the new read's own.
                            if let Some(other_first) = self.get(other).duplicate {
                                if let Some(own_first) = self.get(id).duplicate {
                                    let last = self.chain_end(other_first);
                                    self.get_mut(last).duplicate = Some(own_first);
                                }

                                self.get_mut(id).duplicate = Some(other_first);
                            }

                            self.get_mut(other).duplicate = Some(id);
                            self.get_mut(id).original = Some(other);
                        }

                        self.get_mut(other).duplicate_count += 1;
                        (other, id)
                    };

                    self.mark_duplicate(original, duplicate, read_group, DuplicateKind::Pair)?;
                    self.stats[read_group].duplicate += 1;
                }
            }
        } else {
            let barcode = self.barcode(id);
            let single_key = make_single_key(&self.get(id).record, read_group_key, barcode);

            self.stats[read_group].single += 1;
            self.get_mut(id).pos = single_key.this_coord;

            match self.existing(&self.single_hash, &single_key) {
                None => {
                    self.single_hash.insert(single_key, id);
                    self.get_mut(id).single_key = Some(single_key);
                }
                Some(other) => {
                    let other_is_paired = has_mate(&self.get(other).record);

                    let new_wins = !other_is_paired
                        && read_score(self.get(id).record.as_ref())
                            > read_score(self.get(other).record.as_ref());

                    let (original, duplicate) = if new_wins {
                        let other_count = self.get(other).duplicate_count;
                        self.get_mut(id).duplicate_count += other_count;

                        if chain {
                            self.get_mut(id).duplicate = Some(other);
                            self.get_mut(other).original = Some(id);
                        }

                        self.single_hash.insert(single_key, id);
                        (id, other)
                    } else {
                        // Against one of a pair, a single read is always the duplicate.
                        if chain {
                            if let Some(other_first) = self.get(other).duplicate {
                                self.get_mut(id).duplicate = Some(other_first);
                            }

                            self.get_mut(other).duplicate = Some(id);
                            self.get_mut(id).original = Some(other);
                        }

                        self.get_mut(other).duplicate_count += 1;
                        (other, id)
                    };

                    self.mark_duplicate(original, duplicate, read_group, DuplicateKind::Single)?;
                    self.stats[read_group].single_duplicate += 1;
                }
            }
        }

        Ok(())
    }

    /// The read a key maps to, if that read is still buffered.
    fn existing(&self, hash: &HashMap<Key, u64>, key: &Key) -> Option<u64> {
        hash.get(key).copied().filter(|&id| self.read(id).is_some())
    }

    /// `mark_duplicates`: flags `duplicate`, tags it (`-t`, `-d`), and records it for `-S`.
    fn mark_duplicate(
        &mut self,
        original: u64,
        duplicate: u64,
        read_group: usize,
        kind: DuplicateKind,
    ) -> anyhow::Result<()> {
        let options = self.options;
        let incoming_warnings = self.warnings;
        let original_name = name(self.get(original).record.as_ref()).to_vec();
        let duplicate_name = name(self.get(duplicate).record.as_ref()).to_vec();

        let is_optical = options.optical_distance > 0
            && is_optical_duplicate(
                &self.coordinates,
                &original_name,
                &duplicate_name,
                options.optical_distance,
                &mut self.warnings,
            );

        let mut tag_result = Ok(());

        self.edit(duplicate, |raw| {
            let flags = flags_of(raw) | Flags::DUPLICATE;
            set_flags(raw, flags);

            if options.tag {
                tag_result = bam_raw::update_string_field(raw, ORIGINAL_NAME, &original_name);
            }

            if options.optical_distance > 0 {
                let duplicate_type = if is_optical { OPTICAL } else { LIBRARY };
                // As in samtools, a dt field of another type is left as it is.
                let _ = bam_raw::update_string_field(raw, DUPLICATE_TYPE, duplicate_type);
            }

            Ok(())
        })?;

        tag_result.context("[markdup] error, unable to append 'do' tag.")?;

        if is_optical {
            let stats = &mut self.stats[read_group];

            match kind {
                DuplicateKind::Pair => stats.optical += 1,
                DuplicateKind::Single => stats.single_optical += 1,
            }
        }

        if self.warnings == MAX_WARNINGS && incoming_warnings != self.warnings {
            eprintln!(
                "[markdup] warning, {} decipher read name warnings.  New warnings will not be \
                 reported.",
                self.warnings
            );
        }

        if options.supplementary && has_other_records(self.get(duplicate).record.as_ref())? {
            self.duplicates
                .entry(duplicate_name)
                .or_insert_with(|| DuplicateName {
                    original: options.tag.then_some(original_name),
                    optical: is_optical,
                    read_group,
                });
        }

        Ok(())
    }

    /// `optical_retag`: changes a duplicate's type to optical.
    fn retag_optical(&mut self, id: u64, read_group: usize) -> anyhow::Result<()> {
        self.edit(id, |raw| {
            bam_raw::update_string_field(raw, DUPLICATE_TYPE, OPTICAL)
        })
        .context("[markdup] error, unable to update 'dt' tag.")?;

        let record = &self.get(id).record;
        let paired = has_mate(record);
        let marks_others = self.options.supplementary && has_other_records(record.as_ref())?;
        let read_name = name(record.as_ref()).to_vec();

        if paired {
            self.stats[read_group].optical += 1;
        } else {
            self.stats[read_group].single_optical += 1;
        }

        if marks_others {
            match self.duplicates.get_mut(&read_name) {
                Some(entry) => entry.optical = true,
                None => bail!(
                    "[markdup] error, duplicate name {} not found in hash.",
                    String::from_utf8_lossy(&read_name)
                ),
            }
        }

        Ok(())
    }

    /// `find_duplicate_chains`: checks the chain of the read that `id` duplicates (or is).
    fn check_chain(&mut self, id: u64) -> anyhow::Result<()> {
        let mut root = id;

        while let Some(original) = self.get(root).original {
            root = original;
        }

        let read_group = self.get(root).read_group;

        // A kept read whose name has no coordinates leaves the chain unchecked for optical
        // duplicates.
        if let Some(checks) = self.check_chain_against_original(root, read_group)?
            && self.options.optical_distance > 0
        {
            self.check_duplicate_chain(checks, read_group)?;
        }

        Ok(())
    }

    /// `check_chain_against_original`: points `do` tags at the kept read, retags duplicates
    /// optically close to it, and collects the chain for [`Self::check_duplicate_chain`]. Returns
    /// `None` if the kept read's name has no coordinates.
    fn check_chain_against_original(
        &mut self,
        root: u64,
        read_group: usize,
    ) -> anyhow::Result<Option<Vec<Check>>> {
        let options = self.options;
        let original_name = name(self.get(root).record.as_ref()).to_vec();

        let parsed = (options.optical_distance > 0)
            .then(|| self.coordinates.parse(&original_name, &mut self.warnings));
        let unmatched = matches!(parsed, Some(Err(CoordinateError::Unmatched)));
        let original_coordinates = parsed.and_then(Result::ok);

        let coordinates_failed = options.optical_distance > 0 && original_coordinates.is_none();

        let mut checks = Vec::new();
        let mut current = self.get(root).duplicate;

        while let Some(id) = current {
            self.get_mut(id).dup_checked = true;

            let raw = self.get(id).record.as_ref();
            let paired = has_mate(&self.get(id).record);
            let read_name = name(raw).to_vec();

            let mut check = Check {
                id,
                name: read_name.clone(),
                prefix: Vec::new(),
                x: -1,
                y: -1,
                has_coordinates: false,
                optical: false,
                paired,
                qc_fail: flags_of(raw).is_qc_fail(),
                score: 0,
                mate_score: 0,
            };

            if options.tag {
                let needs_update = match bam_raw::string_field(raw, ORIGINAL_NAME) {
                    Ok(Some(name)) => name != original_name,
                    Ok(None) => false,
                    Err(_) => bail!(
                        "[markdup] error, 'do' tag has wrong type for read {}.",
                        String::from_utf8_lossy(&read_name)
                    ),
                };

                if needs_update {
                    self.edit(id, |raw| {
                        bam_raw::update_string_field(raw, ORIGINAL_NAME, &original_name)
                    })
                    .context("[markdup] error, unable to update 'do' tag.")?;
                }
            }

            if let Some(original) = &original_coordinates {
                let raw = self.get(id).record.as_ref();
                check.optical = matches!(bam_raw::string_field(raw, DUPLICATE_TYPE), Ok(Some(t)) if t == OPTICAL);

                let is_optical = optical_duplicate_partial(
                    &self.coordinates,
                    &original_name,
                    original,
                    &read_name,
                    &mut check,
                    options.optical_distance,
                    &mut self.warnings,
                );

                if !check.optical && is_optical {
                    self.retag_optical(id, read_group)?;
                    check.optical = true;
                }

                let record = &self.get(id).record;
                check.score = read_score(record.as_ref());

                if paired {
                    check.mate_score = mate_score(record)?;
                }
            }

            current = self.get(id).duplicate;
            checks.push(check);
        }

        self.get_mut(root).dup_checked = true;

        if unmatched {
            bail!(
                "read name {} does not match --read-coords",
                String::from_utf8_lossy(&original_name)
            );
        }

        Ok((!coordinates_failed).then_some(checks))
    }

    /// `check_duplicate_chain`: checks the duplicates of one kept read against each other, and
    /// retags the worse of each optically close pair.
    fn check_duplicate_chain(
        &mut self,
        mut checks: Vec<Check>,
        read_group: usize,
    ) -> anyhow::Result<()> {
        let distance = self.options.optical_distance;

        // samtools compares duplicates whose names have no coordinates using whatever its check
        // list held before (uninitialized fields), so leaving them out is the closest match.
        checks.retain(|check| check.has_coordinates);

        // `chain_sort`, which compares x as an `int` difference. The sort is stable, as glibc's
        // `qsort` is when it can allocate.
        checks.sort_by(|a, b| {
            a.prefix
                .len()
                .cmp(&b.prefix.len())
                .then_with(|| a.prefix.cmp(&b.prefix))
                .then_with(|| ((a.x - b.x) as i32).cmp(&0))
        });

        let len = checks.len();
        let mut curr = 0;

        while curr + 1 < len {
            // As in samtools, only a same-length prefix that differs ends the group.
            let mut end = curr + 1;

            while end < len
                && !(checks[curr].prefix.len() == checks[end].prefix.len()
                    && checks[curr].prefix != checks[end].prefix)
            {
                end += 1;
            }

            while curr < end {
                let mut count = curr + 1;

                while count < end && checks[count].x - checks[curr].x <= distance {
                    let (current, chk) = (&checks[curr], &checks[count]);

                    if (current.optical && chk.optical) || (current.y - chk.y).abs() > distance {
                        count += 1;
                        continue;
                    }

                    let chk_is_duplicate = if current.paired != chk.paired {
                        // A single read is the duplicate of a pair.
                        !chk.paired
                    } else {
                        let (current_score, mut chk_score) = if current.qc_fail != chk.qc_fail {
                            if current.qc_fail { (0, 1) } else { (1, 0) }
                        } else if current.paired {
                            (
                                current.score + current.mate_score,
                                chk.score + chk.mate_score,
                            )
                        } else {
                            (current.score, chk.score)
                        };

                        if current_score == chk_score {
                            if chk.name < current.name {
                                chk_score += 1;
                            } else {
                                chk_score -= 1;
                            }
                        }

                        current_score > chk_score
                    };

                    let retag = if chk_is_duplicate { count } else { curr };

                    if !checks[retag].optical {
                        self.retag_optical(checks[retag].id, read_group)?;
                        checks[retag].optical = true;
                    }

                    count += 1;
                }

                curr += 1;
            }
        }

        Ok(())
    }

    /// Writes reads that can no longer gain duplicates: those behind the window.
    fn flush(
        &mut self,
        prev_tid: i32,
        prev_coord: i64,
        header: &sam::Header,
        out: &mut Output,
    ) -> anyhow::Result<()> {
        while let Some(front) = self.reads.front() {
            let in_window = front.pos + self.options.max_length > prev_coord
                && front.tid == prev_tid
                && (prev_tid != -1 || prev_coord != -1);

            if in_window {
                break;
            }

            self.write_front(header, out)?;
        }

        Ok(())
    }

    fn flush_all(&mut self, header: &sam::Header, out: &mut Output) -> anyhow::Result<()> {
        while !self.reads.is_empty() {
            self.write_front(header, out)?;
        }

        Ok(())
    }

    fn write_front(&mut self, header: &sam::Header, out: &mut Output) -> anyhow::Result<()> {
        let front = self.reads.front().expect("front exists");

        if self.options.check_chain
            && !front.dup_checked
            && (front.original.is_some() || front.duplicate.is_some())
        {
            self.check_chain(self.front_id)
                .context("[markdup] error, duplicate checking failed.")?;
        }

        let read = self.reads.pop_front().expect("front exists");
        self.front_id += 1;

        let is_duplicate = flags_of(read.record.as_ref()).is_duplicate();

        if !self.options.remove || !is_duplicate {
            let mut raw = Vec::from(read.record);

            if self.options.duplicate_count && !is_duplicate {
                let mut stripped = Vec::new();
                bam_raw::strip_fields(&raw, &mut stripped, |tag| tag == *b"dc")?;
                append_integer_tag(&mut stripped, DUPLICATE_COUNT, read.duplicate_count);
                raw = stripped;
            }

            // With -S, names are moved in the second pass, after other records are matched to
            // them. (samtools 1.24 moves names here only at the end of the input, a bug.)
            if let Some(umi) = &self.umi
                && !self.options.supplementary
            {
                umi.move_to_tag(&mut raw)?;
            }

            out.write(header, &raw)?;

            self.stats[read.read_group].writing += 1;
        }

        if let Some(key) = read.pair_key {
            self.pair_hash.remove(&key);
        }

        if let Some(key) = read.single_key {
            self.single_hash.remove(&key);
        }

        Ok(())
    }

    /// The `-S` second pass: marks the supplementary, secondary, and unmapped records of primary
    /// duplicates, and writes every record from the temporary file.
    fn mark_non_primary(
        &mut self,
        mut file: File,
        threads: usize,
        header: &sam::Header,
        writer: &mut AlignmentWriter,
    ) -> anyhow::Result<()> {
        use std::io::Seek as _;

        file.rewind()?;

        let mut reader = BgzfReader::new(file, threads);
        let mut raw = Vec::new();
        let mut buf = [0; size_of::<u32>()];

        while !reader.fill_buf()?.is_empty() {
            reader.read_exact(&mut buf)?;
            let block_size = usize::try_from(u32::from_le_bytes(buf))?;
            raw.resize(block_size, 0);
            reader.read_exact(&mut raw)?;

            let flags = flags_of(&raw);

            if flags.intersects(Flags::SUPPLEMENTARY | Flags::UNMAPPED | Flags::SECONDARY)
                && let Some(entry) = self.duplicates.get(name(&raw))
            {
                let (original, optical, read_group) =
                    (entry.original.clone(), entry.optical, entry.read_group);

                set_flags(&mut raw, flags | Flags::DUPLICATE);
                self.stats[read_group].non_primary_duplicate += 1;

                if self.options.tag
                    && let Some(original) = original
                {
                    bam_raw::update_string_field(&mut raw, ORIGINAL_NAME, &original)
                        .context("[markdup] error, unable to append supplementary 'do' tag.")?;
                }

                if self.options.optical_distance > 0 {
                    let duplicate_type = if optical { OPTICAL } else { LIBRARY };
                    let _ = bam_raw::update_string_field(&mut raw, DUPLICATE_TYPE, duplicate_type);

                    if optical {
                        self.stats[read_group].non_primary_optical += 1;
                    }
                }
            }

            let is_duplicate = flags_of(&raw).is_duplicate();

            if !self.options.remove || !is_duplicate {
                if self.options.duplicate_count && is_duplicate {
                    let mut stripped = Vec::new();
                    bam_raw::strip_fields(&raw, &mut stripped, |tag| tag == *b"dc")?;
                    raw = stripped;
                }

                if let Some(umi) = &self.umi {
                    umi.move_to_tag(&mut raw)?;
                }

                writer.write_raw(header, &raw)?;
            }
        }

        Ok(())
    }
}

/// Whether a duplicate may have records that `-S` marks: it has an `SA` or `XA` tag, or an
/// unmapped mate.
fn has_other_records(record: &[u8]) -> io::Result<bool> {
    Ok(flags_of(record).is_mate_unmapped()
        || bam_raw::find_field(record, *b"SA")?.is_some()
        || bam_raw::find_field(record, *b"XA")?.is_some())
}

/// Counts a coordinate warning, printing the first few.
fn warn(warnings: &mut u64, message: std::fmt::Arguments<'_>) {
    *warnings += 1;

    if *warnings <= MAX_WARNINGS {
        eprintln!("[markdup] warning, {message}");
    }
}

/// `get_coordinates_colons`: the x and y coordinates from an Illumina-style read name
/// (`machine:run:flowcell:lane:tile:x:y`, with an optional `:UMI`, or older forms with 3 or 4
/// colons).
fn colon_coordinates(name: &[u8], warnings: &mut u64) -> Option<Coordinates> {
    let mut separators = 0;
    let (mut x_start, mut y_start) = (0, 0);

    for (i, &b) in name.iter().enumerate() {
        if b == b':' {
            separators += 1;

            match separators {
                2 => x_start = i + 1,
                3 => y_start = i + 1,
                // HiSeq-style names
                4 => {
                    x_start = y_start;
                    y_start = i + 1;
                }
                // Newer Illumina names
                5 => x_start = i + 1,
                6 => y_start = i + 1,
                _ => {}
            }
        }
    }

    let display = String::from_utf8_lossy(name);

    if !matches!(separators, 3 | 4 | 6 | 7) {
        warn(
            warnings,
            format_args!("cannot decipher read name {display} for optical duplicate marking."),
        );
        return None;
    }

    let Some(x) = strtol(&name[x_start..]) else {
        warn(
            warnings,
            format_args!("cannot decipher x coordinate in {display} ."),
        );
        return None;
    };

    let Some(y) = strtol(&name[y_start..]) else {
        warn(
            warnings,
            format_args!("cannot decipher y coordinate in {display} ."),
        );
        return None;
    };

    Some(Coordinates {
        prefix: 0..x_start,
        x,
        y,
    })
}

/// `--read-coords` and `--coords-order`.
fn coordinate_regex(pattern: &str, order: &str) -> anyhow::Result<CoordinateParser> {
    let (t, x, y) = if order.starts_with("txy") || order.starts_with("tyx") {
        (Some(1), 2, 3)
    } else if order.starts_with("xyt") || order.starts_with("yxt") {
        (Some(3), 1, 2)
    } else if order.starts_with("xty") || order.starts_with("ytx") {
        (Some(2), 1, 3)
    } else if order.starts_with("xy") || order.starts_with("yx") {
        (None, 1, 2)
    } else {
        bail!("[markdup] error,  could not recognise regex coordinate order \"{order}\".");
    };

    let regex = compile_regex(pattern).context("[markdup] error, regex fail")?;
    Ok(CoordinateParser::Regex { regex, x, y, t })
}

/// A regex over the bytes of read names.
fn compile_regex(pattern: &str) -> Result<Regex, regex::Error> {
    RegexBuilder::new(pattern).unicode(false).build()
}

impl CoordinateParser {
    /// `get_coordinates`: the coordinates in a read name, counting a warning if they cannot be
    /// read. A name that a `--read-coords` regex does not match is not warned about.
    fn parse(&self, name: &[u8], warnings: &mut u64) -> Result<Coordinates, CoordinateError> {
        let Self::Regex { regex, x, y, t } = self else {
            return colon_coordinates(name, warnings).ok_or(CoordinateError::Unreadable);
        };

        let groups = || {
            let captures = regex.captures(name)?;
            let prefix = match t {
                Some(t) => captures.get(*t)?.range(),
                None => 0..0,
            };
            Some((captures.get(*x)?, captures.get(*y)?, prefix))
        };

        let (x, y, prefix) = groups().ok_or(CoordinateError::Unmatched)?;

        let display = String::from_utf8_lossy(name);
        let mut coordinate = |axis: &str, m: regex::bytes::Match<'_>| {
            // samtools copies each coordinate into a 255-byte buffer.
            if m.len() > 254 {
                warn(
                    warnings,
                    format_args!(
                        "{axis} coordinate string longer than allowed qname length in {display} ({} long).",
                        m.len()
                    ),
                );
                return None;
            }

            let n = strtol(m.as_bytes());

            if n.is_none() {
                warn(
                    warnings,
                    format_args!(
                        "cannot decipher {axis} coordinate in {display} ({}).",
                        String::from_utf8_lossy(m.as_bytes())
                    ),
                );
            }

            n
        };

        let x = coordinate("x", x).ok_or(CoordinateError::Unreadable)?;
        let y = coordinate("y", y).ok_or(CoordinateError::Unreadable)?;

        Ok(Coordinates { prefix, x, y })
    }
}

impl Barcodes {
    /// The barcode's hash for a read's keys (0 if it has none), counting a warning if it cannot
    /// be read.
    fn hash(&self, record: &[u8], warnings: &mut u64) -> i32 {
        let incoming_warnings = *warnings;

        let barcode = match self {
            Self::None => None,
            Self::Tag(tag) => match bam_raw::find_field(record, *tag) {
                Ok(Some(range)) if matches!(record[range.start + 2], b'Z' | b'H') => {
                    Some(&record[range.start + 3..range.end - 1])
                }
                Ok(None) => None,
                _ => {
                    warn(
                        warnings,
                        format_args!(
                            "{} tag wrong type. Aux tag needs to be a string type.",
                            String::from_utf8_lossy(tag).trim_end_matches('\0')
                        ),
                    );
                    None
                }
            },
            Self::Name(regex) => {
                let qname = name(record);

                match regex.captures(qname) {
                    Some(captures) => {
                        if captures.get(1).is_none() {
                            warn(
                                warnings,
                                format_args!(
                                    "barcode regex unable to match substring on {}.",
                                    String::from_utf8_lossy(qname)
                                ),
                            );
                        }

                        captures.get(1).map(|m| m.as_bytes())
                    }
                    None => {
                        warn(
                            warnings,
                            format_args!(
                                "barcode regex match error \"No match\" on {}.",
                                String::from_utf8_lossy(qname)
                            ),
                        );
                        None
                    }
                }
            }
        };

        if *warnings == MAX_WARNINGS && incoming_warnings != *warnings {
            eprintln!(
                "[markdup] warning, {} barcode read warnings.  New warnings will not be reported.",
                *warnings
            );
        }

        barcode.map_or(0, |b| do_hash(b) as i32)
    }
}

/// samtools markdup's `do_hash` (Jenkins' one-at-a-time hash).
fn do_hash(key: &[u8]) -> u32 {
    let mut hash = 0u32;

    for &b in key {
        hash = hash.wrapping_add(u32::from(b));
        hash = hash.wrapping_add(hash << 10);
        hash ^= hash >> 6;
    }

    hash = hash.wrapping_add(hash << 3);
    hash ^= hash >> 11;
    hash.wrapping_add(hash << 15)
}

impl UmiMover {
    /// `move_umi_to_tag`: moves the UMI (the regex's first group) from the read name to an RX
    /// tag, removing one separator next to it.
    fn move_to_tag(&self, record: &mut Vec<u8>) -> anyhow::Result<()> {
        let qname = name(record);

        let Some(umi) = self.regex.captures(qname).and_then(|c| c.get(1)) else {
            return Ok(());
        };

        let Range { mut start, mut end } = umi.range();
        let umi = umi.as_bytes().to_vec();

        if start > 0 && qname[start - 1] == self.separator {
            start -= 1;
        } else if end < qname.len() && qname[end] == self.separator {
            end += 1;
        }

        let mut new_name = Vec::with_capacity(qname.len() - (end - start));
        new_name.extend_from_slice(&qname[..start]);
        new_name.extend_from_slice(&qname[end..]);

        bam_raw::update_string_field(record, UMI, &umi)
            .context("[markdup] error, unable to add RX tag.")?;

        // l_read_name, then the NUL-terminated name at 32.
        let old_len = usize::from(record[8]);
        record.splice(32..32 + old_len - 1, new_name.iter().copied());
        record[8] = u8::try_from(new_name.len() + 1).expect("a shorter name fits");

        Ok(())
    }
}

/// C's `strtol(s, &end, 10)`, or `None` if it reads no digits. Out-of-range values saturate.
fn strtol(s: &[u8]) -> Option<i64> {
    let s = s.trim_ascii_start();

    let (negative, digits) = match s.first() {
        Some(b'-') => (true, &s[1..]),
        Some(b'+') => (false, &s[1..]),
        _ => (false, s),
    };

    let len = digits.iter().take_while(|b| b.is_ascii_digit()).count();

    if len == 0 {
        return None;
    }

    let n = digits[..len].iter().fold(0i64, |n, &d| {
        n.saturating_mul(10).saturating_add(i64::from(d - b'0'))
    });

    Some(if negative { n.saturating_neg() } else { n })
}

/// Whether two coordinates are within `distance` of each other on both axes, given matching
/// prefixes.
fn is_close(original: &Coordinates, x: i64, y: i64, distance: i64) -> bool {
    (original.x - x).abs() <= distance && (original.y - y).abs() <= distance
}

/// `is_optical_duplicate`
fn is_optical_duplicate(
    parser: &CoordinateParser,
    original: &[u8],
    duplicate: &[u8],
    distance: i64,
    warnings: &mut u64,
) -> bool {
    let Ok(o) = parser.parse(original, warnings) else {
        return false;
    };

    let Ok(d) = parser.parse(duplicate, warnings) else {
        return false;
    };

    original[o.prefix.clone()] == duplicate[d.prefix] && is_close(&o, d.x, d.y, distance)
}

/// `optical_duplicate_partial`: like [`is_optical_duplicate`] with the original's coordinates
/// already parsed, also recording the duplicate's coordinates in `check`.
fn optical_duplicate_partial(
    parser: &CoordinateParser,
    original_name: &[u8],
    original: &Coordinates,
    duplicate: &[u8],
    check: &mut Check,
    distance: i64,
    warnings: &mut u64,
) -> bool {
    let Ok(d) = parser.parse(duplicate, warnings) else {
        return false;
    };

    let is_optical = original_name[original.prefix.clone()] == duplicate[d.prefix.clone()]
        && is_close(original, d.x, d.y, distance);

    check.x = d.x;
    check.y = d.y;
    check.has_coordinates = true;
    check.prefix = duplicate[d.prefix].to_vec();

    is_optical
}

fn i32_field(record: &[u8], range: std::ops::Range<usize>) -> i32 {
    le_u32(&record[range]).cast_signed()
}

fn flags_of(record: &[u8]) -> Flags {
    Flags::from_bits_retain(le_u16(&record[field::FLAGS]))
}

fn set_flags(record: &mut [u8], flags: Flags) {
    record[field::FLAGS].copy_from_slice(&flags.bits().to_le_bytes());
}

fn name(record: &[u8]) -> &[u8] {
    crate::sort_order::name(record)
}

/// `calc_score`: the sum of base qualities of at least 15 (missing qualities count as 255).
fn read_score(record: &[u8]) -> i64 {
    let Ok(data) = bam_raw::data_range(record) else {
        return 0;
    };

    let base_count = usize::try_from(le_u32(&record[field::SEQUENCE_LENGTH])).unwrap_or(0);
    let quality_scores = &record[data.start - base_count..data.start];

    quality_scores
        .iter()
        .filter(|&&q| q >= MIN_SCORED_QUALITY)
        .map(|&q| i64::from(q))
        .sum()
}

/// `get_mate_score`, from the `ms` tag.
fn mate_score(record: &bam::Record) -> anyhow::Result<i64> {
    let data = record.data();

    match sam::alignment::record::Data::get(&data, &MATE_SCORE).transpose()? {
        Some(value) => value
            .as_int()
            .context("[markdup] error, ms tag is not an integer"),
        None => {
            bail!("[markdup] error, no ms score tag. Please run samtools fixmate on file first.")
        }
    }
}

fn read_group_of(record: &bam::Record, read_groups: &HashMap<Vec<u8>, usize>) -> io::Result<usize> {
    let data = record.data();

    Ok(
        match sam::alignment::record::Data::get(&data, &Tag::READ_GROUP).transpose()? {
            Some(Value::String(id)) => read_groups.get(id.as_ref() as &[u8]).copied().unwrap_or(0),
            _ => 0,
        },
    )
}

/// `has_mate`: the read's mate is present and mapped.
fn has_mate(record: &bam::Record) -> bool {
    let raw = record.as_ref();
    let flags = flags_of(raw);
    let mtid = i32_field(raw, field::MATE_REFERENCE_SEQUENCE_ID);
    let mpos = i32_field(raw, 24..28);

    flags.is_segmented() && !flags.is_mate_unmapped() && !(mtid == -1 && mpos == -1)
}

/// The CIGAR as `(kind, len)` pairs.
fn cigar_ops(record: &bam::Record) -> io::Result<Vec<(Kind, i64)>> {
    record
        .cigar()
        .iter()
        .map(|result| result.map(|op| (op.kind(), i64::try_from(op.len()).unwrap_or(i64::MAX))))
        .collect()
}

/// HTSlib's `bam_endpos`.
fn end_pos(record: &bam::Record, ops: &[(Kind, i64)]) -> i64 {
    let raw = record.as_ref();
    let pos = i64::from(i32_field(raw, field::POSITION));

    let span: i64 = if flags_of(raw).is_unmapped() {
        0
    } else {
        ops.iter()
            .filter(|(kind, _)| kind.consumes_reference())
            .map(|(_, len)| len)
            .sum()
    };

    pos + span.max(1)
}

fn is_clip(kind: Kind) -> bool {
    matches!(kind, Kind::SoftClip | Kind::HardClip)
}

/// `unclipped_start(b, 1)`
fn unclipped_start(record: &bam::Record, ops: &[(Kind, i64)]) -> i64 {
    let pos = i64::from(i32_field(record.as_ref(), field::POSITION));
    let clipped: i64 = ops
        .iter()
        .take_while(|(kind, _)| is_clip(*kind))
        .map(|(_, len)| len)
        .sum();

    pos - clipped + 1
}

/// `unclipped_end(b, 1)`
fn unclipped_end(record: &bam::Record, ops: &[(Kind, i64)]) -> i64 {
    let clipped: i64 = ops
        .iter()
        .rev()
        .take_while(|(kind, _)| is_clip(*kind))
        .map(|(_, len)| len)
        .sum();

    end_pos(record, ops) + clipped
}

/// Parses a CIGAR string as samtools' `unclipped_other_*` do: numbers via `strtol`, and an
/// operation without a number counts 1.
fn parse_cigar_string(s: &[u8]) -> Vec<(u8, i64)> {
    let mut ops = Vec::new();
    let mut i = 0;

    while i < s.len() && s[i] != b'*' {
        let start = i;

        while i < s.len() && s[i].is_ascii_digit() {
            i += 1;
        }

        let n = if i > start {
            std::str::from_utf8(&s[start..i])
                .ok()
                .and_then(|t| t.parse().ok())
                .unwrap_or(i64::MAX)
        } else {
            1
        };

        if i < s.len() {
            ops.push((s[i], n));
            i += 1;
        }
    }

    ops
}

/// `unclipped_other_start(mpos, MC, 1)`
fn unclipped_other_start(mpos: i64, cigar: &[(u8, i64)]) -> i64 {
    let mut clipped = 0;

    for &(op, n) in cigar {
        match op {
            b'S' | b'H' => clipped += n,
            _ => break,
        }
    }

    mpos - clipped + 1
}

/// `unclipped_other_end(mpos, MC, 1)`
fn unclipped_other_end(mpos: i64, cigar: &[(u8, i64)]) -> i64 {
    let mut reference_position = 0;
    let mut skip = true;

    for &(op, n) in cigar {
        match op {
            b'M' | b'D' | b'N' | b'=' | b'X' => {
                reference_position += n;
                skip = false;
            }
            b'S' | b'H' if !skip => reference_position += n,
            _ => {}
        }
    }

    mpos + reference_position
}

/// `make_single_key`
fn make_single_key(record: &bam::Record, read_group: i32, barcode: i32) -> Key {
    let ops = cigar_ops(record).unwrap_or_default();
    let raw = record.as_ref();
    let this_ref = i32_field(raw, field::REFERENCE_SEQUENCE_ID) + 1;

    let (this_coord, orientation) = if flags_of(raw).is_reverse_complemented() {
        (unclipped_end(record, &ops), O_RR)
    } else {
        (unclipped_start(record, &ops), O_FF)
    };

    Key {
        this_coord,
        this_ref,
        barcode,
        read_group,
        single: true,
        orientation,
        ..Key::default()
    }
}

/// `make_pair_key`
fn make_pair_key(
    record: &bam::Record,
    mode: Mode,
    read_group: i32,
    barcode: i32,
) -> anyhow::Result<Key> {
    let ops = cigar_ops(record)?;
    let raw = record.as_ref();
    let flags = flags_of(raw);

    let this_ref = i32_field(raw, field::REFERENCE_SEQUENCE_ID) + 1;
    let other_ref = i32_field(raw, field::MATE_REFERENCE_SEQUENCE_ID) + 1;
    let mpos = i64::from(i32_field(raw, 24..28));
    let pos = i64::from(i32_field(raw, field::POSITION));

    let mut this_coord = unclipped_start(record, &ops);
    let this_end = unclipped_end(record, &ops);

    let data = record.data();
    let mate_cigar = match sam::alignment::record::Data::get(&data, &MATE_CIGAR).transpose()? {
        Some(Value::String(s)) => parse_cigar_string(s.as_ref()),
        Some(_) => bail!(
            "[markdup] error, MC tag wrong type. Please use the MC tag provided by samtools fixmate."
        ),
        None => bail!("[markdup] error, no MC tag. Please run samtools fixmate on file first."),
    };

    let other_end = unclipped_other_end(mpos, &mate_cigar);
    let mut other_coord = unclipped_other_start(mpos, &mate_cigar);

    let is_rev = flags.is_reverse_complemented();
    let is_mrev = flags.is_mate_reverse_complemented();
    let is_read1 = flags.is_first_segment();

    let leftmost;
    let orientation;

    match mode {
        Mode::T => {
            leftmost = if this_ref != other_ref {
                this_ref < other_ref
            } else if is_rev == is_mrev {
                if is_rev {
                    this_end <= other_end
                } else {
                    this_coord <= other_coord
                }
            } else if is_rev {
                this_end <= other_coord
            } else {
                this_coord <= other_end
            };

            if leftmost {
                if is_rev == is_mrev {
                    other_coord = other_end;
                    orientation = match (is_rev, is_read1) {
                        (false, true) | (true, false) => O_FF,
                        (false, false) | (true, true) => O_RR,
                    };
                } else if !is_rev {
                    orientation = O_FR;
                    other_coord = other_end;
                } else {
                    orientation = O_RF;
                    this_coord = this_end;
                }
            } else if is_rev == is_mrev {
                this_coord = this_end;
                orientation = match (is_rev, is_read1) {
                    (false, true) | (true, false) => O_RR,
                    (false, false) | (true, true) => O_FF,
                };
            } else if !is_rev {
                orientation = O_RF;
                other_coord = other_end;
            } else {
                orientation = O_FR;
                this_coord = this_end;
            }
        }
        Mode::S => {
            let difference = if this_ref != other_ref {
                i64::from(this_ref - other_ref)
            } else if is_rev == is_mrev {
                if is_rev {
                    this_end - other_end
                } else {
                    this_coord - other_coord
                }
            } else if is_rev {
                this_end - other_coord
            } else {
                this_coord - other_end
            };

            leftmost = match difference.cmp(&0) {
                std::cmp::Ordering::Less => true,
                std::cmp::Ordering::Greater => false,
                std::cmp::Ordering::Equal => {
                    if pos == mpos {
                        is_read1
                    } else {
                        pos < mpos
                    }
                }
            };

            orientation = match (leftmost, is_rev == is_mrev, is_rev) {
                (true, true, false) => O_FF,
                (true, true, true) => O_RR,
                (true, false, false) => O_FR,
                (true, false, true) => O_RF,
                (false, true, false) => O_RR,
                (false, true, true) => O_FF,
                (false, false, false) => O_RF,
                (false, false, true) => O_FR,
            };

            this_coord = if is_rev {
                unclipped_end(record, &ops)
            } else {
                unclipped_start(record, &ops)
            };

            other_coord = if is_mrev {
                unclipped_other_end(mpos, &mate_cigar)
            } else {
                unclipped_other_start(mpos, &mate_cigar)
            };
        }
    }

    Ok(Key {
        this_coord,
        other_coord,
        this_ref,
        other_ref,
        barcode,
        read_group,
        single: false,
        leftmost: if leftmost { R_LE } else { R_RI },
        orientation,
    })
}

/// Appends an integer tag with the smallest type that holds it, as HTSlib does.
fn append_integer_tag(record: &mut Vec<u8>, tag: Tag, n: i64) {
    let tag: [u8; 2] = tag.into();
    record.extend_from_slice(&tag);

    if let Ok(n) = u8::try_from(n) {
        record.push(b'C');
        record.push(n);
    } else if let Ok(n) = i8::try_from(n) {
        record.push(b'c');
        record.extend_from_slice(&n.to_le_bytes());
    } else if let Ok(n) = u16::try_from(n) {
        record.push(b'S');
        record.extend_from_slice(&n.to_le_bytes());
    } else if let Ok(n) = i16::try_from(n) {
        record.push(b's');
        record.extend_from_slice(&n.to_le_bytes());
    } else if let Ok(n) = u32::try_from(n) {
        record.push(b'I');
        record.extend_from_slice(&n.to_le_bytes());
    } else {
        record.push(b'i');
        let n = i32::try_from(n).unwrap_or(i32::MAX);
        record.extend_from_slice(&n.to_le_bytes());
    }
}

/// The Picard library size estimate samtools reports.
fn estimate_library_size(paired_reads: i64, paired_duplicate_reads: i64, optical: i64) -> i64 {
    let coverage = |x: f64, c: f64, n: f64| c / x - 1.0 + (-n / x).exp();

    let non_optical_pairs = (paired_reads - optical) / 2;
    let unique_pairs = (paired_reads - paired_duplicate_reads) / 2;
    let duplicate_pairs = (paired_duplicate_reads - optical) / 2;

    if non_optical_pairs == 0
        || duplicate_pairs == 0
        || unique_pairs == 0
        || non_optical_pairs <= duplicate_pairs
    {
        eprintln!(
            "[markdup] warning, unable to calculate estimated library size. Read pairs \
             {non_optical_pairs} should be greater than duplicate pairs {duplicate_pairs}, which \
             should both be non zero."
        );
        return 0;
    }

    let (u, n) = (unique_pairs as f64, non_optical_pairs as f64);
    let (mut m, mut big_m) = (1.0, 100.0);

    if coverage(m * u, u, n) < 0.0 {
        eprintln!("[markdup] warning, unable to calculate estimated library size.");
        return 0;
    }

    while coverage(big_m * u, u, n) > 0.0 {
        big_m *= 10.0;
    }

    for _ in 0..40 {
        let r = (m + big_m) / 2.0;
        let value = coverage(r * u, u, n);

        if value > 0.0 {
            m = r;
        } else if value < 0.0 {
            big_m = r;
        } else {
            break;
        }
    }

    #[expect(
        clippy::cast_possible_truncation,
        reason = "samtools truncates the estimate to an unsigned long the same way"
    )]
    let size = (u * (m + big_m) / 2.0) as i64;
    size
}

/// `write_stats`. The library size is estimated (possibly warning) before the title is written.
fn write_stats<W: Write>(out: &mut W, title: &str, name: &str, stats: &Stats) -> io::Result<()> {
    let library_size = estimate_library_size(stats.pair, stats.duplicate, stats.optical);

    writeln!(out, "{title}{name}")?;

    for (label, value) in stat_fields(stats, library_size) {
        writeln!(out, "{label}: {value}")?;
    }

    Ok(())
}

/// The reported statistics, in samtools' order.
fn stat_fields(stats: &Stats, library_size: i64) -> [(&'static str, i64); 15] {
    let primary_total = stats.single_duplicate + stats.duplicate;

    [
        ("READ", stats.reading),
        ("WRITTEN", stats.writing),
        ("EXCLUDED", stats.excluded),
        ("EXAMINED", stats.examined),
        ("PAIRED", stats.pair),
        ("SINGLE", stats.single),
        ("DUPLICATE PAIR", stats.duplicate),
        ("DUPLICATE SINGLE", stats.single_duplicate),
        ("DUPLICATE PAIR OPTICAL", stats.optical),
        ("DUPLICATE SINGLE OPTICAL", stats.single_optical),
        ("DUPLICATE NON PRIMARY", stats.non_primary_duplicate),
        ("DUPLICATE NON PRIMARY OPTICAL", stats.non_primary_optical),
        ("DUPLICATE PRIMARY TOTAL", primary_total),
        (
            "DUPLICATE TOTAL",
            primary_total + stats.non_primary_duplicate,
        ),
        ("ESTIMATED_LIBRARY_SIZE", library_size),
    ]
}

/// `write_json_stats`: one block of fields at `indent`, without a trailing newline.
fn write_json_stats<W: Write>(
    out: &mut W,
    indent: &str,
    group_name: Option<&str>,
    stats: &Stats,
) -> io::Result<()> {
    let library_size = estimate_library_size(stats.pair, stats.duplicate, stats.optical);

    if let Some(name) = group_name {
        writeln!(out, "{indent}\"READ GROUP\": \"{name}\",")?;
    }

    let fields = stat_fields(stats, library_size);

    for (i, (label, value)) in fields.iter().enumerate() {
        let separator = if i + 1 < fields.len() { ",\n" } else { "" };
        write!(out, "{indent}\"{label}\": {value}{separator}")?;
    }

    Ok(())
}

/// The `--json` report. As in samtools, names are written without escaping.
fn write_json_report<W: Write>(
    out: &mut W,
    command_line: &str,
    total: &Stats,
    use_read_groups: bool,
    group_stats: &[Stats],
    group_names: &[String],
) -> io::Result<()> {
    const INDENT_1: &str = "    ";
    const INDENT_2: &str = "        ";
    const INDENT_3: &str = "            ";

    writeln!(out, "{{")?;
    writeln!(out, "{INDENT_1}\"COMMAND\": \"{command_line}\",")?;
    write_json_stats(out, INDENT_1, None, total)?;
    writeln!(out, "{}", if use_read_groups { "," } else { "" })?;

    if use_read_groups {
        writeln!(out, "{INDENT_1}\"READ GROUPS\": [")?;

        if group_stats[0].reading > 0 {
            writeln!(out, "{INDENT_2}{{")?;
            write_json_stats(out, INDENT_3, Some("ungrouped"), &group_stats[0])?;
            writeln!(out)?;
            writeln!(out, "{INDENT_2}}},")?;
        }

        for (i, (name, stats)) in group_names.iter().zip(&group_stats[1..]).enumerate() {
            writeln!(out, "{INDENT_2}{{")?;
            write_json_stats(out, INDENT_3, Some(name), stats)?;
            writeln!(out)?;
            let separator = if i + 1 < group_names.len() { "," } else { "" };
            writeln!(out, "{INDENT_2}}}{separator}")?;
        }

        writeln!(out, "{INDENT_1}]")?;
    }

    writeln!(out, "}}")
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_coordinates() {
        let mut warnings = 0;

        let c = colon_coordinates(b"M1:7:FC1:1:1101:1234:5678", &mut warnings).expect("parsed");
        assert_eq!((c.prefix, c.x, c.y), (0..16, 1234, 5678));

        let c =
            colon_coordinates(b"M1:7:FC1:1:1101:1234:5678:ACGT", &mut warnings).expect("parsed");
        assert_eq!((c.prefix, c.x, c.y), (0..16, 1234, 5678));

        // HiSeq style: machine:lane:tile:x:y
        let c = colon_coordinates(b"HS:1:1101:12:34", &mut warnings).expect("parsed");
        assert_eq!((c.prefix, c.x, c.y), (0..10, 12, 34));

        assert_eq!(warnings, 0);
        assert!(colon_coordinates(b"read1", &mut warnings).is_none());
        assert!(colon_coordinates(b"a:b:c:x:y", &mut warnings).is_none());
        assert_eq!(warnings, 2);
    }

    #[test]
    fn test_coordinate_regex() -> anyhow::Result<()> {
        let mut warnings = 0;

        let parser = coordinate_regex(r"^(\w+)_(\d+)_(\d+)$", "txy")?;
        let c = parser
            .parse(b"tile7_120_340", &mut warnings)
            .expect("parsed");
        assert_eq!((c.prefix, c.x, c.y), (0..5, 120, 340));

        let parser = coordinate_regex(r"_(\d+)_(\d+)$", "yx")?;
        let c = parser.parse(b"a_1_2", &mut warnings).expect("parsed");
        assert_eq!((c.prefix, c.x, c.y), (0..0, 1, 2));

        // A name the regex does not match is skipped silently.
        assert!(matches!(
            parser.parse(b"a", &mut warnings),
            Err(CoordinateError::Unmatched)
        ));
        assert_eq!(warnings, 0);

        assert!(coordinate_regex("(x)", "abc").is_err());
        Ok(())
    }

    #[test]
    fn test_do_hash() {
        assert_eq!(do_hash(b""), 0);
        assert_eq!(do_hash(b"a"), 0xca2e_9442);
    }

    #[test]
    fn test_is_optical_duplicate() {
        let mut warnings = 0;
        let original = b"M:7:FC:1:1101:1000:2000";

        assert!(is_optical_duplicate(
            &CoordinateParser::Colons,
            original,
            b"M:7:FC:1:1101:1100:1900",
            100,
            &mut warnings
        ));
        assert!(!is_optical_duplicate(
            &CoordinateParser::Colons,
            original,
            b"M:7:FC:1:1101:1101:2000",
            100,
            &mut warnings
        ));
        assert!(!is_optical_duplicate(
            &CoordinateParser::Colons,
            original,
            b"M:7:FC:1:1102:1000:2000",
            100,
            &mut warnings
        ));
        assert_eq!(warnings, 0);
    }

    #[test]
    fn test_strtol() {
        assert_eq!(strtol(b"123:4"), Some(123));
        assert_eq!(strtol(b" -5x"), Some(-5));
        assert_eq!(strtol(b"x1"), None);
        assert_eq!(strtol(b"99999999999999999999"), Some(i64::MAX));
    }

    #[test]
    fn test_unclipped_other_positions() {
        let cigar = parse_cigar_string(b"5S10M2D3M4S2H");
        assert_eq!(unclipped_other_start(100, &cigar), 96);
        assert_eq!(unclipped_other_end(100, &cigar), 100 + 10 + 2 + 3 + 4 + 2);
        assert!(parse_cigar_string(b"*").is_empty());
    }
}
