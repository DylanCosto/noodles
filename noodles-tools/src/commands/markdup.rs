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
    path::{Path, PathBuf},
};

mod keys;
mod names;
mod optical;
mod record;
mod report;

use anyhow::{Context, bail};
use noodles_bam as bam;
use noodles_sam::{
    self as sam,
    alignment::record::{Flags, data::field::Tag},
    header::record::value::map::header::{sort_order::QUERY_NAME, tag::SORT_ORDER},
};

use crate::{
    alignment::{AlignmentReader, RawEncoder, ReadOptions},
    bam_raw::{self, field},
    bgzf_stream::{BgzfReader, BgzfWriter},
    output::{self, AlignmentWriter, OutputFormat, WriteOptions},
    program,
};

use self::{
    keys::{make_pair_key, make_single_key},
    names::{
        Barcodes, CoordinateError, CoordinateParser, UmiMover, compile_regex, coordinate_regex,
    },
    optical::{is_optical_duplicate, optical_duplicate_partial},
    record::{
        append_integer_tag, flags_of, has_mate, has_other_records, i32_field, mate_score, name,
        read_group_of, read_score, set_flags,
    },
    report::{write_json_report, write_stats},
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
    let mut encoder = RawEncoder::new();
    let mut cleared = Vec::new();

    reader.for_each(&header.clone(), |rec| {
        let mut raw: Vec<u8> = rec.encode_raw(&header, &mut encoder)?.to_vec();

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

/// Counts a coordinate warning, printing the first few.
fn warn(warnings: &mut u64, message: std::fmt::Arguments<'_>) {
    *warnings += 1;

    if *warnings <= MAX_WARNINGS {
        eprintln!("[markdup] warning, {message}");
    }
}
