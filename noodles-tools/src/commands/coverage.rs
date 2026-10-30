//! `ntools coverage`: per-reference coverage summary or histogram (samtools coverage).

mod histogram;

use std::{
    fs::File,
    io::{self, BufWriter, Write},
    path::{Path, PathBuf},
    thread,
};

use anyhow::{Context, bail};
use crossbeam_channel::{Receiver, Sender};
use noodles_sam::{self as sam, alignment::record::Flags};

use crate::{
    alignment::{AlignmentIndex, AlignmentReader, ReadOptions, Rec},
    depth::{Cell, CountOptions, DEFAULT_EXCLUDE_FLAGS, DepthWindow, cigar_query_length, format_g},
    file_list, flags, hts_region, output,
};

use self::histogram::{Blocks, Plot};
use crate::global_args::{IgnoredWriteIndex, InputFormatOptions, Verbosity};

/// samtools' default `-d`.
const DEFAULT_MAX_DEPTH: i64 = 1_000_000;

/// HTSlib's pileup limit when `-d` is negative (`bam_plp_init`).
const PILEUP_DEFAULT_MAX_COUNT: usize = 8000;

/// Histogram bins without `-w` or a terminal width.
const DEFAULT_BINS: i64 = 50;

/// Positions sent per message from a reader thread.
const BATCH_LEN: usize = 1 << 16;

#[derive(Debug, clap::Args)]
#[command(disable_help_flag = true)]
pub struct Args {
    /// Coordinate-sorted SAM, BAM, or CRAM files.
    #[arg(required_unless_present = "bam_list")]
    inputs: Vec<String>,

    /// Read the input file names from this file, one per line, instead of the arguments.
    #[arg(short = 'b', long = "bam-list")]
    bam_list: Option<PathBuf>,

    /// Exclude reads shorter than this.
    #[arg(
        short = 'l',
        long = "min-read-len",
        default_value_t = 0,
        allow_negative_numbers = true
    )]
    min_read_length: i64,

    /// Exclude reads with mapping quality below this value.
    #[arg(
        short = 'q',
        long = "min-MQ",
        alias = "min-mq",
        default_value_t = 0,
        allow_negative_numbers = true
    )]
    min_mapping_quality: i64,

    /// Exclude bases with base quality below this value.
    #[arg(short = 'Q', long = "min-BQ", alias = "min-bq", default_value_t = 0)]
    min_base_quality: u8,

    /// Only include reads with at least one of these flags.
    #[arg(long = "rf", alias = "incl-flags", value_parser = flags::parse)]
    require_flags: Option<Flags>,

    /// Exclude reads with any of these flags [UNMAP,SECONDARY,QCFAIL,DUP].
    #[arg(long = "ff", alias = "excl-flags", value_parser = flags::parse)]
    exclude_flags: Option<Flags>,

    /// Maximum reads piled up from one start position, per input [1000000]; 0 for no limit.
    #[arg(
        short = 'd',
        long = "depth",
        default_value_t = DEFAULT_MAX_DEPTH,
        allow_negative_numbers = true
    )]
    max_depth: i64,

    /// Ignore positions with a depth below this.
    #[arg(long = "min-depth", default_value_t = 1, allow_negative_numbers = true)]
    min_depth: i64,

    /// Show a histogram of covered bases instead of the table.
    #[arg(short = 'm', long = "histogram")]
    histogram: bool,

    /// Plot the depth instead of the table.
    #[arg(short = 'D', long = "plot-depth")]
    plot_depth: bool,

    /// Draw the histogram with ASCII characters only.
    #[arg(short = 'A', long = "ascii")]
    ascii: bool,

    /// Number of histogram bins [terminal width - 40].
    #[arg(short = 'w', long = "n-bins", allow_negative_numbers = true)]
    n_bins: Option<i64>,

    /// Write output to this file instead of stdout.
    #[arg(short = 'o', long = "output")]
    output: Option<PathBuf>,

    /// Do not print the header line of the table.
    #[arg(short = 'H', long = "no-header")]
    no_header: bool,

    /// Only summarize this region.
    #[arg(short = 'r', long = "region")]
    region: Option<String>,

    /// Reference FASTA (indexed) for CRAM input.
    #[arg(long = "reference")]
    reference: Option<PathBuf>,

    /// Additional decompression threads, per input.
    #[arg(short = '@', long = "threads", default_value_t = 0)]
    threads: usize,

    /// Print help.
    #[arg(short = 'h', long, action = clap::ArgAction::Help)]
    help: Option<bool>,

    #[command(flatten)]
    input_format_options: InputFormatOptions,

    #[command(flatten)]
    write_index: IgnoredWriteIndex,

    #[command(flatten)]
    verbosity: Verbosity,
}

/// Whether the table starts with a line of column names (`-H` omits it).
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
enum ColumnNames {
    Print,
    Omit,
}

/// What is written.
#[derive(Clone, Copy, Debug)]
enum Mode {
    Table,
    Histogram {
        plot: Plot,
        blocks: Blocks,
        n_bins: i64,
    },
}

impl Mode {
    fn new(args: &Args) -> Self {
        if !(args.histogram || args.plot_depth || args.ascii || args.n_bins.is_some()) {
            return Self::Table;
        }

        let plot = if args.plot_depth {
            Plot::Depth
        } else {
            Plot::Breadth
        };

        let blocks = if args.ascii {
            Blocks::Ascii
        } else {
            Blocks::Utf8
        };

        // As in samtools, the terminal decides the width unless -w or -o is given.
        let fixed = args.n_bins.is_some() || args.output.is_some();
        let mut n_bins = args.n_bins.unwrap_or(DEFAULT_BINS);

        if n_bins <= 0 || !fixed {
            let columns = histogram::terminal_columns();
            n_bins = if columns > 60 { columns - 40 } else { 40 };
        }

        Self::Histogram {
            plot,
            blocks,
            n_bins,
        }
    }
}

/// Read filters, as samtools coverage applies them while reading.
#[derive(Clone, Copy, Debug)]
struct ReadFilter {
    exclude: Flags,
    /// Reads need at least one of these (when not empty).
    require_any: Flags,
    min_mapping_quality: i64,
    min_read_length: i64,
}

impl ReadFilter {
    fn passes(&self, record: &Rec<'_>, flags: Flags, mapping_quality: u8) -> io::Result<bool> {
        if flags.intersects(self.exclude)
            || (!self.require_any.is_empty() && !flags.intersects(self.require_any))
            || i64::from(mapping_quality) < self.min_mapping_quality
        {
            return Ok(false);
        }

        if self.min_read_length > 0 {
            let len = i64::try_from(cigar_query_length(record)?).unwrap_or(i64::MAX);

            if len < self.min_read_length {
                return Ok(false);
            }
        }

        Ok(true)
    }
}

/// The region of `-r`: a reference sequence ID and a 0-based, half-open interval.
#[derive(Clone, Copy, Debug)]
struct Region {
    id: usize,
    start: i64,
    end: i64,
}

/// What every input is scanned with.
struct Scan {
    filter: ReadFilter,
    count_options: CountOptions,
    /// HTSlib's pileup `maxcnt`, if it can be reached.
    max_count: Option<usize>,
    region: Option<String>,
    reference_count: usize,
}

/// Per-reference read counts from one input.
#[derive(Clone, Copy, Debug, Default)]
struct ReadCounts {
    /// Every read, before filtering (samtools' `n_reads`).
    reads: u64,
    /// Reads that pass the filters.
    selected: u64,
    mapping_quality_sum: u64,
}

pub fn run(args: Args) -> anyhow::Result<()> {
    let sources = match &args.bam_list {
        Some(list) => file_list::read(list)?,
        None => args.inputs.clone(),
    };

    let mut options = ReadOptions {
        threads: args.threads,
        reference: args.reference.clone(),
    };
    args.input_format_options.configure(&mut options)?;

    let mut inputs = Vec::with_capacity(sources.len());

    for src in &sources {
        let (mut reader, path) = AlignmentReader::open(src, &options)?;
        args.input_format_options.apply(&mut reader)?;
        // As in samtools (CRAM_OPT_DECODE_MD set to 0); coverage never reads MD or NM.
        reader.keep_stored_md_nm();
        let header = reader
            .read_header()
            .with_context(|| format!("Could not read header for \"{src}\""))?;
        inputs.push((reader, path, header));
    }

    let header = inputs[0].2.clone();

    // samtools reports nothing for the unmapped reads ("*").
    let region = match args.region.as_deref() {
        Some("*") => return Ok(()),
        Some(raw) => Some(parse_region(&header, raw)?),
        None => None,
    };

    let mode = Mode::new(&args);

    let max_count = match args.max_depth {
        0 => None,
        n if n < 0 => Some(PILEUP_DEFAULT_MAX_COUNT),
        n => Some(usize::try_from(n).unwrap_or(usize::MAX)),
    };

    let scan = Scan {
        filter: ReadFilter {
            exclude: args.exclude_flags.unwrap_or(DEFAULT_EXCLUDE_FLAGS),
            require_any: args.require_flags.unwrap_or_else(Flags::empty),
            min_mapping_quality: args.min_mapping_quality,
            min_read_length: args.min_read_length,
        },
        count_options: CountOptions {
            min_base_quality: args.min_base_quality,
            include_deletions: false,
            track_quality: true,
            include_skips: false,
        },
        max_count,
        region: args.region.clone(),
        reference_count: header.reference_sequences().len(),
    };

    let mut summary = Summary::new(&header, region, mode, args.min_depth);

    let counts = if let [(reader, path, input_header)] = inputs.as_mut_slice() {
        scan_input(
            reader,
            path.as_deref(),
            input_header,
            &scan,
            &mut |id, pos, cell| {
                summary.add(id, pos, &cell);
                Ok(())
            },
        )
        .map(|counts| vec![counts])?
    } else {
        scan_inputs(&mut inputs, &scan, &mut summary)?
    };

    for counts in counts {
        summary.add_read_counts(&counts);
    }

    let sink: output::Sink = match args.output.as_deref().filter(|p| *p != Path::new("-")) {
        Some(path) => Box::new(
            File::create(path)
                .with_context(|| format!("Cannot open \"{}\" for writing", path.display()))?,
        ),
        None => output::raw_stdout()?,
    };

    let mut out = BufWriter::new(sink);

    match mode {
        Mode::Table => {
            let column_names = if args.no_header {
                ColumnNames::Omit
            } else {
                ColumnNames::Print
            };
            summary.write_table(&mut out, &header, column_names)?;
        }
        Mode::Histogram { plot, blocks, .. } => {
            summary.write_histograms(&mut out, &header, plot, blocks)?;
        }
    }

    out.flush()?;

    if summary.missing_quality {
        eprintln!(
            "samtools coverage: Warning:  Missing quality values in alignments.  Mean base \
             quality calculated only on available values."
        );
    }

    Ok(())
}

fn parse_region(header: &sam::Header, raw: &str) -> anyhow::Result<Region> {
    let reference_sequences = header.reference_sequences();

    let Some((id, start, end)) =
        hts_region::parse_region(raw, |name| reference_sequences.get_index_of(name))
    else {
        bail!(
            "Failed to parse region \"{raw}\". Check the region format or region name presence \
             in the file"
        );
    };

    let end = if end == hts_region::POS_MAX {
        reference_sequences.get_index(id).map_or(0, |(_, rs)| {
            i64::try_from(rs.length().get()).unwrap_or(i64::MAX)
        })
    } else {
        end
    };

    Ok(Region { id, start, end })
}

/// Scans several inputs in parallel, summing their depths at each position.
fn scan_inputs(
    inputs: &mut [(AlignmentReader, Option<PathBuf>, sam::Header)],
    scan: &Scan,
    summary: &mut Summary,
) -> anyhow::Result<Vec<Vec<ReadCounts>>> {
    thread::scope(|scope| {
        let mut streams = Vec::with_capacity(inputs.len());
        let mut handles = Vec::with_capacity(inputs.len());

        for (reader, path, header) in inputs.iter_mut() {
            let (tx, rx) = crossbeam_channel::bounded::<Vec<Point>>(4);

            handles
                .push(scope.spawn(move || send_depths(reader, path.as_deref(), header, scan, &tx)));
            streams.push(Stream::new(rx));
        }

        while let Some((id, pos, cell)) = next_position(&mut streams) {
            summary.add(id, pos, &cell);
        }

        drop(streams);

        handles
            .into_iter()
            .map(|handle| {
                handle
                    .join()
                    .map_err(|_| anyhow::anyhow!("coverage thread panicked"))?
            })
            .collect()
    })
}

/// `(reference sequence ID, 1-based position, cell)`.
type Point = (usize, usize, Cell);

/// Scans one input, sending its depths in batches.
fn send_depths(
    reader: &mut AlignmentReader,
    path: Option<&Path>,
    header: &sam::Header,
    scan: &Scan,
    tx: &Sender<Vec<Point>>,
) -> anyhow::Result<Vec<ReadCounts>> {
    let mut batch = Vec::with_capacity(BATCH_LEN);

    let counts = scan_input(reader, path, header, scan, &mut |id, pos, cell| {
        batch.push((id, pos, cell));

        if batch.len() >= BATCH_LEN {
            let full = std::mem::replace(&mut batch, Vec::with_capacity(BATCH_LEN));
            tx.send(full)
                .map_err(|_| io::Error::new(io::ErrorKind::BrokenPipe, "coverage stopped"))?;
        }

        Ok(())
    })?;

    if !batch.is_empty() {
        let _ = tx.send(batch);
    }

    Ok(counts)
}

/// One input's ordered depths.
struct Stream {
    rx: Receiver<Vec<Point>>,
    batch: Vec<Point>,
    i: usize,
}

impl Stream {
    fn new(rx: Receiver<Vec<Point>>) -> Self {
        Self {
            rx,
            batch: Vec::new(),
            i: 0,
        }
    }

    fn peek(&mut self) -> Option<&Point> {
        while self.i >= self.batch.len() {
            self.batch = self.rx.recv().ok()?;
            self.i = 0;
        }

        Some(&self.batch[self.i])
    }
}

/// Pops the smallest position across streams, with the inputs' cells there summed.
fn next_position(streams: &mut [Stream]) -> Option<Point> {
    let next = streams
        .iter_mut()
        .filter_map(|s| s.peek().map(|&(id, pos, _)| (id, pos)))
        .min()?;

    let mut sum = Cell::default();

    for stream in streams.iter_mut() {
        if let Some(&(id, pos, cell)) = stream.peek()
            && (id, pos) == next
        {
            sum.depth += cell.depth;
            sum.quality_sum += cell.quality_sum;
            sum.missing_quality += cell.missing_quality;
            stream.i += 1;
        }
    }

    Some((next.0, next.1, sum))
}

/// Reads one input, passing every position a read piles up on (in order) to `emit`, and returns
/// its read counts by reference sequence.
fn scan_input<F>(
    reader: &mut AlignmentReader,
    path: Option<&Path>,
    header: &sam::Header,
    scan: &Scan,
    emit: &mut F,
) -> anyhow::Result<Vec<ReadCounts>>
where
    F: FnMut(usize, usize, Cell) -> io::Result<()>,
{
    let mut counts = vec![ReadCounts::default(); scan.reference_count];
    let reference_count = header.reference_sequences().len();

    let mut window = match scan.max_count {
        Some(max_count) => DepthWindow::with_max_count(scan.count_options, max_count),
        None => DepthWindow::new(scan.count_options),
    };
    let mut last_id = 0;

    let mut visit = |record: Rec<'_>| -> anyhow::Result<()> {
        let flags = record.flags()?;
        let id = record
            .reference_sequence_id(header)?
            .filter(|&id| id < reference_count);
        if let Some(c) = id.and_then(|id| counts.get_mut(id)) {
            c.reads += 1;
        }

        let mapping_quality = record.mapping_quality()?.unwrap_or(u8::MAX);

        if !scan.filter.passes(&record, flags, mapping_quality)? {
            return Ok(());
        }

        let Some(id) = id else {
            return Ok(());
        };

        if let Some(c) = counts.get_mut(id) {
            c.selected += 1;
            c.mapping_quality_sum += u64::from(mapping_quality);
        }

        // The pileup leaves out unmapped reads.
        if flags.is_unmapped() {
            return Ok(());
        }

        if id < last_id {
            bail!("the input is not sorted (chromosomes out of order)");
        }

        last_id = id;

        window.add(&record, header, emit)?;

        Ok(())
    };

    match &scan.region {
        Some(raw) => {
            let path = path.context("coverage -r needs an indexed file, not stdin")?;
            let index = AlignmentIndex::read(reader.format(), path, None)
                .with_context(|| format!("Failed to load index for \"{}\"", path.display()))?;
            let names = header.reference_sequences();
            let Some((id, start, end)) =
                hts_region::parse_region(raw, |name| names.get_index_of(name))
            else {
                bail!("Failed to parse region \"{raw}\" in \"{}\"", path.display());
            };
            let start = usize::try_from(start).unwrap_or(0);
            let end = (end != hts_region::POS_MAX).then(|| usize::try_from(end).unwrap_or(0));

            reader.for_each_in_interval(&index, header, id, start, end, &mut visit)?;
        }
        None => reader.for_each(header, &mut visit)?,
    }

    window.finish(emit)?;

    Ok(counts)
}

/// One reference sequence's results.
#[derive(Clone, Debug, Default)]
struct Tally {
    counts: ReadCounts,
    /// Whether any read piled up here; samtools reports these sequences first.
    covered: bool,
    /// The reported 0-based, half-open range.
    start: i64,
    end: i64,
    covered_bases: u64,
    depth_sum: u64,
    quality_sum: u64,
    quality_count: u64,
    /// Histogram bins (allocated once covered) and the bases in each.
    bins: Vec<u64>,
    bin_width: i64,
}

impl Tally {
    fn len(&self) -> f64 {
        // As samtools computes it: `(double) end - beg`.
        self.end as f64 - self.start as f64
    }

    /// Sizes the histogram as samtools does: up to `n_bins` bins of equal width, leaving out
    /// the remainder at the end.
    fn allocate_bins(&mut self, n_bins: i64) {
        let len = self.end - self.start;
        let n = n_bins.min(len);
        self.bin_width = len / n.max(1);
        self.bins = vec![0; usize::try_from(n).unwrap_or(0)];
    }

    fn mean_base_quality(&self) -> f64 {
        if self.quality_count > 0 {
            self.quality_sum as f64 / self.quality_count as f64
        } else {
            0.0
        }
    }
}

/// Accumulates positions and read counts per reference sequence.
struct Summary {
    tallies: Vec<Tally>,
    region: Option<Region>,
    min_depth: u64,
    /// The histogram's bin count and what it counts, if drawn.
    histogram: Option<(i64, Plot)>,
    /// Whether a counted base had no quality.
    missing_quality: bool,
}

impl Summary {
    fn new(header: &sam::Header, region: Option<Region>, mode: Mode, min_depth: i64) -> Self {
        let mut tallies: Vec<Tally> = header
            .reference_sequences()
            .values()
            .map(|rs| Tally {
                end: i64::try_from(rs.length().get()).unwrap_or(i64::MAX),
                ..Tally::default()
            })
            .collect();

        if let Some(region) = region
            && let Some(tally) = tallies.get_mut(region.id)
        {
            tally.start = region.start;
            tally.end = region.end;
        }

        let histogram = match mode {
            Mode::Table => None,
            Mode::Histogram { plot, n_bins, .. } => Some((n_bins, plot)),
        };

        // A region's histogram is drawn even without reads.
        if let (Some(region), Some((n_bins, _))) = (region, histogram)
            && let Some(tally) = tallies.get_mut(region.id)
        {
            tally.allocate_bins(n_bins);
        }

        Self {
            tallies,
            region,
            // samtools ignores a --min-depth below 1.
            min_depth: u64::try_from(min_depth).unwrap_or(0).max(1),
            histogram,
            missing_quality: false,
        }
    }

    /// Adds a position some read piles up on (1-based), with the inputs' cells summed.
    #[inline]
    fn add(&mut self, id: usize, position: usize, cell: &Cell) {
        let Some(tally) = self.tallies.get_mut(id) else {
            return;
        };

        if !tally.covered {
            tally.covered = true;

            if let (None, Some((n_bins, _))) = (self.region, self.histogram) {
                tally.allocate_bins(n_bins);
            }
        }

        let position = position as i64 - 1;

        if position < tally.start || position >= tally.end {
            return;
        }

        self.missing_quality |= cell.missing_quality > 0;

        let depth = u64::from(cell.depth);
        // `min_depth` is at least 1, as a covered base needs a depth.
        let counted = depth >= self.min_depth;

        if counted {
            tally.covered_bases += 1;
            tally.depth_sum += depth;
            tally.quality_sum += cell.quality_sum;
            // Without deletions and skips, every counted base either has a quality or not.
            tally.quality_count += u64::from(cell.depth - cell.missing_quality);
        }

        if let Some((_, plot)) = self.histogram {
            let bin = usize::try_from((position - tally.start) / tally.bin_width)
                .ok()
                .and_then(|bin| tally.bins.get_mut(bin));

            match (plot, bin) {
                (Plot::Depth, Some(n)) => *n += depth,
                (Plot::Breadth, Some(n)) if counted => *n += 1,
                _ => {}
            }
        }
    }

    fn add_read_counts(&mut self, counts: &[ReadCounts]) {
        for (tally, c) in self.tallies.iter_mut().zip(counts) {
            tally.counts.reads += c.reads;
            tally.counts.selected += c.selected;
            tally.counts.mapping_quality_sum += c.mapping_quality_sum;
        }
    }

    /// The reference sequences reported, in samtools' order: those with reads piled up, then
    /// (without a region) the others.
    fn reported(&self) -> Vec<usize> {
        if let Some(region) = self.region {
            return vec![region.id];
        }

        let covered = (0..self.tallies.len()).filter(|&id| self.tallies[id].covered);
        let uncovered = (0..self.tallies.len()).filter(|&id| !self.tallies[id].covered);
        covered.chain(uncovered).collect()
    }

    fn write_table<W: Write>(
        &self,
        out: &mut W,
        header: &sam::Header,
        column_names: ColumnNames,
    ) -> io::Result<()> {
        let names: Vec<_> = header.reference_sequences().keys().collect();
        // As in samtools, the column names are written before the first row only.
        let mut print_header = column_names == ColumnNames::Print;

        for id in self.reported() {
            let tally = &self.tallies[id];

            if print_header {
                writeln!(
                    out,
                    "#rname\tstartpos\tendpos\tnumreads\tcovbases\tcoverage\tmeandepth\tmeanbaseq\tmeanmapq"
                )?;
                print_header = false;
            }

            let len = tally.len();
            let c = tally.counts;
            let mean_mapping_quality = if c.selected > 0 {
                c.mapping_quality_sum as f64 / c.selected as f64
            } else {
                0.0
            };

            out.write_all(names[id])?;
            writeln!(
                out,
                "\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
                tally.start + 1,
                tally.end,
                c.selected,
                tally.covered_bases,
                format_g(100.0 * tally.covered_bases as f64 / len, 6),
                format_g(tally.depth_sum as f64 / len, 6),
                format_g(tally.mean_base_quality(), 3),
                format_g(mean_mapping_quality, 3),
            )?;
        }

        Ok(())
    }

    fn write_histograms<W: Write>(
        &self,
        out: &mut W,
        header: &sam::Header,
        plot: Plot,
        blocks: Blocks,
    ) -> io::Result<()> {
        let ids: Vec<usize> = match self.region {
            Some(region) => vec![region.id],
            None => (0..self.tallies.len())
                .filter(|&id| self.tallies[id].covered)
                .collect(),
        };

        for (i, &id) in ids.iter().enumerate() {
            if i > 0 {
                writeln!(out)?;
            }

            let Some((name, reference_sequence)) = header.reference_sequences().get_index(id)
            else {
                continue;
            };

            histogram::write(
                out,
                name,
                reference_sequence.length().get(),
                &self.tallies[id],
                plot,
                blocks,
            )?;
        }

        Ok(())
    }
}
