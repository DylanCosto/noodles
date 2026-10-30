//! `ntools sort`: sort alignments by coordinate, read name, or tag (samtools sort).

use std::path::{Path, PathBuf};

use anyhow::{Context, bail};
use noodles_sam as sam;

use crate::{
    alignment::{AlignmentReader, RawEncoder, ReadOptions, Rec},
    commands::index::{IndexKind, index_file},
    output::{self, AlignmentWriter, OutputFormat, WriteOptions},
    program,
    runs::{self, Chunk, SpillFile},
    sort_order::{Order, TagTieBreak},
    source::STDIN,
    template_coordinate::{self, HexIds, Libraries},
};

/// samtools' default per-thread memory limit.
const DEFAULT_MEMORY_PER_THREAD: &str = "768M";

#[derive(Debug, clap::Args)]
pub struct Args {
    /// Input SAM, BAM, or CRAM file (`-` for stdin).
    #[arg(default_value = STDIN)]
    input: String,

    /// Write output to this file instead of stdout (`FILE##idx##INDEX` also names the index for
    /// --write-index).
    #[arg(short = 'o', long = "output")]
    output: Option<PathBuf>,

    /// Output format [default: from the -o extension, else BAM].
    #[arg(short = 'O', long = "output-fmt", ignore_case = true)]
    output_format: Option<OutputFormat>,

    /// BGZF compression level for BAM output.
    #[arg(short = 'l', default_value_t = 6)]
    level: u8,

    /// Output uncompressed BAM.
    #[arg(short = 'u')]
    uncompressed: bool,

    /// Maximum memory per thread before spilling to temporary files (K/M/G suffixes).
    #[arg(short = 'm', default_value = DEFAULT_MEMORY_PER_THREAD, value_parser = parse_memory)]
    memory_per_thread: usize,

    /// Sort by read name, comparing numbers naturally.
    #[arg(short = 'n', overrides_with = "template_coordinate")]
    name_natural: bool,

    /// Sort by read name, bytewise (also with -n).
    #[arg(short = 'N', overrides_with = "template_coordinate")]
    name_ascii: bool,

    /// Sort by template-coordinate: the unclipped 5' positions of both ends of the template
    /// (the mate's from its MC tag, as added by fixmate -m), library, CB, MI, then read name.
    /// The last of -n, -N, and this option applies; -t overrides it.
    #[arg(
        long = "template-coordinate",
        overrides_with_all = ["name_natural", "name_ascii"]
    )]
    template_coordinate: bool,

    /// Sort by the value of TAG (records without it first), then by position, or by read name
    /// with -n/-N.
    #[arg(short = 't', value_name = "TAG", value_parser = parse_tag)]
    tag: Option<[u8; 2]>,

    /// Do not use the reverse strand. Accepted for compatibility: in samtools it only affects
    /// minimiser sorting (-M), which ntools does not support.
    #[arg(short = 'R')]
    forward_only: bool,

    /// Squash homopolymers. Accepted for compatibility: in samtools it only affects minimiser
    /// sorting (-M), which ntools does not support.
    #[arg(short = 'H')]
    squash_homopolymers: bool,

    /// Prefix for temporary files [default: next to the output, or in the current directory].
    #[arg(short = 'T')]
    temp_prefix: Option<PathBuf>,

    /// Reference FASTA (indexed) for CRAM input or output.
    #[arg(long = "reference")]
    reference: Option<PathBuf>,

    /// Additional threads for sorting, compression, and decompression.
    #[arg(short = '@', long = "threads", default_value_t = 0)]
    threads: usize,

    /// Index the output (coordinate order only): CSI for BAM, or the index named with
    /// `-o FILE##idx##INDEX` (BAI if it ends in .bai); CRAI for CRAM.
    #[arg(long = "write-index")]
    write_index: bool,

    /// Do not add an @PG line to the header.
    #[arg(long = "no-PG")]
    no_pg: bool,
}

pub fn run(args: Args) -> anyhow::Result<()> {
    let order = order(&args);

    let (output, index) = match args.output.as_deref() {
        Some(path) => {
            let (output, index) = split_index_path(path);
            (Some(output), index)
        }
        None => (None, None),
    };

    let format = args
        .output_format
        .or_else(|| output.as_deref().and_then(OutputFormat::from_path))
        .unwrap_or(OutputFormat::Bam);

    let index_request = if !args.write_index {
        None
    } else if order == Order::Coordinate {
        Some(IndexRequest::new(output.as_deref(), index, format)?)
    } else {
        eprintln!(
            "[W::bam_sort] Ignoring --write-index as it only works for position sorted files."
        );
        None
    };

    let read_options = ReadOptions {
        threads: args.threads,
        reference: args.reference.clone(),
    };

    let (mut reader, _) = AlignmentReader::open(&args.input, &read_options)?;
    let mut header = reader.read_header()?;

    let memory_limit = args.memory_per_thread.saturating_mul(args.threads + 1);
    let temp_dir = create_temp_dir(args.temp_prefix.as_deref(), output.as_deref())?;

    let sorted = SortedRecords::collect(
        &mut reader,
        &header,
        order,
        memory_limit,
        args.threads,
        temp_dir.path(),
        |_| Ok(true),
    )?;

    set_sort_order(&mut header, order);

    if !args.no_pg {
        program::add(&mut header)?;
    }

    let write_options = WriteOptions {
        format,
        level: output::parse_level(if args.uncompressed { 0 } else { args.level })?,
        threads: args.threads,
        reference: args.reference,
    };

    let mut writer = AlignmentWriter::create(output.as_deref(), &write_options, &header)?;
    writer.write_header(&header)?;
    sorted.write_to(&mut writer, &header)?;
    writer.finish(&header)?;

    if let Some(request) = index_request {
        request.write(args.threads)?;
    }

    Ok(())
}

/// The sort order, as samtools picks it: the last of `-n`/`-N` (bytewise if `-N` is given at
/// all) and `--template-coordinate`, except that `-t` sorts by tag, then by name if `-n`/`-N`
/// came last.
fn order(args: &Args) -> Order {
    let then = name_tie_break(args.name_natural, args.name_ascii);

    match (args.tag, then) {
        (Some(tag), then) => Order::Tag { tag, then },
        (None, TagTieBreak::NameNatural) => Order::NameNatural,
        (None, TagTieBreak::NameAscii) => Order::NameAscii,
        (None, TagTieBreak::Coordinate) if args.template_coordinate => Order::TemplateCoordinate,
        (None, TagTieBreak::Coordinate) => Order::Coordinate,
    }
}

/// The name order of `-n` and `-N` (bytewise if `-N` is given at all), or `Coordinate` for
/// neither.
pub fn name_tie_break(name_natural: bool, name_ascii: bool) -> TagTieBreak {
    if name_ascii {
        TagTieBreak::NameAscii
    } else if name_natural {
        TagTieBreak::NameNatural
    } else {
        TagTieBreak::Coordinate
    }
}

/// The delimiter between an output file name and its index file name (`HTS_IDX_DELIM`).
const INDEX_DELIMITER: &str = "##idx##";

/// Splits `FILE##idx##INDEX` into the file and the index; other paths have no index.
pub fn split_index_path(path: &Path) -> (PathBuf, Option<PathBuf>) {
    match path.to_str().and_then(|s| s.split_once(INDEX_DELIMITER)) {
        Some((file, index)) => (PathBuf::from(file), Some(PathBuf::from(index))),
        None => (path.to_path_buf(), None),
    }
}

/// An index to build for an output file once it is written (`--write-index`, `auto_index`).
pub struct IndexRequest {
    output: PathBuf,
    index: Option<PathBuf>,
    kind: IndexKind,
}

impl IndexRequest {
    /// Checks that `output` (`None` or `-` for stdout) can be indexed. The index goes to
    /// `index`, or next to the output: CSI for BAM unless `index` ends in `.bai`, and CRAI for
    /// CRAM.
    pub fn new(
        output: Option<&Path>,
        index: Option<PathBuf>,
        format: OutputFormat,
    ) -> anyhow::Result<Self> {
        let Some(output) = output.filter(|p| *p != Path::new("-")) else {
            bail!("--write-index needs an output file");
        };

        let is_bai = index
            .as_deref()
            .is_some_and(|p| p.as_os_str().as_encoded_bytes().ends_with(b".bai"));

        let kind = match format {
            OutputFormat::Bam if is_bai => IndexKind::Default,
            OutputFormat::Bam => IndexKind::Csi,
            OutputFormat::Cram => IndexKind::Default,
            OutputFormat::Sam => bail!("--write-index needs BAM or CRAM output"),
        };

        Ok(Self {
            output: output.to_path_buf(),
            index,
            kind,
        })
    }

    /// Indexes the written output.
    pub fn write(self, threads: usize) -> anyhow::Result<()> {
        index_file(&self.output, self.index, self.kind, threads)
    }
}

/// Parses a tag name; samtools uses its first two characters.
pub fn parse_tag(s: &str) -> anyhow::Result<[u8; 2]> {
    match *s.as_bytes() {
        [a, b, ..] => Ok([a, b]),
        _ => bail!("invalid tag: {s}"),
    }
}

/// Records read and sorted in memory, with any sorted runs spilled to temporary files.
pub struct SortedRecords<'d> {
    chunk: Chunk,
    spill_paths: Vec<PathBuf>,
    order: Order,
    threads: usize,
    temp_dir: &'d Path,
}

impl<'d> SortedRecords<'d> {
    /// Reads every record that `keep` accepts, spilling sorted runs to `temp_dir` whenever the
    /// buffered records reach `memory_limit` bytes.
    pub fn collect<K>(
        reader: &mut AlignmentReader,
        header: &sam::Header,
        order: Order,
        memory_limit: usize,
        threads: usize,
        temp_dir: &'d Path,
        mut keep: K,
    ) -> anyhow::Result<Self>
    where
        K: FnMut(&Rec<'_>) -> anyhow::Result<bool>,
    {
        let mut sorted = Self {
            chunk: Chunk::new(order),
            spill_paths: Vec::new(),
            order,
            threads,
            temp_dir,
        };

        let mut encoder = RawEncoder::new();

        // Template-coordinate keys depend on the header, so records carry them.
        let libraries =
            (order == Order::TemplateCoordinate).then(|| Libraries::from_header(header));
        let mut keyed = Vec::new();

        reader.for_each(header, |record| {
            if !keep(&record)? {
                return Ok(());
            }

            // Non-BAM records are encoded once so sorting and spilling work on raw bytes.
            let raw = record.encode_raw(header, &mut encoder)?;

            match &libraries {
                None => sorted.chunk.push(raw),
                Some(libraries) => {
                    keyed.clear();
                    let library = libraries.get(raw, HexIds::Ignore)?;
                    template_coordinate::push_keyed(&mut keyed, raw, library)?;
                    sorted.chunk.push(&keyed);
                }
            }

            if sorted.chunk.size() >= memory_limit {
                sorted.spill()?;
            }

            Ok(())
        })?;

        Ok(sorted)
    }

    fn spill(&mut self) -> anyhow::Result<()> {
        let index = self.spill_paths.len();
        let path = runs::spill(
            &mut self.chunk,
            self.order,
            self.temp_dir,
            index,
            self.threads,
        )?;
        self.spill_paths.push(path);
        Ok(())
    }

    /// Writes the records in order. `header` must have the input's reference sequences.
    pub fn write_to(
        mut self,
        writer: &mut AlignmentWriter,
        header: &sam::Header,
    ) -> anyhow::Result<()> {
        let order = self.order;

        if self.spill_paths.is_empty() {
            self.chunk.drain_sorted(order, self.threads, |record| {
                writer.write_raw(header, unkeyed(order, record))
            })?;

            return Ok(());
        }

        if !self.chunk.is_empty() {
            self.spill()?;
        }

        // Spill files are decompressed on the calling thread; there may be many of them.
        let mut sources = self
            .spill_paths
            .iter()
            .map(|path| SpillFile::open(path, 0))
            .collect::<Result<Vec<_>, _>>()
            .context("opening temporary files")?;

        runs::merge(&mut sources, order, |record| {
            writer.write_raw(header, unkeyed(order, record))?;
            Ok(())
        })
    }
}

/// A sorted record without the key that [`Order::TemplateCoordinate`] appends.
pub fn unkeyed(order: Order, record: &[u8]) -> &[u8] {
    match order {
        Order::TemplateCoordinate => template_coordinate::record(record),
        _ => record,
    }
}

/// Sets `@HD SO` (and `SS`, and `GO` for template-coordinate) for `order` as samtools does,
/// creating `@HD` if needed, and removes `GO` otherwise.
pub fn set_sort_order(header: &mut sam::Header, order: Order) {
    let (sort_order, subsort_order) = order.header_values();
    let group_order = (order == Order::TemplateCoordinate).then_some("query");
    set_sort_fields(header, sort_order, subsort_order.as_deref(), group_order);
}

/// Sets `@HD SO`, `SS`, and `GO` (removing `SS` and `GO` when they are `None`), creating `@HD`
/// if needed (`sam_hdr_update_sort` and `set_sort_order` in samtools).
pub fn set_sort_fields(
    header: &mut sam::Header,
    sort_order: &str,
    subsort_order: Option<&str>,
    group_order: Option<&str>,
) {
    use noodles_sam::header::record::value::{
        Map,
        map::header::tag::{GROUP_ORDER, SORT_ORDER, SUBSORT_ORDER},
    };

    let hd = header.header_mut().get_or_insert_with(Map::default);
    let fields = hd.other_fields_mut();
    fields.insert(SORT_ORDER, sort_order.into());

    for (tag, value) in [(GROUP_ORDER, group_order), (SUBSORT_ORDER, subsort_order)] {
        match value {
            Some(value) => {
                fields.insert(tag, value.into());
            }
            None => {
                fields.shift_remove(&tag);
            }
        }
    }
}

pub fn create_temp_dir(
    prefix: Option<&Path>,
    output: Option<&Path>,
) -> anyhow::Result<tempfile::TempDir> {
    let (dir, name) = match prefix {
        Some(prefix) => (
            prefix.parent().filter(|p| !p.as_os_str().is_empty()),
            prefix
                .file_name()
                .and_then(|s| s.to_str())
                .unwrap_or("ntools-sort"),
        ),
        None => (
            output
                .and_then(Path::parent)
                .filter(|p| !p.as_os_str().is_empty()),
            "ntools-sort",
        ),
    };

    let dir = dir.unwrap_or_else(|| Path::new("."));

    tempfile::Builder::new()
        .prefix(&format!("{name}."))
        .tempdir_in(dir)
        .with_context(|| format!("creating a temporary directory in {}", dir.display()))
}

/// Parses a memory size with an optional K, M, or G suffix (powers of 1024), as samtools does.
pub fn parse_memory(s: &str) -> anyhow::Result<usize> {
    let s = s.trim();
    let (digits, multiplier) = match s.char_indices().last() {
        Some((i, 'k' | 'K')) => (&s[..i], 1 << 10),
        Some((i, 'm' | 'M')) => (&s[..i], 1 << 20),
        Some((i, 'g' | 'G')) => (&s[..i], 1 << 30),
        _ => (s, 1),
    };

    let n: usize = digits
        .parse()
        .with_context(|| format!("invalid memory size: {s}"))?;

    if n == 0 {
        bail!("memory size must be positive: {s}");
    }

    n.checked_mul(multiplier)
        .with_context(|| format!("memory size overflows: {s}"))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_parse_memory() -> anyhow::Result<()> {
        assert_eq!(parse_memory("768M")?, 768 << 20);
        assert_eq!(parse_memory("2g")?, 2 << 30);
        assert_eq!(parse_memory("100")?, 100);
        assert!(parse_memory("0").is_err());
        assert!(parse_memory("1X").is_err());
        Ok(())
    }
}
