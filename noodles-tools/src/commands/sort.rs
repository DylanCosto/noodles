//! `ntools sort`: sort alignments by coordinate or read name (samtools sort).

use std::path::{Path, PathBuf};

use anyhow::{Context, bail};
use noodles_bam as bam;
use noodles_sam::{self as sam, alignment::io::Write as _};

use crate::{
    alignment::{AlignmentReader, ReadOptions, Rec},
    output::{self, AlignmentWriter, OutputFormat, WriteOptions},
    program,
    runs::{self, Chunk, SpillFile},
    sort_order::Order,
    source::STDIN,
};

/// samtools' default per-thread memory limit.
const DEFAULT_MEMORY_PER_THREAD: &str = "768M";

#[derive(Debug, clap::Args)]
pub struct Args {
    /// Input SAM, BAM, or CRAM file (`-` for stdin).
    #[arg(default_value = STDIN)]
    input: String,

    /// Write output to this file instead of stdout.
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
    #[arg(short = 'n', conflicts_with = "name_ascii")]
    name_natural: bool,

    /// Sort by read name, bytewise.
    #[arg(short = 'N')]
    name_ascii: bool,

    /// Prefix for temporary files [default: next to the output, or in the current directory].
    #[arg(short = 'T')]
    temp_prefix: Option<PathBuf>,

    /// Reference FASTA (indexed) for CRAM input or output.
    #[arg(long = "reference")]
    reference: Option<PathBuf>,

    /// Additional threads for sorting, compression, and decompression.
    #[arg(short = '@', long = "threads", default_value_t = 0)]
    threads: usize,

    /// Do not add an @PG line to the header.
    #[arg(long = "no-PG")]
    no_pg: bool,
}

pub fn run(args: Args) -> anyhow::Result<()> {
    let order = if args.name_natural {
        Order::NameNatural
    } else if args.name_ascii {
        Order::NameAscii
    } else {
        Order::Coordinate
    };

    let format = args
        .output_format
        .or_else(|| args.output.as_deref().and_then(OutputFormat::from_path))
        .unwrap_or(OutputFormat::Bam);

    let read_options = ReadOptions {
        threads: args.threads,
        reference: args.reference.clone(),
    };

    let (mut reader, _) = AlignmentReader::open(&args.input, &read_options)?;
    let mut header = reader.read_header()?;

    let memory_limit = args.memory_per_thread.saturating_mul(args.threads + 1);
    let temp_dir = create_temp_dir(args.temp_prefix.as_deref(), args.output.as_deref())?;

    let sorted = SortedRecords::collect(
        &mut reader,
        &header,
        order,
        memory_limit,
        args.threads,
        temp_dir.path(),
        |_| Ok(true),
    )?;

    set_sort_order(&mut header, order)?;

    if !args.no_pg {
        program::add(&mut header)?;
    }

    let write_options = WriteOptions {
        format,
        level: output::parse_level(if args.uncompressed { 0 } else { args.level })?,
        threads: args.threads,
        reference: args.reference,
    };

    let mut writer = AlignmentWriter::create(args.output.as_deref(), &write_options, &header)?;
    writer.write_header(&header)?;
    sorted.write_to(&mut writer, &header)?;
    writer.finish(&header)?;

    Ok(())
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

        let mut encoder = bam::io::Writer::from(Vec::new());

        reader.for_each(header, |record| {
            if !keep(&record)? {
                return Ok(());
            }

            match record {
                Rec::Bam(r) => sorted.chunk.push(r.as_ref()),
                _ => {
                    // Non-BAM records are encoded once so sorting and spilling work on raw bytes.
                    encoder.get_mut().clear();
                    encoder.write_alignment_record(header, record.as_dyn())?;
                    sorted.chunk.push(&encoder.get_ref()[size_of::<u32>()..]);
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
        if self.spill_paths.is_empty() {
            self.chunk
                .drain_sorted(self.order, self.threads, |record| {
                    writer.write_raw(header, record)
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

        runs::merge(&mut sources, self.order, |record| {
            writer.write_raw(header, record)?;
            Ok(())
        })
    }
}

/// Sets `@HD SO` (and `SS` for name orders) as samtools does, creating `@HD` if needed.
pub fn set_sort_order(header: &mut sam::Header, order: Order) -> anyhow::Result<()> {
    use noodles_sam::header::record::value::{
        Map,
        map::header::tag::{SORT_ORDER, SUBSORT_ORDER},
    };

    let (sort_order, subsort_order) = order.header_values();

    let hd = header.header_mut().get_or_insert_with(Map::default);
    let fields = hd.other_fields_mut();
    fields.insert(SORT_ORDER, sort_order.into());

    match subsort_order {
        Some(ss) => {
            fields.insert(SUBSORT_ORDER, ss.into());
        }
        None => {
            fields.shift_remove(&SUBSORT_ORDER);
        }
    }

    Ok(())
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
