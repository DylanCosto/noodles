//! `ntools faidx`: index a FASTA file or extract regions from it (samtools faidx).

use std::{
    collections::HashMap,
    fs::{self, File},
    io::{self, BufWriter, Write},
    num::NonZero,
    path::{Path, PathBuf},
};

use anyhow::Context;
use noodles_core::Region;
use noodles_fasta::{self as fasta, fai};

use crate::sequence::reverse_complement;

#[derive(Debug, clap::Args)]
pub struct Args {
    /// FASTA file.
    input: PathBuf,

    /// Regions to extract (`chr1`, `chr1:100-200`). With none, builds the `.fai` index.
    regions: Vec<String>,

    /// Read regions from this file, one per line.
    #[arg(short = 'r', long = "region-file")]
    region_file: Option<PathBuf>,

    /// Write output to this file instead of stdout.
    #[arg(short = 'o', long = "output")]
    output: Option<PathBuf>,

    /// Line width for extracted sequences [default: the input's line width, as in samtools].
    #[arg(short = 'n', long = "length")]
    line_width: Option<NonZero<usize>>,

    /// Reverse complement the extracted sequences.
    #[arg(short = 'i', long = "reverse-complement")]
    reverse_complement: bool,
}

pub fn run(args: Args) -> anyhow::Result<()> {
    let mut regions = args.regions;

    if let Some(path) = &args.region_file {
        let text =
            fs::read_to_string(path).with_context(|| format!("reading {}", path.display()))?;
        regions.extend(
            text.lines()
                .map(str::trim)
                .filter(|l| !l.is_empty())
                .map(String::from),
        );
    }

    if regions.is_empty() {
        return build_index(&args.input);
    }

    let index = fasta::fs::read_associated_index(&args.input).with_context(|| {
        format!(
            "reading the index of {} (is it indexed?)",
            args.input.display()
        )
    })?;

    let records: HashMap<&[u8], &fai::Record> = index
        .as_ref()
        .iter()
        .map(|record| (record.name().as_ref(), record))
        .collect();

    let file =
        File::open(&args.input).with_context(|| format!("opening {}", args.input.display()))?;

    let sink: Box<dyn Write> = match &args.output {
        Some(path) => {
            Box::new(File::create(path).with_context(|| format!("creating {}", path.display()))?)
        }
        None => Box::new(io::stdout().lock()),
    };
    let mut writer = BufWriter::new(sink);

    let mut raw = Vec::new();
    let mut sequence = Vec::new();

    for raw_region in &regions {
        let region: Region = raw_region
            .parse()
            .with_context(|| format!("invalid region: {raw_region}"))?;

        let name: &[u8] = region.name().as_ref();
        let record = *records
            .get(name)
            .with_context(|| format!("{} is not in the index", region.name()))?;

        let length = usize::try_from(record.length()).context("sequence length overflows usize")?;
        let interval = region.interval();
        let start = interval.start().map_or(1, usize::from);
        let requested_end = interval.end().map(usize::from);

        if requested_end.is_some_and(|end| end > length) {
            eprintln!("[faidx] Truncated sequence: {raw_region}");
        }

        let end = requested_end.unwrap_or(length).min(length);

        read_bases(&file, record, start, end, &mut raw, &mut sequence)
            .with_context(|| format!("retrieving {raw_region}"))?;

        if args.reverse_complement {
            reverse_complement(&mut sequence);
            writeln!(writer, ">{raw_region}/rc")?;
        } else {
            writeln!(writer, ">{raw_region}")?;
        }

        let line_width = match args.line_width {
            Some(n) => n.get(),
            None => usize::try_from(record.line_base_count().get())
                .context("line width overflows usize")?,
        };

        for line in sequence.chunks(line_width) {
            writer.write_all(line)?;
            writer.write_all(b"\n")?;
        }
    }

    writer.flush()?;

    Ok(())
}

/// Reads the bases of the 1-based closed interval `[start, end]` into `sequence`, reading the
/// whole span at once and dropping line endings, as samtools does.
fn read_bases(
    file: &File,
    record: &fai::Record,
    start: usize,
    end: usize,
    raw: &mut Vec<u8>,
    sequence: &mut Vec<u8>,
) -> anyhow::Result<()> {
    sequence.clear();

    if start > end || start == 0 {
        return Ok(());
    }

    let line_bases = record.line_base_count().get();
    let line_width = record.line_width().get();

    // The file offset of the base at a 0-based position.
    let offset = |position: u64| {
        record.position() + position / line_bases * line_width + position % line_bases
    };

    let first = offset(u64::try_from(start - 1)?);
    let last = offset(u64::try_from(end - 1)?);
    let len = usize::try_from(last - first + 1)?;

    raw.resize(len, 0);
    read_exact_at(file, raw, first)?;

    // Copy whole lines between line endings rather than filtering byte by byte.
    for line in raw.split(|&b| b == b'\n') {
        sequence.extend_from_slice(line.strip_suffix(b"\r").unwrap_or(line));
    }

    Ok(())
}

#[cfg(unix)]
fn read_exact_at(file: &File, buf: &mut [u8], offset: u64) -> io::Result<()> {
    use std::os::unix::fs::FileExt;
    file.read_exact_at(buf, offset)
}

#[cfg(not(unix))]
fn read_exact_at(file: &File, buf: &mut [u8], offset: u64) -> io::Result<()> {
    use std::io::{Read, Seek, SeekFrom};

    let mut file = file;
    file.seek(SeekFrom::Start(offset))?;
    file.read_exact(buf)
}

fn build_index(src: &Path) -> anyhow::Result<()> {
    let index = fasta::fs::index(src).with_context(|| format!("indexing {}", src.display()))?;

    let mut dst = src.as_os_str().to_owned();
    dst.push(".fai");
    let dst = PathBuf::from(dst);

    let file = File::create(&dst).with_context(|| format!("creating {}", dst.display()))?;
    let mut writer = fai::io::Writer::new(BufWriter::new(file));
    writer.write_index(&index)?;
    writer.get_mut().flush()?;

    Ok(())
}
