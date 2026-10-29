//! `ntools faidx`: index a FASTA file or extract regions from it (samtools faidx).

use std::{
    fs::{self, File},
    io::{self, BufWriter, Write},
    num::NonZero,
    path::{Path, PathBuf},
};

use anyhow::Context;
use noodles_core::Region;
use noodles_fasta::{self as fasta, fai};

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

    let mut reader = fasta::io::indexed_reader::Builder::default()
        .build_from_path(&args.input)
        .with_context(|| format!("opening {} (is it indexed?)", args.input.display()))?;

    let sink: Box<dyn Write> = match &args.output {
        Some(path) => {
            Box::new(File::create(path).with_context(|| format!("creating {}", path.display()))?)
        }
        None => Box::new(io::stdout().lock()),
    };
    let mut writer = BufWriter::new(sink);

    let mut sequence = Vec::new();

    for raw_region in &regions {
        let region: Region = raw_region
            .parse()
            .with_context(|| format!("invalid region: {raw_region}"))?;

        let record = reader
            .query(&region)
            .with_context(|| format!("retrieving {raw_region}"))?;

        sequence.clear();
        sequence.extend_from_slice(record.sequence().as_ref());

        if args.reverse_complement {
            reverse_complement(&mut sequence);
            writeln!(writer, ">{raw_region}/rc")?;
        } else {
            writeln!(writer, ">{raw_region}")?;
        }

        let index_record = index_record(reader.index(), &region)?;

        if region.interval().end().is_some_and(|end| {
            u64::try_from(usize::from(end)).is_ok_and(|end| end > index_record.length())
        }) {
            eprintln!("[faidx] Truncated sequence: {raw_region}");
        }

        let line_width = match args.line_width {
            Some(n) => n.get(),
            None => usize::try_from(index_record.line_base_count().get())
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

fn index_record<'i>(index: &'i fai::Index, region: &Region) -> anyhow::Result<&'i fai::Record> {
    index
        .as_ref()
        .iter()
        .find(|record| record.name() == region.name())
        .with_context(|| format!("{} is not in the index", region.name()))
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

/// Reverse complements IUPAC nucleotides in place, preserving case, as samtools does.
fn reverse_complement(sequence: &mut [u8]) {
    sequence.reverse();

    for base in sequence {
        *base = complement(*base);
    }
}

fn complement(base: u8) -> u8 {
    let c = match base.to_ascii_uppercase() {
        b'A' => b'T',
        b'C' => b'G',
        b'G' => b'C',
        b'T' | b'U' => b'A',
        b'R' => b'Y',
        b'Y' => b'R',
        b'K' => b'M',
        b'M' => b'K',
        b'B' => b'V',
        b'V' => b'B',
        b'D' => b'H',
        b'H' => b'D',
        other => other,
    };

    if base.is_ascii_lowercase() {
        c.to_ascii_lowercase()
    } else {
        c
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_reverse_complement() {
        let mut s = b"ACGTNacgtnRY".to_vec();
        reverse_complement(&mut s);
        assert_eq!(s, b"RYnacgtNACGT");
    }
}
