//! `ntools samples`: list the samples of alignment files (samtools samples).

use std::{
    fs::File,
    io::{self, BufRead, BufWriter, Write},
    path::{Path, PathBuf},
    process::ExitCode,
};

use anyhow::{Context, bail};
use noodles_fasta as fasta;
use noodles_sam as sam;

use crate::{
    alignment::{AlignmentIndex, AlignmentReader, ReadOptions},
    commands::faidx,
    khash::StringSet,
};

/// Whether an input has an index (`-i`).
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
enum IndexStatus {
    Found,
    Missing,
}

#[derive(Debug, clap::Args)]
#[command(disable_help_flag = true)]
pub struct Args {
    /// Alignment files [default: paths read from stdin, one per line].
    inputs: Vec<String>,

    /// Print a header line before the results.
    #[arg(short = 'h')]
    header: bool,

    /// Report whether each file is indexed.
    #[arg(short = 'i')]
    test_index: bool,

    /// Read group tag holding the sample.
    #[arg(short = 'T', default_value = "SM")]
    tag: String,

    /// Output file [default: stdout].
    #[arg(short = 'o')]
    output: Option<PathBuf>,

    /// An indexed FASTA file to match files' sequences against; repeatable.
    #[arg(short = 'f')]
    references: Vec<PathBuf>,

    /// A file listing indexed FASTA files, one per line.
    #[arg(short = 'F')]
    reference_list: Option<PathBuf>,

    /// Inputs are followed by their index files (or, on stdin, are `path<TAB>index` lines).
    #[arg(short = 'X')]
    index_files: bool,

    /// Print help.
    #[arg(short = '?', long = "help", action = clap::ArgAction::Help)]
    help: Option<bool>,
}

/// An indexed FASTA file's sequence names and lengths.
struct Reference {
    path: String,
    sequences: Vec<(Vec<u8>, u64)>,
}

pub fn run(args: Args) -> anyhow::Result<ExitCode> {
    let tag: [u8; 2] = args
        .tag
        .as_bytes()
        .try_into()
        .map_err(|_| anyhow::anyhow!("[samples] Length of tag \"{}\" is not 2.", args.tag))?;

    // samtools searches the references given last first.
    let mut references = Vec::new();

    for path in &args.references {
        references.push(load_reference(&path.to_string_lossy())?);
    }

    if let Some(list) = &args.reference_list {
        let file = File::open(list)
            .with_context(|| format!("[samples] Cannot open \"{}\"", list.display()))?;

        for line in io::BufReader::new(file).lines() {
            references.push(load_reference(&line?)?);
        }
    }

    references.reverse();

    let sink: Box<dyn Write> = match &args.output {
        Some(path) => Box::new(File::create(path).with_context(|| {
            format!("[samples] Cannot open \"{}\" for writing", path.display())
        })?),
        None => Box::new(io::stdout().lock()),
    };
    let mut out = BufWriter::new(sink);

    if args.header {
        out.write_all(b"#")?;
        out.write_all(&tag)?;
        out.write_all(b"\tPATH")?;

        if args.test_index {
            out.write_all(b"\tINDEX")?;
        }

        if !references.is_empty() {
            out.write_all(b"\tREFERENCE")?;
        }

        out.write_all(b"\n")?;
    }

    let lister = Lister {
        tag,
        test_index: args.test_index,
        references: &references,
    };

    let mut inputs: Vec<(String, Option<String>)> = Vec::new();

    if args.inputs.is_empty() {
        for line in io::stdin().lock().lines() {
            let line = line?;

            if args.index_files {
                match line.split_once('\t') {
                    Some((path, index)) if !index.is_empty() => {
                        inputs.push((path.into(), Some(index.into())));
                    }
                    _ => {
                        bail!("[samples] Expected path-to-bam(tab)path-to-index but got \"{line}\"")
                    }
                }
            } else {
                inputs.push((line, None));
            }
        }
    } else if args.index_files {
        if !args.inputs.len().is_multiple_of(2) {
            bail!(
                "[samples] Odd number of filenames detected! Each BAM file should have an index file"
            );
        }

        let (paths, indexes) = args.inputs.split_at(args.inputs.len() / 2);
        inputs.extend(
            paths
                .iter()
                .zip(indexes)
                .map(|(path, index)| (path.clone(), Some(index.clone()))),
        );
    } else {
        inputs.extend(args.inputs.iter().map(|path| (path.clone(), None)));
    }

    for (path, index) in &inputs {
        lister.print_samples(&mut out, path, index.as_deref())?;
    }

    out.flush()?;

    Ok(ExitCode::SUCCESS)
}

struct Lister<'a> {
    tag: [u8; 2],
    test_index: bool,
    references: &'a [Reference],
}

impl Lister<'_> {
    /// `print_samples`: one line per distinct sample in a file's read groups (`.` if none), in
    /// samtools' (hash) order.
    fn print_samples<W: Write>(
        &self,
        out: &mut W,
        path: &str,
        index: Option<&str>,
    ) -> anyhow::Result<()> {
        let options = ReadOptions::default();
        let (mut reader, local) = AlignmentReader::open(path, &options)
            .with_context(|| format!("[samples] Failed to open \"{path}\" for reading"))?;
        let header = reader
            .read_header()
            .with_context(|| format!("[samples] Failed to read the header from \"{path}\""))?;

        let index_status = self.test_index.then(|| {
            let found = local.as_deref().is_some_and(|local| {
                AlignmentIndex::read(reader.format(), local, index.map(Path::new)).is_ok()
            });

            if found {
                IndexStatus::Found
            } else {
                IndexStatus::Missing
            }
        });

        let mut samples = StringSet::new();

        for (id, read_group) in header.read_groups() {
            let value = if self.tag == *b"ID" {
                Some(id.as_slice())
            } else {
                read_group
                    .other_fields()
                    .get(&self.tag)
                    .map(|v| v.as_slice())
            };

            if let Some(value) = value {
                samples.insert(value);
            }
        }

        let reference = self.reference_for(&header);

        if samples.is_empty() {
            self.print_sample(out, b".", path, index_status, reference)?;
        } else {
            for sample in samples.iter() {
                self.print_sample(out, sample, path, index_status, reference)?;
            }
        }

        Ok(())
    }

    fn print_sample<W: Write>(
        &self,
        out: &mut W,
        sample: &[u8],
        path: &str,
        index_status: Option<IndexStatus>,
        reference: Option<&Reference>,
    ) -> io::Result<()> {
        out.write_all(sample)?;
        write!(out, "\t{path}")?;

        match index_status {
            Some(IndexStatus::Found) => out.write_all(b"\tY")?,
            Some(IndexStatus::Missing) => out.write_all(b"\tN")?,
            None => {}
        }

        if !self.references.is_empty() {
            write!(out, "\t{}", reference.map_or(".", |r| r.path.as_str()))?;
        }

        writeln!(out)
    }

    /// The first reference with the header's sequence names and lengths, in order.
    fn reference_for(&self, header: &sam::Header) -> Option<&Reference> {
        let sequences = header.reference_sequences();

        self.references.iter().find(|reference| {
            reference.sequences.len() == sequences.len()
                && reference
                    .sequences
                    .iter()
                    .zip(sequences)
                    .all(|((name, len), (n, rs))| {
                        name.as_slice() == n.as_slice() && *len == rs.length().get() as u64
                    })
        })
    }
}

/// `fai_load`: reads a FASTA file's index, building it if it is missing.
fn load_reference(path: &str) -> anyhow::Result<Reference> {
    let src = Path::new(path);

    let index = match fasta::fs::read_associated_index(src) {
        Ok(index) => index,
        Err(_) => {
            faidx::build_index(src)
                .with_context(|| format!("[samples] Cannot load index from \"{path}\""))?;
            fasta::fs::read_associated_index(src)
                .with_context(|| format!("[samples] Cannot load index from \"{path}\""))?
        }
    };

    Ok(Reference {
        path: path.into(),
        sequences: index
            .as_ref()
            .iter()
            .map(|record| (record.name().to_vec(), record.length()))
            .collect(),
    })
}
