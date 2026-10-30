//! `ntools dict`: write a SAM sequence dictionary for a FASTA file (samtools dict).

use std::{
    collections::HashSet,
    fs::File,
    io::{self, BufRead, BufReader, BufWriter, Read, Write},
    path::{Path, PathBuf},
};

use anyhow::Context;
use flate2::bufread::MultiGzDecoder;
use md5::{Digest, Md5};

#[derive(Debug, clap::Args)]
#[command(disable_help_flag = true)]
pub struct Args {
    /// FASTA file.
    input: PathBuf,

    /// Assembly identifier (AS).
    #[arg(short = 'a', long = "assembly")]
    assembly: Option<String>,

    /// Species (SP).
    #[arg(short = 's', long = "species")]
    species: Option<String>,

    /// URI (UR) [default: file:// and the input's absolute path].
    #[arg(short = 'u', long = "uri")]
    uri: Option<String>,

    /// Add an AN tag naming the sequence with `chr` added or removed.
    #[arg(short = 'A', long = "alias", visible_alias = "alternative-name")]
    alias: bool,

    /// Add AH:* to the sequences named in this alternate locus (`.alt`) file.
    #[arg(short = 'l', long = "alt")]
    alt: Option<PathBuf>,

    /// Do not print the @HD line.
    #[arg(short = 'H', long = "no-header")]
    no_header: bool,

    /// Write output to this file instead of stdout.
    #[arg(short = 'o', long = "output")]
    output: Option<PathBuf>,

    /// Print help.
    #[arg(long, action = clap::ArgAction::Help)]
    help: Option<bool>,
}

pub fn run(args: Args) -> anyhow::Result<()> {
    let alternates = args.alt.as_deref().map(read_alt_names).transpose()?;

    let uri = match &args.uri {
        Some(uri) => uri.clone(),
        None => {
            // As samtools' realpath(), symbolic links are resolved.
            let path = std::fs::canonicalize(&args.input)
                .with_context(|| format!("resolving {}", args.input.display()))?;
            format!("file://{}", path.display())
        }
    };

    let file =
        File::open(&args.input).with_context(|| format!("opening {}", args.input.display()))?;
    let reader = BufReader::with_capacity(1 << 17, file);

    let sink: Box<dyn Write> = match args.output.as_deref().filter(|p| *p != Path::new("-")) {
        Some(path) => {
            Box::new(File::create(path).with_context(|| format!("creating {}", path.display()))?)
        }
        None => Box::new(io::stdout().lock()),
    };
    let mut out = BufWriter::new(sink);

    if !args.no_header {
        writeln!(out, "@HD\tVN:1.0\tSO:unsorted")?;
    }

    let mut current: Option<Sequence> = None;
    let mut line = Vec::new();
    let mut reader = reader;

    loop {
        line.clear();

        if reader.read_until(b'\n', &mut line)? == 0 {
            break;
        }

        let text = line.strip_suffix(b"\n").unwrap_or(&line);
        let text = text.strip_suffix(b"\r").unwrap_or(text);

        if let Some(definition) = text.strip_prefix(b">") {
            if let Some(sequence) = current.take() {
                sequence.write(&mut out, &uri, &args, alternates.as_ref())?;
            }

            // As in samtools, the name ends at the first whitespace.
            let name_end = definition
                .iter()
                .position(u8::is_ascii_whitespace)
                .unwrap_or(definition.len());

            current = Some(Sequence::new(&definition[..name_end]));
        } else if let Some(sequence) = &mut current {
            sequence.add(text);
        }
    }

    if let Some(sequence) = current {
        sequence.write(&mut out, &uri, &args, alternates.as_ref())?;
    }

    out.flush()?;

    Ok(())
}

/// A sequence being read: its name, length, and the MD5 of its uppercased bases.
struct Sequence {
    name: Vec<u8>,
    length: usize,
    hasher: Md5,
    upper: Vec<u8>,
}

impl Sequence {
    fn new(name: &[u8]) -> Self {
        Self {
            name: name.to_vec(),
            length: 0,
            hasher: Md5::new(),
            upper: Vec::new(),
        }
    }

    fn add(&mut self, line: &[u8]) {
        self.upper.clear();
        self.upper.extend(
            line.iter()
                .filter(|b| !b.is_ascii_whitespace())
                .map(u8::to_ascii_uppercase),
        );

        self.length += self.upper.len();
        self.hasher.update(&self.upper);
    }

    fn write<W: Write>(
        self,
        out: &mut W,
        uri: &str,
        args: &Args,
        alternates: Option<&HashSet<Vec<u8>>>,
    ) -> io::Result<()> {
        let digest = self.hasher.finalize();

        out.write_all(b"@SQ\tSN:")?;
        out.write_all(&self.name)?;
        write!(out, "\tLN:{}\tM5:", self.length)?;

        for byte in digest {
            write!(out, "{byte:02x}")?;
        }

        if alternates.is_some_and(|names| names.contains(&self.name)) {
            out.write_all(b"\tAH:*")?;
        }

        if args.alias {
            write_alias(out, &self.name)?;
        }

        write!(out, "\tUR:{uri}")?;

        if let Some(assembly) = &args.assembly {
            write!(out, "\tAS:{assembly}")?;
        }

        if let Some(species) = &args.species {
            write!(out, "\tSP:{species}")?;
        }

        writeln!(out)
    }
}

/// Writes the AN tag of `-A`: the name without its `chr` prefix, or with one added, and for
/// mitochondria both the `M` and `MT` forms.
fn write_alias<W: Write>(out: &mut W, name: &[u8]) -> io::Result<()> {
    let base = match name.strip_prefix(b"chr") {
        Some(base) => {
            out.write_all(b"\tAN:")?;
            out.write_all(base)?;
            base
        }
        None => {
            out.write_all(b"\tAN:chr")?;
            out.write_all(name)?;
            name
        }
    };

    match base {
        b"M" => out.write_all(b",chrMT,MT"),
        b"MT" => out.write_all(b",chrM,M"),
        _ => Ok(()),
    }
}

/// Reads the sequence names of an alternate locus file: the first field of each line that is
/// not a header (`@`) line. The file may be gzip-compressed.
fn read_alt_names(path: &Path) -> anyhow::Result<HashSet<Vec<u8>>> {
    const GZIP_MAGIC_NUMBER: [u8; 2] = [0x1f, 0x8b];

    let file = File::open(path).with_context(|| format!("Cannot open {}", path.display()))?;
    let mut reader = BufReader::new(file);

    let mut text = Vec::new();

    if reader.fill_buf()?.starts_with(&GZIP_MAGIC_NUMBER) {
        MultiGzDecoder::new(reader).read_to_end(&mut text)?;
    } else {
        reader.read_to_end(&mut text)?;
    }

    Ok(text
        .split(|&b| b == b'\n')
        .map(|line| line.strip_suffix(b"\r").unwrap_or(line))
        .filter(|line| !line.is_empty() && !line.starts_with(b"@"))
        .map(|line| line.split(|&b| b == b'\t').next().unwrap_or(line).to_vec())
        .collect())
}
