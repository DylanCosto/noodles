//! `ntools dict`: write a SAM sequence dictionary for a FASTA file (samtools dict).

use std::{
    fs::File,
    io::{self, BufRead, BufReader, BufWriter, Write},
    path::{Path, PathBuf},
};

use anyhow::Context;
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
    let uri = match &args.uri {
        Some(uri) => uri.clone(),
        None => {
            let path = std::path::absolute(&args.input)
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
                sequence.write(&mut out, &uri, &args)?;
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
        sequence.write(&mut out, &uri, &args)?;
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

    fn write<W: Write>(self, out: &mut W, uri: &str, args: &Args) -> io::Result<()> {
        let digest = self.hasher.finalize();

        out.write_all(b"@SQ\tSN:")?;
        out.write_all(&self.name)?;
        write!(out, "\tLN:{}\tM5:", self.length)?;

        for byte in digest {
            write!(out, "{byte:02x}")?;
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
