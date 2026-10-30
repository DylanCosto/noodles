//! `ntools quickcheck`: check that alignment files look intact (samtools quickcheck).

use std::{
    fs::File,
    io::{self, Read, Seek, SeekFrom},
    path::{Path, PathBuf},
    process::ExitCode,
};

use crate::{
    alignment::{AlignmentReader, ReadOptions},
    bam_blocks,
    cram_file::{EOF_2, EOF_3},
    source::{Compression, Format, Input},
};

#[derive(Debug, clap::Args)]
pub struct Args {
    /// Files to check.
    #[arg(required = true)]
    inputs: Vec<PathBuf>,

    /// Print the names of files that fail to stdout; repeat for progress on stderr.
    #[arg(short = 'v', action = clap::ArgAction::Count)]
    verbose: u8,

    /// Do not print messages about failures to stderr (unless -vv).
    #[arg(short = 'q')]
    quiet: bool,

    /// Do not require @SQ lines (for unmapped files).
    #[arg(short = 'u')]
    unmapped: bool,
}

/// The ways a file can fail, as samtools reports them: the exit status is the bitwise OR of the
/// failures of all files.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
enum Failure {
    Open,
    NotSequence,
    BadHeader,
    NoEofBlock,
    Close,
}

impl Failure {
    fn bit(self) -> u8 {
        match self {
            Self::Open => 2,
            Self::NotSequence => 4,
            Self::BadHeader => 8,
            Self::NoEofBlock => 16,
            Self::Close => 32,
        }
    }
}

/// Reports failures and progress for one file.
struct Report<'a> {
    args: &'a Args,
    name: String,
    state: u8,
}

impl Report<'_> {
    /// Records a failure, with samtools' message unless `-q` (and not `-vv`).
    fn fail(&mut self, failure: Failure, message: &str) {
        self.state |= failure.bit();

        if !self.args.quiet || self.args.verbose >= 2 {
            eprintln!("{} {message}", self.name);
        }
    }

    /// Prints a progress message with `-vvv`.
    fn progress(&self, message: &str) {
        if self.args.verbose >= 3 {
            eprintln!("{} {message}", self.name);
        }
    }
}

pub fn run(args: Args) -> anyhow::Result<ExitCode> {
    if args.verbose >= 2 {
        eprintln!("verbosity set to {}", args.verbose);
    }

    let mut status = 0;

    for path in &args.inputs {
        let mut report = Report {
            args: &args,
            name: path.display().to_string(),
            state: 0,
        };

        if args.verbose >= 3 {
            eprintln!("checking {}", report.name);
        }

        check(path, &mut report);

        if report.state != 0 && args.verbose >= 1 {
            println!("{}", report.name);
        }

        status |= report.state;
    }

    Ok(ExitCode::from(status))
}

/// `main_quickcheck`'s checks of one file: that it opens, is alignment data, has a header with
/// reference sequences (unless `-u`), and ends with an EOF block where its format has one.
fn check(path: &Path, report: &mut Report<'_>) {
    let Some(src) = path.to_str() else {
        report.fail(Failure::Open, "could not be opened for reading.");
        return;
    };

    let Ok(input) = Input::open(src) else {
        report.fail(Failure::Open, "could not be opened for reading.");
        return;
    };

    if report.args.verbose >= 3 {
        eprintln!("opened {}", report.name);
    }

    if is_empty(path) {
        // HTSlib does not recognize an empty file as any format.
        report.fail(Failure::NotSequence, "was not identified as sequence data.");
        return;
    }

    report.progress("is sequence data");

    let (format, compression) = (input.format, input.compression);

    match read_reference_sequence_count(input) {
        Ok(0) if !report.args.unmapped => {
            report.fail(Failure::BadHeader, "had no targets in header.");
        }
        Ok(n) => report.progress(&format!("has {n} targets in header.")),
        Err(e) => {
            report.fail(
                Failure::BadHeader,
                "caused an error whilst reading its header.",
            );

            // HTSlib's BGZF reader keeps the error of a truncated block, so closing the file
            // fails too.
            if compression == Compression::Bgzf && is_truncation(&e) {
                report.state |= Failure::Close.bit();
            }
        }
    }

    match eof_state(path, format, compression) {
        Ok(EofState::Present) => report.progress("has good EOF block."),
        Ok(EofState::Missing) => report.fail(
            Failure::NoEofBlock,
            "was missing EOF block when one should be present.",
        ),
        Ok(EofState::NotApplicable) => report
            .progress("cannot be checked for EOF block because its filetype does not contain one."),
        Err(_) => report.fail(
            Failure::NoEofBlock,
            "caused an error whilst checking for EOF block.",
        ),
    }
}

fn is_empty(path: &Path) -> bool {
    path.metadata().is_ok_and(|metadata| metadata.len() == 0)
}

fn read_reference_sequence_count(input: Input) -> anyhow::Result<usize> {
    let mut reader = AlignmentReader::from_input(input, &ReadOptions::default())?;
    let header = reader.read_header()?;
    Ok(header.reference_sequences().len())
}

/// Whether an error comes from data that ends early.
fn is_truncation(e: &anyhow::Error) -> bool {
    e.chain()
        .filter_map(|cause| cause.downcast_ref::<io::Error>())
        .any(|e| e.kind() == io::ErrorKind::UnexpectedEof)
}

#[derive(Debug, Eq, PartialEq)]
enum EofState {
    Present,
    Missing,
    /// The format has no EOF marker (uncompressed SAM, CRAM before 2.1).
    NotApplicable,
}

/// `hts_check_EOF`
fn eof_state(path: &Path, format: Format, compression: Compression) -> io::Result<EofState> {
    let mut file = File::open(path)?;
    let len = file.metadata()?.len();

    let present = match (format, compression) {
        (Format::Cram, _) => {
            let mut definition = [0; 6];
            file.read_exact(&mut definition)?;

            let eof: &[u8] = match (definition[4], definition[5]) {
                (3.., _) => &EOF_3,
                (2, 1..) => &EOF_2,
                _ => return Ok(EofState::NotApplicable),
            };

            ends_with_eof_container(&mut file, len, eof)?
        }
        (_, Compression::Bgzf) => bam_blocks::ends_with_eof(&mut file, len)?,
        (_, Compression::None) => return Ok(EofState::NotApplicable),
    };

    Ok(if present {
        EofState::Present
    } else {
        EofState::Missing
    })
}

/// `cram_check_EOF`: whether the file ends with the EOF container `eof`.
fn ends_with_eof_container(file: &mut File, len: u64, eof: &[u8]) -> io::Result<bool> {
    let eof_len = eof.len() as u64;

    if len < eof_len {
        return Ok(false);
    }

    let mut buf = vec![0; eof.len()];
    file.seek(SeekFrom::Start(len - eof_len))?;
    file.read_exact(&mut buf)?;

    Ok(matches_eof_container(&buf, eof))
}

/// Like HTSlib, compares only the low 4 bits of byte 8 (the last ITF-8 byte of the reference
/// sequence ID), which early Java and C implementations wrote differently.
fn matches_eof_container(buf: &[u8], eof: &[u8]) -> bool {
    const MASKED: usize = 8;

    buf.len() == eof.len()
        && buf.iter().zip(eof).enumerate().all(|(i, (&a, &b))| {
            if i == MASKED {
                a & 0x0f == b & 0x0f
            } else {
                a == b
            }
        })
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_failure_bits() {
        let all = [
            Failure::Open,
            Failure::NotSequence,
            Failure::BadHeader,
            Failure::NoEofBlock,
            Failure::Close,
        ]
        .iter()
        .fold(0, |bits, failure| bits | failure.bit());

        assert_eq!(all, 2 | 4 | 8 | 16 | 32);
    }

    #[test]
    fn test_matches_eof_container() {
        // HTSlib writes the CRAM 2.1 EOF container's reference sequence ID as ff ff ff ff 0f.
        let mut htslib_eof = EOF_2;
        htslib_eof[8] = 0x0f;
        assert!(matches_eof_container(&htslib_eof, &EOF_2));
        assert!(matches_eof_container(&EOF_3, &EOF_3));

        let mut other = EOF_3;
        other[20] ^= 1;
        assert!(!matches_eof_container(&other, &EOF_3));
        assert!(!matches_eof_container(&EOF_3[1..], &EOF_3));
    }
}
