//! samtools' shared long options (`sam_opts.c`), for the commands whose samtools counterparts
//! accept them.

use std::path::{Path, PathBuf};

use anyhow::{Context, bail};
use noodles_cram::io::reader::MdNm;

use crate::{
    alignment::{AlignmentReader, ReadOptions},
    commands::index::{IndexKind, index_file},
    output::OutputFormat,
};

/// `--verbosity`, HTSlib's logging level. It changes only diagnostics, so it is accepted and
/// ignored.
#[derive(Debug, clap::Args)]
pub struct Verbosity {
    /// Ignored (HTSlib's logging level).
    #[arg(long = "verbosity", value_name = "INT", allow_negative_numbers = true)]
    verbosity: Option<i32>,
}

/// `--write-index` for commands that, like their samtools counterparts, write no output to
/// index: samtools accepts and ignores it.
#[derive(Debug, clap::Args)]
pub struct IgnoredWriteIndex {
    /// Ignored, as in samtools: there is no alignment output to index.
    #[arg(long = "write-index")]
    write_index: bool,
}

/// `--input-fmt-option OPT[=VAL]`, HTSlib's input format options.
#[derive(Debug, Default, clap::Args)]
pub struct InputFormatOptions {
    /// Input format option: `decode_md=0` keeps CRAM MD/NM tags as stored, `reference=FILE` is
    /// --reference, and `nthreads=N` is ignored.
    #[arg(long = "input-fmt-option", value_name = "OPT[=VAL]")]
    input_format_options: Vec<String>,
}

impl InputFormatOptions {
    /// Options given by a command's own argument.
    pub fn new(input_format_options: Vec<String>) -> Self {
        Self {
            input_format_options,
        }
    }

    /// These options, after the ones a `--input-fmt FORMAT,OPT[=VAL]...` value adds; the format
    /// itself is detected, so its name is ignored.
    pub fn with_input_format(&self, input_format: Option<&str>) -> Self {
        let mut input_format_options: Vec<String> = input_format
            .map(|s| s.split(',').skip(1).map(str::to_owned).collect())
            .unwrap_or_default();

        input_format_options.extend(self.input_format_options.iter().cloned());

        Self {
            input_format_options,
        }
    }

    /// Checks the options, for commands that read no alignments they would change.
    pub fn validate(&self) -> anyhow::Result<()> {
        self.parse().map(|_| ())
    }

    /// Takes `reference=FILE` as the reference when `--reference` was not given.
    pub fn configure(&self, read_options: &mut ReadOptions) -> anyhow::Result<()> {
        for option in self.parse()? {
            if let InputOption::Reference(path) = option
                && read_options.reference.is_none()
            {
                read_options.reference = Some(path);
            }
        }

        Ok(())
    }

    /// Applies the options that change how an open reader decodes (`decode_md`).
    pub fn apply(&self, reader: &mut AlignmentReader) -> anyhow::Result<()> {
        for option in self.parse()? {
            if let InputOption::DecodeMd(md_nm) = option {
                reader.set_md_nm(md_nm);
            }
        }

        Ok(())
    }

    fn parse(&self) -> anyhow::Result<Vec<InputOption>> {
        self.input_format_options
            .iter()
            .map(|s| InputOption::parse(s))
            .collect()
    }
}

/// A supported HTSlib input format option.
enum InputOption {
    DecodeMd(MdNm),
    Reference(PathBuf),
    /// A performance option, ignored.
    Ignored,
}

impl InputOption {
    fn parse(s: &str) -> anyhow::Result<Self> {
        let (name, value) = s.split_once('=').unwrap_or((s, ""));

        match name {
            "decode_md" => {
                // HTSlib reads a missing value as 1 and the value with `strtol`.
                let n: i64 = if value.is_empty() {
                    1
                } else {
                    value
                        .parse()
                        .with_context(|| format!("invalid input format option: {s}"))?
                };

                Ok(Self::DecodeMd(if n == 0 {
                    MdNm::AsStored
                } else {
                    MdNm::Regenerate
                }))
            }
            "reference" if !value.is_empty() => Ok(Self::Reference(PathBuf::from(value))),
            "nthreads" => Ok(Self::Ignored),
            _ => bail!("unsupported input format option: {s}"),
        }
    }
}

/// `--output-fmt-option OPT[=VAL]`, HTSlib's output format options.
#[derive(Debug, Default, clap::Args)]
pub struct OutputFormatOptions {
    /// Output format option: `level=N` (BGZF compression level); `nthreads=N` is ignored.
    #[arg(long = "output-fmt-option", value_name = "OPT=VAL")]
    output_format_options: Vec<String>,
}

impl OutputFormatOptions {
    /// The compression level from `level=N`, or `default`.
    pub fn level(&self, default: u8) -> anyhow::Result<u8> {
        let mut level = default;

        for option in &self.output_format_options {
            match option.split_once('=') {
                Some(("level", value)) => {
                    level = value
                        .parse()
                        .with_context(|| format!("invalid compression level: {value}"))?;
                }
                Some(("nthreads", _)) => {}
                _ => bail!("unsupported output format option: {option}"),
            }
        }

        Ok(level)
    }
}

/// Indexes an alignment output file as samtools' `--write-index` does (`auto_index`): CSI for
/// BAM, CRAI for CRAM.
pub fn write_index(
    output: Option<&Path>,
    format: OutputFormat,
    threads: usize,
) -> anyhow::Result<()> {
    let Some(output) = output.filter(|p| *p != Path::new("-")) else {
        bail!("--write-index needs an output file");
    };

    let kind = match format {
        OutputFormat::Bam => IndexKind::Csi,
        OutputFormat::Cram => IndexKind::Default,
        OutputFormat::Sam => bail!("--write-index needs BAM or CRAM output"),
    };

    index_file(output, None, kind, threads)
}
