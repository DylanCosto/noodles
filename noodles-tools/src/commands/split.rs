//! `ntools split`: write one file per read group or tag value (samtools split).

use std::{
    collections::HashMap,
    path::{Path, PathBuf},
};

use anyhow::{Context, bail};
use bstr::{BString, ByteSlice};
use noodles_sam::{
    self as sam,
    alignment::record::data::field::{Tag, Value},
};

use crate::{
    alignment::{AlignmentReader, ReadOptions},
    output::{self, AlignmentWriter, OutputFormat, WriteOptions},
    program,
};

/// samtools' default output name when splitting by read group: basename, index, extension.
const DEFAULT_FORMAT: &str = "%*_%#.%.";

/// samtools' default output name when splitting by tag: basename, tag value, extension.
const DEFAULT_TAG_FORMAT: &str = "%*_%!.%.";

#[derive(Debug, clap::Args)]
pub struct Args {
    /// Input SAM, BAM, or CRAM file.
    input: String,

    /// Write reads without a (known) read group or tag value to this file.
    #[arg(short = 'u')]
    unaccounted: Option<PathBuf>,

    /// Output name format: %* input basename, %# group index, %! group ID or tag value,
    /// %. extension, %% [default: %*_%#.%., or %*_%!.%. with -d].
    #[arg(short = 'f')]
    format: Option<String>,

    /// Zero-pad %# to this many digits.
    #[arg(short = 'p', default_value_t = 0)]
    pad: usize,

    /// Split by the value of this tag instead of by read group.
    #[arg(short = 'd')]
    tag: Option<String>,

    /// Maximum number of outputs when splitting by tag.
    #[arg(short = 'M', long = "max-split", default_value_t = 100)]
    max_split: usize,

    /// Output format.
    #[arg(
        short = 'O',
        long = "output-fmt",
        ignore_case = true,
        default_value = "bam"
    )]
    output_format: OutputFormat,

    /// Print the name of each output file to stderr.
    #[arg(short = 'v')]
    verbose: bool,

    /// Reference FASTA (indexed) for CRAM input or output.
    #[arg(long = "reference")]
    reference: Option<PathBuf>,

    /// Additional threads for compression and decompression.
    #[arg(short = '@', long = "threads", default_value_t = 0)]
    threads: usize,

    /// Do not add an @PG line to the headers.
    #[arg(long = "no-PG")]
    no_pg: bool,
}

/// One output file.
struct Group {
    writer: AlignmentWriter,
    header: sam::Header,
}

pub fn run(args: Args) -> anyhow::Result<()> {
    let options = ReadOptions {
        threads: args.threads,
        reference: args.reference.clone(),
    };

    let (mut reader, _) = AlignmentReader::open(&args.input, &options)?;
    let mut header = reader.read_header()?;

    if !args.no_pg {
        program::add(&mut header)?;
    }

    let write_options = WriteOptions {
        format: args.output_format,
        level: output::parse_level(6)?,
        threads: args.threads,
        reference: args.reference.clone(),
    };

    let tag = match &args.tag {
        Some(raw) => {
            let tag: [u8; 2] = raw
                .as_bytes()
                .try_into()
                .with_context(|| format!("invalid tag: {raw}"))?;
            Tag::from(tag)
        }
        None => Tag::READ_GROUP,
    };

    let splitting_by_read_group = args.tag.is_none();
    let format = args.format.clone().unwrap_or_else(|| {
        if splitting_by_read_group {
            DEFAULT_FORMAT.to_owned()
        } else {
            DEFAULT_TAG_FORMAT.to_owned()
        }
    });
    let pad = args.pad;
    let basename = basename(&args.input);
    let extension = extension(args.output_format);

    let mut groups: Vec<Group> = Vec::new();
    let mut ids: HashMap<BString, usize> = HashMap::new();

    if splitting_by_read_group {
        for (i, id) in header.read_groups().keys().enumerate() {
            let group_header = read_group_header(&header, id);
            let path = output_path(&format, &basename, i, pad, id.as_ref(), extension)?;
            groups.push(open_group(
                &path,
                group_header,
                &write_options,
                args.verbose,
            )?);
            ids.insert(id.clone(), i);
        }
    }

    let mut unaccounted = match &args.unaccounted {
        Some(path) => {
            let mut writer = AlignmentWriter::create(Some(path), &write_options, &header)?;
            writer.write_header(&header)?;
            Some(writer)
        }
        None => None,
    };

    reader.for_each(&header, |record| {
        let value = match record.as_dyn().data().get(&tag).transpose()? {
            Some(Value::String(s)) => Some(BString::from(s.as_bytes())),
            Some(Value::Character(c)) => Some(BString::from(vec![c])),
            Some(Value::Int8(n)) => Some(n.to_string().into()),
            Some(Value::UInt8(n)) => Some(n.to_string().into()),
            Some(Value::Int16(n)) => Some(n.to_string().into()),
            Some(Value::UInt16(n)) => Some(n.to_string().into()),
            Some(Value::Int32(n)) => Some(n.to_string().into()),
            Some(Value::UInt32(n)) => Some(n.to_string().into()),
            _ => None,
        };

        let index = match value {
            Some(value) => match ids.get(&value) {
                Some(&i) => Some(i),
                None if !splitting_by_read_group => {
                    if groups.len() >= args.max_split {
                        bail!("more than {} tag values; raise -M", args.max_split);
                    }

                    let i = groups.len();
                    let path = output_path(&format, &basename, i, pad, &value, extension)?;
                    groups.push(open_group(
                        &path,
                        header.clone(),
                        &write_options,
                        args.verbose,
                    )?);
                    ids.insert(value, i);
                    Some(i)
                }
                None => None,
            },
            None => None,
        };

        match (index, &mut unaccounted) {
            (Some(i), _) => {
                let group = &mut groups[i];
                group.writer.write(&group.header, &record)?;
            }
            (None, Some(writer)) => writer.write(&header, &record)?,
            // samtools' message and exit status.
            (None, None) => bail!(
                "Read \"{}\" has no {} tag.",
                record
                    .as_dyn()
                    .name()
                    .map_or_else(|| "*".into(), |name| name.to_string()),
                String::from_utf8_lossy(tag.as_ref())
            ),
        }

        Ok(())
    })?;

    for group in groups {
        group.writer.finish(&group.header)?;
    }

    if let Some(writer) = unaccounted {
        writer.finish(&header)?;
    }

    Ok(())
}

fn open_group(
    path: &Path,
    header: sam::Header,
    options: &WriteOptions,
    verbose: bool,
) -> anyhow::Result<Group> {
    if verbose {
        eprintln!("{}", path.display());
    }

    let mut writer = AlignmentWriter::create(Some(path), options, &header)?;
    writer.write_header(&header)?;
    Ok(Group { writer, header })
}

/// The header with only the given read group's `@RG` line.
fn read_group_header(header: &sam::Header, id: &BString) -> sam::Header {
    let mut header = header.clone();
    header.read_groups_mut().retain(|rg_id, _| rg_id == id);
    header
}

/// Expands samtools' output name format.
fn output_path(
    format: &str,
    basename: &str,
    index: usize,
    pad: usize,
    id: &[u8],
    extension: &str,
) -> anyhow::Result<PathBuf> {
    let mut name = String::new();
    let mut chars = format.chars();

    while let Some(c) = chars.next() {
        if c != '%' {
            name.push(c);
            continue;
        }

        match chars.next() {
            Some('%') => name.push('%'),
            Some('*') => name.push_str(basename),
            Some('#') => name.push_str(&format!("{index:0pad$}")),
            Some('!') => name.push_str(&String::from_utf8_lossy(id)),
            Some('.') => name.push_str(extension),
            Some(other) => bail!("invalid format specifier: %{other}"),
            None => bail!("format ends with %"),
        }
    }

    Ok(PathBuf::from(name))
}

/// The input name without directories or its last extension, as `%*` expands.
fn basename(input: &str) -> String {
    let file_name = Path::new(input)
        .file_name()
        .map_or_else(|| input.to_owned(), |s| s.to_string_lossy().into_owned());

    match file_name.rfind('.') {
        Some(i) if i > 0 => file_name[..i].to_owned(),
        _ => file_name,
    }
}

fn extension(format: OutputFormat) -> &'static str {
    match format {
        OutputFormat::Sam => "sam",
        OutputFormat::Bam => "bam",
        OutputFormat::Cram => "cram",
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_output_path() -> anyhow::Result<()> {
        assert_eq!(
            output_path(DEFAULT_FORMAT, "in", 2, 0, b"rg2", "bam")?,
            PathBuf::from("in_2.bam")
        );
        assert_eq!(
            output_path(DEFAULT_FORMAT, "in", 2, 3, b"rg2", "bam")?,
            PathBuf::from("in_002.bam")
        );
        assert_eq!(
            output_path("x/%!-%%.%.", "in", 0, 0, b"A", "sam")?,
            PathBuf::from("x/A-%.sam")
        );
        assert!(output_path("%q", "in", 0, 0, b"A", "sam").is_err());
        Ok(())
    }

    #[test]
    fn test_basename() {
        assert_eq!(basename("/data/sample.sorted.bam"), "sample.sorted");
        assert_eq!(basename("-"), "-");
    }
}
