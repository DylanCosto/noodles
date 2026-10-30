//! `ntools view`: convert, filter, count, and query alignments (samtools view).

use std::{
    io::{self, Write},
    path::{Path, PathBuf},
};

use anyhow::{Context, bail};
use noodles_bam as bam;
use noodles_cram as cram;
use noodles_fasta as fasta;
use noodles_sam::{
    self as sam,
    alignment::{io::Write as _, record::Flags},
};

use crate::{
    alignment::{AlignmentReader, Queries, ReadOptions, Rec, encode_cram_record},
    bam_raw,
    bam_sam::{self, Names},
    commands::index::{IndexKind, index_file},
    filter::{self, Filter, Subsample, TagFilter},
    flags,
    output::{self, AlignmentWriter, OutputFormat, WriteOptions},
    parallel_cram, program,
    regions::RegionSet,
    source::Format,
};

#[derive(Debug, clap::Args)]
#[command(disable_help_flag = true)]
pub struct Args {
    /// Input SAM, BAM, or CRAM file (`-` for stdin).
    input: String,

    /// Regions to query (`chr1`, `chr1:100-200`, or `*` for unplaced unmapped reads); requires an
    /// index.
    regions: Vec<String>,

    /// Write output to this file instead of stdout.
    #[arg(short = 'o', long = "output")]
    output: Option<PathBuf>,

    /// Output format.
    #[arg(short = 'O', long = "output-fmt", ignore_case = true)]
    output_format: Option<OutputFormat>,

    /// Output BAM.
    #[arg(short = 'b', long = "bam")]
    bam: bool,

    /// Output CRAM (requires -T).
    #[arg(short = 'C', long = "cram")]
    cram: bool,

    /// Output uncompressed BAM.
    #[arg(short = 'u', long = "uncompressed")]
    uncompressed: bool,

    /// Use fast BAM compression (level 1).
    #[arg(short = '1', long = "fast")]
    fast: bool,

    /// BGZF compression level for BAM output.
    #[arg(short = 'l', default_value_t = 6)]
    level: u8,

    /// Include the header in SAM output.
    #[arg(short = 'h', long = "with-header")]
    with_header: bool,

    /// Print the header only.
    #[arg(short = 'H', long = "header-only")]
    header_only: bool,

    /// Print only the count of matching records.
    #[arg(short = 'c', long = "count")]
    count: bool,

    /// Write records that fail the filters to this file.
    #[arg(short = 'U', long = "unoutput")]
    unselected_output: Option<PathBuf>,

    /// Only include reads overlapping this BED file.
    #[arg(short = 'L', long = "target-file")]
    bed: Option<PathBuf>,

    /// Merge overlapping regions and output each record once, in file order.
    #[arg(short = 'M', long = "use-index")]
    merge_regions: bool,

    /// Index file to use for region queries.
    #[arg(short = 'X', long = "customized-index")]
    index: Option<PathBuf>,

    /// Only include reads with all of these flags set.
    #[arg(short = 'f', long = "require-flags", value_parser = flags::parse)]
    require_flags: Option<Flags>,

    /// Exclude reads with any of these flags set.
    #[arg(short = 'F', long = "excl-flags", value_parser = flags::parse)]
    exclude_flags: Option<Flags>,

    /// Exclude reads with all of these flags set.
    #[arg(short = 'G', value_parser = flags::parse)]
    exclude_all_flags: Option<Flags>,

    /// Exclude reads with mapping quality below this value.
    #[arg(short = 'q', long = "min-MQ", default_value_t = 0)]
    min_mapping_quality: u8,

    /// Subsample: keep a FRAC fraction of templates, seeded by INT (`INT.FRAC`).
    #[arg(short = 's', long = "subsample", value_parser = Subsample::parse)]
    subsample: Option<Subsample>,

    /// Only include reads whose names are listed in this file.
    #[arg(short = 'N', long = "qname-file")]
    qname_file: Option<PathBuf>,

    /// Only include reads with this tag, and optionally this value (`TAG[:VALUE]`).
    #[arg(short = 'd', long = "tag")]
    tags: Vec<String>,

    /// Remove these tags from output records (comma-separated; may be repeated).
    #[arg(short = 'x', long = "remove-tag", value_delimiter = ',')]
    remove_tags: Vec<String>,

    /// Index the output file (CSI for BAM, CRAI for CRAM).
    #[arg(long = "write-index")]
    write_index: bool,

    /// Reference FASTA (indexed) for CRAM input or output.
    #[arg(short = 'T', long = "reference")]
    reference: Option<PathBuf>,

    /// Additional threads for decoding, formatting, and compression.
    #[arg(short = '@', long = "threads", default_value_t = 0)]
    threads: usize,

    /// Do not add an @PG line to the header.
    #[arg(long = "no-PG")]
    no_pg: bool,

    /// Print help.
    #[arg(long, action = clap::ArgAction::Help)]
    help: Option<bool>,
}

/// How records are reached.
enum Traversal {
    /// Every record, in file order.
    All,
    /// samtools regions, one indexed query each.
    Regions(Vec<String>),
    /// Merged intervals through the index, each record once.
    Set(RegionSet),
}

pub fn run(args: Args) -> anyhow::Result<()> {
    let format = output_format(&args)?;

    let read_options = ReadOptions {
        threads: args.threads,
        reference: args.reference.clone(),
    };

    let (mut reader, path) = AlignmentReader::open(&args.input, &read_options)?;
    let mut header = reader.read_header()?;

    let (traversal, mut filter) = plan(&args, &header, &reader, path.as_deref())?;
    filter.require = args.require_flags.unwrap_or_else(Flags::empty);
    filter.exclude = args.exclude_flags.unwrap_or_else(Flags::empty);
    filter.exclude_all = args.exclude_all_flags;
    filter.min_mapping_quality = args.min_mapping_quality;
    filter.subsample = args.subsample;
    filter.tag = TagFilter::from_args(&args.tags)?;
    filter.names = args
        .qname_file
        .as_deref()
        .map(filter::read_names)
        .transpose()?;
    let filter = filter;

    let remove_tags = parse_tags(&args.remove_tags)?;
    let decoding = if matches!(traversal, Traversal::All)
        && args.threads > 0
        && remove_tags.is_empty()
        && args.unselected_output.is_none()
    {
        Decoding::Parallel
    } else {
        Decoding::Serial
    };

    if args.count {
        let n = count(
            &mut reader,
            path.as_deref(),
            &args,
            &header,
            &traversal,
            &filter,
            decoding,
        )?;
        writeln!(io::stdout().lock(), "{n}")?;
        return Ok(());
    }

    if !args.no_pg {
        program::add(&mut header)?;
    }

    let level = if args.uncompressed {
        0
    } else if args.fast {
        1
    } else {
        args.level
    };

    let write_options = WriteOptions {
        format,
        level: output::parse_level(level)?,
        threads: args.threads,
        reference: args.reference.clone(),
    };

    let mut writer = AlignmentWriter::create(args.output.as_deref(), &write_options, &header)?;
    let mut unselected = args
        .unselected_output
        .as_deref()
        .map(|dst| AlignmentWriter::create(Some(dst), &write_options, &header))
        .transpose()?;

    let writes_header = format != OutputFormat::Sam || args.with_header || args.header_only;

    if writes_header {
        writer.write_header(&header)?;
    }

    if let Some(unselected) = &mut unselected
        && writes_header
    {
        unselected.write_header(&header)?;
    }

    if !args.header_only {
        let parallel_kind = writer.preformatted_kind();

        if let (AlignmentReader::Cram(cram_reader, repository, _), Decoding::Parallel, Some(kind)) =
            (&mut reader, decoding, parallel_kind)
        {
            write_cram_in_parallel(
                cram_reader,
                repository,
                &header,
                &filter,
                kind,
                args.threads,
                &mut writer,
            )?;
        } else {
            let mut stripped = Vec::new();
            let mut encoder = bam::io::Writer::from(Vec::new());

            visit(
                &mut reader,
                path.as_deref(),
                args.index.as_deref(),
                &header,
                &traversal,
                |rec| {
                    if !filter.matches(&rec, &header)? {
                        if let Some(unselected) = &mut unselected {
                            unselected.write(&header, &rec)?;
                        }
                        return Ok(());
                    }

                    if remove_tags.is_empty() {
                        writer.write(&header, &rec)?;
                    } else {
                        let raw = raw_bytes(&rec, &header, &mut encoder)?;
                        bam_raw::strip_fields(raw, &mut stripped, |tag| {
                            remove_tags.contains(&tag)
                        })?;
                        writer.write_raw(&header, &stripped)?;
                    }

                    Ok(())
                },
            )?;
        }
    }

    writer.finish(&header)?;

    if let Some(unselected) = unselected {
        unselected.finish(&header)?;
    }

    if args.write_index {
        let Some(output) = args.output.as_deref().filter(|p| *p != Path::new("-")) else {
            bail!("--write-index needs an output file (-o)");
        };

        let kind = match format {
            OutputFormat::Bam => IndexKind::Csi,
            OutputFormat::Cram => IndexKind::Default,
            OutputFormat::Sam => bail!("--write-index needs BAM or CRAM output"),
        };

        index_file(output, None, kind, args.threads)?;
    }

    Ok(())
}

/// Chooses the traversal and any region filter: `-M` merges the regions given; `-L` restricts to
/// a BED file, through the index when one exists and no regions are given.
fn plan(
    args: &Args,
    header: &sam::Header,
    reader: &AlignmentReader,
    path: Option<&Path>,
) -> anyhow::Result<(Traversal, Filter)> {
    let mut filter = Filter::default();

    let bed = match &args.bed {
        Some(bed) => {
            let mut set = RegionSet::new(header);
            set.add_bed(header, bed)?;
            set.merge();
            Some(set)
        }
        None => None,
    };

    if !args.regions.is_empty() {
        let traversal = if args.merge_regions {
            if args.regions.iter().any(|r| r == "*") {
                bail!("-M cannot be combined with the unmapped region *");
            }

            let mut set = RegionSet::new(header);
            for region in &args.regions {
                set.add_region(header, region)?;
            }
            set.merge();
            Traversal::Set(set)
        } else {
            Traversal::Regions(args.regions.clone())
        };

        filter.regions = bed;
        return Ok((traversal, filter));
    }

    let Some(bed) = bed else {
        return Ok((Traversal::All, filter));
    };

    let indexed = path.is_some()
        && matches!(reader.format(), Format::Bam | Format::Cram)
        && (args.index.is_some() || path.is_some_and(has_index));

    if indexed {
        Ok((Traversal::Set(bed), filter))
    } else {
        filter.regions = Some(bed);
        Ok((Traversal::All, filter))
    }
}

fn has_index(path: &Path) -> bool {
    ["bai", "csi", "crai"].iter().any(|ext| {
        let mut s = path.as_os_str().to_owned();
        s.push(".");
        s.push(ext);
        Path::new(&s).exists()
    }) || path.with_extension("bai").exists()
}

/// Whether CRAM containers are decoded on worker threads (whole-file reads with `-@`).
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
enum Decoding {
    Serial,
    Parallel,
}

fn count(
    reader: &mut AlignmentReader,
    path: Option<&Path>,
    args: &Args,
    header: &sam::Header,
    traversal: &Traversal,
    filter: &Filter,
    decoding: Decoding,
) -> anyhow::Result<u64> {
    let mut n = 0u64;

    if let (AlignmentReader::Cram(cram_reader, repository, _), Decoding::Parallel) =
        (&mut *reader, decoding)
    {
        parallel_cram::decode(
            cram_reader,
            header,
            repository,
            args.threads,
            |count: &mut u64, record| {
                if filter.matches(&Rec::Dyn(record), header)? {
                    *count += 1;
                }
                Ok(())
            },
            |count| {
                n += count;
                Ok(())
            },
        )?;

        return Ok(n);
    }

    visit(
        reader,
        path,
        args.index.as_deref(),
        header,
        traversal,
        |rec| {
            if filter.matches(&rec, header)? {
                n += 1;
            }
            Ok(())
        },
    )?;

    Ok(n)
}

/// Returns a record's raw BAM encoding, encoding non-BAM records with `encoder`.
fn raw_bytes<'a>(
    rec: &'a Rec<'_>,
    header: &sam::Header,
    encoder: &'a mut bam::io::Writer<Vec<u8>>,
) -> io::Result<&'a [u8]> {
    match rec {
        Rec::Bam(record) => Ok(record.as_ref()),
        _ => {
            encoder.get_mut().clear();
            encoder.write_alignment_record(header, rec.as_dyn())?;
            Ok(&encoder.get_ref()[size_of::<u32>()..])
        }
    }
}

fn parse_tags(raw: &[String]) -> anyhow::Result<Vec<[u8; 2]>> {
    raw.iter()
        .filter(|s| !s.is_empty())
        .map(|s| {
            s.as_bytes()
                .try_into()
                .with_context(|| format!("invalid tag: {s}"))
        })
        .collect()
}

/// SAM text for the records of one CRAM container that pass the filter, and a buffer for
/// encoding records.
#[derive(Default)]
struct SamChunk {
    text: Vec<u8>,
    record: bam::Record,
}

/// `block_size`-prefixed BAM records for one CRAM container that pass the filter, and a buffer
/// for encoding records.
#[derive(Default)]
struct BamChunk {
    records: Vec<u8>,
    record: bam::Record,
}

/// Decodes CRAM containers on `threads` workers, formatting records there as SAM or BAM.
fn write_cram_in_parallel<R>(
    reader: &mut cram::io::Reader<R>,
    repository: &fasta::Repository,
    header: &sam::Header,
    filter: &Filter,
    kind: OutputFormat,
    threads: usize,
    writer: &mut AlignmentWriter,
) -> anyhow::Result<()>
where
    R: io::Read + Send,
{
    // Records are encoded as BAM on the workers, so SAM output uses the raw BAM formatter.
    match kind {
        OutputFormat::Sam => {
            let names = Names::new(header);

            parallel_cram::decode(
                reader,
                header,
                repository,
                threads,
                |chunk: &mut SamChunk, record| {
                    encode_cram_record(header, record, &mut chunk.record)?;

                    if filter.matches(&Rec::Bam(&chunk.record), header)?
                        && !bam_sam::format_record(&names, chunk.record.as_ref(), &mut chunk.text)
                    {
                        let mut writer = sam::io::Writer::new(&mut chunk.text);
                        writer.write_alignment_record(header, &chunk.record)?;
                    }

                    Ok(())
                },
                |chunk| Ok(writer.write_sam_text(chunk.text)?),
            )
        }
        OutputFormat::Bam => parallel_cram::decode(
            reader,
            header,
            repository,
            threads,
            |chunk: &mut BamChunk, record| {
                encode_cram_record(header, record, &mut chunk.record)?;

                if filter.matches(&Rec::Bam(&chunk.record), header)? {
                    let raw = chunk.record.as_ref();
                    let block_size = u32::try_from(raw.len())
                        .map_err(|e| io::Error::new(io::ErrorKind::InvalidInput, e))?;
                    chunk.records.extend(block_size.to_le_bytes());
                    chunk.records.extend_from_slice(raw);
                }

                Ok(())
            },
            |chunk| Ok(writer.write_bam_records(&chunk.records)?),
        ),
        OutputFormat::Cram => bail!("CRAM-to-CRAM output is written serially"),
    }
}

fn visit<F>(
    reader: &mut AlignmentReader,
    path: Option<&Path>,
    index: Option<&Path>,
    header: &sam::Header,
    traversal: &Traversal,
    f: F,
) -> anyhow::Result<()>
where
    F: FnMut(Rec<'_>) -> anyhow::Result<()>,
{
    let queries = match traversal {
        Traversal::All => return reader.for_each(header, f),
        Traversal::Regions(regions) => Queries::Regions(regions),
        Traversal::Set(set) => Queries::Set(set),
    };

    let Some(path) = path else {
        bail!("region queries need an indexed file, not stdin");
    };

    reader.for_each_in_queries(path, index, header, &queries, f)
}

fn output_format(args: &Args) -> anyhow::Result<OutputFormat> {
    let explicit = [
        args.output_format,
        (args.bam || args.uncompressed).then_some(OutputFormat::Bam),
        args.cram.then_some(OutputFormat::Cram),
    ];

    let mut chosen = explicit.into_iter().flatten();

    match (chosen.next(), chosen.next()) {
        (Some(a), Some(b)) if a != b => bail!("conflicting output formats: {a:?} and {b:?}"),
        (Some(format), _) => Ok(format),
        (None, _) => Ok(args
            .output
            .as_deref()
            .and_then(OutputFormat::from_path)
            .unwrap_or(OutputFormat::Sam)),
    }
}
