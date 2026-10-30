//! `ntools view`: convert, filter, count, and query alignments (samtools view).

use std::{
    collections::HashSet,
    fs,
    io::{self, Write},
    mem,
    num::NonZeroUsize,
    path::{Path, PathBuf},
};

use anyhow::{Context, bail};
use bstr::BString;
use noodles_bam as bam;
use noodles_cram as cram;
use noodles_fasta as fasta;
use noodles_sam::{
    self as sam,
    alignment::{io::Write as _, record::Flags},
    header::record::value::{Map, map::ReferenceSequence},
};

use crate::{
    alignment::{AlignmentReader, Queries, RawEncoder, ReadOptions, Rec, encode_cram_record},
    bam_raw,
    bam_sam::{self, Names},
    commands::index::{IndexKind, index_file},
    expr::Expression,
    filter::{self, Filter, ReadGroupFilter, ReadGroupPresence, Subsample, TagFilter},
    flags,
    output::{self, AlignmentWriter, OutputFormat, WriteOptions},
    parallel_cram, program,
    regions::{RegionSet, record_interval},
    source::Format,
};

/// The BGZF compression level htslib uses by default (zlib's).
const DEFAULT_LEVEL: u8 = 6;

/// The region name samtools uses for unplaced, unmapped records.
const UNMAPPED_REGION: &str = "*";

/// The notes of samtools' long help (`-?`), for the options ntools supports.
const NOTES: &str = "\
Notes:

1. The input format (SAM, BAM, or CRAM) is detected automatically.

2. The file given with -t lists a reference name and length on each line, separated by a tab
   (the first two columns of a .fai file). It is used for SAM input without @SQ lines.

3. SAM->BAM conversion:  ntools view -b -o out.bam in.sam

4. BAM->SAM conversion:  ntools view -h in.bam

5. A region is `chr1', `chr2:1,000', or `chr3:1000-2,000'. Regions need a sorted, indexed BAM
   or CRAM file.

6. -M/--use-index computes overlaps with -L BED regions and region arguments through the index.
   This is faster, omits duplicates, and outputs the reads in file order.

7. -e takes a filter expression, e.g. 'mapq >= 30 && [NM] < 5 && qname =~ \"^read\"'. Fields:
   qname, flag (and flag.paired, flag.proper_pair, flag.unmap, flag.munmap, flag.reverse,
   flag.mreverse, flag.read1, flag.read2, flag.secondary, flag.qcfail, flag.dup,
   flag.supplementary), rname, refid, pos, endpos, mapq, cigar, ncigar, rnext, mrname, mrefid,
   pnext, mpos, tlen, seq, qual, qlen, rlen, sclen, hclen, library, and tags as [XX].
   Functions: length, min, max, avg, exists, default, exp, log, pow, sqrt. Operators, by
   precedence: ! ~ (unary + -), * / %, + -, &, ^, |, < <= > >=, == != =~ !~, && ||. Regular
   expressions are POSIX extended.
";

#[derive(Debug, clap::Args)]
#[command(disable_help_flag = true, after_long_help = NOTES)]
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

    /// Output format option; only `level=N`, the BGZF compression level, is supported.
    #[arg(long = "output-fmt-option", value_name = "OPT=VAL")]
    output_format_options: Vec<String>,

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
    #[arg(short = 'U', long = "unoutput", visible_alias = "output-unselected")]
    unselected_output: Option<PathBuf>,

    /// Set the unmapped flag on reads that fail the filters (clearing their CIGAR, mapping
    /// quality, and template length) and write them to the output.
    #[arg(short = 'p', long = "unmap")]
    unmap: bool,

    /// Retrieve complete pairs: also output the mates, wherever they are, of paired reads in
    /// the regions or -L BED file. Requires an index.
    #[arg(short = 'P', long = "fetch-pairs")]
    fetch_pairs: bool,

    /// Reference names and lengths (the first two columns of a `.fai` file), for SAM input
    /// without `@SQ` lines.
    #[arg(short = 't', long = "fai-reference")]
    fai_reference: Option<PathBuf>,

    /// Only include reads overlapping this BED file.
    #[arg(short = 'L', long = "target-file", visible_alias = "targets-file")]
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
    #[arg(short = 'F', long = "excl-flags", visible_alias = "exclude-flags", value_parser = flags::parse)]
    exclude_flags: Option<Flags>,

    /// Exclude reads with all of these flags set.
    #[arg(short = 'G', value_parser = flags::parse)]
    exclude_all_flags: Option<Flags>,

    /// Only include reads in this read group, or in no read group (may be repeated).
    #[arg(short = 'r', long = "read-group", visible_alias = "readgroup")]
    read_groups: Vec<String>,

    /// Only include reads in a read group listed in this file (whitespace-separated IDs), or in
    /// no read group; with `^FILE`, exclude the listed read groups instead.
    #[arg(
        short = 'R',
        long = "read-group-file",
        visible_alias = "readgroup-file",
        value_name = "[^]FILE"
    )]
    read_group_files: Vec<String>,

    /// Exclude reads without a read group.
    #[arg(
        short = 'n',
        long = "exclude-no-read-group",
        visible_aliases = ["excl-no-read-group", "excl-no-readgroup", "exclude-no-readgroup"]
    )]
    exclude_no_read_group: bool,

    /// Only include reads whose read group belongs to this library (`LB`).
    #[arg(short = 'l', long = "library")]
    library: Option<String>,

    /// Exclude reads with mapping quality below this value.
    #[arg(short = 'q', long = "min-MQ", visible_aliases = ["min-mq", "min-mapq"], default_value_t = 0)]
    min_mapping_quality: u8,

    /// Only include reads covering at least this many query bases (by CIGAR, with hard clips).
    #[arg(
        short = 'm',
        long = "min-qlen",
        default_value_t = 0,
        allow_negative_numbers = true
    )]
    min_query_length: i64,

    /// Only include reads matching this filter expression (see the notes in the long help).
    #[arg(
        short = 'e',
        long = "expr",
        visible_alias = "expression",
        allow_hyphen_values = true
    )]
    expression: Option<String>,

    /// Subsample: keep a FRAC fraction of templates, seeded by INT (`INT.FRAC`).
    #[arg(short = 's', long = "subsample", value_parser = Subsample::parse)]
    subsample: Option<Subsample>,

    /// Only include reads whose names are listed in this file.
    #[arg(short = 'N', long = "qname-file")]
    qname_file: Option<PathBuf>,

    /// Only include reads with this tag, and optionally this value (`TAG[:VALUE]`).
    #[arg(short = 'd', long = "tag")]
    tags: Vec<String>,

    /// Only include reads with this tag whose value is listed in FILE (`TAG:FILE`).
    #[arg(short = 'D', long = "tag-file", value_name = "TAG:FILE")]
    tag_files: Vec<String>,

    /// Remove these tags from output records (comma-separated; may be repeated).
    #[arg(short = 'x', long = "remove-tag", value_delimiter = ',')]
    remove_tags: Vec<String>,

    /// Collapse the backward (`B`) CIGAR operation.
    #[arg(short = 'B', long = "remove-B")]
    remove_backward_ops: bool,

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

    /// Ignored (the input format is detected automatically).
    #[arg(short = 'S')]
    sam_input: bool,

    /// Print long help, including notes on regions and filter expressions.
    #[arg(short = '?', long = "help", action = clap::ArgAction::HelpLong)]
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
    /// The records of these intervals and the mates of their paired reads (`-P`).
    Pairs(RegionSet),
}

pub fn run(args: Args) -> anyhow::Result<()> {
    if args.unmap && args.unselected_output.is_some() {
        bail!("Options --unoutput and --unmap are mutually exclusive.");
    }

    if args.count && args.fetch_pairs {
        bail!("The options -P and -c cannot be combined");
    }

    let format = output_format(&args)?;
    let reference = reference(&args);

    let read_options = ReadOptions {
        threads: args.threads,
        reference: reference.clone(),
    };

    let (mut reader, path) = AlignmentReader::open(&args.input, &read_options)?;
    let mut header = reader.read_header()?;

    if reader.format() == Format::Sam
        && header.reference_sequences().is_empty()
        && let Some(fai) = fai_path(&args)
    {
        add_fai_reference_sequences(&mut header, &fai)?;
    }

    let (traversal, mut filter) = plan(&args, &header, &reader, path.as_deref())?;
    configure_filter(&mut filter, &args, &mut header)?;
    let filter = filter;

    let remove_tags = parse_tags(&args.remove_tags)?;
    let edits = Edits::new(&args);

    let decoding = if matches!(traversal, Traversal::All)
        && args.threads > 0
        && remove_tags.is_empty()
        && args.unselected_output.is_none()
        && edits.is_none()
    {
        Decoding::Parallel
    } else {
        Decoding::Serial
    };

    if args.count {
        let mut selector = Selector::new(&filter, edits);
        let n = count(
            &mut reader,
            path.as_deref(),
            &args,
            &header,
            &traversal,
            &mut selector,
            decoding,
        )?;
        writeln!(io::stdout().lock(), "{n}")?;
        return Ok(());
    }

    if !args.no_pg {
        program::add(&mut header)?;
    }

    let write_options = write_options(&args, format, reference)?;

    let mut writer = AlignmentWriter::create(args.output.as_deref(), &write_options, &header)?;

    let mut rejected = if args.unmap {
        Rejected::Unmapped
    } else if let Some(dst) = args.unselected_output.as_deref() {
        Rejected::Written(AlignmentWriter::create(Some(dst), &write_options, &header)?)
    } else {
        Rejected::Dropped
    };

    let writes_header = format != OutputFormat::Sam || args.with_header || args.header_only;

    if writes_header {
        writer.write_header(&header)?;
    }

    if let Rejected::Written(unselected) = &mut rejected
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
            let mut sink = Sink {
                selector: Selector::new(&filter, edits),
                writer: &mut writer,
                rejected: &mut rejected,
                remove_tags: &remove_tags,
                encoder: RawEncoder::new(),
                stripped: Vec::new(),
            };

            if let Traversal::Pairs(regions) = &traversal {
                fetch_pairs(
                    &mut reader,
                    path.as_deref(),
                    args.index.as_deref(),
                    &header,
                    regions,
                    &mut sink,
                )?;
            } else {
                visit(
                    &mut reader,
                    path.as_deref(),
                    args.index.as_deref(),
                    &header,
                    &traversal,
                    |rec| sink.emit(rec, &header),
                )?;
            }
        }
    }

    writer.finish(&header)?;

    if let Rejected::Written(unselected) = rejected {
        unselected.finish(&header)?;
    }

    if args.write_index {
        write_output_index(&args, format)?;
    }

    Ok(())
}

/// Sets the record filter's options (`-f`, `-F`, `-G`, `-q`, `-m`, `-r`/`-R`, `-n`, `-l`, `-s`,
/// `-d`/`-D`, `-N`, `-e`). As in samtools, `-r`/`-R` also remove the `@RG` lines of read groups
/// whose records are dropped.
fn configure_filter(
    filter: &mut Filter,
    args: &Args,
    header: &mut sam::Header,
) -> anyhow::Result<()> {
    filter.read_group_ids = ReadGroupFilter::from_args(&args.read_groups, &args.read_group_files)?;

    if let Some(read_groups) = &filter.read_group_ids {
        read_groups.retain_header_lines(header);
    }

    filter.read_group_presence = if args.exclude_no_read_group {
        ReadGroupPresence::Required
    } else {
        ReadGroupPresence::Optional
    };

    filter.require = args.require_flags.unwrap_or_else(Flags::empty);
    filter.exclude = args.exclude_flags.unwrap_or_else(Flags::empty);
    filter.exclude_all = args.exclude_all_flags;
    filter.min_mapping_quality = args.min_mapping_quality;
    filter.min_query_length = args.min_query_length;
    filter.read_groups = args
        .library
        .as_deref()
        .map(|library| filter::library_read_groups(header, library));
    filter.subsample = args.subsample;
    filter.tag = TagFilter::from_args(&args.tags, &args.tag_files)?;
    filter.names = args
        .qname_file
        .as_deref()
        .map(filter::read_names)
        .transpose()?;
    filter.expression = args
        .expression
        .as_deref()
        .map(Expression::parse)
        .transpose()?;

    Ok(())
}

/// The CRAM reference: `-T`, or else the FASTA of a `-t` index (`ref.fa` for `ref.fa.fai`), as
/// HTSlib loads it.
fn reference(args: &Args) -> Option<PathBuf> {
    args.reference.clone().or_else(|| {
        let fai = args.fai_reference.as_deref()?;
        let fasta = fai.to_str()?.strip_suffix(".fai")?;
        Path::new(fasta).is_file().then(|| PathBuf::from(fasta))
    })
}

/// The reference index giving `@SQ` lines to SAM input without them: `-t`, or the `.fai` of the
/// `-T` FASTA when it exists.
fn fai_path(args: &Args) -> Option<PathBuf> {
    args.fai_reference.clone().or_else(|| {
        let mut fai = args.reference.clone()?.into_os_string();
        fai.push(".fai");
        let fai = PathBuf::from(fai);
        fai.is_file().then_some(fai)
    })
}

/// Adds an `@SQ` line for each line of a reference index (`NAME<TAB>LENGTH...`), as HTSlib does
/// when reading SAM without `@SQ` lines. Lines without a tab are skipped.
fn add_fai_reference_sequences(header: &mut sam::Header, path: &Path) -> anyhow::Result<()> {
    let text = fs::read(path).with_context(|| format!("reading {}", path.display()))?;

    for line in text.split(|&b| b == b'\n') {
        let Some(tab) = memchr::memchr(b'\t', line) else {
            continue;
        };

        let (name, rest) = (&line[..tab], &line[tab + 1..]);
        let digits: Vec<u8> = rest
            .iter()
            .copied()
            .skip_while(u8::is_ascii_whitespace)
            .take_while(u8::is_ascii_digit)
            .collect();

        let length = std::str::from_utf8(&digits)
            .ok()
            .and_then(|s| s.parse::<usize>().ok())
            .and_then(NonZeroUsize::new)
            .with_context(|| {
                format!(
                    "{}: invalid length for {}",
                    path.display(),
                    String::from_utf8_lossy(name)
                )
            })?;

        header
            .reference_sequences_mut()
            .insert(BString::from(name), Map::<ReferenceSequence>::new(length));
    }

    Ok(())
}

/// The output options, with the compression level from `-u`, `-1`, or `--output-fmt-option`.
fn write_options(
    args: &Args,
    format: OutputFormat,
    reference: Option<PathBuf>,
) -> anyhow::Result<WriteOptions> {
    let mut level = DEFAULT_LEVEL;

    for option in &args.output_format_options {
        match option.split_once('=') {
            Some(("level", value)) => {
                level = value
                    .parse()
                    .with_context(|| format!("invalid compression level: {value}"))?;
            }
            _ => bail!("unsupported output format option: {option}"),
        }
    }

    if args.uncompressed {
        level = 0;
    } else if args.fast {
        level = 1;
    }

    Ok(WriteOptions {
        format,
        level: output::parse_level(level)?,
        threads: args.threads,
        reference,
    })
}

/// Indexes the output file (`--write-index`).
fn write_output_index(args: &Args, format: OutputFormat) -> anyhow::Result<()> {
    let Some(output) = args.output.as_deref().filter(|p| *p != Path::new("-")) else {
        bail!("--write-index needs an output file (-o)");
    };

    let kind = match format {
        OutputFormat::Bam => IndexKind::Csi,
        OutputFormat::Cram => IndexKind::Default,
        OutputFormat::Sam => bail!("--write-index needs BAM or CRAM output"),
    };

    index_file(output, None, kind, args.threads)?;

    Ok(())
}

/// Chooses the traversal and any region filter: `-M` merges the regions given; `-L` restricts to
/// a BED file, through the index when one exists and no regions are given; `-P` reads the
/// regions (intersected with `-L`) through the index, then their mates.
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

    let region_set = || -> anyhow::Result<RegionSet> {
        if args.regions.iter().any(|r| r == UNMAPPED_REGION) {
            bail!("-M and -P cannot be combined with the unmapped region *");
        }

        let mut set = RegionSet::new(header);
        for region in &args.regions {
            set.add_region(header, region)?;
        }
        set.merge();
        Ok(set)
    };

    if args.fetch_pairs {
        let set = match (args.regions.is_empty(), bed) {
            (false, Some(bed)) => region_set()?.intersection(&bed),
            (false, None) => region_set()?,
            (true, Some(bed)) => bed,
            (true, None) => bail!("No regions or BED file have been provided. Aborting."),
        };

        return Ok((Traversal::Pairs(set), filter));
    }

    if !args.regions.is_empty() {
        let traversal = if args.merge_regions {
            Traversal::Set(region_set()?)
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

/// Record edits samtools makes around filtering (`-B` after the expression, `-p` for rejected
/// records), which work on an encoded copy of each record.
struct Edits {
    remove_backward_ops: bool,
}

impl Edits {
    /// Returns `None` when no record is edited, so records are filtered as they are read.
    fn new(args: &Args) -> Option<Self> {
        // `-p` is ignored when counting.
        let unmaps = args.unmap && !args.count;

        (args.remove_backward_ops || unmaps).then_some(Self {
            remove_backward_ops: args.remove_backward_ops,
        })
    }
}

/// Applies the filter, and any edits, to records.
struct Selector<'f> {
    filter: &'f Filter,
    edits: Option<Edits>,
    encoder: RawEncoder,
    /// The edited copy of the current record.
    record: bam::Record,
}

impl<'f> Selector<'f> {
    fn new(filter: &'f Filter, edits: Option<Edits>) -> Self {
        Self {
            filter,
            edits,
            encoder: RawEncoder::new(),
            record: bam::Record::default(),
        }
    }

    /// Returns whether a record is selected, and the record as edited (`process_aln`): tested
    /// against the expression, then stripped of `B` operations, then tested against the other
    /// filters.
    fn select<'a>(&'a mut self, rec: Rec<'a>, header: &sam::Header) -> io::Result<(bool, Rec<'a>)> {
        let Some(edits) = &self.edits else {
            let selected = self.filter.matches(&rec, header)?;
            return Ok((selected, rec));
        };

        let mut buf = Vec::from(mem::take(&mut self.record));
        buf.clear();
        buf.extend_from_slice(rec.encode_raw(header, &mut self.encoder)?);
        self.record = bam::Record::try_from(buf)?;

        if !self
            .filter
            .matches_expression(&Rec::Bam(&self.record), header)?
        {
            return Ok((false, Rec::Bam(&self.record)));
        }

        if edits.remove_backward_ops {
            let mut buf = Vec::from(mem::take(&mut self.record));
            bam_raw::remove_backward_ops(&mut buf)?;
            self.record = bam::Record::try_from(buf)?;
        }

        let selected = self
            .filter
            .matches_fields(&Rec::Bam(&self.record), header)?;
        Ok((selected, Rec::Bam(&self.record)))
    }

    /// Unmaps the record last selected (`-p`); only valid with edits.
    fn unmap(&mut self) -> io::Result<Rec<'_>> {
        let mut buf = Vec::from(mem::take(&mut self.record));
        bam_raw::unmap(&mut buf)?;
        self.record = bam::Record::try_from(buf)?;
        Ok(Rec::Bam(&self.record))
    }
}

/// Where records that fail the filters go.
enum Rejected {
    Dropped,
    /// Unmapped and written to the output (`-p`).
    Unmapped,
    /// Written as they are to another file (`-U`).
    Written(AlignmentWriter),
}

/// Filters records and writes them.
struct Sink<'a, 'f> {
    selector: Selector<'f>,
    writer: &'a mut AlignmentWriter,
    rejected: &'a mut Rejected,
    remove_tags: &'a [[u8; 2]],
    encoder: RawEncoder,
    stripped: Vec<u8>,
}

impl Sink<'_, '_> {
    fn emit(&mut self, rec: Rec<'_>, header: &sam::Header) -> anyhow::Result<()> {
        let (selected, rec) = self.selector.select(rec, header)?;

        if selected {
            if self.remove_tags.is_empty() {
                self.writer.write(header, &rec)?;
            } else {
                let raw = rec.encode_raw(header, &mut self.encoder)?;
                bam_raw::strip_fields(raw, &mut self.stripped, |tag| {
                    self.remove_tags.contains(&tag)
                })?;
                self.writer.write_raw(header, &self.stripped)?;
            }

            return Ok(());
        }

        match self.rejected {
            Rejected::Dropped => {}
            Rejected::Written(unselected) => unselected.write(header, &rec)?,
            Rejected::Unmapped => {
                let rec = self.selector.unmap()?;
                self.writer.write(header, &rec)?;
            }
        }

        Ok(())
    }
}

fn count(
    reader: &mut AlignmentReader,
    path: Option<&Path>,
    args: &Args,
    header: &sam::Header,
    traversal: &Traversal,
    selector: &mut Selector<'_>,
    decoding: Decoding,
) -> anyhow::Result<u64> {
    let mut n = 0u64;

    if let (AlignmentReader::Cram(cram_reader, repository, _), Decoding::Parallel) =
        (&mut *reader, decoding)
    {
        let filter = selector.filter;

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
            if selector.select(rec, header)?.0 {
                n += 1;
            }
            Ok(())
        },
    )?;

    Ok(n)
}

/// `-P` (`fetch_pairs_collect_mates`): reads the regions once to collect the names and mate
/// positions of the paired reads that pass the filters and whose mates lie outside the regions,
/// then reads the regions and the mate positions (and the unplaced reads, if a mate is
/// unplaced), writing the records that overlap the regions or have a collected name.
fn fetch_pairs(
    reader: &mut AlignmentReader,
    path: Option<&Path>,
    index: Option<&Path>,
    header: &sam::Header,
    regions: &RegionSet,
    sink: &mut Sink<'_, '_>,
) -> anyhow::Result<()> {
    let Some(path) = path else {
        bail!("-P needs an indexed file, not stdin");
    };

    let mut names: HashSet<BString> = HashSet::new();
    let mut mates = regions.clone();
    let mut unplaced_mates = false;

    reader.for_each_in_queries(path, index, header, &Queries::Set(regions), |rec| {
        if !rec.flags()?.is_segmented() {
            return Ok(());
        }

        let mate_id = rec.mate_reference_sequence_id(header)?;
        let mate_start = rec.as_dyn().mate_alignment_start().transpose()?;

        if let (Some(id), Some(start)) = (mate_id, mate_start)
            && mate_in_regions(regions, id, usize::from(start))
        {
            return Ok(());
        }

        let name = record_name(&rec);

        if !sink.selector.select(rec, header)?.0 {
            return Ok(());
        }

        names.insert(name);

        match (mate_id, mate_start) {
            (None, _) => unplaced_mates = true,
            (Some(id), Some(start)) => {
                let start = usize::from(start);
                mates.add_interval(id, start, start);
            }
            (Some(_), None) => {}
        }

        Ok(())
    })?;

    mates.merge();

    let mut emit = |rec: Rec<'_>| -> anyhow::Result<()> {
        let in_regions = record_interval(rec.as_dyn(), header)?
            .is_some_and(|(id, start, end)| regions.first_overlap(id, start, end).is_some());

        if in_regions || names.contains(&record_name(&rec)) {
            sink.emit(rec, header)?;
        }

        Ok(())
    };

    reader.for_each_in_queries(path, index, header, &Queries::Set(&mates), &mut emit)?;

    if unplaced_mates {
        let unmapped = [UNMAPPED_REGION.to_owned()];
        reader.for_each_in_queries(path, index, header, &Queries::Regions(&unmapped), emit)?;
    }

    Ok(())
}

/// Whether a mate at 1-based `position` is in the regions as samtools tests it: as the empty
/// interval at the position (`bed_overlap(mpos, mpos)`), which overlaps an interval only if the
/// interval also covers the base before it.
fn mate_in_regions(regions: &RegionSet, reference_sequence_id: usize, position: usize) -> bool {
    regions
        .first_overlap(reference_sequence_id, position, position)
        .is_some_and(|i| regions.intervals(reference_sequence_id)[i].0 < position)
}

fn record_name(rec: &Rec<'_>) -> BString {
    rec.as_dyn()
        .name()
        .map_or_else(|| BString::from("*"), |name| name.to_owned())
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
        Traversal::Set(set) | Traversal::Pairs(set) => Queries::Set(set),
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
