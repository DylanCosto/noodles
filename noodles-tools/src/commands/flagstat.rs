//! `ntools flagstat`: counts of records by flag category (samtools flagstat).

use std::{
    fmt,
    io::{self, Write},
    path::PathBuf,
};

use crate::alignment::{AlignmentReader, ReadOptions, Rec};
use crate::global_args::{IgnoredWriteIndex, InputFormatOptions, Verbosity};

/// Mate-on-another-reference records at or above this mapping quality are counted separately.
const MIN_HQ_MAPPING_QUALITY: u8 = 5;

#[derive(Clone, Copy, Debug, Eq, PartialEq, clap::ValueEnum)]
enum OutputFormat {
    Default,
    Tsv,
    Json,
}

#[derive(Debug, clap::Args)]
pub struct Args {
    /// Input SAM, BAM, or CRAM file (`-` for stdin).
    input: String,

    /// Output format.
    #[arg(short = 'O', long = "output-fmt", value_enum, default_value_t = OutputFormat::Default, ignore_case = true)]
    output_format: OutputFormat,

    /// Reference FASTA (indexed) for CRAM input.
    #[arg(short = 'T', long = "reference")]
    reference: Option<PathBuf>,

    /// Additional threads for BGZF decompression.
    #[arg(short = '@', long = "threads", default_value_t = 0)]
    threads: usize,

    #[command(flatten)]
    input_format_options: InputFormatOptions,

    #[command(flatten)]
    write_index: IgnoredWriteIndex,

    #[command(flatten)]
    verbosity: Verbosity,
}

#[derive(Debug, Default)]
struct Counts {
    total: u64,
    primary: u64,
    secondary: u64,
    supplementary: u64,
    duplicate: u64,
    primary_duplicate: u64,
    mapped: u64,
    primary_mapped: u64,
    paired: u64,
    read_1: u64,
    read_2: u64,
    proper_pair: u64,
    mate_mapped: u64,
    singleton: u64,
    mate_on_other_reference: u64,
    mate_on_other_reference_hq: u64,
}

impl Counts {
    fn add(&mut self, record: &Rec<'_>, header: &noodles_sam::Header) -> io::Result<()> {
        let flags = record.flags()?;

        self.total += 1;

        if !flags.is_unmapped() {
            self.mapped += 1;
        }

        if flags.is_duplicate() {
            self.duplicate += 1;
        }

        if flags.is_secondary() {
            self.secondary += 1;
            return Ok(());
        } else if flags.is_supplementary() {
            self.supplementary += 1;
            return Ok(());
        }

        self.primary += 1;

        if !flags.is_unmapped() {
            self.primary_mapped += 1;
        }

        if flags.is_duplicate() {
            self.primary_duplicate += 1;
        }

        if !flags.is_segmented() {
            return Ok(());
        }

        self.paired += 1;

        if flags.is_properly_segmented() && !flags.is_unmapped() {
            self.proper_pair += 1;
        }

        if flags.is_first_segment() {
            self.read_1 += 1;
        }

        if flags.is_last_segment() {
            self.read_2 += 1;
        }

        if flags.is_mate_unmapped() && !flags.is_unmapped() {
            self.singleton += 1;
        }

        if !flags.is_unmapped() && !flags.is_mate_unmapped() {
            self.mate_mapped += 1;

            if record.mate_reference_sequence_id(header)? != record.reference_sequence_id(header)? {
                self.mate_on_other_reference += 1;

                // samtools compares the raw byte, so a missing quality (255) counts as high.
                let mapq = record.mapping_quality()?.unwrap_or(u8::MAX);
                if mapq >= MIN_HQ_MAPPING_QUALITY {
                    self.mate_on_other_reference_hq += 1;
                }
            }
        }

        Ok(())
    }
}

/// One output statistic: its labels, its count, and the count it is a percentage of, if any.
struct Stat {
    text: &'static str,
    tsv: &'static str,
    json: &'static str,
    count: fn(&Counts) -> u64,
    percent_of: Option<fn(&Counts) -> u64>,
}

const STATS: [Stat; 16] = [
    Stat {
        text: "in total (QC-passed reads + QC-failed reads)",
        tsv: "total (QC-passed reads + QC-failed reads)",
        json: "total",
        count: |c| c.total,
        percent_of: None,
    },
    stat("primary", |c| c.primary),
    stat("secondary", |c| c.secondary),
    stat("supplementary", |c| c.supplementary),
    stat("duplicates", |c| c.duplicate),
    stat("primary duplicates", |c| c.primary_duplicate),
    Stat {
        percent_of: Some(|c| c.total),
        ..stat("mapped", |c| c.mapped)
    },
    Stat {
        percent_of: Some(|c| c.primary),
        ..stat("primary mapped", |c| c.primary_mapped)
    },
    stat("paired in sequencing", |c| c.paired),
    stat("read1", |c| c.read_1),
    stat("read2", |c| c.read_2),
    Stat {
        percent_of: Some(|c| c.paired),
        ..stat("properly paired", |c| c.proper_pair)
    },
    stat("with itself and mate mapped", |c| c.mate_mapped),
    Stat {
        percent_of: Some(|c| c.paired),
        ..stat("singletons", |c| c.singleton)
    },
    stat("with mate mapped to a different chr", |c| {
        c.mate_on_other_reference
    }),
    Stat {
        json: "with mate mapped to a different chr (mapQ >= 5)",
        ..stat("with mate mapped to a different chr (mapQ>=5)", |c| {
            c.mate_on_other_reference_hq
        })
    },
];

const fn stat(label: &'static str, count: fn(&Counts) -> u64) -> Stat {
    Stat {
        text: label,
        tsv: label,
        json: label,
        count,
        percent_of: None,
    }
}

/// A percentage formatted as samtools does: two decimals, or `N/A` for a zero denominator.
struct Percent(u64, u64);

impl Percent {
    fn value(&self) -> Option<f64> {
        #[expect(
            clippy::cast_precision_loss,
            reason = "read counts are far below 2^52, and samtools also divides as floating point"
        )]
        (self.1 != 0).then(|| self.0 as f64 / self.1 as f64 * 100.0)
    }
}

impl fmt::Display for Percent {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self.value() {
            Some(p) => write!(f, "{p:.2}%"),
            None => f.write_str("N/A"),
        }
    }
}

pub fn run(args: Args) -> anyhow::Result<()> {
    let mut options = ReadOptions {
        threads: args.threads,
        reference: args.reference,
    };
    args.input_format_options.configure(&mut options)?;

    let (mut reader, _) = AlignmentReader::open(&args.input, &options)?;
    args.input_format_options.apply(&mut reader)?;
    let header = reader.read_header()?;

    let mut pass = Counts::default();
    let mut fail = Counts::default();

    reader.for_each(&header, |record| {
        let counts = if record.flags()?.is_qc_fail() {
            &mut fail
        } else {
            &mut pass
        };

        counts.add(&record, &header)?;
        Ok(())
    })?;

    let mut stdout = io::BufWriter::new(io::stdout().lock());

    match args.output_format {
        OutputFormat::Default => write_default(&mut stdout, &pass, &fail)?,
        OutputFormat::Tsv => write_tsv(&mut stdout, &pass, &fail)?,
        OutputFormat::Json => write_json(&mut stdout, &pass, &fail)?,
    }

    stdout.flush()?;

    Ok(())
}

fn write_default<W: Write>(w: &mut W, pass: &Counts, fail: &Counts) -> io::Result<()> {
    for s in &STATS {
        write!(w, "{} + {} {}", (s.count)(pass), (s.count)(fail), s.text)?;

        if let Some(of) = s.percent_of {
            write!(
                w,
                " ({} : {})",
                Percent((s.count)(pass), of(pass)),
                Percent((s.count)(fail), of(fail))
            )?;
        }

        writeln!(w)?;
    }

    Ok(())
}

fn write_tsv<W: Write>(w: &mut W, pass: &Counts, fail: &Counts) -> io::Result<()> {
    for s in &STATS {
        writeln!(w, "{}\t{}\t{}", (s.count)(pass), (s.count)(fail), s.tsv)?;

        if let Some(of) = s.percent_of {
            writeln!(
                w,
                "{}\t{}\t{} %",
                Percent((s.count)(pass), of(pass)),
                Percent((s.count)(fail), of(fail)),
                s.tsv
            )?;
        }
    }

    Ok(())
}

fn write_json<W: Write>(w: &mut W, pass: &Counts, fail: &Counts) -> io::Result<()> {
    // The spacing (including trailing spaces) matches samtools byte for byte.
    writeln!(w, "{{")?;

    for (i, (name, counts)) in [("QC-passed reads", pass), ("QC-failed reads", fail)]
        .into_iter()
        .enumerate()
    {
        writeln!(w, " \"{name}\": {{ ")?;

        let mut fields = Vec::new();

        for s in &STATS {
            fields.push(format!("  \"{}\": {}", s.json, (s.count)(counts)));

            if let Some(of) = s.percent_of {
                let value = match Percent((s.count)(counts), of(counts)).value() {
                    Some(p) => format!("{p:.2}"),
                    None => "null".to_owned(),
                };
                fields.push(format!("  \"{} %\": {value}", s.json));
            }
        }

        writeln!(w, "{} ", fields.join(", \n"))?;
        writeln!(w, " }}{}", if i == 0 { "," } else { "" })?;
    }

    writeln!(w, "}}")
}
