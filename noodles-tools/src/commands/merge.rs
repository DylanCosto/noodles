//! `ntools merge`: merge sorted BAM files (samtools merge).

use std::{io::BufReader, path::PathBuf};

use anyhow::{Context, bail};
use noodles_bam as bam;
use noodles_sam::{self as sam, header::ReferenceSequences};

use crate::{
    bgzf_stream::BgzfReader,
    commands::sort::set_sort_order,
    output::{self, AlignmentWriter, OutputFormat, WriteOptions},
    program,
    runs::{self, RecordSource},
    sort_order::{Order, remap_reference_sequence_ids},
    source::{Format, Input, Stream},
};

#[derive(Debug, clap::Args)]
pub struct Args {
    /// Output file, then sorted input BAM files (or give the output with -o).
    #[arg(required = true)]
    files: Vec<PathBuf>,

    /// Write output to this file; all positional arguments are then inputs.
    #[arg(short = 'o', long = "output")]
    output: Option<PathBuf>,

    /// Overwrite the output file if it exists.
    #[arg(short = 'f', long = "force")]
    force: bool,

    /// Output format [default: from the output extension, else BAM].
    #[arg(short = 'O', long = "output-fmt", ignore_case = true)]
    output_format: Option<OutputFormat>,

    /// BGZF compression level for BAM output.
    #[arg(short = 'l', default_value_t = 6)]
    level: u8,

    /// Output uncompressed BAM.
    #[arg(short = 'u')]
    uncompressed: bool,

    /// Inputs are sorted by read name, compared naturally.
    #[arg(short = 'n', conflicts_with = "name_ascii")]
    name_natural: bool,

    /// Inputs are sorted by read name, bytewise.
    #[arg(short = 'N')]
    name_ascii: bool,

    /// Reference FASTA (indexed) for CRAM output.
    #[arg(long = "reference")]
    reference: Option<PathBuf>,

    /// Additional threads for compression and decompression.
    #[arg(short = '@', long = "threads", default_value_t = 0)]
    threads: usize,

    /// Do not add an @PG line to the header.
    #[arg(long = "no-PG")]
    no_pg: bool,
}

/// A sorted BAM input, with reference sequence IDs translated to the merged header's.
struct BamInput {
    reader: BgzfReader<BufReader<Stream>>,
    reference_sequence_id_map: Option<Vec<i32>>,
}

impl RecordSource for BamInput {
    fn next_record(&mut self, buf: &mut Vec<u8>) -> std::io::Result<bool> {
        let has_record = runs::read_raw_record(&mut self.reader, buf)?;

        if has_record && let Some(map) = &self.reference_sequence_id_map {
            remap_reference_sequence_ids(buf, map);
        }

        Ok(has_record)
    }
}

pub fn run(args: Args) -> anyhow::Result<()> {
    let order = if args.name_natural {
        Order::NameNatural
    } else if args.name_ascii {
        Order::NameAscii
    } else {
        Order::Coordinate
    };

    let mut files = args.files.into_iter();
    let output = match args.output {
        Some(output) => output,
        None => files.next().context("missing output file")?,
    };
    let inputs: Vec<PathBuf> = files.collect();

    if inputs.is_empty() {
        bail!("no input files");
    }

    if !args.force && output.exists() {
        bail!("{} exists; use -f to overwrite it", output.display());
    }

    let threads_per_input = args.threads / inputs.len();
    let mut header: Option<sam::Header> = None;
    let mut sources = Vec::with_capacity(inputs.len());
    let mut coordinate_order_lost = false;

    for path in &inputs {
        let src = path.to_str().context("input path is not valid UTF-8")?;
        let input = Input::open(src)?;

        if input.format != Format::Bam {
            bail!("{src}: only BAM inputs can be merged");
        }

        let mut reader = bam::io::Reader::from(BgzfReader::new(input.reader, threads_per_input));
        let input_header = reader
            .read_header()
            .with_context(|| format!("reading header of {src}"))?;

        let map = match &mut header {
            None => {
                header = Some(input_header);
                None
            }
            Some(merged) => merge_header(merged, &input_header)
                .with_context(|| format!("merging header of {src}"))?,
        };

        // If the merged header orders this input's references differently, its records are no
        // longer coordinate-sorted after translation. samtools warns and marks the output unsorted.
        if order == Order::Coordinate
            && let Some(map) = &map
            && !map.is_sorted_by(|a, b| a < b)
        {
            eprintln!("[merge] Order of targets in file {src} caused coordinate sort to be lost");
            coordinate_order_lost = true;
        }

        sources.push(BamInput {
            reader: reader.into_inner(),
            reference_sequence_id_map: map,
        });
    }

    let mut header = header.context("no input headers")?;
    set_sort_order(&mut header, order)?;

    if coordinate_order_lost {
        set_unsorted(&mut header);
    }

    if !args.no_pg {
        program::add(&mut header)?;
    }

    let format = args
        .output_format
        .or_else(|| OutputFormat::from_path(&output))
        .unwrap_or(OutputFormat::Bam);

    let write_options = WriteOptions {
        format,
        level: output::parse_level(if args.uncompressed { 0 } else { args.level })?,
        threads: args.threads,
        reference: args.reference,
    };

    let mut writer = AlignmentWriter::create(Some(&output), &write_options)?;
    writer.write_header(&header)?;

    runs::merge(&mut sources, order, |record| {
        writer.write_raw(&header, record)?;
        Ok(())
    })?;

    writer.finish(&header)?;

    Ok(())
}

fn set_unsorted(header: &mut sam::Header) {
    use noodles_sam::header::record::value::map::header::{sort_order::UNSORTED, tag::SORT_ORDER};

    if let Some(hd) = header.header_mut() {
        hd.other_fields_mut().insert(SORT_ORDER, UNSORTED.into());
    }
}

/// Adds `other`'s reference sequences, read groups, programs, and comments that `merged` lacks.
///
/// Read groups and programs whose IDs are already present are kept from the earlier input and
/// records are not rewritten, as with `samtools merge -c -p`.
///
/// Returns the map from `other`'s reference sequence IDs to the merged IDs, or `None` when they
/// are the same.
fn merge_header(merged: &mut sam::Header, other: &sam::Header) -> anyhow::Result<Option<Vec<i32>>> {
    let map = merge_reference_sequences(
        merged.reference_sequences_mut(),
        other.reference_sequences(),
    )?;

    for (id, read_group) in other.read_groups() {
        merged
            .read_groups_mut()
            .entry(id.clone())
            .or_insert_with(|| read_group.clone());
    }

    for (id, program) in other.programs().as_ref() {
        if !merged.programs().as_ref().contains_key(id) {
            merged
                .programs_mut()
                .as_mut()
                .insert(id.clone(), program.clone());
        }
    }

    for comment in other.comments() {
        if !merged.comments().contains(comment) {
            merged.comments_mut().push(comment.clone());
        }
    }

    let identity = map
        .iter()
        .enumerate()
        .all(|(i, &id)| usize::try_from(id).is_ok_and(|id| id == i));

    Ok((!identity).then_some(map))
}

fn merge_reference_sequences(
    merged: &mut ReferenceSequences,
    other: &ReferenceSequences,
) -> anyhow::Result<Vec<i32>> {
    other
        .iter()
        .map(|(name, reference_sequence)| {
            let id = match merged.get_full(name) {
                Some((id, _, existing)) => {
                    if existing.length() != reference_sequence.length() {
                        bail!(
                            "reference sequence {name} has length {} here but {} in an earlier input",
                            reference_sequence.length(),
                            existing.length()
                        );
                    }

                    id
                }
                None => {
                    merged.insert(name.clone(), reference_sequence.clone());
                    merged.len() - 1
                }
            };

            i32::try_from(id).context("too many reference sequences")
        })
        .collect()
}
