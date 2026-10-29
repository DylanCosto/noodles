//! `ntools ampliconclip`: clip primer regions from reads (samtools ampliconclip).
//!
//! This follows `bam_clip` in samtools 1.24: a read is matched against the BED regions (primers)
//! at its 5' end (or both ends), allowing `--tolerance` bases, and the longest overlapping
//! region is soft or hard clipped from it. Reads clipped to nothing are unmapped (`--unmap-len`).
//!
//! samtools quirks are kept: the `OA` tag's first field is the read name rather than the
//! reference name, and its MAPQ is not comma-separated from the CIGAR when there is no `NM`
//! tag; trimming the 3' end of a read ending in a deletion drops CIGAR operations. Unlike
//! samtools, clipped records get a recalculated BAM bin.

use std::{
    fs::File,
    io::{self, BufRead, BufReader, Write},
    path::PathBuf,
};

use anyhow::{Context, bail};
use noodles_bam as bam;
use noodles_sam::{self as sam, alignment::io::Write as _};

use crate::{
    alignment::{AlignmentReader, ReadOptions, Rec},
    bam_raw::{self, le_u16, le_u32},
    output::{self, AlignmentWriter, OutputFormat, WriteOptions},
    program,
};

const CIGAR_SOFT_CLIP: u32 = 4;
const CIGAR_HARD_CLIP: u32 = 5;

const FLAG_UNMAPPED: u16 = 0x4;
const FLAG_REVERSE: u16 = 0x10;
const FLAG_QC_FAIL: u16 = 0x200;

#[derive(Debug, clap::Args)]
pub struct Args {
    /// Input file.
    input: String,

    /// BED file of regions (such as amplicon primers) to clip.
    #[arg(short = 'b')]
    bed: PathBuf,

    /// Output file.
    #[arg(short = 'o', default_value = "-")]
    output: PathBuf,

    /// Write statistics to this file [default: stderr].
    #[arg(short = 'f')]
    stats_file: Option<PathBuf>,

    /// Output uncompressed BAM.
    #[arg(short = 'u')]
    uncompressed: bool,

    /// Output format [default: from the output extension, else BAM].
    #[arg(short = 'O', long = "output-fmt", ignore_case = true)]
    output_format: Option<OutputFormat>,

    /// Soft clip primers (the default).
    #[arg(long = "soft-clip", overrides_with = "hard_clip")]
    soft_clip: bool,

    /// Hard clip primers.
    #[arg(long = "hard-clip", overrides_with = "soft_clip")]
    hard_clip: bool,

    /// Clip both the 5' and 3' ends.
    #[arg(long = "both-ends")]
    both_ends: bool,

    /// Match only regions on the read's strand (from the BED file's sixth column).
    #[arg(long = "strand")]
    strand: bool,

    /// Output only clipped reads.
    #[arg(long = "clipped")]
    clipped: bool,

    /// Mark unclipped, mapped reads as QC failures.
    #[arg(long = "fail")]
    fail: bool,

    /// Do not output reads this size or shorter.
    #[arg(long = "filter-len", default_value = "-1", value_parser = parse_atoi, allow_negative_numbers = true)]
    filter_len: i64,

    /// Mark reads this size or shorter as QC failures.
    #[arg(long = "fail-len", default_value = "-1", value_parser = parse_atoi, allow_negative_numbers = true)]
    fail_len: i64,

    /// Unmap reads this size or shorter.
    #[arg(long = "unmap-len", default_value = "0", value_parser = parse_atoi, allow_negative_numbers = true)]
    unmap_len: i64,

    /// Do not output excluded reads (unmapped, QC failures, or on references with no regions).
    #[arg(long = "no-excluded")]
    no_excluded: bool,

    /// Write filtered reads to this file.
    #[arg(long = "rejects-file")]
    rejects_file: Option<PathBuf>,

    /// Write the number of reads clipped by each region to this file.
    #[arg(long = "primer-counts")]
    primer_counts: Option<PathBuf>,

    /// Add an OA tag with the original alignment to clipped reads.
    #[arg(long = "original")]
    original: bool,

    /// Keep the NM and MD tags of clipped reads.
    #[arg(long = "keep-tag")]
    keep_tag: bool,

    /// Match regions within this many bases.
    #[arg(long = "tolerance", default_value = "5", value_parser = parse_atoi, allow_negative_numbers = true)]
    tolerance: i64,

    /// Reference FASTA (indexed) for CRAM input or output.
    #[arg(long = "reference")]
    reference: Option<PathBuf>,

    /// Additional threads for compression and decompression.
    #[arg(short = '@', long = "threads", default_value_t = 0)]
    threads: usize,

    /// Do not add an @PG line to the header.
    #[arg(long = "no-PG")]
    no_pg: bool,
}

/// `clipping_type`
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
enum Clipping {
    Soft,
    Hard,
}

/// `bed_entry_t`
#[derive(Debug)]
struct Region {
    left: i64,
    right: i64,
    name: Option<String>,
    score: Option<String>,
    reverse: bool,
    reads: i64,
}

/// `bed_entry_list_t`: one reference's regions, sorted by end.
#[derive(Debug, Default)]
struct Regions {
    regions: Vec<Region>,
    longest: i64,
}

#[derive(Debug, Default)]
struct Counts {
    reads: i64,
    forward: i64,
    reverse: i64,
    both: i64,
    not_clipped: i64,
    excluded: i64,
    filtered: i64,
    failed: i64,
    written: i64,
}

pub fn run(args: Args) -> anyhow::Result<()> {
    let clipping = if args.hard_clip {
        Clipping::Hard
    } else {
        Clipping::Soft
    };

    let tolerance = if args.tolerance < 0 {
        eprintln!(
            "[ampliconclip] warning: invalid tolerance of {}, resetting tolerance to default of 5.",
            args.tolerance
        );
        5
    } else {
        args.tolerance
    };

    let (names, mut beds) = read_bed(&args.bed, args.strand)
        .context("[ampliconclip] error: unable to load bed file.")?;

    let read_options = ReadOptions {
        threads: args.threads,
        reference: args.reference.clone(),
    };

    let (mut reader, _) = AlignmentReader::open(&args.input, &read_options)?;
    let mut header = reader
        .read_header()
        .context("[ampliconclip] error: could not read header")?;

    // Changing positions can break coordinate order.
    set_unknown_order(&mut header);

    if !args.no_pg {
        program::add(&mut header)?;
    }

    let level = output::parse_level(if args.uncompressed { 0 } else { 6 })?;
    let write_options = |path: &PathBuf| WriteOptions {
        format: args
            .output_format
            .or_else(|| OutputFormat::from_path(path))
            .unwrap_or(OutputFormat::Bam),
        level,
        threads: args.threads,
        reference: args.reference.clone(),
    };

    let mut writer =
        AlignmentWriter::create(Some(&args.output), &write_options(&args.output), &header)?;
    writer.write_header(&header)?;

    let mut rejects = args
        .rejects_file
        .as_ref()
        .map(|path| AlignmentWriter::create(Some(path), &write_options(path), &header))
        .transpose()?;

    if let Some(rejects) = &mut rejects {
        rejects.write_header(&header)?;
    }

    let names_by_id: Vec<String> = header
        .reference_sequences()
        .keys()
        .map(|name| name.to_string())
        .collect();

    let options = ClipOptions {
        clipping,
        tolerance,
        use_strand: args.strand,
        original: args.original,
        delete_tags: !args.keep_tag,
    };

    let mut counts = Counts::default();
    let mut encoder = bam::io::Writer::from(Vec::new());
    let mut last_tid = None;
    let mut sites: Option<usize> = None;

    reader.for_each(&header.clone(), |rec| {
        let mut raw: Vec<u8> = match &rec {
            Rec::Bam(r) => r.as_ref().to_vec(),
            _ => {
                encoder.get_mut().clear();
                encoder.write_alignment_record(&header, rec.as_dyn())?;
                encoder.get_ref()[size_of::<u32>()..].to_vec()
            }
        };

        counts.reads += 1;

        let tid = le_u32(&raw[0..4]).cast_signed();

        if last_tid != Some(tid) {
            last_tid = Some(tid);
            sites = usize::try_from(tid)
                .ok()
                .and_then(|i| names_by_id.get(i))
                .and_then(|name| names.iter().position(|n| n == name));
        }

        let flags = le_u16(&raw[14..16]);
        let mut filter = false;

        match sites {
            Some(i) if flags & (FLAG_UNMAPPED | FLAG_QC_FAIL) == 0 => {
                let regions = &mut beds[i];
                let clipped = if args.both_ends {
                    clip_both_ends(&mut raw, regions, &options, &mut counts)?
                } else {
                    clip_five_prime(&mut raw, regions, &options, &mut counts)?
                };

                if !clipped && args.fail {
                    add_flags(&mut raw, FLAG_QC_FAIL);
                }

                let length = active_query_len(&raw);

                if args.fail_len >= 0 && length <= args.fail_len {
                    add_flags(&mut raw, FLAG_QC_FAIL);
                }

                if args.filter_len >= 0 && length <= args.filter_len {
                    filter = true;
                }

                if args.unmap_len >= 0 && length <= args.unmap_len {
                    raw = unmap(&raw);
                }

                if le_u16(&raw[14..16]) & FLAG_QC_FAIL != 0 {
                    counts.failed += 1;
                }

                if args.clipped && !clipped {
                    filter = true;
                }
            }
            _ => {
                counts.excluded += 1;

                if args.no_excluded {
                    filter = true;
                }
            }
        }

        if filter {
            if let Some(rejects) = &mut rejects {
                rejects.write_raw(&header, &raw)?;
            }

            counts.filtered += 1;
        } else {
            writer.write_raw(&header, &raw)?;
            counts.written += 1;
        }

        Ok(())
    })?;

    writer.finish(&header)?;

    if let Some(rejects) = rejects {
        rejects.finish(&header)?;
    }

    let command_line = std::env::args().collect::<Vec<_>>().join(" ");
    let stats = format!(
        "COMMAND: {command_line}\nTOTAL READS: {}\nTOTAL CLIPPED: {}\nFORWARD CLIPPED: {}\n\
         REVERSE CLIPPED: {}\nBOTH CLIPPED: {}\nNOT CLIPPED: {}\nEXCLUDED: {}\nFILTERED: {}\n\
         FAILED: {}\nWRITTEN: {}\n",
        counts.reads,
        counts.forward + counts.reverse,
        counts.forward,
        counts.reverse,
        counts.both,
        counts.not_clipped,
        counts.excluded,
        counts.filtered,
        counts.failed,
        counts.written,
    );

    match &args.stats_file {
        Some(path) => match File::create(path) {
            Ok(mut file) => file.write_all(stats.as_bytes())?,
            Err(_) => {
                eprintln!(
                    "[ampliconclip] warning: cannot write stats to {}.",
                    path.display()
                );
                io::stderr().write_all(stats.as_bytes())?;
            }
        },
        None => io::stderr().write_all(stats.as_bytes())?,
    }

    if let Some(path) = &args.primer_counts {
        let mut text = String::from("#CHR\tLEFT\tRIGHT\tNAME\tSCORE\tSTRAND\tNUM_CLIPPED\n");

        for (name, regions) in names.iter().zip(&beds) {
            for region in &regions.regions {
                let strand = match (args.strand, region.reverse) {
                    (false, _) => ".",
                    (true, true) => "-",
                    (true, false) => "+",
                };

                // glibc prints NULL strings as "(null)".
                text.push_str(&format!(
                    "{name}\t{}\t{}\t{}\t{}\t{strand}\t{}\n",
                    region.left,
                    region.right,
                    region.name.as_deref().unwrap_or("(null)"),
                    region.score.as_deref().unwrap_or("(null)"),
                    region.reads
                ));
            }
        }

        match File::create(path) {
            Ok(mut file) => file.write_all(text.as_bytes())?,
            Err(_) => {
                eprintln!(
                    "[ampliconclip] warning: cannot write count summary to {}.",
                    path.display()
                );
                io::stderr().write_all(text.as_bytes())?;
            }
        }
    }

    Ok(())
}

struct ClipOptions {
    clipping: Clipping,
    tolerance: i64,
    use_strand: bool,
    original: bool,
    delete_tags: bool,
}

/// Clips the 5' end (the end for reverse reads); returns whether the read was clipped.
fn clip_five_prime(
    raw: &mut Vec<u8>,
    regions: &mut Regions,
    options: &ClipOptions,
    counts: &mut Counts,
) -> anyhow::Result<bool> {
    let original = options.original.then(|| original_alignment(raw));
    let reverse = le_u16(&raw[14..16]) & FLAG_REVERSE != 0;

    let pos = if reverse {
        end_position(raw)
    } else {
        position(raw)
    };

    let size = matching_clip_site(regions, pos, reverse, options);

    if size == 0 {
        counts.not_clipped += 1;
        return Ok(false);
    }

    if reverse {
        *raw = trim_right(raw, size as u32, options.clipping);
        counts.reverse += 1;
    } else {
        *raw = trim_left(raw, size as u32, options.clipping);
        counts.forward += 1;
    }

    finish_clipped(raw, original.as_deref(), options)?;
    Ok(true)
}

/// `--both-ends`: clips the left, then the right end.
fn clip_both_ends(
    raw: &mut Vec<u8>,
    regions: &mut Regions,
    options: &ClipOptions,
    counts: &mut Counts,
) -> anyhow::Result<bool> {
    let original = options.original.then(|| original_alignment(raw));

    let size = matching_clip_site(regions, position(raw), false, options);
    let left = size > 0;

    if left {
        *raw = trim_left(raw, size as u32, options.clipping);
        counts.forward += 1;
    }

    let size = matching_clip_site(regions, end_position(raw), true, options);
    let right = size > 0;

    if right {
        *raw = trim_right(raw, size as u32, options.clipping);
        counts.reverse += 1;
    }

    if left || right {
        finish_clipped(raw, original.as_deref(), options)?;
    }

    if left && right {
        counts.both += 1;
    } else if !left && !right {
        counts.not_clipped += 1;
    }

    Ok(left || right)
}

/// Adds the OA tag and removes NM and MD, as asked, from a clipped read.
fn finish_clipped(
    raw: &mut Vec<u8>,
    original: Option<&[u8]>,
    options: &ClipOptions,
) -> anyhow::Result<()> {
    if let Some(original) = original {
        bam_raw::update_string_field(raw, *b"OA", original)?;
    }

    if options.delete_tags {
        for tag in [*b"NM", *b"MD"] {
            if let Some(range) = bam_raw::find_field(raw, tag)? {
                raw.drain(range);
            }
        }
    }

    Ok(())
}

/// `matching_clip_site`: the length to clip for the longest region the read's end falls in
/// (within the tolerance), or 0; counts the read against that region.
fn matching_clip_site(
    regions: &mut Regions,
    pos: i64,
    reverse: bool,
    options: &ClipOptions,
) -> i64 {
    let tol = options.tolerance;
    let sites = &mut regions.regions;
    let longest = regions.longest;

    // samtools keeps this in an `int`.
    let pos_tol = if reverse {
        if pos > tol { pos - tol } else { 0 }
    } else {
        pos
    };

    let (mut l, mut r) = (0, sites.len());
    let mut mid = sites.len() / 2;

    while r - l > 1 {
        if sites[mid].right <= pos_tol {
            l = mid;
        } else {
            r = mid;
        }

        mid = (l + r) / 2;
    }

    let mut size = 0;
    let mut used = None;

    for (i, site) in sites.iter().enumerate().skip(l) {
        if options.use_strand && reverse != site.reverse {
            continue;
        }

        let (mod_left, mod_right) = if reverse {
            (site.left, site.right + tol)
        } else {
            ((site.left - tol).max(0), site.right)
        };

        if pos + longest + tol < mod_right {
            break;
        }

        if pos >= mod_left && pos <= mod_right {
            let candidate = if reverse {
                pos - site.left
            } else {
                site.right - pos
            };

            if size < candidate {
                size = candidate;
                used = Some(i);
            }
        }
    }

    if let Some(i) = used {
        sites[i].reads += 1;
    }

    size
}

/// The parts of an encoded BAM record that clipping changes.
struct Parts<'a> {
    raw: &'a [u8],
    cigar: Vec<u32>,
    l_seq: usize,
    seq: &'a [u8],
    qual: &'a [u8],
    aux: &'a [u8],
}

impl<'a> Parts<'a> {
    fn new(raw: &'a [u8]) -> Self {
        let l_name = usize::from(raw[8]);
        let n_cigar = usize::from(le_u16(&raw[12..14]));
        let l_seq = le_u32(&raw[16..20]) as usize;
        let cigar_start = 32 + l_name;
        let seq_start = cigar_start + 4 * n_cigar;
        let qual_start = seq_start + l_seq.div_ceil(2);
        let aux_start = qual_start + l_seq;

        Self {
            raw,
            cigar: raw[cigar_start..seq_start]
                .as_chunks::<4>()
                .0
                .iter()
                .map(|b| u32::from_le_bytes(*b))
                .collect(),
            l_seq,
            seq: &raw[seq_start..qual_start],
            qual: &raw[qual_start..aux_start],
            aux: &raw[aux_start..],
        }
    }

    /// Encodes a record with this one's name and tags, and a new position, CIGAR, sequence (as
    /// packed bases), and qualities.
    fn encode(&self, pos: i64, cigar: &[u32], l_seq: usize, seq: &[u8], qual: &[u8]) -> Vec<u8> {
        let l_name = usize::from(self.raw[8]);
        let mut dst = Vec::with_capacity(self.raw.len());
        dst.extend_from_slice(&self.raw[..32]);
        dst[4..8].copy_from_slice(&(pos as i32).to_le_bytes());
        dst[12..14].copy_from_slice(&(cigar.len() as u16).to_le_bytes());
        dst[16..20].copy_from_slice(&(l_seq as u32).to_le_bytes());
        dst.extend_from_slice(&self.raw[32..32 + l_name]);

        for op in cigar {
            dst.extend_from_slice(&op.to_le_bytes());
        }

        dst.extend_from_slice(seq);
        dst.extend_from_slice(qual);
        dst.extend_from_slice(self.aux);
        set_bin(&mut dst);
        dst
    }
}

/// `bam_cigar_type`: bit 1 for query-consuming, bit 2 for reference-consuming operations.
fn cigar_type(op: u32) -> u32 {
    (0x3c1a7 >> ((op & 0xf) << 1)) & 3
}

fn cigar_op(len: u32, op: u32) -> u32 {
    (len << 4) | op
}

/// `bam_trim_left`: clips `bases` reference bases from the start of the alignment.
fn trim_left(raw: &[u8], bases: u32, clipping: Clipping) -> Vec<u8> {
    let parts = Parts::new(raw);
    let cigar = &parts.cigar;
    let mut ref_remove = bases;
    let mut qry_removed = 0;
    let mut hardclip = 0;
    let mut new_pos = position(raw);
    let mut i = 0;

    while i < cigar.len() {
        let (op, len) = (cigar[i] & 0xf, cigar[i] >> 4);

        if op == CIGAR_HARD_CLIP {
            hardclip += len;
        } else {
            let ty = cigar_type(op);

            if ty & 2 != 0 {
                if len <= ref_remove {
                    ref_remove -= len;
                } else {
                    break;
                }

                new_pos += i64::from(len);
            }

            if ty & 1 != 0 {
                qry_removed += len;
            }
        }

        i += 1;
    }

    if i < cigar.len() {
        let ty = cigar_type(cigar[i] & 0xf);

        if ty & 2 != 0 {
            new_pos += i64::from(ref_remove);
        }

        if ty & 1 != 0 {
            qry_removed += ref_remove;
        }
    } else {
        if clipping == Clipping::Hard {
            return parts.encode(position(raw), &[], 0, &[], &[]);
        }

        qry_removed = parts.l_seq as u32;
    }

    let mut new_cigar = Vec::with_capacity(cigar.len() + 2);

    match clipping {
        Clipping::Hard if hardclip + qry_removed > 0 => {
            new_cigar.push(cigar_op(hardclip + qry_removed, CIGAR_HARD_CLIP));
        }
        Clipping::Hard => {}
        Clipping::Soft => {
            if hardclip > 0 {
                new_cigar.push(cigar_op(hardclip, CIGAR_HARD_CLIP));
            }

            if qry_removed > 0 {
                new_cigar.push(cigar_op(qry_removed, CIGAR_SOFT_CLIP));
            }
        }
    }

    if i < cigar.len() && cigar[i] >> 4 > ref_remove {
        new_cigar.push(cigar_op((cigar[i] >> 4) - ref_remove, cigar[i] & 0xf));
        new_cigar.extend_from_slice(&cigar[i + 1..]);
    }

    // Soft clipping keeps every base.
    let removed = match clipping {
        Clipping::Soft => 0,
        Clipping::Hard => (qry_removed as usize).min(parts.l_seq),
    };

    let l_seq = parts.l_seq - removed;
    let seq = shifted_bases(parts.seq, removed, l_seq);
    parts.encode(new_pos, &new_cigar, l_seq, &seq, &parts.qual[removed..])
}

/// `bam_trim_right`: clips `bases` reference bases from the end of the alignment.
fn trim_right(raw: &[u8], bases: u32, clipping: Clipping) -> Vec<u8> {
    let parts = Parts::new(raw);
    let cigar = &parts.cigar;
    let mut ref_remove = bases;
    let mut qry_removed = 0;
    let mut hardclip = 0;
    let mut i = cigar.len() as isize - 1;

    while i >= 0 {
        let (op, len) = (cigar[i as usize] & 0xf, cigar[i as usize] >> 4);

        if op == CIGAR_HARD_CLIP {
            hardclip += len;
        } else {
            let ty = cigar_type(op);

            if ty & 2 != 0 {
                if len <= ref_remove {
                    ref_remove -= len;
                } else {
                    break;
                }
            }

            if ty & 1 != 0 {
                qry_removed += len;
            }
        }

        i -= 1;
    }

    let mut j: isize;

    if i >= 0 {
        if cigar_type(cigar[i as usize] & 0xf) & 1 != 0 {
            qry_removed += ref_remove;
        }

        j = i;

        if qry_removed > 0 {
            j += 1;
        }

        if hardclip > 0 && (clipping == Clipping::Soft || qry_removed == 0) {
            j += 1;
        }
    } else {
        if clipping == Clipping::Hard {
            return parts.encode(position(raw), &[], 0, &[], &[]);
        }

        qry_removed = parts.l_seq as u32;
        j = 0;

        if hardclip > 0 && clipping == Clipping::Soft {
            j += 1;
        }
    }

    // Written by index from the end, as samtools does: a CIGAR ending in a deletion that is
    // clipped without removing bases loses its first operations.
    let mut new_cigar = vec![0; cigar.len() + 2];
    let mut new_n = 0;

    match clipping {
        Clipping::Hard => {
            if hardclip + qry_removed > 0 {
                new_cigar[j as usize] = cigar_op(hardclip + qry_removed, CIGAR_HARD_CLIP);
                new_n += 1;
            }
        }
        Clipping::Soft => {
            if hardclip > 0 {
                new_cigar[j as usize] = cigar_op(hardclip, CIGAR_HARD_CLIP);
                new_n += 1;

                if qry_removed > 0 {
                    j -= 1;
                }
            }

            if qry_removed > 0 {
                new_cigar[j as usize] = cigar_op(qry_removed, CIGAR_SOFT_CLIP);
                new_n += 1;
            }
        }
    }

    if j > 0 {
        j -= 1;
        let op = cigar[i as usize];
        new_cigar[j as usize] = cigar_op((op >> 4) - ref_remove, op & 0xf);
        new_n += 1;
    }

    while j > 0 {
        j -= 1;
        i -= 1;
        new_cigar[j as usize] = cigar[i as usize];
        new_n += 1;
    }

    new_cigar.truncate(new_n);

    let removed = match clipping {
        Clipping::Soft => 0,
        Clipping::Hard => (qry_removed as usize).min(parts.l_seq),
    };

    let l_seq = parts.l_seq - removed;
    let seq = shifted_bases(parts.seq, 0, l_seq);
    parts.encode(position(raw), &new_cigar, l_seq, &seq, &parts.qual[..l_seq])
}

/// Packed bases `start..start + len` of a packed sequence.
fn shifted_bases(seq: &[u8], start: usize, len: usize) -> Vec<u8> {
    // Whole bytes are copied, as samtools does, leaving an odd length's last low half as it was.
    if start.is_multiple_of(2) {
        return seq[start / 2..start / 2 + len.div_ceil(2)].to_vec();
    }

    let base = |i: usize| (seq[i / 2] >> (4 * (1 - i % 2))) & 0xf;
    let mut bases = vec![0; len.div_ceil(2)];

    for k in 0..len {
        bases[k / 2] |= base(start + k) << (4 * (1 - k % 2));
    }

    bases
}

/// `active_query_len`: the query bases aligned (not soft clipped).
fn active_query_len(raw: &[u8]) -> i64 {
    Parts::new(raw)
        .cigar
        .iter()
        .filter(|&&op| cigar_type(op & 0xf) & 1 != 0 && op & 0xf != CIGAR_SOFT_CLIP)
        .map(|&op| i64::from(op >> 4))
        .sum()
}

/// `bam_set1` with the unmapped flag: no CIGAR and MAPQ 0, with the same sequence and tags.
fn unmap(raw: &[u8]) -> Vec<u8> {
    let parts = Parts::new(raw);
    let mut dst = parts.encode(position(raw), &[], parts.l_seq, parts.seq, parts.qual);
    dst[9] = 0;
    add_flags(&mut dst, FLAG_UNMAPPED);
    set_bin(&mut dst);
    dst
}

/// `tag_original_data`: the OA value (an existing one, then this alignment), with samtools'
/// read name in the reference name's place.
fn original_alignment(raw: &[u8]) -> Vec<u8> {
    let parts = Parts::new(raw);
    let mut oa = Vec::new();

    if let Ok(Some(range)) = bam_raw::find_field(raw, *b"OA") {
        if raw[range.start + 2] == b'Z' {
            oa.extend_from_slice(&raw[range.start + 3..range.end - 1]);
        } else {
            oa.extend_from_slice(b"(null)");
        }
    }

    let flags = le_u16(&raw[14..16]);
    let name = &raw[32..32 + usize::from(raw[8]) - 1];
    oa.extend_from_slice(name);
    let strand = if flags & FLAG_REVERSE != 0 { '-' } else { '+' };
    oa.extend_from_slice(format!(",{},{strand},", position(raw) + 1).as_bytes());

    for &op in &parts.cigar {
        oa.extend_from_slice((op >> 4).to_string().as_bytes());
        oa.push(b"MIDNSHP=XBBBBBBB"[(op & 0xf) as usize]);
    }

    let mapq = raw[9];

    match integer_tag(raw, *b"NM") {
        Some(nm) => oa.extend_from_slice(format!(",{mapq},{nm};").as_bytes()),
        None => oa.extend_from_slice(format!("{mapq},;").as_bytes()),
    }

    oa
}

/// htslib's `bam_aux2i` of a tag: `Some(0)` for a non-integer field.
fn integer_tag(raw: &[u8], tag: [u8; 2]) -> Option<i64> {
    let range = bam_raw::find_field(raw, tag).ok()??;
    let v = &raw[range.start + 3..range.end];

    Some(match raw[range.start + 2] {
        b'c' => i64::from(v[0].cast_signed()),
        b'C' => i64::from(v[0]),
        b's' => i64::from(i16::from_le_bytes([v[0], v[1]])),
        b'S' => i64::from(u16::from_le_bytes([v[0], v[1]])),
        b'i' => i64::from(i32::from_le_bytes([v[0], v[1], v[2], v[3]])),
        b'I' => i64::from(u32::from_le_bytes([v[0], v[1], v[2], v[3]])),
        _ => 0,
    })
}

fn position(raw: &[u8]) -> i64 {
    i64::from(le_u32(&raw[4..8]).cast_signed())
}

/// `bam_endpos`: the end of the alignment (one past the start if it covers no reference).
fn end_position(raw: &[u8]) -> i64 {
    let flags = le_u16(&raw[14..16]);
    let len = if flags & FLAG_UNMAPPED != 0 {
        0
    } else {
        reference_len(&Parts::new(raw).cigar)
    };

    position(raw) + len.max(1)
}

fn reference_len(cigar: &[u32]) -> i64 {
    cigar
        .iter()
        .filter(|&&op| cigar_type(op & 0xf) & 2 != 0)
        .map(|&op| i64::from(op >> 4))
        .sum()
}

fn add_flags(raw: &mut [u8], flags: u16) {
    let flags = le_u16(&raw[14..16]) | flags;
    raw[14..16].copy_from_slice(&flags.to_le_bytes());
}

/// Sets the BAM bin from the position and CIGAR (`hts_reg2bin`).
fn set_bin(raw: &mut [u8]) {
    let pos = position(raw);
    let flags = le_u16(&raw[14..16]);
    let len = if flags & FLAG_UNMAPPED != 0 {
        0
    } else {
        reference_len(&Parts::new(raw).cigar)
    };

    let bin = reg2bin(pos, pos + len.max(1));
    raw[10..12].copy_from_slice(&bin.to_le_bytes());
}

/// `hts_reg2bin(beg, end, 14, 5)`: the first bin of each level is (8^level - 1) / 7.
fn reg2bin(beg: i64, end: i64) -> u16 {
    let end = end - 1;

    let bin = if beg >> 14 == end >> 14 {
        4681 + (beg >> 14)
    } else if beg >> 17 == end >> 17 {
        585 + (beg >> 17)
    } else if beg >> 20 == end >> 20 {
        73 + (beg >> 20)
    } else if beg >> 23 == end >> 23 {
        9 + (beg >> 23)
    } else if beg >> 26 == end >> 26 {
        1 + (beg >> 26)
    } else {
        0
    };

    bin as u16
}

fn set_unknown_order(header: &mut sam::Header) {
    use noodles_sam::header::record::value::map::header::tag::SORT_ORDER;

    if let Some(hd) = header.header_mut()
        && hd
            .other_fields()
            .get(&SORT_ORDER)
            .is_some_and(|so| so == "coordinate")
    {
        hd.other_fields_mut().insert(SORT_ORDER, "unknown".into());
    }
}

/// `load_bed_file_multi_ref`: regions by reference, in the order references first appear, each
/// sorted by end.
fn read_bed(path: &PathBuf, use_strand: bool) -> anyhow::Result<(Vec<String>, Vec<Regions>)> {
    let file = File::open(path)
        .with_context(|| format!("[amplicon] unable to open file {}.", path.display()))?;

    let mut names: Vec<String> = Vec::new();
    let mut beds: Vec<Regions> = Vec::new();

    for (n, line) in BufReader::new(file).split(b'\n').enumerate() {
        let mut line = line?;

        if line.last() == Some(&b'\r') {
            line.pop();
        }

        if line.is_empty()
            || line[0] == b'#'
            || line.starts_with(b"track ")
            || line.starts_with(b"browser ")
        {
            continue;
        }

        let line = String::from_utf8_lossy(&line);
        let (columns, fields) = scan_bed_line(&line);
        let needed = if use_strand { 6 } else { 3 };

        if columns < needed {
            bail!(
                "[amplicon] error: invalid bed file format in line {} of {}.\n\
                 Parsed {columns} columns, but need at least {needed}\n\
                 (N.B. ref/chrom name limited to 1023 characters.)",
                n + 1,
                path.display()
            );
        }

        let reverse = if use_strand {
            match fields.strand {
                Some('+') => false,
                Some('-') => true,
                c => bail!(
                    "[amplicon] error: bad strand value in line {}, expecting '+' or '-', found '{}'.",
                    n + 1,
                    c.unwrap_or(' ')
                ),
            }
        } else {
            false
        };

        let i = match names.iter().position(|name| *name == fields.reference) {
            Some(i) => i,
            None => {
                names.push(fields.reference.clone());
                beds.push(Regions::default());
                names.len() - 1
            }
        };

        let regions = &mut beds[i];
        regions.longest = regions.longest.max(fields.right - fields.left);
        regions.regions.push(Region {
            left: fields.left,
            right: fields.right,
            name: (columns >= 4).then(|| fields.name.clone()),
            score: (columns >= 5).then(|| fields.score.clone()),
            reverse,
            reads: 0,
        });
    }

    if names.is_empty() {
        bail!("no regions in {}", path.display());
    }

    for regions in &mut beds {
        regions.regions.sort_by_key(|region| region.right);
    }

    Ok((names, beds))
}

#[derive(Debug, Default)]
struct BedFields {
    reference: String,
    left: i64,
    right: i64,
    name: String,
    score: String,
    strand: Option<char>,
}

/// `sscanf(line, "%1023s %ld %ld %1023s %1023s %c", ...)`: the number of fields read (-1 if
/// the line is blank) and their values.
fn scan_bed_line(line: &str) -> (i32, BedFields) {
    let mut fields = BedFields::default();
    let mut rest = line;

    let token = |rest: &mut &str| -> Option<String> {
        let s = rest.trim_start();
        let end = s
            .find(|c: char| c.is_ascii_whitespace())
            .unwrap_or(s.len())
            .min(1023);
        let end = (0..=end)
            .rev()
            .find(|&i| s.is_char_boundary(i))
            .unwrap_or(0);

        if end == 0 {
            *rest = s;
            return None;
        }

        *rest = &s[end..];
        Some(s[..end].to_string())
    };

    let integer = |rest: &mut &str| -> Option<i64> {
        let s = rest.trim_start();
        let sign_len = usize::from(s.starts_with(['+', '-']));
        let digits = s[sign_len..]
            .find(|c: char| !c.is_ascii_digit())
            .unwrap_or(s.len() - sign_len);

        if digits == 0 {
            *rest = s;
            return None;
        }

        let end = sign_len + digits;
        *rest = &s[end..];
        s[..end].parse().ok()
    };

    let Some(reference) = token(&mut rest) else {
        return (-1, fields);
    };
    fields.reference = reference;

    let Some(left) = integer(&mut rest) else {
        return (1, fields);
    };
    fields.left = left;

    let Some(right) = integer(&mut rest) else {
        return (2, fields);
    };
    fields.right = right;

    let Some(name) = token(&mut rest) else {
        return (3, fields);
    };
    fields.name = name;

    let Some(score) = token(&mut rest) else {
        return (4, fields);
    };
    fields.score = score;

    match rest.trim_start().chars().next() {
        Some(c) => {
            fields.strand = Some(c);
            (6, fields)
        }
        None => (5, fields),
    }
}

/// C's `atoi`: an optionally signed decimal prefix, 0 if there is none.
fn parse_atoi(s: &str) -> Result<i64, String> {
    let s = s.trim_start();
    let sign_len = usize::from(s.starts_with(['+', '-']));
    let digits = s[sign_len..]
        .find(|c: char| !c.is_ascii_digit())
        .unwrap_or(s.len() - sign_len);

    Ok(s[..sign_len + digits].parse::<i64>().unwrap_or(0) as i32 as i64)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_scan_bed_line() {
        let (n, f) = scan_bed_line("chr1\t10\t30\tp1\t60\t-");
        assert_eq!(n, 6);
        assert_eq!((f.left, f.right, f.strand), (10, 30, Some('-')));

        assert_eq!(scan_bed_line("chr1 10 30").0, 3);
        assert_eq!(scan_bed_line("chr1 x 30").0, 1);
        assert_eq!(scan_bed_line("   ").0, -1);
    }

    #[test]
    fn test_shifted_bases() {
        // ACGT
        let seq = [0x12, 0x48];
        assert_eq!(shifted_bases(&seq, 1, 3), [0x24, 0x80]);
        assert_eq!(shifted_bases(&seq, 2, 2), [0x48]);
    }

    #[test]
    fn test_reg2bin() {
        assert_eq!(reg2bin(-1, 0), 4680);
        assert_eq!(reg2bin(0, 1), 4681);
    }
}
