//! `ntools reference`: rebuild a reference from alignments (samtools reference).
//!
//! This follows samtools 1.24's `reference.c`. By default each reference sequence is filled with
//! `N`s and then patched from every record's sequence, CIGAR, and `MD` tag. With `-e`, a CRAM
//! file's embedded reference blocks are copied out instead. Unlike samtools, reference skips
//! (`N` CIGAR operations) advance the reference position, so spliced reads land in the right
//! place, and a region with no data is labelled with its range, as one with data is.

use std::{
    fs::File,
    io::{self, BufReader, BufWriter, Read, SeekFrom, Write},
    path::PathBuf,
};

use anyhow::{Context, bail};
use noodles_cram as cram;
use noodles_sam as sam;

use crate::global_args::{IgnoredWriteIndex, Verbosity};
use crate::{
    alignment::{AlignmentIndex, AlignmentReader, RawEncoder, ReadOptions, Rec},
    bam_raw::{BASES, Record},
    cram_raw::read_slice_header,
    hts_region::{self, POS_MAX},
};

#[derive(Debug, clap::Args)]
pub struct Args {
    /// Input file [default: stdin].
    input: Option<String>,

    /// Write the reference to this file.
    #[arg(short = 'o', long = "output")]
    output: Option<PathBuf>,

    /// Do not report coverage on stderr.
    #[arg(short = 'q', long = "quiet")]
    quiet: bool,

    /// Copy the reference embedded in a CRAM file instead of using MD tags.
    #[arg(short = 'e', long = "embedded")]
    embedded: bool,

    /// Only output this region (requires an index).
    #[arg(short = 'r', long = "region")]
    region: Option<String>,

    /// Additional decompression threads.
    #[arg(short = '@', long = "threads", default_value_t = 0)]
    threads: usize,

    /// FASTA reference for decoding CRAM input.
    #[arg(long = "reference")]
    reference: Option<PathBuf>,

    #[command(flatten)]
    write_index: IgnoredWriteIndex,

    #[command(flatten)]
    verbosity: Verbosity,
}

/// The FASTA line width.
const LINE_WIDTH: usize = 60;

// BAM CIGAR operation codes.
const CIGAR_MATCH: u32 = 0;
const CIGAR_DELETION: u32 = 2;
const CIGAR_SKIP: u32 = 3;
const CIGAR_SOFT_CLIP: u32 = 4;
const CIGAR_EQUAL: u32 = 7;
const CIGAR_DIFF: u32 = 8;

pub fn run(args: Args) -> anyhow::Result<()> {
    let src = args.input.as_deref().unwrap_or("-");

    let mut out: Box<dyn Write> = match &args.output {
        Some(path) => Box::new(BufWriter::new(
            File::create(path).with_context(|| format!("creating {}", path.display()))?,
        )),
        None => Box::new(BufWriter::new(io::stdout().lock())),
    };

    let mut dumper = Dumper {
        out: &mut out,
        verbose: !args.quiet,
    };

    if args.embedded {
        embedded_reference(src, args.region.as_deref(), &mut dumper)
    } else {
        md_reference(&args, src, &mut dumper)
    }
    .context("[reference]")?;

    out.flush()?;

    Ok(())
}

/// A region to report: a reference sequence ID and a 0-based, half-open interval.
#[derive(Clone, Copy)]
struct Region {
    reference_sequence_id: usize,
    start: usize,
    /// The end, or `None` for the end of the reference sequence.
    end: Option<usize>,
}

fn parse_region(header: &sam::Header, raw: &str) -> anyhow::Result<Region> {
    let reference_sequences = header.reference_sequences();

    let (reference_sequence_id, start, end) =
        hts_region::parse_region(raw, |name| reference_sequences.get_index_of(name))
            .with_context(|| format!("failed to parse region '{raw}'"))?;

    Ok(Region {
        reference_sequence_id,
        start: usize::try_from(start)?,
        end: (end < POS_MAX).then(|| usize::try_from(end)).transpose()?,
    })
}

/// Writes reference sequences as FASTA.
struct Dumper<'a> {
    out: &'a mut dyn Write,
    verbose: bool,
}

impl Dumper<'_> {
    /// Writes `sequence` (a whole reference sequence), or the region's part of it (htslib's
    /// `dump_ref`).
    fn dump(
        &mut self,
        header: &sam::Header,
        reference_sequence_id: usize,
        sequence: &[u8],
        region: Option<Region>,
    ) -> io::Result<()> {
        let name = header
            .reference_sequences()
            .get_index(reference_sequence_id)
            .map(|(name, _)| name.as_ref())
            .unwrap_or_default();

        let mut sequence = sequence;

        match region {
            Some(Region { start, end, .. })
                if start > 0 || end.is_some_and(|e| e < sequence.len()) =>
            {
                let end = end.unwrap_or(sequence.len());
                self.out.write_all(b">")?;
                self.out.write_all(name)?;
                writeln!(self.out, ":{}-{end}", start + 1)?;
                sequence = sequence
                    .get(start..end.min(sequence.len()))
                    .unwrap_or_default();
            }
            _ => {
                self.out.write_all(b">")?;
                self.out.write_all(name)?;
                writeln!(self.out)?;
            }
        }

        for line in sequence.chunks(LINE_WIDTH) {
            self.out.write_all(line)?;
            writeln!(self.out)?;
        }

        if self.verbose {
            let n = sequence.iter().filter(|&&b| b == b'N').count();
            let coverage = 100.0 - n as f64 * 100.0 / sequence.len() as f64;
            eprintln!(
                "Dump ref {reference_sequence_id} len {}, coverage {coverage:.2}%",
                sequence.len()
            );
        }

        Ok(())
    }

    /// Writes an all-`N` region, for a region query with no data.
    fn dump_empty(&mut self, header: &sam::Header, region: Region) -> io::Result<()> {
        let len = reference_sequence_len(header, region.reference_sequence_id);
        let sequence = vec![b'N'; len];
        self.dump(
            header,
            region.reference_sequence_id,
            &sequence,
            Some(region),
        )
    }
}

fn reference_sequence_len(header: &sam::Header, reference_sequence_id: usize) -> usize {
    header
        .reference_sequences()
        .get_index(reference_sequence_id)
        .map(|(_, rs)| rs.length().get())
        .unwrap_or_default()
}

/// The reference sequence being rebuilt.
struct Current {
    /// The reference sequence ID, or `None` for unplaced records.
    reference_sequence_id: Option<usize>,
    sequence: Vec<u8>,
}

impl Current {
    fn new(header: &sam::Header, reference_sequence_id: Option<usize>) -> Self {
        let len = reference_sequence_id.map_or(0, |id| reference_sequence_len(header, id));

        Self {
            reference_sequence_id,
            sequence: vec![b'N'; len],
        }
    }

    fn dump(
        &self,
        dumper: &mut Dumper<'_>,
        header: &sam::Header,
        region: Option<Region>,
    ) -> io::Result<()> {
        match self.reference_sequence_id {
            Some(id) => dumper.dump(header, id, &self.sequence, region),
            None => Ok(()),
        }
    }
}

/// Rebuilds references from records' sequences and MD tags (samtools' `MD2ref`).
fn md_reference(args: &Args, src: &str, dumper: &mut Dumper<'_>) -> anyhow::Result<()> {
    let options = ReadOptions {
        threads: args.threads,
        reference: args.reference.clone(),
    };

    let (mut reader, path) = AlignmentReader::open(src, &options)
        .with_context(|| format!("failed to open file '{src}'"))?;
    let header = reader.read_header()?;

    let region = args
        .region
        .as_deref()
        .map(|raw| parse_region(&header, raw))
        .transpose()?;

    let mut current: Option<Current> = None;
    let mut encoder = RawEncoder::new();
    let mut cigar = Vec::new();

    let mut visit = |rec: Rec<'_>| -> anyhow::Result<()> {
        let record = rec.encode_raw(&header, &mut encoder)?;

        let record = Record::new(record)?;
        let tid = record.reference_sequence_id();
        let reference_sequence_id = usize::try_from(tid).ok();

        if current
            .as_ref()
            .is_none_or(|c| c.reference_sequence_id != reference_sequence_id)
        {
            if let Some(c) = &current {
                c.dump(dumper, &header, region)?;
            }

            current = Some(Current::new(&header, reference_sequence_id));
        }

        let sequence = &mut current.as_mut().expect("current is set above").sequence;
        patch_reference(&record, &mut cigar, sequence)
    };

    match region {
        Some(r) => {
            let path = path.context("region queries need an indexed file, not stdin")?;
            let index = AlignmentIndex::read(reader.format(), &path, None)
                .context("Failed to load the index")?;
            reader.for_each_in_interval(
                &index,
                &header,
                r.reference_sequence_id,
                r.start,
                r.end,
                &mut visit,
            )?;
        }
        None => reader.for_each(&header, &mut visit)?,
    }

    match (&current, region) {
        (Some(c), _) if c.reference_sequence_id.is_some() => c.dump(dumper, &header, region)?,
        (_, Some(r)) => dumper.dump_empty(&header, r)?,
        _ => {}
    }

    Ok(())
}

/// Copies a record's bases into `reference` where its MD tag says they match, and the MD tag's
/// bases where they differ (samtools' `build_ref`). Records without an `MD:Z` tag are skipped.
fn patch_reference(
    record: &Record<'_>,
    cigar: &mut Vec<u32>,
    reference: &mut [u8],
) -> anyhow::Result<()> {
    let Some(md) = record.find_field(*b"MD")? else {
        return Ok(());
    };

    let raw = record.as_bytes();

    if raw[md.start + 2] != b'Z' {
        return Ok(());
    }

    let md = &raw[md.start + 3..md.end - 1];
    let l_seq = record.sequence_len();

    record.read_cigar(cigar)?;

    let pos = record.position();
    let mut ops = CigarWalker::new(cigar);
    let mut iseq = 0;
    let mut iref = i64::from(pos);

    let mut set = |iref: i64, base: u8| {
        if let Some(b) = usize::try_from(iref)
            .ok()
            .and_then(|i| reference.get_mut(i))
        {
            *b = base;
        }
    };

    let mut i = 0;

    while iseq < l_seq && i < md.len() {
        if md[i].is_ascii_digit() {
            let digits = md[i..].iter().take_while(|b| b.is_ascii_digit()).count();
            let mut len = md[i..i + digits].iter().fold(0usize, |n, &d| {
                n.saturating_mul(10).saturating_add(usize::from(d - b'0'))
            });
            i += digits;

            while iseq < l_seq && len > 0 {
                let op = ops.next(&mut iseq, &mut iref)?;

                if op != CIGAR_MATCH && op != CIGAR_EQUAL {
                    bail!("MD:Z and CIGAR are incompatible");
                }

                if iseq >= l_seq {
                    bail!("MD:Z and CIGAR are incompatible");
                }

                set(iref, BASES[usize::from(record.base(iseq))]);
                iseq += 1;
                iref += 1;
                len -= 1;
            }
        } else if md[i] == b'^' {
            i += 1;

            while i < md.len() && md[i].is_ascii_alphabetic() {
                if ops.next(&mut iseq, &mut iref)? != CIGAR_DELETION {
                    bail!("MD:Z and CIGAR are incompatible");
                }

                set(iref, md[i]);
                i += 1;
                iref += 1;
            }
        } else {
            let op = ops.next(&mut iseq, &mut iref)?;

            if op != CIGAR_MATCH && op != CIGAR_DIFF {
                bail!("MD:Z and CIGAR are incompatible");
            }

            set(iref, md[i]);
            i += 1;
            iref += 1;
            iseq += 1;
        }
    }

    Ok(())
}

/// Steps through CIGAR operations one base at a time, skipping insertions, clips, padding, and
/// reference skips (samtools' `next_cigar_op`).
struct CigarWalker<'a> {
    cigar: &'a [u32],
    op: u32,
    len: u32,
}

impl<'a> CigarWalker<'a> {
    fn new(cigar: &'a [u32]) -> Self {
        Self {
            cigar,
            op: 0,
            len: 0,
        }
    }

    /// Returns the next M, D, =, or X operation base, advancing `iseq` over skipped query bases
    /// and `iref` over skipped reference bases.
    fn next(&mut self, iseq: &mut usize, iref: &mut i64) -> anyhow::Result<u32> {
        loop {
            while self.len == 0 {
                let (&raw, rest) = self
                    .cigar
                    .split_first()
                    .context("MD:Z is longer than the CIGAR")?;
                self.cigar = rest;
                self.op = raw & 0xf;
                self.len = raw >> 4;
            }

            match self.op {
                CIGAR_MATCH | CIGAR_DELETION | CIGAR_EQUAL | CIGAR_DIFF => {
                    self.len -= 1;
                    return Ok(self.op);
                }
                // Insertions and soft clips consume the query.
                1 | CIGAR_SOFT_CLIP => *iseq += self.len as usize,
                // samtools does not advance the reference here, misplacing spliced reads.
                CIGAR_SKIP => *iref += i64::from(self.len),
                _ => {}
            }

            self.len = 0;
        }
    }
}

/// Copies the embedded reference blocks of a CRAM file (samtools' `cram2ref`).
fn embedded_reference(
    src: &str,
    region: Option<&str>,
    dumper: &mut Dumper<'_>,
) -> anyhow::Result<()> {
    if src == "-" {
        if region.is_some() {
            bail!("region queries need an indexed file, not stdin");
        }

        let mut reader = cram::io::Reader::new(io::stdin().lock());
        let header = read_cram_header(&mut reader, src)?;
        return copy_embedded(&mut reader, &header, None, dumper);
    }

    let file = File::open(src).with_context(|| format!("failed to open file '{src}'"))?;
    let mut reader = cram::io::Reader::new(BufReader::new(file));
    let header = read_cram_header(&mut reader, src)?;

    let Some(region) = region.map(|raw| parse_region(&header, raw)).transpose()? else {
        return copy_embedded(&mut reader, &header, None, dumper);
    };

    let index = cram::fs::read_associated_index(src).context("Failed to load the index")?;

    // Start at the first container with a slice of the region's reference sequence that ends
    // after the region start.
    let offset = index
        .iter()
        .filter(|record| record.reference_sequence_id() == Some(region.reference_sequence_id))
        .filter(|record| {
            record
                .alignment_start()
                .is_none_or(|start| start.get() + record.alignment_span() > region.start)
        })
        .map(|record| record.offset())
        .min();

    match offset {
        Some(offset) => {
            reader.seek(SeekFrom::Start(offset))?;
            copy_embedded(&mut reader, &header, Some(region), dumper)
        }
        None => Ok(dumper.dump_empty(&header, region)?),
    }
}

fn read_cram_header<R: Read>(
    reader: &mut cram::io::Reader<R>,
    src: &str,
) -> anyhow::Result<sam::Header> {
    reader.read_header().map_err(|e| {
        if e.kind() == io::ErrorKind::InvalidData {
            anyhow::anyhow!("{src} is not a CRAM file: {e}")
        } else {
            e.into()
        }
    })
}

/// Rebuilds references from the embedded reference blocks of the containers from the reader's
/// position on.
fn copy_embedded<R: Read>(
    reader: &mut cram::io::Reader<R>,
    header: &sam::Header,
    region: Option<Region>,
    dumper: &mut Dumper<'_>,
) -> anyhow::Result<()> {
    let mut current: Option<Current> = None;
    let mut container = cram::io::reader::Container::default();

    'containers: while reader.read_container(&mut container)? != 0 {
        let mut src = container.src();
        cram::io::reader::read_block(&mut src)?; // compression header

        for _ in container.header().landmarks() {
            let block = cram::io::reader::read_block(&mut src)?;
            let slice_header = read_slice_header(&mut &block.decode()?[..])?;
            let tid = slice_header.reference_sequence_id;
            let reference_sequence_id = usize::try_from(tid).ok();

            if let Some(r) = region {
                if reference_sequence_id != Some(r.reference_sequence_id) {
                    // Slices of other reference sequences could precede the region's in a
                    // multi-reference container.
                    if current.is_some() {
                        break 'containers;
                    }

                    for _ in 0..slice_header.block_count {
                        cram::io::reader::read_block(&mut src)?;
                    }

                    continue;
                }

                if r.end
                    .is_some_and(|end| i64::from(slice_header.alignment_start) > end as i64)
                {
                    break 'containers;
                }
            }

            if slice_header.embedded_reference.is_none() && tid != -1 {
                bail!("CRAM file has slice without embedded reference");
            }

            if current
                .as_ref()
                .is_none_or(|c| c.reference_sequence_id != reference_sequence_id)
            {
                if let Some(c) = &current {
                    c.dump(dumper, header, region)?;
                }

                current = Some(Current::new(header, reference_sequence_id));
            }

            let sequence = &mut current.as_mut().expect("current is set above").sequence;

            for _ in 0..slice_header.block_count {
                let block = cram::io::reader::read_block(&mut src)?;

                if block.content_type == cram::container::block::ContentType::ExternalData
                    && Some(block.content_id) == slice_header.embedded_reference
                {
                    let data = block.decode()?;
                    let start = usize::try_from(slice_header.alignment_start)
                        .unwrap_or(0)
                        .saturating_sub(1);
                    let end = (start + data.len()).min(sequence.len());

                    if end > start {
                        sequence[start..end].copy_from_slice(&data[..end - start]);
                    }
                }
            }
        }
    }

    match (&current, region) {
        (Some(c), _) if c.reference_sequence_id.is_some() => c.dump(dumper, header, region)?,
        (_, Some(r)) => dumper.dump_empty(header, r)?,
        _ => {}
    }

    Ok(())
}

#[cfg(test)]
mod tests {
    use noodles_bam as bam;
    use noodles_sam::alignment::io::Write as _;

    use super::*;

    fn encode(line: &str) -> anyhow::Result<Vec<u8>> {
        let header: sam::Header = "@SQ\tSN:chr\tLN:20\n".parse()?;
        let record = sam::Record::try_from(line.as_bytes())?;
        let mut writer = bam::io::Writer::from(Vec::new());
        writer.write_alignment_record(&header, &record)?;
        Ok(writer.get_ref()[size_of::<u32>()..].to_vec())
    }

    fn patch(line: &str) -> anyhow::Result<String> {
        let record = encode(line)?;
        let mut reference = vec![b'N'; 20];
        patch_reference(&Record::new(&record)?, &mut Vec::new(), &mut reference)?;
        Ok(String::from_utf8(reference)?)
    }

    #[test]
    fn test_patch_reference() -> anyhow::Result<()> {
        // A mismatch, an insertion, a deletion, and a soft clip.
        assert_eq!(
            patch("r\t0\tchr\t3\t60\t1S3M1I2M2D2M\t*\t0\t0\tTACTGGAAC\t*\tMD:Z:1G3^TT2")?,
            "NNAGTGATTACNNNNNNNNN"
        );

        // Bases after a reference skip land after the skipped bases.
        assert_eq!(
            patch("r\t0\tchr\t2\t60\t2M3N2M\t*\t0\t0\tACGT\t*\tMD:Z:4")?,
            "NACNNNGTNNNNNNNNNNNN"
        );

        // Records without an MD tag are skipped, and an MD tag longer than the CIGAR is an error.
        assert_eq!(
            patch("r\t0\tchr\t2\t60\t2M\t*\t0\t0\tAC\t*")?,
            "N".repeat(20)
        );
        assert!(patch("r\t0\tchr\t2\t60\t1M1I\t*\t0\t0\tAC\t*\tMD:Z:2").is_err());

        Ok(())
    }
}
