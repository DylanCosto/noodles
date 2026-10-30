//! `ntools checksum`: order-independent checksums of alignment data (samtools checksum).
//!
//! This follows samtools 1.24's `bam_checksum.c` (after biobambam's `bamseqchksum`): for each
//! record, CRC32s of flag + sequence, + name, + qualities, + selected tags (and optionally
//! position, CIGAR, and mate fields) are combined multiplicatively modulo 2^31 - 1, overall,
//! for records without a read group, and per read group. `-m` merges such reports.

use std::{
    collections::HashMap,
    fs::File,
    io::{self, BufRead, BufReader, Write},
    path::PathBuf,
};

use anyhow::{Context, bail};

use crate::{
    alignment::{AlignmentReader, RawEncoder, ReadOptions},
    bam_raw::{BASES, Record, field},
    commands::flags,
    sanitize::{Sanitize, sanitize_raw},
    sequence::Strand,
};

/// The checksums' modulus.
const PRIME: u64 = (1 << 31) - 1;

const FLAG_QC_FAIL: u16 = 0x200;

const SEQ_NT16_COMPLEMENT_STR: &[u8; 16] = b"=TGKCYSBAWRDMHVN";

#[derive(Debug, clap::Args)]
pub struct Args {
    /// Input files (or, with -m, checksum reports) [default: stdin].
    inputs: Vec<String>,

    /// Skip records with any of these flags.
    #[arg(short = 'F', long = "exclude-flags", default_value = "0x900")]
    exclude_flags: String,

    /// Skip records without all of these flags.
    #[arg(short = 'f', long = "require-flags", default_value = "0")]
    require_flags: String,

    /// Flags to include in checksums.
    #[arg(short = 'b', long = "flag-mask", default_value = "0x0c1")]
    flag_mask: String,

    /// Do not reverse complement sequences of reverse-strand records.
    #[arg(short = 'c', long = "no-rev-comp")]
    no_rev_comp: bool,

    /// Tags to checksum (`XX,YY`), or `*,XX,YY` for all but those.
    #[arg(short = 't', long = "tags", default_value = "BC,FI,QT,RT,TC")]
    tags: String,

    /// Order-specific checksums (repeat to count within each group).
    #[arg(short = 'O', long = "in-order", action = clap::ArgAction::Count)]
    in_order: u8,

    /// Also checksum the reference and position.
    #[arg(short = 'P', long = "check-pos")]
    check_pos: bool,

    /// Also checksum the mapping quality and CIGAR.
    #[arg(short = 'C', long = "check-cigar")]
    check_cigar: bool,

    /// Also checksum the mate reference, position, and template length.
    #[arg(short = 'M', long = "check-mate")]
    check_mate: bool,

    /// Fix records first (as fixmate -z).
    #[arg(short = 'z', long = "sanitize")]
    sanitize: Option<String>,

    /// Stop after this many records (0 for all).
    #[arg(short = 'N', long = "count", default_value = "0")]
    count: String,

    /// Write the report to this file.
    #[arg(short = 'o', long = "output")]
    output: Option<PathBuf>,

    /// Also show QC pass and fail lines.
    #[arg(short = 'q', long = "show-qc")]
    show_qc: bool,

    /// Show lines with zero counts.
    #[arg(short = 'v', long = "verbose", action = clap::ArgAction::Count)]
    verbose: u8,

    /// Checksum everything: -PCMOc -b 0xfff -f0 -F0 -z all,cigarx -t *,cF,MD,NM.
    #[arg(short = 'a', long = "all")]
    all: bool,

    /// Tab-delimited output.
    #[arg(short = 'T', long = "tabs")]
    tabs: bool,

    /// Merge checksum reports instead of reading alignments.
    #[arg(short = 'm', long = "merge")]
    merge: bool,

    /// Report in bamseqchksum's format.
    #[arg(short = 'B', long = "bamseqchksum")]
    bamseqchksum: bool,

    /// Reference FASTA (indexed) for CRAM input.
    #[arg(long = "reference")]
    reference: Option<PathBuf>,

    /// Additional threads for decompression.
    #[arg(short = '@', long = "threads", default_value_t = 0)]
    threads: usize,
}

/// `opts`
#[derive(Clone, Debug)]
struct Options {
    require_flags: u16,
    exclude_flags: u16,
    flag_mask: i32,
    rev_comp: bool,
    in_order: u8,
    sanitize: Sanitize,
    check_pos: bool,
    check_cigar: bool,
    check_mate: bool,
    tag_str: String,
    tags: Vec<Vec<u8>>,
    records: i64,
    verbose: bool,
    show_pass: bool,
    show_fail: bool,
    show_combine: bool,
    tabs: bool,
    compat: bool,
}

/// `crcs_t`: one record's CRCs.
#[derive(Clone, Copy, Debug)]
struct Crcs {
    seq: u32,
    name: u32,
    qual: u32,
    aux: u32,
    pos: u32,
    cigar: u32,
    mate: u32,
}

/// `sums_t`: the combined hashes of a group, for all, QC-passed, and QC-failed records.
#[derive(Clone, Copy, Debug)]
struct Sums {
    seq: [u64; 3],
    name: [u64; 3],
    qual: [u64; 3],
    aux: [u64; 3],
    pos: [u64; 3],
    cigar: [u64; 3],
    mate: [u64; 3],
    count: [u64; 3],
}

impl Default for Sums {
    fn default() -> Self {
        Self {
            seq: [1; 3],
            name: [1; 3],
            qual: [1; 3],
            aux: [1; 3],
            pos: [1; 3],
            cigar: [1; 3],
            mate: [1; 3],
            count: [0; 3],
        }
    }
}

impl Sums {
    /// `sums_update_row`
    fn update_row(&mut self, row: usize, c: &Crcs, count_crc: u32, n: u64) {
        self.seq[row] = update_hash(self.seq[row], count_crc ^ c.seq);
        self.name[row] = update_hash(self.name[row], count_crc ^ c.name);
        self.qual[row] = update_hash(self.qual[row], count_crc ^ c.qual);
        self.aux[row] = update_hash(self.aux[row], count_crc ^ c.aux);
        self.pos[row] = update_hash(self.pos[row], count_crc ^ c.pos);
        self.cigar[row] = update_hash(self.cigar[row], count_crc ^ c.cigar);
        self.mate[row] = update_hash(self.mate[row], count_crc ^ c.mate);
        self.count[row] += n;
    }

    /// `sums_update`: the all row, and the pass or fail row.
    fn update(&mut self, qc: Qc, c: &Crcs, options: &Options, count: u64) {
        let count_crc = match options.in_order {
            0 => 0,
            1 => crc32(0, &count.to_le_bytes()),
            _ => crc32(0, &self.count[0].to_le_bytes()),
        };

        self.update_row(0, c, count_crc, 1);

        if options.show_pass && qc == Qc::Pass {
            self.update_row(1, c, count_crc, 1);
        }

        if options.show_fail && qc == Qc::Fail {
            self.update_row(2, c, count_crc, 1);
        }
    }
}

/// Whether a read passed QC (its QC-fail flag is unset).
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
enum Qc {
    Pass,
    Fail,
}

/// Which tags `-t` selects.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
enum TagSelection {
    /// The listed tags.
    Listed,
    /// Every tag but those listed (a list starting with `*`).
    AllExceptListed,
}

/// `update_hash`: multiplication modulo 2^31 - 1, with 0 (and the modulus) replaced by 1.
fn update_hash(hash: u64, crc: u32) -> u64 {
    let mut crc = u64::from(crc) & PRIME;

    if crc == 0 || crc == PRIME {
        crc = 1;
    }

    (hash * crc) % PRIME
}

/// zlib's `crc32(crc, buf, len)`.
fn crc32(crc: u32, data: &[u8]) -> u32 {
    let mut hasher = crc32fast::Hasher::new_with_initial(crc);
    hasher.update(data);
    hasher.finalize()
}

pub fn run(args: Args) -> anyhow::Result<()> {
    let flag = |s: &str| {
        flags::parse(s)
            .filter(|&n| n >= 0)
            .with_context(|| format!("[checksum] could not parse flag {s}"))
    };

    let mut options = Options {
        require_flags: flag(&args.require_flags)? as u16,
        exclude_flags: flag(&args.exclude_flags)? as u16,
        flag_mask: flag(&args.flag_mask)?,
        rev_comp: !args.no_rev_comp,
        in_order: args.in_order,
        sanitize: args
            .sanitize
            .as_deref()
            .map(Sanitize::parse)
            .transpose()?
            .unwrap_or(Sanitize(0)),
        check_pos: args.check_pos,
        check_cigar: args.check_cigar,
        check_mate: args.check_mate,
        tag_str: args.tags.clone(),
        tags: Vec::new(),
        records: strtoll(&args.count),
        verbose: args.verbose > 0,
        show_pass: args.show_qc || args.bamseqchksum,
        show_fail: args.show_qc,
        show_combine: true,
        tabs: args.tabs,
        compat: args.bamseqchksum,
    };

    if args.all {
        options.require_flags = 0;
        options.exclude_flags = 0;
        options.flag_mask = -1;
        options.rev_comp = false;
        options.in_order = 1;
        options.check_pos = true;
        options.check_cigar = true;
        options.check_mate = true;
        options.sanitize = Sanitize(Sanitize::ALL | Sanitize::CIGARX);
        options.tag_str = "*,cF,MD,NM".into();
    }

    match parse_tags(&options.tag_str) {
        Some(tags) => options.tags = tags,
        // samtools reports the error and carries on with no tags.
        None => eprintln!("Bad tag string.  Should be XX,YY,... syntax"),
    }

    let mut out: Box<dyn Write> = match &args.output {
        Some(path) => Box::new(io::BufWriter::new(
            File::create(path).with_context(|| format!("{}", path.display()))?,
        )),
        None => Box::new(io::BufWriter::new(io::stdout().lock())),
    };

    let result = if args.merge {
        merge(&mut options, &args.inputs, &mut out)
    } else {
        let inputs = if args.inputs.is_empty() {
            vec!["-".to_string()]
        } else {
            args.inputs.clone()
        };

        let read_options = ReadOptions {
            threads: args.threads,
            reference: args.reference.clone(),
        };

        inputs.iter().try_for_each(|input| {
            let report = checksum(&options, &read_options, input)?;
            report.write(input, &options, &mut out)?;
            Ok(())
        })
    };

    out.flush()?;
    result.context("[checksum] Failed to process data")
}

/// The checksums of one input: overall, without a read group, and by read group.
#[derive(Default)]
struct Report {
    all: Sums,
    no_read_group: Sums,
    read_groups: HashMap<Vec<u8>, Sums>,
}

fn checksum(options: &Options, read_options: &ReadOptions, src: &str) -> anyhow::Result<Report> {
    let tag_index = tag_index(options)?;

    let (mut reader, _) = AlignmentReader::open(src, read_options)
        .with_context(|| format!("[checksum] Cannot open input file \"{src}\""))?;
    let header = reader.read_header()?;
    let lengths: Vec<usize> = header
        .reference_sequences()
        .values()
        .map(|rs| rs.length().get())
        .collect();

    let mut report = Report::default();
    let mut remaining = options.records;
    let mut encoder = RawEncoder::new();
    let mut record = Vec::new();
    let mut buffers = Buffers::default();
    let selection = if options.tag_str.starts_with('*') {
        TagSelection::AllExceptListed
    } else {
        TagSelection::Listed
    };

    reader.for_each_while(&header, |rec| {
        record.clear();

        record.extend_from_slice(rec.encode_raw(&header, &mut encoder)?);

        let flags = Record::new(&record)?.flags();

        if flags & options.exclude_flags != 0
            || flags & options.require_flags != options.require_flags
        {
            return Ok(true);
        }

        if options.sanitize.0 != 0 {
            sanitize_raw(&mut record, &lengths, options.sanitize)?;
        }

        let r = Record::new(&record)?;
        let flags = r.flags();
        let (crcs, read_group) = read_crcs(&r, options, selection, &tag_index, &mut buffers)?;

        let qc = if flags & FLAG_QC_FAIL != 0 {
            Qc::Fail
        } else {
            Qc::Pass
        };

        report.add(qc, read_group, &crcs, options);

        if remaining != 0 {
            remaining -= 1;

            if remaining == 0 {
                return Ok(false);
            }
        }

        Ok(true)
    })?;

    Ok(report)
}

/// The tag positions in the requested list (`tag_keep`), checked as samtools does.
fn tag_index(options: &Options) -> anyhow::Result<HashMap<[u8; 2], usize>> {
    let mut tag_index: HashMap<[u8; 2], usize> = HashMap::new();

    for (i, tag) in options.tags.iter().enumerate() {
        if tag[0] == b'*' {
            continue;
        }

        if !(is_tag_char(tag[0]) && is_tag_char(tag[1])) {
            bail!(
                "[checksum] Illegal tag ID '{}'",
                String::from_utf8_lossy(&tag[..2])
            );
        }

        tag_index.insert([tag[0], tag[1]], i + 1);
    }

    Ok(tag_index)
}

/// Scratch buffers for [`read_crcs`].
#[derive(Default)]
struct Buffers {
    seq: Vec<u8>,
    qual: Vec<u8>,
    aux: Vec<u8>,
}

/// One read's checksums, and its read group (the last `RG:Z` value among the selected tags).
fn read_crcs(
    r: &Record<'_>,
    options: &Options,
    selection: TagSelection,
    tag_index: &HashMap<[u8; 2], usize>,
    buffers: &mut Buffers,
) -> anyhow::Result<(Crcs, Option<Vec<u8>>)> {
    let flags = r.flags();
    let masked = (i32::from(flags) & options.flag_mask) as u8;

    // The strand to undo: reverse-strand reads are reverse complemented unless `-c`.
    let strand = if options.rev_comp {
        Strand::of_flags(flags)
    } else {
        Strand::Forward
    };
    fill_seq_qual(r, strand, &mut buffers.seq, &mut buffers.qual);

    let crc = crc32(0, &[masked]);
    let crc_seq = crc32(crc, &buffers.seq);

    // samtools includes the name's NUL terminator.
    let name = &r.as_bytes()[field::NAME..=field::NAME + r.name().len()];
    let crc = crc32(0, name);
    let crc = crc32(crc, &[masked]);
    let crc_name = crc32(crc, &buffers.seq);

    let crc_qual = crc32(crc_seq, &buffers.qual);

    let read_group = select_tags(
        r.data(),
        selection,
        tag_index,
        options.tags.len(),
        &mut buffers.aux,
    )?;
    let crc_aux = crc32(crc_seq, &buffers.aux);

    let mut crcs = Crcs {
        seq: crc_seq,
        name: crc_name,
        qual: crc_qual,
        aux: crc_aux,
        pos: 0,
        cigar: 0,
        mate: 0,
    };

    if options.check_pos {
        let mut buf = Vec::with_capacity(12);
        buf.extend_from_slice(&r.reference_sequence_id().to_le_bytes());
        buf.extend_from_slice(&i64::from(r.position()).to_le_bytes());
        crcs.pos = crc32(crc_seq, &buf);
    }

    if options.check_mate {
        let mut buf = Vec::with_capacity(20);
        buf.extend_from_slice(&r.mate_reference_sequence_id().to_le_bytes());
        buf.extend_from_slice(&i64::from(r.mate_position()).to_le_bytes());
        buf.extend_from_slice(&i64::from(r.template_length()).to_le_bytes());
        crcs.mate = crc32(crc_seq, &buf);
    }

    if options.check_cigar {
        let crc = crc32(crc_seq, &u32::from(r.mapping_quality()).to_le_bytes());
        crcs.cigar = crc32(crc, r.raw_cigar());
    }

    Ok((crcs, read_group))
}

impl Report {
    /// Adds a read's checksums overall, and to its read group or the ungrouped sums.
    fn add(&mut self, qc: Qc, read_group: Option<Vec<u8>>, crcs: &Crcs, options: &Options) {
        let count = match read_group {
            Some(read_group) => {
                let sums = self.read_groups.entry(read_group).or_default();
                let count = sums.count[0];
                sums.update(qc, crcs, options, count);
                count
            }
            None => {
                let count = self.no_read_group.count[0];
                self.no_read_group.update(qc, crcs, options, count);
                count
            }
        };

        self.all.update(qc, crcs, options, count);
    }
}

/// `fill_seq_qual`: the bases (reverse complemented if asked) and qualities + 33.
fn fill_seq_qual(r: &Record<'_>, strand: Strand, seq: &mut Vec<u8>, qual: &mut Vec<u8>) {
    let len = r.sequence_len();
    let packed = r.packed_sequence();
    let quals = r.quality_scores();
    let code = |i: usize| usize::from((packed[i / 2] >> ((!i & 1) << 2)) & 0x0f);
    let q = |i: usize| quals[i].wrapping_add(33);

    seq.clear();
    qual.clear();

    if strand == Strand::Reverse {
        seq.extend((0..len).rev().map(|i| SEQ_NT16_COMPLEMENT_STR[code(i)]));
        qual.extend((0..len).rev().map(q));
    } else {
        seq.extend((0..len).map(|i| BASES[code(i)]));
        qual.extend((0..len).map(q));
    }
}

fn is_tag_char(c: u8) -> bool {
    (b'0'..=b'z').contains(&c)
}

/// `hash_aux`: the selected tags, canonicalised, into `dst` (in the requested order, or with
/// `*` in tag order); returns the read group (the last `RG:Z` value).
fn select_tags(
    data: &[u8],
    selection: TagSelection,
    tag_index: &HashMap<[u8; 2], usize>,
    n_tags: usize,
    dst: &mut Vec<u8>,
) -> io::Result<Option<Vec<u8>>> {
    let mut read_group = None;
    let mut fields: Vec<(&[u8], usize)> = Vec::new();
    let mut selected: Vec<Option<&[u8]>> = vec![None; n_tags];
    let mut rest = data;
    let mut n = 0;

    while !rest.is_empty() {
        let len = field_len(rest)?;
        let field = &rest[..len];
        rest = &rest[len..];

        if field[..3] == *b"RGZ" {
            let value = &field[3..len - 1];
            read_group = Some(value.to_vec());
        }

        if !(is_tag_char(field[0]) && is_tag_char(field[1])) {
            continue;
        }

        let index = tag_index.get(&[field[0], field[1]]).copied();

        if selection == TagSelection::AllExceptListed {
            if index.is_none() {
                fields.push((field, n));
                n += 1;
            }
        } else if let Some(i) = index {
            selected[i - 1] = Some(field);
        }
    }

    dst.clear();

    if selection == TagSelection::AllExceptListed {
        fields.sort_by_key(|&(field, n)| (field[0], field[1], n));

        for (field, _) in fields {
            push_canonical(dst, field);
        }
    } else {
        for field in selected.into_iter().flatten() {
            push_canonical(dst, field);
        }
    }

    Ok(read_group)
}

/// `canonical_tag`: integer tags in the smallest type (unsigned if not negative).
fn push_canonical(dst: &mut Vec<u8>, field: &[u8]) {
    let value: i64 = match field[2] {
        b'c' => i64::from(field[3].cast_signed()),
        b'C' => i64::from(field[3]),
        b's' => i64::from(i16::from_le_bytes([field[3], field[4]])),
        b'S' => i64::from(u16::from_le_bytes([field[3], field[4]])),
        b'i' => i64::from(i32::from_le_bytes([field[3], field[4], field[5], field[6]])),
        b'I' => i64::from(u32::from_le_bytes([field[3], field[4], field[5], field[6]])),
        _ => {
            dst.extend_from_slice(field);
            return;
        }
    };

    let code = if value >= 0 {
        if value <= 255 {
            b'C'
        } else if value <= 65535 {
            b'S'
        } else {
            b'I'
        }
    } else if value >= -128 {
        b'c'
    } else if value >= -32768 {
        b's'
    } else {
        b'i'
    };

    if code == field[2] {
        dst.extend_from_slice(field);
        return;
    }

    dst.extend_from_slice(&[field[0], field[1], code]);

    let bytes = value.to_le_bytes();

    match code {
        b'C' | b'c' => dst.push(bytes[0]),
        b'S' | b's' => dst.extend_from_slice(&bytes[..2]),
        _ => dst.extend_from_slice(&bytes[..4]),
    }
}

/// The length of the auxiliary field at the start of `src`.
fn field_len(src: &[u8]) -> io::Result<usize> {
    let invalid = || io::Error::new(io::ErrorKind::InvalidData, "invalid auxiliary field");
    let ty = *src.get(2).ok_or_else(invalid)?;
    let size = |t: u8| match t {
        b'A' | b'c' | b'C' => Some(1),
        b's' | b'S' => Some(2),
        b'i' | b'I' | b'f' => Some(4),
        b'd' => Some(8),
        _ => None,
    };

    let len = match ty {
        b'Z' | b'H' => 3 + src[3..].iter().position(|&b| b == 0).ok_or_else(invalid)? + 1,
        b'B' => {
            let subtype = *src.get(3).ok_or_else(invalid)?;
            let n = src.get(4..8).ok_or_else(invalid)?;
            let n = u32::from_le_bytes([n[0], n[1], n[2], n[3]]) as usize;
            8 + n * size(subtype).ok_or_else(invalid)?
        }
        t => 3 + size(t).ok_or_else(invalid)?,
    };

    if len > src.len() {
        return Err(invalid());
    }

    Ok(len)
}

/// `parse_tags`: tags separated by commas, each two characters or `*`.
fn parse_tags(s: &str) -> Option<Vec<Vec<u8>>> {
    let tags: Vec<Vec<u8>> = s.split(',').map(|t| t.as_bytes().to_vec()).collect();

    if s.is_empty() {
        return Some(Vec::new());
    }

    tags.iter()
        .all(|t| t.len() == 2 || t.as_slice() == b"*")
        .then_some(tags)
}

/// C's `strtoll(s, NULL, 0)`.
fn strtoll(s: &str) -> i64 {
    let s = s.trim_start();
    let (negative, s) = match s.strip_prefix('-') {
        Some(rest) => (true, rest),
        None => (false, s.strip_prefix('+').unwrap_or(s)),
    };

    let (radix, digits) = if let Some(hex) = s.strip_prefix("0x").or_else(|| s.strip_prefix("0X")) {
        (16, hex)
    } else if s.starts_with('0') {
        (8, s)
    } else {
        (10, s)
    };

    let n = digits
        .chars()
        .map_while(|c| c.to_digit(radix))
        .fold(0i64, |n, d| {
            n.saturating_mul(i64::from(radix))
                .saturating_add(i64::from(d))
        });

    if negative { -n } else { n }
}

impl Report {
    /// `checksum_report`
    fn write<W: Write>(&self, name: &str, options: &Options, out: &mut W) -> io::Result<()> {
        let mut read_groups: Vec<(&Vec<u8>, &Sums)> = self.read_groups.iter().collect();
        read_groups.sort_by(|a, b| a.0.cmp(b.0));

        if options.compat {
            let options = Options {
                tabs: true,
                show_pass: true,
                verbose: true,
                show_combine: false,
                ..options.clone()
            };

            out.write_all(
                b"###\tset\tcount\t\tb_seq\tname_b_seq\tb_seq_qual\tb_seq_tags(BC,FI,QT,RT,TC)\n",
            )?;
            write_sums(out, &options, &self.all, b"all")?;
            write_sums(out, &options, &self.no_read_group, b"")?;

            for (read_group, sums) in read_groups {
                write_sums(out, &options, sums, read_group)?;
            }

            return Ok(());
        }

        let sep = |tab: &'static str, spaces: &'static str| if options.tabs { tab } else { spaces };

        writeln!(out, "# Checksum 1.0 for file:{}{name}", sep("\t", " "))?;
        writeln!(
            out,
            "# Aux tags:{}{}",
            sep("\t", "          "),
            options.tag_str
        )?;
        writeln!(
            out,
            "# BAM flags:{}{}",
            sep("\t", "         "),
            flags::names(options.flag_mask)
        )?;

        if options.tabs {
            write!(out, "\n# Group\tQC\tcount\tflag+seq\t+name\t+qual\t+aux")?;
        } else {
            write!(
                out,
                "\n# Group    QC          count  flag+seq  +name     +qual     +aux    "
            )?;
        }

        if options.check_pos {
            out.write_all(sep("\t+chr/pos", "  +chr/pos").as_bytes())?;
        }

        if options.check_cigar {
            out.write_all(sep("\t+cigar", "  +cigar  ").as_bytes())?;
        }

        if options.check_mate {
            out.write_all(sep("\t+mate", "  +mate   ").as_bytes())?;
        }

        out.write_all(sep("\tcombined\n", "  combined\n").as_bytes())?;

        write_sums(out, options, &self.all, b"all")?;

        if options.verbose || self.no_read_group.count[0] + self.no_read_group.count[1] != 0 {
            write_sums(out, options, &self.no_read_group, b"-")?;
        }

        for (read_group, sums) in read_groups {
            write_sums(out, options, sums, read_group)?;
        }

        Ok(())
    }
}

/// `sums_report`: one group's lines.
fn write_sums<W: Write>(out: &mut W, options: &Options, sums: &Sums, set: &[u8]) -> io::Result<()> {
    const QC: [&str; 3] = ["all", "pass", "fail"];

    for (r, qc) in QC.iter().enumerate() {
        if (r == 1 && !options.show_pass) || (r == 2 && !options.show_fail) {
            continue;
        }

        if !options.verbose && sums.count[r] == 0 {
            continue;
        }

        if options.tabs {
            out.write_all(set)?;
            write!(
                out,
                "\t{}\t{}\t{}{:x}\t{:x}\t{:x}\t{:x}",
                qc,
                sums.count[r],
                if options.compat { "\t" } else { "" },
                sums.seq[r],
                sums.name[r],
                sums.qual[r],
                sums.aux[r]
            )?;

            for (check, value) in [
                (options.check_pos, sums.pos[r]),
                (options.check_cigar, sums.cigar[r]),
                (options.check_mate, sums.mate[r]),
            ] {
                if check {
                    write!(out, "\t{value:x}")?;
                }
            }
        } else {
            // `%-10s`
            out.write_all(set)?;
            out.write_all(&b"          "[..10usize.saturating_sub(set.len())])?;
            write!(
                out,
                " {:<4} {:>12}  {:08x}  {:08x}  {:08x}  {:08x}",
                qc, sums.count[r], sums.seq[r], sums.name[r], sums.qual[r], sums.aux[r]
            )?;

            for (check, value) in [
                (options.check_pos, sums.pos[r]),
                (options.check_cigar, sums.cigar[r]),
                (options.check_mate, sums.mate[r]),
            ] {
                if check {
                    write!(out, "  {value:08x}")?;
                }
            }
        }

        // samtools combines the sequence hash twice and omits the qualities.
        let mut combined = 1;
        combined = update_hash(combined, (sums.count[r] >> 32) as u32);
        combined = update_hash(combined, sums.count[r] as u32);
        combined = update_hash(combined, sums.seq[r] as u32);
        combined = update_hash(combined, sums.name[r] as u32);
        combined = update_hash(combined, sums.seq[r] as u32);
        combined = update_hash(combined, sums.aux[r] as u32);

        for (check, value) in [
            (options.check_pos, sums.pos[r]),
            (options.check_cigar, sums.cigar[r]),
            (options.check_mate, sums.mate[r]),
        ] {
            if check {
                combined = update_hash(combined, value as u32);
            }
        }

        if options.show_combine {
            if options.tabs {
                writeln!(out, "\t{combined:x}")?;
            } else {
                writeln!(out, "  {combined:08x}")?;
            }
        } else {
            writeln!(out)?;
        }
    }

    Ok(())
}

/// Which column of a checksum report a field is.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
enum Column {
    Group,
    Qc,
    Count,
    Seq,
    Name,
    Qual,
    Aux,
    Pos,
    Cigar,
    Mate,
    Combined,
}

/// `combine`: merges checksum reports, as if their inputs had been merged.
fn merge<W: Write>(options: &mut Options, inputs: &[String], out: &mut W) -> anyhow::Result<()> {
    let mut report = Report::default();
    let mut tag_str: Option<String> = None;

    for input in inputs {
        parse_report(options, &mut tag_str, input, &mut report)
            .with_context(|| format!("Failed to parse checksum file '{input}'"))?;
    }

    options.tag_str = tag_str.unwrap_or_default();
    report.write("merge", options, out)?;

    Ok(())
}

/// `sums_parse`: adds one report's groups to `report`.
fn parse_report(
    options: &mut Options,
    tag_str: &mut Option<String>,
    src: &str,
    report: &mut Report,
) -> anyhow::Result<()> {
    let file = File::open(src).with_context(|| src.to_string())?;
    let mut columns: Vec<Column> = Vec::new();
    let mut have_header = false;
    let mut bambam = false;
    let mut crcs = Crcs {
        seq: 1,
        name: 1,
        qual: 1,
        aux: 1,
        pos: 1,
        cigar: 1,
        mate: 1,
    };

    for line in BufReader::new(file).lines() {
        let line = line?;

        if let Some(rest) = line.strip_prefix("# Checksum") {
            let mut parts = rest.trim_start().split('.');
            let major = parts.next().and_then(|s| s.parse::<i32>().ok());
            let minor = parts
                .next()
                .and_then(|s| s.split(|c: char| !c.is_ascii_digit()).next())
                .and_then(|s| s.parse::<i32>().ok());

            if let (Some(major), Some(minor)) = (major, minor)
                && (major != 1 || minor != 0)
            {
                bail!("Unsupported checksum output version");
            }

            continue;
        }

        if line.starts_with("# Group") {
            if have_header {
                bail!("Invalid header, repeat of header line");
            }

            have_header = true;

            columns = samtools_columns(&line, options)?;

            continue;
        }

        if let Some(rest) = line.strip_prefix("# Aux tags:") {
            let tags = rest.trim_start();

            if !tags.is_empty() {
                set_tags(tag_str, tags)?;
            }

            continue;
        }

        if let Some(rest) = line.strip_prefix("# BAM flags:") {
            let value = rest.trim_start();

            if !value.is_empty() {
                options.flag_mask = flags::parse(value).unwrap_or(-1);
            }

            continue;
        }

        if line.starts_with("###\t") {
            if have_header {
                bail!("Invalid header, repeat of header line");
            }

            have_header = true;
            bambam = true;

            columns = bambam_columns(&line, tag_str)?;

            continue;
        }

        if line.is_empty() || line.starts_with('#') {
            continue;
        }

        if line.starts_with("all ") || line.starts_with("all\t") {
            continue;
        }

        add_row(report, &line, &columns, bambam, &mut crcs)?;
    }

    Ok(())
}

/// The columns of a samtools report, from its `# Group` header line, noting which optional
/// checksums it has.
fn samtools_columns(line: &str, options: &mut Options) -> anyhow::Result<Vec<Column>> {
    let mut columns = Vec::new();

    for token in line[2..].split_ascii_whitespace() {
        columns.push(match token {
            "Group" => Column::Group,
            "QC" => Column::Qc,
            "count" => Column::Count,
            "flag+seq" => Column::Seq,
            "+name" => Column::Name,
            "+qual" => Column::Qual,
            "+aux" => Column::Aux,
            "+chr/pos" => {
                options.check_pos = true;
                Column::Pos
            }
            "+cigar" => {
                options.check_cigar = true;
                Column::Cigar
            }
            "+mate" => {
                options.check_mate = true;
                Column::Mate
            }
            "combined" => Column::Combined,
            t => bail!("Unrecognised header token '{t}'"),
        });
    }

    Ok(columns)
}

/// The columns of a biobambam `bamseqchksum` report, from its `###` header line.
fn bambam_columns(line: &str, tag_str: &mut Option<String>) -> anyhow::Result<Vec<Column>> {
    let mut columns = Vec::new();

    for token in line.split_ascii_whitespace() {
        columns.push(match token {
            "###" => Column::Group,
            "set" => Column::Qc,
            "count" => Column::Count,
            "b_seq" => Column::Seq,
            "name_b_seq" => Column::Name,
            "b_seq_qual" => Column::Qual,
            t if t.starts_with("b_seq_tags(") => {
                let tags = &t[11..];
                let tags = tags
                    .strip_suffix(')')
                    .unwrap_or(&tags[..tags.len().saturating_sub(1)]);
                set_tags(tag_str, tags)?;
                Column::Aux
            }
            t => bail!("Unrecognised header token '{t}'"),
        });
    }

    Ok(columns)
}

/// Records the aux tags of a report, which must match those of the others.
fn set_tags(tag_str: &mut Option<String>, tags: &str) -> anyhow::Result<()> {
    match tag_str {
        None => *tag_str = Some(tags.to_string()),
        Some(t) if t != tags => bail!("Aux tag mismatch, {t} - {tags}"),
        Some(_) => {}
    }

    Ok(())
}

/// Adds one data line of a report to `report`. `crcs` keeps the values of columns the report
/// lacks.
fn add_row(
    report: &mut Report,
    line: &str,
    columns: &[Column],
    bambam: bool,
    crcs: &mut Crcs,
) -> anyhow::Result<()> {
    let mut fields: Vec<&str> = Vec::new();

    if bambam && line.starts_with('\t') {
        fields.push("");
    }

    fields.extend(line.split_ascii_whitespace().take(11 - fields.len()));

    let min_fields = if bambam { 7 } else { 8 };

    if fields.len() < min_fields || fields.len() != columns.len() {
        bail!("Incorrect number of columns in line: {line}");
    }

    let mut qc = 0;
    let mut count = 0;
    let hex = |s: &str| u32::from_str_radix(s, 16).unwrap_or(0);

    for (&column, &field) in columns.iter().zip(&fields) {
        match column {
            Column::Qc => {
                qc = match field {
                    "all" => 0,
                    "pass" => 1,
                    "fail" => 2,
                    _ => bail!("invalid QC column: {field}"),
                };
            }
            Column::Count => count = field.parse().unwrap_or(0),
            Column::Seq => crcs.seq = hex(field),
            Column::Name => crcs.name = hex(field),
            Column::Qual => crcs.qual = hex(field),
            Column::Aux => crcs.aux = hex(field),
            Column::Pos => crcs.pos = hex(field),
            Column::Cigar => crcs.cigar = hex(field),
            Column::Mate => crcs.mate = hex(field),
            Column::Group | Column::Combined => {}
        }
    }

    let group = fields[0];

    if group == "-" || group.is_empty() {
        report.no_read_group.update_row(qc, crcs, 0, count);
    } else {
        report
            .read_groups
            .entry(group.as_bytes().to_vec())
            .or_default()
            .update_row(qc, crcs, 0, count);
    }

    report.all.update_row(qc, crcs, 0, count);

    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_update_hash() {
        assert_eq!(update_hash(1, 0), 1);
        assert_eq!(update_hash(5, 3), 15);
        assert_eq!(update_hash(2, PRIME as u32), 2);
    }

    #[test]
    fn test_push_canonical() {
        let mut dst = Vec::new();
        push_canonical(&mut dst, b"XAi\x05\0\0\0");
        assert_eq!(dst, b"XAC\x05");
    }

    #[test]
    fn test_parse_tags() {
        assert_eq!(
            parse_tags("BC,QT"),
            Some(vec![b"BC".to_vec(), b"QT".to_vec()])
        );
        assert_eq!(parse_tags("*,MD").map(|t| t.len()), Some(2));
        assert_eq!(parse_tags("BCD"), None);
    }
}
