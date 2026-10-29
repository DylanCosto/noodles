//! `ntools consensus`: a consensus sequence from aligned reads (samtools consensus).
//!
//! A port of `bam_consensus.c` and `consensus_pileup.c` from samtools 1.24. The output matches
//! samtools byte for byte. Like samtools, `-@` with an indexed input splits each reference into
//! blocks (`-Z`) that are piled up independently on the worker threads; this can give slightly
//! different output from the serial path at block boundaries, so both paths are kept, and each
//! matches samtools run the same way.

mod caller;
mod options;
mod pileup;
mod region;
mod tables;

use std::{
    collections::HashMap,
    fs::File,
    io::{self, Write},
    path::Path,
    process::ExitCode,
    sync::{Arc, Mutex},
    thread,
};

use anyhow::{Context, bail};
use noodles_bam as bam;
use noodles_sam::{self as sam, alignment::io::Write as _};

use self::{
    caller::{Caller, ConsProbs, Mode},
    options::{Format, Options, Parsed, USAGE},
    pileup::{Column, Entry, Handler, POS_MAX, Pileup},
};
use crate::{
    alignment::{AlignmentIndex, AlignmentReader, ReadOptions, Rec},
    reference::IndexedFasta,
};

/// Arguments, parsed as samtools parses them (see [`options`]).
#[derive(Debug, clap::Args)]
pub struct Args {
    /// Options and the input file, as for `samtools consensus`.
    #[arg(trailing_var_arg = true, allow_hyphen_values = true, num_args = 0..)]
    args: Vec<String>,
}

/// The output buffer size above which it is written out.
const FLUSH_SIZE: usize = 1 << 20;

/// Runs `consensus`: prints the usage for no input file (or `--help`), and fails with the usage for
/// invalid options, as samtools does.
pub fn run(args: Args) -> anyhow::Result<ExitCode> {
    if matches!(args.args.as_slice(), [a] if a == "--help") {
        print!("{USAGE}");
        return Ok(ExitCode::SUCCESS);
    }

    let opts = match options::parse(&args.args) {
        Ok(Parsed::Run(opts)) => opts,
        Ok(Parsed::Usage) => {
            print!("{USAGE}");
            return Ok(ExitCode::SUCCESS);
        }
        Err(e) => {
            eprintln!("[consensus] {e:#}");
            eprint!("{USAGE}");
            return Ok(ExitCode::FAILURE);
        }
    };

    consensus(&opts)?;
    Ok(ExitCode::SUCCESS)
}

/// Shared, read-only state.
struct Env<'a> {
    opts: &'a Options,
    caller: Caller,
    names: Vec<Vec<u8>>,
    lengths: Vec<i64>,
    refs: Option<RefStore>,
}

impl Env<'_> {
    /// `sam_hdr_tid2name`
    fn name(&self, tid: i32) -> Option<&[u8]> {
        usize::try_from(tid)
            .ok()
            .and_then(|i| self.names.get(i))
            .map(Vec::as_slice)
    }

    /// `sam_hdr_tid2len`
    fn len(&self, tid: i32) -> i64 {
        usize::try_from(tid)
            .ok()
            .and_then(|i| self.lengths.get(i))
            .copied()
            .unwrap_or(0)
    }
}

/// Whole reference sequences from the `-T` FASTA, shared by the threads.
struct RefStore {
    fasta: IndexedFasta,
    cache: Mutex<Vec<(i32, Arc<Vec<u8>>)>>,
}

impl RefStore {
    /// The number of sequences kept.
    const CAPACITY: usize = 2;

    /// `fai_fetch64` of a whole sequence; `None` if the FASTA does not have it.
    fn fetch(&self, tid: i32, name: &[u8]) -> anyhow::Result<Option<Arc<Vec<u8>>>> {
        let mut cache = self.cache.lock().map_err(|_| anyhow::anyhow!("poisoned"))?;

        if let Some((_, seq)) = cache.iter().find(|(t, _)| *t == tid) {
            return Ok(Some(Arc::clone(seq)));
        }

        let mut seq = Vec::new();

        if !self.fasta.fetch(name, 0, i64::MAX, &mut seq)? {
            return Ok(None);
        }

        let seq = Arc::new(seq);

        if cache.len() >= Self::CAPACITY {
            cache.remove(0);
        }

        cache.push((tid, Arc::clone(&seq)));
        Ok(Some(seq))
    }
}

/// The byte at `i` of a reference sequence, or NUL past its end.
fn ref_byte(seq: &[u8], i: i64) -> u8 {
    usize::try_from(i)
        .ok()
        .and_then(|i| seq.get(i))
        .copied()
        .unwrap_or(0)
}

/// An iterator's region: the reference sequence and the 0-based, half-open interval.
#[derive(Clone, Copy, Debug)]
struct Iter {
    tid: i32,
    beg: i64,
    end: i64,
}

/// Per-job state (`ctx` in C) and the pileup callbacks.
struct Ctx<'e, 'a> {
    env: &'e Env<'a>,
    /// Pileup output (`ks_pileup`) and, in serial mode, all other output.
    out: Vec<u8>,
    /// Where `out` goes in serial mode.
    sink: Option<&'e mut dyn Write>,
    /// Drop `ks_pileup` writes: serial pileup output with `-@` (samtools never writes it).
    drop_ks: bool,
    /// Whether `empty_pileup2` calls from the pileup callback write to `ks_pileup`.
    threaded: bool,
    seq: Vec<u8>,
    qual: Vec<u8>,
    ks_ins_start: i64,
    first_pos: i64,
    last_pos: i64,
    last_tid: i32,
    iter: Option<Iter>,
    /// `c->ref`
    cur_ref: Option<Arc<Vec<u8>>>,
    /// The thread data reference cache (`tdata->ref_tid`, `tdata->ref`).
    tdata: (i32, Option<Arc<Vec<u8>>>),
    /// Scratch space for [`Caller::init_read`].
    bases: Vec<u8>,
}

impl<'e, 'a> Ctx<'e, 'a> {
    fn new(env: &'e Env<'a>, iter: Option<Iter>) -> Self {
        Self {
            env,
            out: Vec::new(),
            sink: None,
            drop_ks: false,
            threaded: env.opts.threads > 0,
            seq: Vec::new(),
            qual: Vec::new(),
            ks_ins_start: -1,
            first_pos: 0,
            last_pos: -1,
            last_tid: -1,
            iter,
            cur_ref: None,
            tdata: (-1, None),
            bases: Vec::new(),
        }
    }

    fn flush(&mut self) -> io::Result<()> {
        if let Some(sink) = &mut self.sink {
            sink.write_all(&self.out)?;
            self.out.clear();
        }

        Ok(())
    }

    fn maybe_flush(&mut self) -> io::Result<()> {
        if self.out.len() >= FLUSH_SIZE {
            self.flush()?;
        }

        Ok(())
    }

    /// `update_ref`: makes `cur_ref` the sequence of `tid`. Returns its length, 0 without a
    /// reference, or -1 if it cannot be fetched.
    fn update_ref(&mut self, tid: i32) -> anyhow::Result<i64> {
        let Some(store) = &self.env.refs else {
            return Ok(0);
        };

        if tid == self.tdata.0
            && let Some(seq) = &self.tdata.1
        {
            self.cur_ref = Some(Arc::clone(seq));
            return Ok(seq.len() as i64);
        }

        self.tdata = (tid, None);
        self.cur_ref = None;

        let Some(name) = self.env.name(tid) else {
            return Ok(-1);
        };

        let Some(seq) = store.fetch(tid, name)? else {
            return Ok(-1);
        };

        let len = seq.len() as i64;
        self.tdata.1 = Some(Arc::clone(&seq));
        self.cur_ref = Some(seq);
        Ok(len)
    }

    /// `empty_pileup2`: pileup lines for uncovered positions `[start, end)`, with the reference
    /// base if there is one. Returns whether the reference was unavailable (with `-T`).
    fn empty_pileup(
        &mut self,
        tid: i32,
        threaded: bool,
        start: i64,
        end: i64,
    ) -> anyhow::Result<bool> {
        let mut err = false;
        let mut rseq = None;

        if self.env.refs.is_some() {
            err = self.update_ref(tid)? <= 0;

            if !err {
                rseq = self.cur_ref.clone();
            }
        }

        if threaded && self.drop_ks {
            return Ok(err);
        }

        let name = self.env.name(tid).unwrap_or(b"(null)");
        let mut i = start;

        while i < end {
            self.out.extend_from_slice(name);
            self.out.push(b'\t');
            push_int(&mut self.out, i + 1);
            self.out.extend_from_slice(b"\t0\t0\t");
            self.out
                .push(rseq.as_deref().map_or(b'N', |seq| ref_byte(seq, i)));
            self.out.extend_from_slice(b"\t0\t*\t*\n");
            i += 1;

            if self.out.len() >= FLUSH_SIZE {
                self.flush()?;
            }
        }

        Ok(err)
    }

    /// `basic_pileup`: one pileup line per column.
    fn basic_pileup(&mut self, column: &Column<'_>) -> anyhow::Result<()> {
        let opts = self.env.opts;
        let (tid, pos, nth) = (column.tid, column.pos, column.nth);

        if !opts.show_ins && nth != 0 {
            return Ok(());
        }

        if let Some(it) = self.iter
            && (it.beg >= pos || it.end < pos)
        {
            return Ok(());
        }

        if opts.all_bases > 0 {
            if tid != self.last_tid && self.last_tid >= -1 {
                if self.last_tid >= 0 {
                    // The remainder of the previous reference.
                    let mut len = self.env.len(self.last_tid);

                    if let Some(it) = self.iter {
                        len = len.min(it.end);
                    }

                    if self.empty_pileup(self.last_tid, self.threaded, self.last_pos, len)? {
                        bail!("failed to fetch the reference");
                    }
                }

                self.last_pos = self.iter.map_or(0, |it| it.beg);
            }

            // Any references between the last and this one.
            if self.iter.is_none() && tid > self.last_tid && opts.all_bases > 1 {
                loop {
                    self.last_tid += 1;

                    if self.last_tid >= tid {
                        break;
                    }

                    let len = self.env.len(self.last_tid);
                    if self.empty_pileup(self.last_tid, false, 0, len)? {
                        bail!("failed to fetch the reference");
                    }
                }
            }

            // Any gaps in this reference.
            let gap = if self.last_pos >= 0 && pos > self.last_pos + 1 {
                Some(self.last_pos)
            } else if self.last_pos < 0 {
                Some(self.iter.map_or(0, |it| it.beg))
            } else {
                None
            };

            if let Some(start) = gap
                && self.empty_pileup(tid, self.threaded, start, pos - 1)?
            {
                bail!("failed to fetch the reference");
            }
        }

        let (cb, cq) = self.env.caller.call(column);

        if !opts.show_del && cb == b'*' {
            return Ok(());
        }

        if !self.drop_ks {
            let out = &mut self.out;
            out.extend_from_slice(self.env.name(tid).unwrap_or(b"(null)"));
            out.push(b'\t');
            push_int(out, i64::from(pos as i32));
            out.push(b'\t');
            push_int(out, i64::from(nth));
            out.push(b'\t');
            push_int(out, i64::from(column.depth));
            out.push(b'\t');
            out.push(cb);
            out.push(b'\t');
            push_int(out, i64::from(cq));
            out.push(b'\t');

            let start = out.len();
            let depth = column.depth as usize;
            out.resize(start + 2 * depth + 2, 0);

            let (bases, quals) = out[start..].split_at_mut(depth + 1);

            for (k, p) in column.reads().enumerate() {
                bases[k] = if p.is_rev {
                    if p.base == b'*' {
                        b'#'
                    } else {
                        p.base.to_ascii_lowercase()
                    }
                } else {
                    p.base
                };

                quals[k] = (p.qual.min(93) + 33) as u8;
            }

            bases[depth] = b'\t';
            quals[depth] = b'\n';
        }

        self.last_pos = pos;
        self.last_tid = tid;

        Ok(())
    }

    /// Appends `n` bases of the current reference from `from`, or `N`s, with qualities.
    fn fill(&mut self, from: i64, n: i64, ref_qual: u8, plain_qual: u8) {
        let n = usize::try_from(n).unwrap_or(0);

        match &self.cur_ref {
            Some(seq) => {
                self.seq
                    .extend((0..n as i64).map(|i| ref_byte(seq, from + i)));
                self.qual.extend(std::iter::repeat_n(ref_qual, n));
            }
            None => {
                self.seq.extend(std::iter::repeat_n(b'N', n));
                self.qual.extend(std::iter::repeat_n(plain_qual, n));
            }
        }
    }

    /// `basic_fasta`: appends each column's consensus base to the sequence.
    fn basic_fasta(&mut self, column: &Column<'_>) -> anyhow::Result<()> {
        let opts = self.env.opts;
        let (tid, pos, nth) = (column.tid, column.pos, column.nth);
        let ref_qual = (opts.ref_qual + 33) as u8;

        if !opts.show_ins && nth != 0 {
            return Ok(());
        }

        if let Some(it) = self.iter
            && (it.beg >= pos || it.end < pos)
        {
            return Ok(());
        }

        if self.first_pos > pos {
            self.first_pos = pos;
        }

        while tid != self.last_tid {
            if self.last_tid != -1 {
                if opts.all_bases > 0 {
                    // Fill in the remainder of the previous reference.
                    let mut n = match self.iter {
                        Some(it) => {
                            self.last_pos = self.last_pos.max(it.beg - 1);
                            it.end
                        }
                        None => POS_MAX,
                    };

                    n = n.min(self.env.len(self.last_tid)) - self.last_pos;

                    if n > 0 {
                        if self.cur_ref.is_some() && self.update_ref(self.last_tid)? < 0 {
                            bail!("failed to fetch the reference");
                        }

                        self.fill(self.last_pos, n, ref_qual, b'!');
                    }
                }

                let name = self.env.name(self.last_tid).unwrap_or_default().to_vec();
                let len = self.seq.len();
                self.dump_fastq(&name, len)?;
            }

            if self.update_ref(tid)? < 0 {
                bail!("failed to fetch the reference");
            }

            self.seq.clear();
            self.qual.clear();

            if self.iter.is_none() && opts.all_bases > 1 {
                self.last_tid += 1;

                if self.last_tid < tid {
                    self.last_pos = 0;
                    continue;
                }
            }

            self.last_tid = tid;
            self.last_pos = match (self.iter, opts.all_bases > 0) {
                (Some(it), true) => it.beg,
                (None, true) => 0,
                (_, false) => pos - 1,
            };
        }

        let (cb, cq) = self.env.caller.call(column);

        if !opts.show_del && cb == b'*' {
            self.last_pos = pos;
            self.last_tid = tid;
            return Ok(());
        }

        if opts.mark_ins && nth != 0 && cb != b'*' {
            self.seq.push(b'_');
            self.qual.push(b'_');
        }

        if pos > self.last_pos && (self.last_pos > 0 || opts.all_bases > 0) {
            if self.update_ref(tid)? < 0 {
                bail!("failed to fetch the reference");
            }

            self.fill(self.last_pos, pos - (self.last_pos + 1), ref_qual, b'!');
        }

        if (nth != 0 && opts.show_ins && cb != b'*')
            || cb != b'*'
            || (pos > self.last_pos && opts.show_del)
        {
            if self.ks_ins_start == -1 {
                self.ks_ins_start = self.seq.len() as i64;
            }

            self.seq.push(cb);
            self.qual.push((cq.min(93).wrapping_add(33)) as u8);
        }

        self.last_pos = pos;
        self.last_tid = tid;

        Ok(())
    }

    /// `dump_fastq` of the first `len` bases of the sequence.
    fn dump_fastq(&mut self, name: &[u8], len: usize) -> io::Result<()> {
        let fmt = self.env.opts.fmt;
        write_fastq(
            &mut self.out,
            fmt,
            self.env.opts.line_len,
            name,
            &self.seq[..len],
            &self.qual,
        )?;
        self.maybe_flush()
    }
}

impl Handler for Ctx<'_, '_> {
    fn init(&mut self, entry: &mut Entry) -> anyhow::Result<bool> {
        if self.env.caller.opts.mode == Mode::Simple {
            return Ok(true);
        }

        Ok(self.env.caller.init_read(entry, &mut self.bases))
    }

    fn column(&mut self, column: &Column<'_>) -> anyhow::Result<()> {
        match self.env.opts.fmt {
            Format::Pileup => self.basic_pileup(column)?,
            Format::Fasta | Format::Fastq => self.basic_fasta(column)?,
        }

        self.maybe_flush()?;
        Ok(())
    }
}

/// `dump_fastq`: a FASTA or FASTQ record, wrapped at `line_len`. Lines stop at a NUL byte, as
/// `%.*s` does.
fn write_fastq(
    out: &mut Vec<u8>,
    fmt: Format,
    line_len: usize,
    name: &[u8],
    seq: &[u8],
    qual: &[u8],
) -> io::Result<()> {
    if seq.is_empty() {
        return Ok(());
    }

    let lines = |data: &[u8], out: &mut Vec<u8>| {
        for chunk in data.chunks(line_len) {
            let end = chunk.iter().position(|&b| b == 0).unwrap_or(chunk.len());
            out.extend_from_slice(&chunk[..end]);
            out.push(b'\n');
        }
    };

    out.push(if fmt == Format::Fastq { b'@' } else { b'>' });
    out.extend_from_slice(name);
    out.push(b'\n');
    lines(seq, out);

    if fmt == Format::Fastq {
        out.extend_from_slice(b"+\n");
        let qual = qual.get(..seq.len()).unwrap_or(qual);
        lines(qual, out);
    }

    Ok(())
}

/// Appends `n` in decimal.
fn push_int(dst: &mut Vec<u8>, n: i64) {
    let mut buf = [0u8; 20];
    let mut i = buf.len();
    let mut m = n.unsigned_abs();

    loop {
        i -= 1;
        buf[i] = b'0' + (m % 10) as u8;
        m /= 10;

        if m == 0 {
            break;
        }
    }

    if n < 0 {
        dst.push(b'-');
    }

    dst.extend_from_slice(&buf[i..]);
}

/// `snprintf(name, 1024, ...)`: a region name, truncated to 1023 bytes.
fn region_name(name: &[u8], range: Option<(i64, i64)>) -> Vec<u8> {
    let mut s = name.to_vec();

    if let Some((beg, end)) = range {
        s.extend_from_slice(format!(":{beg}-{end}").as_bytes());
    }

    s.truncate(1023);
    s
}

/// Reads records and passes them to the pileup as raw BAM records.
struct RecordFeeder {
    encoder: bam::io::Writer<Vec<u8>>,
}

impl RecordFeeder {
    fn new() -> Self {
        Self {
            encoder: bam::io::Writer::from(Vec::new()),
        }
    }

    fn push<H: Handler>(
        &mut self,
        header: &sam::Header,
        rec: Rec<'_>,
        pileup: &mut Pileup,
        handler: &mut H,
    ) -> anyhow::Result<()> {
        match rec {
            Rec::Bam(r) => pileup.push(r.as_ref(), handler),
            _ => {
                self.encoder.get_mut().clear();
                self.encoder.write_alignment_record(header, rec.as_dyn())?;
                pileup.push(&self.encoder.get_ref()[4..], handler)
            }
        }
    }
}

/// Runs the consensus with the parsed options.
fn consensus(opts: &Options) -> anyhow::Result<()> {
    let caller = build_caller(opts);

    let read_options = ReadOptions {
        threads: 0,
        reference: opts.reference.clone(),
    };

    let (mut reader, path) = AlignmentReader::open(&opts.input, &read_options)
        .with_context(|| format!("Cannot open input file \"{}\"", opts.input))?;

    let refs = opts
        .reference
        .as_deref()
        .map(|path| {
            IndexedFasta::open(path).map(|fasta| RefStore {
                fasta,
                cache: Mutex::new(Vec::new()),
            })
        })
        .transpose()?;

    let header = reader
        .read_header()
        .with_context(|| format!("Failed to read header for \"{}\"", opts.input))?;

    let env = Env {
        opts,
        caller,
        names: header
            .reference_sequences()
            .keys()
            .map(|name| name.to_vec())
            .collect(),
        lengths: header
            .reference_sequences()
            .values()
            .map(|rs| rs.length().get() as i64)
            .collect(),
        refs,
    };

    let index = if opts.region.is_some() || opts.regions_file.is_some() {
        let index = path
            .as_deref()
            .and_then(|path| AlignmentIndex::read(reader.format(), path, None).ok());

        Some(
            index
                .with_context(|| format!("Cannot load index for input file \"{}\"", opts.input))?,
        )
    } else if opts.threads > 0 {
        let index = path
            .as_deref()
            .and_then(|path| AlignmentIndex::read(reader.format(), path, None).ok());

        if index.is_none() {
            eprintln!("No index: multi-threading is limited to decompression only");
            // Reopen with decompression threads.
            let threaded = ReadOptions {
                threads: opts.threads,
                ..read_options.clone()
            };
            let (r, _) = AlignmentReader::open(&opts.input, &threaded)?;
            reader = r;
            reader.read_header()?;
        }

        index
    } else {
        None
    };

    let sink: Box<dyn Write> = match &opts.output {
        Some(path) => Box::new(File::create(path).with_context(|| format!("{}", path.display()))?),
        None => Box::new(crate::output::raw_stdout()?),
    };
    let mut sink = io::BufWriter::with_capacity(1 << 20, sink);

    match (opts.threads > 0, &index, &path) {
        (true, Some(index), Some(path)) => {
            run_parallel(&env, path, &header, index, &mut sink)?;
        }
        _ => run_serial(&env, &mut reader, &header, index.as_ref(), &mut sink)?,
    }

    sink.flush()?;
    Ok(())
}

/// The probability tables for the options (the `consensus_init` calls of `main_consensus`).
fn build_caller(opts: &Options) -> Caller {
    let mode = opts.caller.mode;
    let mut recall = ConsProbs::default();
    let mut precise = ConsProbs::default();

    if mode != Mode::Simple {
        if mode == Mode::Precise {
            precise = ConsProbs::new(
                opts.p_het,
                opts.p_indel,
                0.3 * opts.het_scale,
                opts.homopoly_redux,
                &opts.qcal,
                Mode::Precise,
            );
        }

        if mode == Mode::Mixed {
            precise = ConsProbs::new(
                opts.p_het.powf(0.7),
                opts.p_indel.powf(0.7),
                0.3 * opts.het_scale,
                opts.homopoly_redux,
                &opts.qcal,
                Mode::Precise,
            );
        }

        recall = ConsProbs::new(
            opts.p_het,
            opts.p_indel,
            opts.het_scale,
            if mode == Mode::Recall {
                opts.homopoly_redux
            } else {
                0.01
            },
            &opts.qcal,
            Mode::Recall,
        );
    }

    Caller::new(opts.caller.clone(), recall, precise)
}

/// Visits the records of an iterator region, or of the whole file.
fn for_each_record<F>(
    reader: &mut AlignmentReader,
    header: &sam::Header,
    index: Option<&AlignmentIndex>,
    iter: Option<Iter>,
    lengths: &[i64],
    f: F,
) -> anyhow::Result<()>
where
    F: FnMut(Rec<'_>) -> anyhow::Result<()>,
{
    match (iter, index) {
        (Some(it), Some(index)) => {
            let tid = usize::try_from(it.tid)?;
            let len = lengths.get(tid).copied().unwrap_or(0);
            let end = if it.end <= len.max(it.beg) {
                Some(usize::try_from(it.end)?)
            } else {
                None
            };
            reader.for_each_in_interval(index, header, tid, usize::try_from(it.beg)?, end, f)
        }
        _ => reader.for_each(header, f),
    }
}

/// `pileup_loop_serial`: the whole file, a region, or each BED interval, on one thread.
fn run_serial(
    env: &Env<'_>,
    reader: &mut AlignmentReader,
    header: &sam::Header,
    index: Option<&AlignmentIndex>,
    sink: &mut dyn Write,
) -> anyhow::Result<()> {
    let opts = env.opts;
    let mut ctx = Ctx::new(env, None);
    ctx.sink = Some(sink);
    ctx.drop_ks = opts.threads > 0 && opts.fmt == Format::Pileup;

    let result = serial_loop(&mut ctx, reader, header, index);

    // Output before an error is written, as samtools has written it.
    ctx.flush()?;
    result
}

/// The body of [`run_serial`].
fn serial_loop(
    ctx: &mut Ctx<'_, '_>,
    reader: &mut AlignmentReader,
    header: &sam::Header,
    index: Option<&AlignmentIndex>,
) -> anyhow::Result<()> {
    let env = ctx.env;
    let opts = env.opts;
    let n_targets = env.names.len() as i32;
    let mut reglist = None;

    if let Some(path) = &opts.regions_file {
        reglist = Some(region::read_bed(path)?.reglist(&env.names)?.into_iter());
    } else if let Some(reg) = &opts.region {
        let (tid, beg, end) = region::parse_region(reg, &env.names)
            .with_context(|| format!("Failed to parse region \"{reg}\""))?;
        ctx.iter = Some(Iter {
            tid: tid as i32,
            beg,
            end,
        });
    }

    let mut feeder = RecordFeeder::new();

    loop {
        if let Some(list) = &mut reglist {
            ctx.seq.clear();
            ctx.qual.clear();
            ctx.last_tid = -1;
            ctx.last_pos = -1;
            ctx.iter = None;

            let Some((chr, mut start, mut end)) = list.next() else {
                break;
            };

            let len = env.len(chr as i32);

            if start > end || start > len {
                eprintln!(
                    "[consensus] Warning: Invalid region \"{}:{start}-{end}\"",
                    String::from_utf8_lossy(&env.names[chr])
                );

                if list.len() == 0 {
                    break;
                }

                continue;
            }

            start = start.max(0);
            end = end.min(len);

            ctx.iter = Some(Iter {
                tid: chr as i32,
                beg: start,
                end,
            });
            ctx.last_pos = start;
        }

        let mut pileup = Pileup::new(opts.filters);
        for_each_record(reader, header, index, ctx.iter, &env.lengths, |rec| {
            feeder.push(header, rec, &mut pileup, ctx)
        })?;
        pileup.finish(ctx)?;

        if opts.fmt == Format::Pileup {
            finish_serial_pileup(ctx, n_targets)?;
        } else {
            finish_serial_fasta(ctx, n_targets)?;
        }

        ctx.flush()?;

        match &reglist {
            Some(list) if list.len() > 0 => {}
            _ => break,
        }
    }

    Ok(())
}

/// The end of `pileup_loop_serial` for pileup output: `-a` and `-aa` padding.
fn finish_serial_pileup(ctx: &mut Ctx<'_, '_>, n_targets: i32) -> anyhow::Result<()> {
    let opts = ctx.env.opts;

    if opts.all_bases > 0 {
        let tid = ctx.iter.map_or(ctx.last_tid, |it| it.tid);
        let mut len = ctx.env.len(tid);
        let mut pos = ctx.last_pos;

        if let Some(it) = ctx.iter {
            len = len.min(it.end);
            pos = pos.max(it.beg);
        }

        if ctx.empty_pileup(tid, false, pos, len)? {
            bail!("failed to fetch the reference");
        }
    }

    while ctx.iter.is_none() && opts.all_bases > 1 {
        ctx.last_tid += 1;

        if ctx.last_tid >= n_targets {
            break;
        }

        let len = ctx.env.len(ctx.last_tid);
        if ctx.empty_pileup(ctx.last_tid, false, 0, len)? {
            bail!("failed to fetch the reference");
        }
    }

    Ok(())
}

/// The end of `pileup_loop_serial` for FASTA/FASTQ output: padding and the last record.
fn finish_serial_fasta(ctx: &mut Ctx<'_, '_>, n_targets: i32) -> anyhow::Result<()> {
    let opts = ctx.env.opts;
    let ref_qual = (opts.ref_qual + 33) as u8;

    loop {
        if opts.all_bases > 0 {
            let tid = ctx.iter.map_or(ctx.last_tid, |it| it.tid);
            let mut len = ctx.env.len(tid);
            let mut pos = ctx.last_pos;

            if let Some(it) = ctx.iter {
                len = len.min(it.end);
                pos = pos.max(it.beg);
                ctx.last_tid = it.tid;
            }

            if pos < len {
                if ctx.update_ref(ctx.last_tid)? < 0 {
                    bail!("failed to fetch the reference");
                }

                // `while (pos++ < len)` fills from `pos` (a -1 start fills one base).
                let from = pos;
                let n = len - pos;
                ctx.fill(from, n, ref_qual, b'!');
            }
        }

        if ctx.last_tid >= 0 {
            let tid = ctx.iter.map_or(ctx.last_tid, |it| it.tid);
            let len = i64::from(ctx.env.len(tid) as i32);
            let base = ctx.env.name(ctx.last_tid).unwrap_or_default();

            let name = match ctx.iter {
                Some(it) if it.beg > 0 || it.end < len => {
                    region_name(base, Some((it.beg + 1, it.end.min(len))))
                }
                _ => region_name(base, None),
            };

            let n = ctx.seq.len();
            ctx.dump_fastq(&name, n)?;
        }

        if ctx.iter.is_none() && opts.all_bases > 1 {
            ctx.last_tid += 1;

            if ctx.last_tid < n_targets {
                ctx.last_pos = 0;
                ctx.seq.clear();
                ctx.qual.clear();
                continue;
            }
        }

        return Ok(());
    }
}

/// A block of a reference sequence piled up on its own (`pileup_job`).
#[derive(Clone, Copy, Debug)]
struct Job {
    tid: i32,
    start: i64,
    end: i64,
}

/// What a job produces.
struct JobOutput {
    job: Job,
    out: Vec<u8>,
    seq: Vec<u8>,
    qual: Vec<u8>,
    ks_ins_start: i64,
    first_pos: i64,
    last_pos: i64,
}

/// The steps of `pileup_loop_parallel`, in output order.
enum Step {
    /// `-aa` output of a reference sequence with no indexed reads.
    Empty(i32),
    /// The start of a region's output.
    Start,
    /// A job's output, by job number.
    Job(usize),
    /// The end of a region: write the FASTA/FASTQ record with this name.
    End(Vec<u8>),
}

/// Whether `sam_itr_queryi` on a whole reference sequence is `finished` at once: the index has
/// no bins for it.
fn index_is_empty(index: &AlignmentIndex, tid: usize) -> bool {
    match index {
        AlignmentIndex::Bam(bam::Index::Bai(index)) => index
            .reference_sequences()
            .get(tid)
            .is_none_or(|rs| rs.bins().is_empty()),
        AlignmentIndex::Bam(bam::Index::Csi(index)) => index
            .reference_sequences()
            .get(tid)
            .is_none_or(|rs| rs.bins().is_empty()),
        AlignmentIndex::Cram(index) => !index
            .iter()
            .any(|record| record.reference_sequence_id() == Some(tid)),
    }
}

/// Lays out the steps and jobs of `pileup_loop_parallel`.
fn plan_parallel(env: &Env<'_>, index: &AlignmentIndex) -> anyhow::Result<(Vec<Step>, Vec<Job>)> {
    let opts = env.opts;
    let mut steps = Vec::new();
    let mut jobs = Vec::new();

    let mut add_region = |tid: usize, start: i64, end: i64, steps: &mut Vec<Step>| {
        let len = env.lengths[tid];
        let (start, end) = (start.max(0), end.min(len));

        let name = if start <= 0 && end >= len {
            region_name(&env.names[tid], None)
        } else {
            region_name(&env.names[tid], Some((start + 1, end)))
        };

        steps.push(Step::Start);

        let mut sub_start = start;
        while sub_start < end {
            jobs.push(Job {
                tid: tid as i32,
                start: sub_start,
                end: (sub_start + opts.span).min(end),
            });
            steps.push(Step::Job(jobs.len() - 1));
            sub_start += opts.span;
        }

        steps.push(Step::End(name));
    };

    if let Some(reg) = &opts.region {
        let (tid, start, end) = region::parse_region(reg, &env.names)
            .with_context(|| format!("Failed to parse region \"{reg}\""))?;
        add_region(tid, start, end, &mut steps);
    } else if let Some(path) = &opts.regions_file {
        for (tid, start, end) in region::read_bed(path)?.reglist(&env.names)? {
            if start > end || start > env.lengths[tid] {
                eprintln!(
                    "[consensus] Warning: Invalid region \"{}:{start}-{end}\"",
                    String::from_utf8_lossy(&env.names[tid])
                );
                continue;
            }

            add_region(tid, start, end, &mut steps);
        }
    } else {
        for tid in 0..env.names.len() {
            if index_is_empty(index, tid) {
                if opts.all_bases > 1 {
                    steps.push(Step::Empty(tid as i32));
                }
                continue;
            }

            add_region(tid, 0, env.lengths[tid], &mut steps);
        }
    }

    Ok((steps, jobs))
}

/// Runs one job: the pileup of the reads overlapping a block.
fn run_job(
    env: &Env<'_>,
    reader: &mut AlignmentReader,
    header: &sam::Header,
    index: &AlignmentIndex,
    job: Job,
) -> anyhow::Result<JobOutput> {
    let opts = env.opts;
    let iter = Iter {
        tid: job.tid,
        beg: job.start,
        end: job.end,
    };

    let mut ctx = Ctx::new(env, Some(iter));
    ctx.first_pos = POS_MAX;

    let mut pileup = Pileup::new(opts.filters);
    let mut feeder = RecordFeeder::new();

    for_each_record(
        reader,
        header,
        Some(index),
        Some(iter),
        &env.lengths,
        |rec| feeder.push(header, rec, &mut pileup, &mut ctx),
    )?;
    pileup.finish(&mut ctx)?;

    if opts.fmt == Format::Pileup && ctx.last_pos < job.end && opts.all_bases > 0 {
        let beg = job.start.max(ctx.last_pos);
        ctx.empty_pileup(job.tid, true, beg, job.end)?;
    }

    Ok(JobOutput {
        job,
        out: ctx.out,
        seq: ctx.seq,
        qual: ctx.qual,
        ks_ins_start: ctx.ks_ins_start,
        first_pos: ctx.first_pos,
        last_pos: ctx.last_pos,
    })
}

/// A region's FASTA/FASTQ sequence as it is assembled from job outputs.
struct Assembly {
    seq: Vec<u8>,
    qual: Vec<u8>,
    used_start: i64,
    used_end: i64,
}

impl Assembly {
    fn new() -> Self {
        Self {
            seq: Vec::new(),
            qual: Vec::new(),
            used_start: -1,
            used_end: 0,
        }
    }

    /// `ref_or_Ns`: appends reference bases (or `N`s) for `[pos, pos + n)`. Returns `false`
    /// if the reference cannot be fetched.
    fn ref_or_ns(
        &mut self,
        main: &mut Ctx<'_, '_>,
        tid: i32,
        pos: i64,
        n: i64,
    ) -> anyhow::Result<bool> {
        let opts = main.env.opts;
        let fastq = opts.fmt == Format::Fastq;
        let count = usize::try_from(n).unwrap_or(0);

        if main.env.refs.is_some() {
            if main.update_ref(tid)? < 0 {
                return Ok(false);
            }

            let seq = main.cur_ref.as_deref().map_or(&[][..], Vec::as_slice);
            self.seq
                .extend((0..n.max(0)).map(|i| ref_byte(seq, pos + i)));

            if fastq {
                self.qual
                    .extend(std::iter::repeat_n((opts.ref_qual + 33) as u8, count));
            }
        } else {
            self.seq.extend(std::iter::repeat_n(b'N', count));

            if fastq {
                self.qual.extend(std::iter::repeat_n(b'!', count));
            }
        }

        Ok(true)
    }

    /// `append_cons`: appends a job's sequence, padding between and within jobs.
    fn append(&mut self, main: &mut Ctx<'_, '_>, c: &JobOutput) -> anyhow::Result<()> {
        let all_bases = main.env.opts.all_bases > 0;
        let job = c.job;

        if all_bases || self.used_start >= 0 {
            let n = if c.ks_ins_start == -1 {
                job.end - job.start
            } else {
                c.first_pos - job.start - 1
            };

            if n != 0 {
                if !self.ref_or_ns(main, job.tid, job.start, n)? {
                    return Ok(());
                }

                if self.used_start == -1 {
                    self.used_start = 0;
                }

                if all_bases {
                    self.used_end = self.seq.len() as i64;
                }
            }
        }

        let skip = usize::try_from(c.ks_ins_start).unwrap_or(0);
        let seq = c.seq.get(skip..).unwrap_or_default();
        let qual = c.qual.get(skip..).unwrap_or_default();

        if !seq.is_empty() {
            self.seq.extend_from_slice(seq);
            self.qual.extend_from_slice(qual);

            let n = job.end - c.last_pos;
            self.used_end = self.seq.len() as i64;
            self.used_start = 0;

            if n != 0 {
                if !self.ref_or_ns(main, job.tid, c.last_pos, n)? {
                    return Ok(());
                }

                if all_bases {
                    self.used_end = self.seq.len() as i64;
                }
            }
        }

        Ok(())
    }
}

/// `pileup_loop_parallel`: jobs of up to `-Z` bases on the worker threads, with their output
/// assembled in order.
fn run_parallel(
    env: &Env<'_>,
    path: &Path,
    header: &sam::Header,
    index: &AlignmentIndex,
    sink: &mut dyn Write,
) -> anyhow::Result<()> {
    let opts = env.opts;
    let (steps, jobs) = plan_parallel(env, index)?;
    let threads = opts.threads.max(1);
    let window = threads * 2;

    let mut main = Ctx::new(env, None);
    main.sink = Some(sink);

    thread::scope(|scope| -> anyhow::Result<()> {
        let (job_tx, job_rx) = crossbeam_channel::unbounded::<usize>();
        let (result_tx, result_rx) = crossbeam_channel::unbounded();

        for _ in 0..threads.min(jobs.len()) {
            let job_rx = job_rx.clone();
            let result_tx = result_tx.clone();
            let jobs = &jobs;

            scope.spawn(move || {
                let mut reader = None;

                for i in job_rx {
                    let result = (|| -> anyhow::Result<JobOutput> {
                        if reader.is_none() {
                            let read_options = ReadOptions {
                                threads: 0,
                                reference: opts.reference.clone(),
                            };
                            let (mut r, _) = AlignmentReader::open(
                                path.to_str().context("invalid path")?,
                                &read_options,
                            )?;
                            r.read_header()?;
                            reader = Some(r);
                        }

                        let reader = reader.as_mut().expect("opened above");
                        run_job(env, reader, header, index, jobs[i])
                    })();

                    if result_tx.send((i, result)).is_err() {
                        break;
                    }
                }
            });
        }

        drop(result_tx);

        let mut dispatched = 0;
        let mut done: HashMap<usize, JobOutput> = HashMap::new();
        let mut assembly = Assembly::new();

        for step in &steps {
            match step {
                Step::Empty(tid) => write_empty_reference(&mut main, *tid)?,
                Step::Start => assembly = Assembly::new(),
                Step::Job(i) => {
                    while dispatched < jobs.len() && dispatched < i + window {
                        job_tx.send(dispatched)?;
                        dispatched += 1;
                    }

                    let output = loop {
                        if let Some(output) = done.remove(i) {
                            break output;
                        }

                        let (j, result) = result_rx.recv()?;
                        done.insert(j, result?);
                    };

                    if opts.fmt == Format::Pileup {
                        main.out.extend_from_slice(&output.out);
                    } else {
                        assembly.append(&mut main, &output)?;
                    }
                }
                Step::End(name) => {
                    if opts.fmt != Format::Pileup {
                        let len = usize::try_from(assembly.used_end).unwrap_or(0);
                        let seq = assembly.seq.get(..len).unwrap_or(&assembly.seq);
                        write_fastq(
                            &mut main.out,
                            opts.fmt,
                            opts.line_len,
                            name,
                            seq,
                            &assembly.qual,
                        )?;
                    }
                }
            }

            main.maybe_flush()?;
        }

        drop(job_tx);
        main.flush()?;
        Ok(())
    })
}

/// `-aa` output of a reference sequence with no reads, in the threaded code.
fn write_empty_reference(main: &mut Ctx<'_, '_>, tid: i32) -> anyhow::Result<()> {
    let opts = main.env.opts;
    let len = main.env.len(tid);

    if opts.fmt == Format::Pileup {
        if main.empty_pileup(tid, false, 0, len)? {
            bail!("failed to fetch the reference");
        }

        return Ok(());
    }

    let job = JobOutput {
        job: Job {
            tid,
            start: 0,
            end: len,
        },
        out: Vec::new(),
        seq: Vec::new(),
        qual: Vec::new(),
        ks_ins_start: -1,
        first_pos: 0,
        last_pos: 0,
    };

    let mut assembly = Assembly::new();
    assembly.append(main, &job)?;

    let name = main.env.name(tid).unwrap_or_default().to_vec();
    let len = usize::try_from(len).unwrap_or(0).min(assembly.seq.len());
    write_fastq(
        &mut main.out,
        opts.fmt,
        opts.line_len,
        &name,
        &assembly.seq[..len],
        &assembly.qual,
    )?;

    Ok(())
}
