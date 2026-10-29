//! htslib's pileup engine (`bam_plp_*` in `sam.c`, 1.24): reads are buffered in coordinate
//! order and resolved against each reference position in turn, with the same CIGAR resolution,
//! maximum depth rule, and overlapping mate quality adjustment, so that `mpileup` sees the same
//! reads and qualities at every position as samtools does.

use std::{collections::HashMap, io};

use anyhow::bail;

use crate::{
    bam_raw::{self, field, le_u16, le_u32},
    filter::{wang_hash, x31_hash},
};

const CIGAR_MATCH: u32 = 0;
const CIGAR_INSERTION: u32 = 1;
const CIGAR_DELETION: u32 = 2;
const CIGAR_SKIP: u32 = 3;
const CIGAR_SOFT_CLIP: u32 = 4;
const CIGAR_HARD_CLIP: u32 = 5;
const CIGAR_PAD: u32 = 6;
const CIGAR_SEQUENCE_MATCH: u32 = 7;
const CIGAR_SEQUENCE_MISMATCH: u32 = 8;

const FLAG_PAIRED: u16 = 0x1;
const FLAG_PROPER_PAIR: u16 = 0x2;
const FLAG_UNMAPPED: u16 = 0x4;
const FLAG_MATE_UNMAPPED: u16 = 0x8;
const FLAG_REVERSE: u16 = 0x10;

fn is_match(op: u32) -> bool {
    matches!(
        op,
        CIGAR_MATCH | CIGAR_SEQUENCE_MATCH | CIGAR_SEQUENCE_MISMATCH
    )
}

/// An alignment record being piled up: its encoded BAM bytes (whose qualities and tags may be
/// adjusted) and its decoded fixed fields.
#[derive(Clone, Debug)]
pub struct Read {
    pub raw: Vec<u8>,
    pub tid: i32,
    pub pos: i64,
    pub flag: u16,
    pub mtid: i32,
    pub mpos: i64,
    pub isize: i32,
    pub l_seq: usize,
    /// The CIGAR ops, from the CG tag when the record has a `kSmN` placeholder, as htslib reads
    /// them.
    pub cigar: Vec<u32>,
    name_end: usize,
    seq_start: usize,
    qual_start: usize,
}

impl Read {
    pub fn new(raw: Vec<u8>) -> io::Result<Self> {
        let i32_at = |offset: usize| le_u32(&raw[offset..offset + 4]).cast_signed();

        let name_len = usize::from(raw[field::NAME_LENGTH]);
        let n_cigar = usize::from(le_u16(&raw[field::CIGAR_OP_COUNT]));
        let l_seq = usize::try_from(le_u32(&raw[field::SEQUENCE_LENGTH]))
            .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))?;

        let cigar_start = field::NAME + name_len;
        let seq_start = cigar_start + 4 * n_cigar;
        let qual_start = seq_start + l_seq.div_ceil(2);

        let mut cigar: Vec<u32> = raw[cigar_start..seq_start]
            .as_chunks::<4>()
            .0
            .iter()
            .map(|b| u32::from_le_bytes(*b))
            .collect();

        if cigar.len() == 2
            && cigar[0] & 0xf == CIGAR_SOFT_CLIP
            && (cigar[0] >> 4) as usize == l_seq
            && cigar[1] & 0xf == CIGAR_SKIP
            && let Some(range) = bam_raw::find_field(&raw, *b"CG")?
            && raw[range.start + 2] == b'B'
            && raw[range.start + 3] == b'I'
        {
            cigar = raw[range.start + 8..range.end]
                .as_chunks::<4>()
                .0
                .iter()
                .map(|b| u32::from_le_bytes(*b))
                .collect();
        }

        Ok(Self {
            tid: i32_at(0),
            pos: i64::from(i32_at(4)),
            flag: le_u16(&raw[field::FLAGS]),
            mtid: i32_at(20),
            mpos: i64::from(i32_at(24)),
            isize: i32_at(28),
            l_seq,
            cigar,
            name_end: cigar_start - 1,
            seq_start,
            qual_start,
            raw,
        })
    }

    pub fn name(&self) -> &[u8] {
        &self.raw[field::NAME..self.name_end]
    }

    pub fn mapq(&self) -> u8 {
        self.raw[9]
    }

    pub fn set_mapq(&mut self, mapq: u8) {
        self.raw[9] = mapq;
    }

    pub fn is_reverse(&self) -> bool {
        self.flag & FLAG_REVERSE != 0
    }

    /// The 4-bit code of base `i`.
    pub fn base(&self, i: usize) -> u8 {
        (self.raw[self.seq_start + i / 2] >> ((!i & 1) << 2)) & 0xf
    }

    /// Sets the 4-bit code of base `i`.
    pub fn set_base(&mut self, i: usize, code: u8) {
        let byte = &mut self.raw[self.seq_start + i / 2];

        if i & 1 == 0 {
            *byte = (*byte & 0x0f) | (code << 4);
        } else {
            *byte = (*byte & 0xf0) | (code & 0x0f);
        }
    }

    pub fn qual(&self) -> &[u8] {
        &self.raw[self.qual_start..self.qual_start + self.l_seq]
    }

    pub fn qual_mut(&mut self) -> &mut [u8] {
        &mut self.raw[self.qual_start..self.qual_start + self.l_seq]
    }

    /// The byte range of the first aux field with `tag`.
    pub fn find_tag(&self, tag: [u8; 2]) -> Option<(usize, usize)> {
        bam_raw::find_field(&self.raw, tag)
            .ok()
            .flatten()
            .map(|r| (r.start, r.end))
    }

    pub fn remove_tag(&mut self, tag: [u8; 2]) {
        if let Some((start, end)) = self.find_tag(tag) {
            self.raw.drain(start..end);
        }
    }

    /// `bam_aux_append` of a Z field.
    pub fn append_string_tag(&mut self, tag: [u8; 2], value: &[u8]) {
        self.raw.extend_from_slice(&tag);
        self.raw.push(b'Z');
        self.raw.extend_from_slice(value);
        self.raw.push(0);
    }

    /// `bam_cigar2rlen`
    fn reference_length(&self) -> i64 {
        self.cigar
            .iter()
            .filter(|&&op| {
                matches!(
                    op & 0xf,
                    CIGAR_MATCH
                        | CIGAR_DELETION
                        | CIGAR_SKIP
                        | CIGAR_SEQUENCE_MATCH
                        | CIGAR_SEQUENCE_MISMATCH
                )
            })
            .map(|&op| i64::from(op >> 4))
            .sum()
    }
}

/// `cstate_t`: the CIGAR op (`k`) a read is at, and the reference (`x`) and query (`y`) positions
/// where that op starts.
#[derive(Clone, Copy, Debug)]
struct CigarState {
    k: i32,
    x: i64,
    y: i64,
    end: i64,
}

/// `lbnode_t`
struct Node {
    read: Read,
    beg: i64,
    end: i64,
    s: CigarState,
}

/// `bam_pileup1_t`: one read at the current position.
#[derive(Clone, Copy, Debug, Default)]
pub struct Entry {
    node: usize,
    pub qpos: i64,
    pub indel: i32,
    pub is_del: bool,
    pub is_head: bool,
    pub is_tail: bool,
    pub is_refskip: bool,
    pub cigar_ind: i32,
}

/// `bam_plp_t`. Buffered reads live in `slots`, and `order` lists the occupied slots in input
/// order (htslib's linked list), so that dropping a read moves only slot indices.
pub struct Pileup {
    slots: Vec<Option<Node>>,
    free: Vec<u32>,
    order: Vec<u32>,
    tid: i32,
    pos: i64,
    max_tid: i32,
    max_pos: i64,
    is_eof: bool,
    maxcnt: usize,
    /// The last record pushed, which htslib leaves in its spare node; the buffer's head when the
    /// buffer is empty.
    last: Option<(i32, i64)>,
    /// Reads whose mates may still overlap them, by name, to their slots.
    overlaps: Option<HashMap<Vec<u8>, u32>>,
    entries: Vec<Entry>,
    out_tid: i32,
    out_pos: i64,
}

impl Pileup {
    pub fn new(maxcnt: usize, overlaps: bool) -> Self {
        Self {
            slots: Vec::new(),
            free: Vec::new(),
            order: Vec::new(),
            tid: 0,
            pos: 0,
            max_tid: -1,
            max_pos: -1,
            is_eof: false,
            maxcnt,
            last: None,
            overlaps: overlaps.then(HashMap::new),
            entries: Vec::new(),
            out_tid: 0,
            out_pos: 0,
        }
    }

    /// The reads at the position [`Self::next`] last prepared.
    pub fn entries(&self) -> &[Entry] {
        &self.entries
    }

    pub fn is_eof(&self) -> bool {
        self.is_eof
    }

    pub fn position(&self) -> (i32, i64) {
        (self.out_tid, self.out_pos)
    }

    pub fn read(&self, entry: &Entry) -> &Read {
        &self.node(entry.node as u32).read
    }

    fn node(&self, slot: u32) -> &Node {
        self.slots[slot as usize]
            .as_ref()
            .expect("ordered slots are occupied")
    }

    /// `bam_plp_push`; `None` marks the end of input.
    pub fn push(&mut self, read: Option<Read>) -> anyhow::Result<()> {
        let Some(read) = read else {
            self.is_eof = true;
            return Ok(());
        };

        if read.tid < 0 || read.flag & FLAG_UNMAPPED != 0 {
            self.overlap_remove(&read);
            return Ok(());
        }

        // The buffer and htslib's spare node.
        if self.tid == read.tid && self.pos == read.pos && self.order.len() + 1 > self.maxcnt {
            self.overlap_remove(&read);
            return Ok(());
        }

        let beg = read.pos;
        let end = read.pos + read.reference_length();

        if read.tid < self.max_tid {
            bail!("The input is not sorted (chromosomes out of order)");
        }

        if read.tid == self.max_tid && beg < self.max_pos {
            bail!("The input is not sorted (reads out of order)");
        }

        self.max_tid = read.tid;
        self.max_pos = beg;
        self.last = Some((read.tid, beg));

        if end > self.pos || read.tid > self.tid {
            let mut node = Node {
                read,
                beg,
                end,
                s: CigarState {
                    k: -1,
                    x: 0,
                    y: 0,
                    end: end - 1,
                },
            };

            let slot = match self.free.pop() {
                Some(slot) => slot,
                None => {
                    self.slots.push(None);
                    u32::try_from(self.slots.len() - 1).expect("depth fits in u32")
                }
            };

            self.overlap_push(&mut node, slot)?;
            self.slots[slot as usize] = Some(node);
            self.order.push(slot);
        }

        Ok(())
    }

    /// `bam_plp64_next`: prepares the next position with reads, if the buffer allows; returns
    /// whether it did.
    pub fn next(&mut self) -> anyhow::Result<bool> {
        if self.is_eof && self.order.is_empty() {
            return Ok(false);
        }

        while self.is_eof
            || self.max_tid > self.tid
            || (self.max_tid == self.tid && self.max_pos > self.pos)
        {
            self.entries.clear();

            let (tid, pos) = (self.tid, self.pos);
            let mut kept = 0;

            // Drop reads that ended before this position, and pile up the rest that started.
            for i in 0..self.order.len() {
                let slot = self.order[i];
                let node = self.slots[slot as usize]
                    .as_mut()
                    .expect("ordered slots are occupied");

                if node.read.tid < tid || (node.read.tid == tid && node.end <= pos) {
                    let node = self.slots[slot as usize].take().expect("occupied");

                    if let Some(overlaps) = &mut self.overlaps
                        && node.read.flag & FLAG_UNMAPPED == 0
                        && node.read.flag & FLAG_PROPER_PAIR != 0
                    {
                        overlaps.remove(node.read.name());
                    }

                    self.free.push(slot);
                    continue;
                }

                if node.read.tid == tid && node.beg <= pos {
                    let entry = resolve_cigar(&node.read, pos, &mut node.s, slot as usize);
                    self.entries.push(entry);
                }

                self.order[kept] = slot;
                kept += 1;
            }

            self.order.truncate(kept);

            self.out_tid = tid;
            self.out_pos = pos;

            let head = match self.order.first() {
                Some(&slot) => {
                    let node = self.node(slot);

                    if self.tid > node.read.tid {
                        bail!("Unsorted input. Pileup aborts");
                    }

                    Some((node.read.tid, node.beg))
                }
                None => self.last,
            };

            let (head_tid, head_beg) = head.unwrap_or((0, 0));

            if self.tid < head_tid {
                self.tid = head_tid;
                self.pos = head_beg;
            } else if self.pos < head_beg {
                self.pos = head_beg;
            } else {
                self.pos += 1;
            }

            if !self.entries.is_empty() {
                return Ok(true);
            }

            if self.is_eof && self.order.is_empty() {
                break;
            }
        }

        Ok(false)
    }

    fn overlap_remove(&mut self, read: &Read) {
        self.overlap_remove_name(read.flag, read.name());
    }

    fn overlap_remove_name(&mut self, flag: u16, name: &[u8]) {
        let Some(overlaps) = &mut self.overlaps else {
            return;
        };

        if flag & FLAG_UNMAPPED != 0 || flag & FLAG_PROPER_PAIR == 0 {
            return;
        }

        overlaps.remove(name);
    }

    /// `overlap_push`
    fn overlap_push(&mut self, node: &mut Node, slot: u32) -> anyhow::Result<()> {
        let Some(overlaps) = &mut self.overlaps else {
            return Ok(());
        };

        let b = &node.read;

        if b.flag & FLAG_MATE_UNMAPPED != 0 || b.flag & FLAG_PROPER_PAIR == 0 {
            return Ok(());
        }

        if (b.mtid >= 0 && b.tid != b.mtid)
            || (i64::from(b.isize).abs() >= 2 * b.l_seq as i64 && b.mpos >= node.end)
        {
            return Ok(());
        }

        match overlaps.get(b.name()) {
            None => {
                if b.mpos >= b.pos || (b.flag & FLAG_PAIRED != 0 && b.mpos == -1) {
                    overlaps.insert(b.name().to_vec(), slot);
                }
            }
            Some(&mate) => {
                let name = b.name().to_vec();
                overlaps.remove(&name);

                let a = &mut self.slots[mate as usize]
                    .as_mut()
                    .expect("overlap entries refer to buffered reads")
                    .read;

                if !tweak_overlap_quality(a, &mut node.read) {
                    bail!("tweak_overlap_quality failed");
                }
            }
        }

        Ok(())
    }
}

/// `resolve_cigar2`
fn resolve_cigar(b: &Read, pos: i64, s: &mut CigarState, node: usize) -> Entry {
    let cigar = &b.cigar;
    let n_cigar = cigar.len() as i32;
    let op_at = |k: i32| cigar[k as usize] & 0xf;
    let len_at = |k: i32| i64::from(cigar[k as usize] >> 4);

    let mut entry = Entry {
        node,
        ..Entry::default()
    };

    if s.k == -1 {
        entry.qpos = 0;

        if n_cigar == 1 {
            if is_match(op_at(0)) {
                s.k = 0;
                s.x = b.pos;
                s.y = 0;
            }
        } else {
            s.x = b.pos;
            s.y = 0;
            let mut k = 0;

            while k < n_cigar {
                let (op, l) = (op_at(k), len_at(k));

                if is_match(op) || op == CIGAR_DELETION || op == CIGAR_SKIP {
                    break;
                } else if op == CIGAR_INSERTION || op == CIGAR_SOFT_CLIP {
                    s.y += l;
                }

                k += 1;
            }

            s.k = k;
        }
    } else {
        let l = len_at(s.k);

        if pos - s.x >= l {
            let op = op_at(s.k + 1);

            if is_match(op) || op == CIGAR_DELETION || op == CIGAR_SKIP {
                if is_match(op_at(s.k)) {
                    s.y += l;
                }

                s.x += l;
                s.k += 1;
            } else {
                if is_match(op_at(s.k)) {
                    s.y += l;
                }

                s.x += l;
                let mut k = s.k + 1;

                while k < n_cigar {
                    let (op, l) = (op_at(k), len_at(k));

                    if is_match(op) || op == CIGAR_DELETION || op == CIGAR_SKIP {
                        break;
                    } else if op == CIGAR_INSERTION || op == CIGAR_SOFT_CLIP {
                        s.y += l;
                    }

                    k += 1;
                }

                s.k = k;
            }
        }
    }

    let (op, l) = (op_at(s.k), len_at(s.k));

    if s.x + l - 1 == pos && s.k + 1 < n_cigar {
        let op2 = op_at(s.k + 1);
        let l2 = len_at(s.k + 1);

        if op2 == CIGAR_DELETION && op != CIGAR_DELETION {
            // Merge consecutive deletions (1D2D is 3D).
            entry.indel = -(l2 as i32);

            for k in s.k + 2..n_cigar {
                if op_at(k) == CIGAR_DELETION {
                    entry.indel -= len_at(k) as i32;
                } else {
                    break;
                }
            }
        } else if op2 == CIGAR_INSERTION {
            entry.indel = l2 as i32;

            for k in s.k + 2..n_cigar {
                let op2 = op_at(k);

                if op2 == CIGAR_INSERTION {
                    entry.indel += len_at(k) as i32;
                } else if op2 != CIGAR_PAD {
                    break;
                }
            }
        } else if op2 == CIGAR_PAD && s.k + 2 < n_cigar {
            let mut l3 = 0;

            for k in s.k + 2..n_cigar {
                let op2 = op_at(k);

                if op2 == CIGAR_INSERTION {
                    l3 += len_at(k);
                } else if op2 == CIGAR_DELETION || is_match(op2) || op2 == CIGAR_SKIP {
                    break;
                }
            }

            if l3 > 0 {
                entry.indel = l3 as i32;
            }
        }
    }

    if is_match(op) {
        entry.qpos = s.y + (pos - s.x);
    } else if op == CIGAR_DELETION || op == CIGAR_SKIP {
        entry.is_del = true;
        entry.qpos = s.y;
        entry.is_refskip = op == CIGAR_SKIP;
    }

    entry.is_head = pos == b.pos;
    entry.is_tail = pos == s.end;
    entry.cigar_ind = s.k;
    entry
}

/// `bam_plp_insertion`: the inserted bases (and pads, as `*`) after this entry, and the length
/// of a deletion directly after them.
pub fn insertion(read: &Read, entry: &Entry, dst: &mut Vec<u8>) -> (usize, usize) {
    const SEQ_NT16_STR: &[u8; 16] = b"=ACMGRSVTWYHKDBN";

    dst.clear();

    if entry.indel <= 0 {
        return (0, 0);
    }

    let cigar = &read.cigar;
    let mut del_len = 0;
    let mut j = 1i64;
    let mut k = entry.cigar_ind as usize + 1;

    while k < cigar.len() {
        let (op, l) = (cigar[k] & 0xf, (cigar[k] >> 4) as usize);

        match op {
            CIGAR_PAD => dst.extend(std::iter::repeat_n(b'*', l)),
            CIGAR_INSERTION => {
                for _ in 0..l {
                    let i = entry.qpos + j - i64::from(entry.is_del);

                    dst.push(if (i as usize) < read.l_seq && i >= 0 {
                        SEQ_NT16_STR[usize::from(read.base(i as usize))]
                    } else {
                        b'N'
                    });

                    j += 1;
                }
            }
            CIGAR_DELETION => {
                del_len = l;
                break;
            }
            _ => break,
        }

        k += 1;
    }

    (dst.len(), del_len)
}

/// `cigar_iref2iseq_set`: moves to the first match base at reference offset `iref`.
fn iref2iseq_set(
    cigar: &[u32],
    ci: &mut usize,
    icig: &mut i64,
    iseq: &mut i64,
    iref: &mut i64,
) -> Option<bool> {
    let mut pos = *iref;

    if pos < 0 {
        return Some(false);
    }

    *icig = 0;
    *iseq = 0;
    *iref = 0;

    while *ci < cigar.len() {
        let (op, n) = (cigar[*ci] & 0xf, i64::from(cigar[*ci] >> 4));

        match op {
            CIGAR_SOFT_CLIP => {
                *ci += 1;
                *iseq += n;
                *icig = 0;
            }
            CIGAR_HARD_CLIP | CIGAR_PAD => {
                *ci += 1;
                *icig = 0;
            }
            _ if is_match(op) => {
                pos -= n;

                if pos < 0 {
                    *icig = n + pos;
                    *iseq += *icig;
                    *iref += *icig;
                    return Some(true);
                }

                *ci += 1;
                *iseq += n;
                *icig = 0;
                *iref += n;
            }
            CIGAR_INSERTION => {
                *ci += 1;
                *iseq += n;
                *icig = 0;
            }
            CIGAR_DELETION | CIGAR_SKIP => {
                pos -= n;

                if pos < 0 {
                    pos = 0;
                }

                *ci += 1;
                *icig = 0;
                *iref += n;
            }
            _ => return None,
        }
    }

    *iseq = -1;
    Some(false)
}

/// `cigar_iref2iseq_next`: moves to the next match base.
fn iref2iseq_next(
    cigar: &[u32],
    ci: &mut usize,
    icig: &mut i64,
    iseq: &mut i64,
    iref: &mut i64,
) -> Option<bool> {
    while *ci < cigar.len() {
        let (op, n) = (cigar[*ci] & 0xf, i64::from(cigar[*ci] >> 4));

        match op {
            _ if is_match(op) => {
                if *icig >= n - 1 {
                    *icig = -1;
                    *ci += 1;
                    continue;
                }

                *iseq += 1;
                *icig += 1;
                *iref += 1;
                return Some(true);
            }
            CIGAR_DELETION | CIGAR_SKIP => {
                *ci += 1;
                *iref += n;
                *icig = -1;
            }
            CIGAR_INSERTION | CIGAR_SOFT_CLIP => {
                *ci += 1;
                *iseq += n;
                *icig = -1;
            }
            CIGAR_HARD_CLIP | CIGAR_PAD => {
                *ci += 1;
                *icig = -1;
            }
            _ => return None,
        }
    }

    *iseq = -1;
    *iref = -1;
    Some(false)
}

/// A read's CIGAR, sequence, and mutable qualities, borrowed apart.
struct Mate<'a> {
    cigar: &'a [u32],
    l_seq: usize,
    seq: &'a [u8],
    qual: &'a mut [u8],
}

impl<'a> Mate<'a> {
    fn new(read: &'a mut Read) -> Self {
        let (head, qual) = read.raw.split_at_mut(read.qual_start);

        Self {
            cigar: &read.cigar,
            l_seq: read.l_seq,
            seq: &head[read.seq_start..],
            qual: &mut qual[..read.l_seq],
        }
    }

    fn base(&self, i: usize) -> u8 {
        (self.seq[i / 2] >> ((!i & 1) << 2)) & 0xf
    }

    fn qual_mut(&mut self, i: usize) -> &mut u8 {
        &mut self.qual[i]
    }
}

/// `tweak_overlap_quality`: for overlapping mates `a` (left) and `b`, zeroes or combines the
/// qualities of one of each pair of overlapping bases. Returns `false` on a malformed CIGAR.
fn tweak_overlap_quality(a: &mut Read, b: &mut Read) -> bool {
    /// `uint8 * 0.8` truncated back to a byte.
    fn scaled(q: u8) -> u8 {
        (f64::from(q) * 0.8) as u8
    }

    let (a_pos, b_pos) = (a.pos, b.pos);

    // Which read's qualities are kept, chosen by a hash of the name.
    let (amul, bmul): (u8, u8) = if wang_hash(x31_hash(a.name())) & 1 != 0 {
        (1, 0)
    } else {
        (0, 1)
    };

    let mut a = Mate::new(a);
    let mut b = Mate::new(b);
    let (a_cigar, b_cigar) = (a.cigar, b.cigar);
    let (mut a_ci, mut b_ci) = (0usize, 0usize);
    let (mut a_icig, mut a_iseq, mut b_icig, mut b_iseq) = (0i64, 0i64, 0i64, 0i64);

    let mut iref = b_pos;
    let mut a_iref = iref - a_pos;
    let mut b_iref = iref - b_pos;

    let Some(mut a_ret) = iref2iseq_set(a_cigar, &mut a_ci, &mut a_icig, &mut a_iseq, &mut a_iref)
    else {
        return false;
    };

    if !a_ret {
        return true;
    }

    let Some(mut b_ret) = iref2iseq_set(b_cigar, &mut b_ci, &mut b_icig, &mut b_iseq, &mut b_iref)
    else {
        return false;
    };

    if !b_ret {
        return true;
    }

    let a_l = a.l_seq as i64;
    let b_l = b.l_seq as i64;

    loop {
        while a_ret && a_iref >= 0 && a_iref < iref - a_pos {
            match iref2iseq_next(a_cigar, &mut a_ci, &mut a_icig, &mut a_iseq, &mut a_iref) {
                Some(r) => a_ret = r,
                None => return false,
            }
        }

        if !a_ret {
            return true;
        }

        while b_ret && b_iref >= 0 && b_iref < iref - b_pos {
            match iref2iseq_next(b_cigar, &mut b_ci, &mut b_icig, &mut b_iseq, &mut b_iref) {
                Some(r) => b_ret = r,
                None => return false,
            }
        }

        if !b_ret {
            return true;
        }

        if iref < a_iref + a_pos {
            iref = a_iref + a_pos;
        }

        if iref < b_iref + b_pos {
            iref = b_iref + b_pos;
        }

        iref += 1;

        // A deletion in one read: catch the other up, amending its qualities as for a mismatch.
        if a_iref + a_pos != b_iref + b_pos {
            if a_iref + a_pos < b_iref + b_pos
                && b_ci > 0
                && b_cigar[b_ci - 1] & 0xf == CIGAR_DELETION
            {
                loop {
                    let q = a.qual_mut(a_iseq as usize);
                    *q = if amul == 1 { scaled(*q) } else { 0 };

                    match iref2iseq_next(a_cigar, &mut a_ci, &mut a_icig, &mut a_iseq, &mut a_iref)
                    {
                        Some(true) => {}
                        Some(false) => return true,
                        None => return false,
                    }

                    if a_iref + a_pos >= b_iref + b_pos {
                        break;
                    }
                }
            } else if a_ci > 0 && a_cigar[a_ci - 1] & 0xf == CIGAR_DELETION {
                loop {
                    let q = b.qual_mut(b_iseq as usize);
                    *q = if bmul == 1 { scaled(*q) } else { 0 };

                    match iref2iseq_next(b_cigar, &mut b_ci, &mut b_icig, &mut b_iseq, &mut b_iref)
                    {
                        Some(true) => {}
                        Some(false) => return true,
                        None => return false,
                    }

                    if b_iref + b_pos >= a_iref + a_pos {
                        break;
                    }
                }
            } else {
                // Other cases, such as reference skips, are not supported.
                continue;
            }
        }

        if a_iseq > a_l || b_iseq > b_l {
            return false;
        }

        // htslib reads past the qualities here; only a malformed CIGAR gets this far.
        if a_iseq == a_l || b_iseq == b_l {
            return true;
        }

        let (ai, bi) = (a_iseq as usize, b_iseq as usize);
        let (aq, bq) = (*a.qual_mut(ai), *b.qual_mut(bi));

        if a.base(ai) == b.base(bi) {
            // Confident: the sum of the qualities, on one read.
            let qual = (u32::from(aq) + u32::from(bq)).min(200) as u8;
            *a.qual_mut(ai) = amul * qual;
            *b.qual_mut(bi) = bmul * qual;
        } else if aq > bq {
            *a.qual_mut(ai) = scaled(aq);
            *b.qual_mut(bi) = 0;
        } else if aq < bq {
            *b.qual_mut(bi) = scaled(bq);
            *a.qual_mut(ai) = 0;
        } else {
            *a.qual_mut(ai) = (f64::from(amul) * 0.8 * f64::from(aq)) as u8;
            *b.qual_mut(bi) = (f64::from(bmul) * 0.8 * f64::from(bq)) as u8;
        }
    }
}
