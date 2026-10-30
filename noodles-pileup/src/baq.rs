//! Base alignment quality (BAQ) and mapping quality capping, ported from htslib 1.24
//! (`probaln.c` and `realn.c`) with the same floating-point operations in the same order, so
//! that adjusted qualities match bit for bit.

use std::sync::LazyLock;

use crate::Read;

const EI: f64 = 0.25;
const EM: f64 = 0.333_333_333_33;

/// `seq_nt16_table`: IUPAC characters to 4-bit codes.
pub const SEQ_NT16_TABLE: [u8; 256] = {
    let mut table = [15; 256];
    let codes = *b"=ACMGRSVTWYHKDBN";
    let mut i = 0;

    while i < 16 {
        table[codes[i] as usize] = i as u8;
        table[codes[i].to_ascii_lowercase() as usize] = i as u8;
        i += 1;
    }

    table[b'U' as usize] = 8;
    table[b'u' as usize] = 8;

    // htslib also reads 0-3 as A, C, G, and T.
    table[b'0' as usize] = 1;
    table[b'1' as usize] = 2;
    table[b'2' as usize] = 4;
    table[b'3' as usize] = 8;

    table
};

/// `seq_nt16_int`: 4-bit codes to 0-3 for A, C, G, T and 4 otherwise.
const SEQ_NT16_INT: [u8; 16] = [4, 0, 1, 4, 2, 4, 4, 4, 3, 4, 4, 4, 4, 4, 4, 4];

/// `probaln_par_t`: gap open and extension probabilities (single precision) and band width.
#[derive(Clone, Copy, Debug)]
struct Params {
    d: f32,
    e: f32,
    bw: i32,
}

/// `g_qual2prob`
static QUAL_TO_PROB: LazyLock<[f32; 256]> =
    LazyLock::new(|| std::array::from_fn(|i| 10f64.powf(-(i as f64) / 10.0) as f32));

/// `set_u`
fn band_index(bw: i32, i: i32, k: i32) -> usize {
    let x = (i - bw).max(0);
    ((k - x + 1) * 3) as usize
}

/// `probaln_glocal`, with the backward pass and MAP: fills `state` and `q` (each `query.len()`
/// long). Returns `false` on invalid input.
fn probaln_glocal(
    reference: &[u8],
    query: &[u8],
    iqual: &[u8],
    c: &Params,
    state: &mut [i32],
    q: &mut [u8],
) -> bool {
    let l_ref = reference.len() as i32;
    let l_query = query.len() as i32;

    if l_query >= i32::MAX - 2 {
        return false;
    }

    if l_ref == 0 || l_query == 0 {
        return true;
    }

    let mut bw = l_ref.max(l_query);

    if bw > c.bw {
        bw = c.bw;
    }

    if bw < (l_ref - l_query).abs() {
        bw = (l_ref - l_query).abs();
    }

    let bw2 = bw * 2 + 1;
    let i_dim = if bw2 < l_ref {
        bw2 as usize * 3 + 6
    } else {
        l_ref as usize * 3 + 6
    };

    let rows = l_query as usize + 1;
    let mut f = vec![0f64; rows * i_dim];
    let mut b = vec![0f64; rows * i_dim];
    let mut s = vec![0f64; l_query as usize + 2];

    let table = &*QUAL_TO_PROB;
    let qual: Vec<f32> = iqual.iter().map(|&q| table[usize::from(q)]).collect();

    let s_m = 1.0 / (2.0 * f64::from(l_query) + 2.0);
    let s_i = s_m;

    let d = c.d;
    let e = c.e;
    let mut m = [0f64; 9];
    m[0] = f64::from(1.0f32 - d - d) * (1.0 - s_m);
    m[1] = f64::from(d) * (1.0 - s_m);
    m[2] = m[1];
    m[3] = f64::from(1.0f32 - e) * (1.0 - s_i);
    m[4] = f64::from(e) * (1.0 - s_i);
    m[5] = 0.0;
    m[6] = f64::from(1.0f32 - e);
    m[7] = 0.0;
    m[8] = f64::from(e);
    let b_m = f64::from((1.0f32 - d) / l_ref as f32);
    let b_i = f64::from(d / l_ref as f32);

    let rf = |k: i32| reference[(k - 1) as usize];
    let qy = |i: i32| query[(i - 1) as usize];

    // forward: f[0]
    f[band_index(bw, 0, 0)] = 1.0;
    s[0] = 1.0;

    // f[1]
    {
        let end = if l_ref < bw + 1 { l_ref } else { bw + 1 };
        let mut sum = 0.0;

        for k in 1..=end {
            let e = if rf(k) > 3 || query[0] > 3 {
                1.0
            } else if rf(k) == query[0] {
                1.0 - f64::from(qual[0])
            } else {
                f64::from(qual[0]) * EM
            };

            let u = i_dim + band_index(bw, 1, k);
            f[u] = e * b_m;
            f[u + 1] = EI * b_i;
            sum += f[u] + f[u + 1];
        }

        s[1] = sum;
    }

    // f[2..l_query]
    for i in 2..=l_query {
        let qli = f64::from(qual[(i - 1) as usize]);
        let qyi = qy(i);
        let beg = 1.max(i - bw);
        let end = l_ref.min(i + bw);

        let e_table = [qli * EM, 1.0 - qli, 1.0, 1.0];
        let mm = 1.0 / s[(i - 1) as usize];

        let xm = [
            mm * m[0],
            mm * m[3],
            mm * m[6],
            EI * mm * m[1],
            EI * mm * m[4],
        ];

        let (prev, cur) = f.split_at_mut(i as usize * i_dim);
        let fi1 = &prev[(i as usize - 1) * i_dim..];
        let fi = &mut cur[..i_dim];

        let mut xi = band_index(bw, i, beg);
        let mut yi = band_index(bw, i - 1, beg - 1);

        let mut l_x0 = m[2] * fi[xi];
        let mut l_x2 = m[8] * fi[xi + 2];
        let mut sum = 0.0;

        for k in beg..=end {
            let r = rf(k);
            let cond = usize::from(r > 3 || qyi > 3) * 2 + usize::from(r == qyi);
            let z0 = xm[0] * fi1[yi];
            let z1 = xm[1] * fi1[yi + 1];
            let z2 = xm[2] * fi1[yi + 2];
            let z3 = xm[3] * fi1[yi + 3];
            let z4 = xm[4] * fi1[yi + 4];
            fi[xi] = e_table[cond] * (z0 + z1 + z2);
            fi[xi + 1] = z3 + z4;
            fi[xi + 2] = l_x0 + l_x2;
            sum += fi[xi] + fi[xi + 1] + fi[xi + 2];
            l_x0 = m[2] * fi[xi];
            l_x2 = m[8] * fi[xi + 2];
            xi += 3;
            yi += 3;
        }

        s[i as usize] = sum;
    }

    // f[l_query + 1]
    {
        let mm = 1.0 / s[l_query as usize];
        let mut sum = 0.0;
        let row = l_query as usize * i_dim;

        for k in 1..=l_ref {
            let u = band_index(bw, l_query, k);

            if u < 3 || u >= i_dim {
                continue;
            }

            sum += mm * f[row + u] * s_m + mm * f[row + u + 1] * s_i;
        }

        s[l_query as usize + 1] = sum;
    }

    // backward: b[l_query]
    {
        let row = l_query as usize * i_dim;
        let (sl, sl1) = (s[l_query as usize], s[l_query as usize + 1]);

        for k in 1..=l_ref {
            let u = band_index(bw, l_query, k);

            if u < 3 || u >= i_dim {
                continue;
            }

            b[row + u] = s_m / sl / sl1;
            b[row + u + 1] = s_i / sl / sl1;
        }
    }

    // b[l_query - 1..1]
    for i in (1..l_query).rev() {
        let beg = 1.max(i - bw);
        let end = l_ref.min(i + bw);
        let y = if i > 1 { 1.0 } else { 0.0 };
        let qli1 = f64::from(qual[i as usize]);
        let qyi1 = query[i as usize];
        let e_table = [qli1 * EM, 1.0 - qli1, 1.0, 1.0];

        let (cur, next) = b.split_at_mut((i as usize + 1) * i_dim);
        let bi = &mut cur[i as usize * i_dim..];
        let bi1 = &next[..i_dim];

        let mut xi = band_index(bw, i, end) as isize;
        let mut yi = band_index(bw, i + 1, end) as isize;

        let mut xi_5 = bi[(xi + 5) as usize];
        let e1 = EI * m[1];
        let e4 = EI * m[4];
        let n = 1.0 / s[i as usize];

        let mut k = end;

        while k >= beg {
            let (x, yy) = (xi as usize, yi as usize);

            let e = if k >= l_ref {
                0.0
            } else {
                let r = reference[k as usize];
                e_table[usize::from(r > 3 || qyi1 > 3) * 2 + usize::from(r == qyi1)] * bi1[yy + 3]
            };

            bi[x + 1] = e * m[3] + e4 * bi1[yy + 1];
            bi[x] = e * m[0] + e1 * bi1[yy + 1] + m[2] * xi_5;
            bi[x + 2] = (e * m[6] + m[8] * xi_5) * y;
            xi_5 = bi[x + 2];
            bi[x + 1] *= n;
            bi[x] *= n;
            bi[x + 2] *= n;

            k -= 1;
            xi -= 3;
            yi -= 3;
        }
    }

    // MAP
    for i in 1..=l_query {
        let row = i as usize * i_dim;
        let beg = 1.max(i - bw);
        let end = l_ref.min(i + bw);
        let mm = 1.0 / s[i as usize];
        let (mut sum, mut max) = (0.0, 0.0);
        let mut max_k = -1;
        let mut u = band_index(bw, i, beg);

        for k in beg..=end {
            let z1 = mm * f[row + u] * b[row + u];
            let z2 = mm * f[row + u + 1] * b[row + u + 1];
            let which = i32::from(z2 > z1);
            let zm = if which == 1 { z2 } else { z1 };

            if zm > max {
                max = zm;
                max_k = ((k - 1) << 2) | which;
            }

            sum += z1 + z2;
            u += 3;
        }

        max /= sum;
        state[(i - 1) as usize] = max_k;

        let k = c_int(-4.343 * (1.0 - max).ln() + 0.499);
        q[(i - 1) as usize] = if k > 100 { 99 } else { k as u8 };
    }

    true
}

/// Why BAQ was not computed.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
#[non_exhaustive]
pub enum BaqSkipped {
    /// Unmapped, empty, or with no matches or with reference skips (htslib returns -1).
    NotApplicable,
    /// Already adjusted (-3), or tags were converted (0).
    Tagged,
    /// A ZQ tag is malformed (-4).
    Invalid,
}

impl std::error::Error for BaqSkipped {}

impl std::fmt::Display for BaqSkipped {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.write_str(match self {
            Self::NotApplicable => "BAQ does not apply to this read",
            Self::Tagged => "the read already has BAQ tags",
            Self::Invalid => "invalid BAQ input",
        })
    }
}

/// `sam_prob_realn`'s `BAQ_APPLY`, `BAQ_EXTEND`, and `BAQ_REDO` flags.
#[derive(Clone, Copy, Debug, Default)]
pub struct BaqMode {
    /// Cap the qualities by their BAQ and record the change in a ZQ tag, instead of only storing
    /// the BAQ in a BQ tag.
    pub apply: bool,
    /// Use extended BAQ.
    pub extend: bool,
    /// Ignore existing BQ tags.
    pub redo: bool,
}

/// `sam_prob_realn`: computes the BAQ of a read's bases against `reference`, recording it in a
/// BQ tag, or with [`BaqMode::apply`] capping the qualities by it and recording the change in a
/// ZQ tag.
///
/// # Errors
///
/// Returns why the read was left unchanged: it is not aligned or has no bases or matches
/// ([`BaqSkipped::NotApplicable`]), it already carries the requested tag
/// ([`BaqSkipped::Tagged`]), or its ZQ tag is malformed ([`BaqSkipped::Invalid`]).
pub fn prob_realn(read: &mut Read, reference: &[u8], mode: BaqMode) -> Result<(), BaqSkipped> {
    let BaqMode {
        apply,
        extend,
        redo,
    } = mode;
    let l_qseq = read.sequence_len();

    let mut params = Params {
        d: 0.001,
        e: 0.1,
        bw: 10,
    };

    if l_qseq > 1000 {
        params.d = 1e-7;
        params.e = 1e-1;
    }

    if read.flags() & 0x4 != 0 || l_qseq == 0 || read.quality_scores()[0] == 0xff {
        return Err(BaqSkipped::NotApplicable);
    }

    // `realn_check_tag`: a Z field as long as the sequence.
    let is_valid = |read: &Read, tag: [u8; 2]| {
        read.tag_range(tag)
            .is_some_and(|r| read.raw[r.start + 2] == b'Z' && r.end - r.start - 4 == l_qseq)
    };

    let mut bq = read.tag_range(*b"BQ").is_some();
    let fix_bq = bq && !redo && !is_valid(read, *b"BQ");
    let mut zq = read.tag_range(*b"ZQ").is_some();

    if zq && !is_valid(read, *b"ZQ") {
        return Err(BaqSkipped::Invalid);
    }

    if bq && redo {
        read.remove_tag(*b"BQ");
        bq = false;
    }

    if bq && zq {
        read.remove_tag(*b"ZQ");
        zq = false;
    }

    if !zq && fix_bq {
        read.remove_tag(*b"BQ");
        bq = false;
    }

    if bq || zq {
        // A tag of the kind requested means there is nothing to do; one of the other kind is
        // converted.
        if (apply && zq) || (!apply && bq) {
            return Err(BaqSkipped::Tagged);
        }

        if zq {
            // ZQ to BQ
            let start = read.tag_range(*b"ZQ").expect("ZQ is present").start;
            let values: Vec<u8> = read.raw[start + 3..start + 3 + l_qseq].to_vec();

            for (q, &z) in read.quality_scores_mut().iter_mut().zip(&values) {
                *q = (i32::from(*q) + i32::from(z) - 64) as u8;
            }

            read.raw[start] = b'B';
            return Err(BaqSkipped::Tagged);
        }

        let start = read.tag_range(*b"BQ").expect("BQ is present").start;
        let values: Vec<u8> = read.raw[start + 3..start + 3 + l_qseq].to_vec();

        for (q, &b) in read.quality_scores_mut().iter_mut().zip(&values) {
            *q = if u32::from(*q) + 64 < u32::from(b) {
                0
            } else {
                (i32::from(*q) - (i32::from(b) - 64)) as u8
            };
        }

        read.raw[start] = b'Z';
        return Err(BaqSkipped::Tagged);
    }

    // The start and end of the alignment
    let (mut x, mut y) = (read.position(), 0i64);
    let (mut yb, mut ye, mut xb, mut xe) = (-1i64, -1i64, -1i64, -1i64);

    for &op in read.cigar() {
        let (kind, l) = (op & 0xf, i64::from(op >> 4));

        match kind {
            0 | 7 | 8 => {
                if yb < 0 {
                    yb = y;
                }

                if xb < 0 {
                    xb = x;
                }

                ye = y + l;
                xe = x + l;
                x += l;
                y += l;
            }
            4 | 1 => y += l,
            2 => x += l,
            3 => return Err(BaqSkipped::NotApplicable),
            _ => {}
        }
    }

    if xb == -1 {
        return Err(BaqSkipped::NotApplicable);
    }

    let l_qseq_i = l_qseq as i64;

    #[expect(
        clippy::cast_possible_truncation,
        reason = "htslib passes the difference to abs(int)"
    )]
    let diff = (((xe - xb) - (ye - yb)) as i32).abs();
    let mut bw = 7;

    if diff > bw {
        bw = diff + 3;
    }

    params.bw = bw;

    let bw = i64::from(bw);
    xb -= yb + bw / 2;

    if xb < 0 {
        xb = 0;
    }

    xe += l_qseq_i - ye + bw / 2;

    if xe - xb - l_qseq_i > bw {
        xb += (xe - xb - l_qseq_i - bw) / 2;
        xe -= (xe - xb - l_qseq_i - bw) / 2;
    }

    let tseq: Vec<u8> = (0..l_qseq)
        .map(|i| SEQ_NT16_INT[usize::from(read.base(i))])
        .collect();

    let mut tref = Vec::with_capacity(usize::try_from(xe - xb).unwrap_or(0));
    let ref_len = reference.len() as i64;

    for i in xb..xe {
        if i >= ref_len || reference[i as usize] == 0 {
            xe = i;
            break;
        }

        tref.push(SEQ_NT16_INT[usize::from(SEQ_NT16_TABLE[usize::from(reference[i as usize])])]);
    }

    let qual: Vec<u8> = read.quality_scores().to_vec();
    let mut state = vec![0i32; l_qseq];
    let mut q = vec![0u8; l_qseq];

    if xe < xb {
        return Err(BaqSkipped::Invalid);
    }

    let tref = &tref[..(xe - xb) as usize];

    if !probaln_glocal(tref, &tseq, &qual, &params, &mut state, &mut q) {
        return Err(BaqSkipped::Invalid);
    }

    // htslib starts from the qualities, so bases outside matches get no adjustment.
    let mut bq = qual.clone();

    if extend {
        extended_baq(read.cigar(), read.position(), xb, &state, &q, &mut bq);

        for (b, &q) in bq.iter_mut().zip(&qual) {
            *b = 64 + q.saturating_sub(*b);
        }
    } else {
        capped_baq(read.cigar(), read.position(), xb, &state, &q, &mut bq);

        for (b, &q) in bq.iter_mut().zip(&qual) {
            *b = q.wrapping_sub(*b).wrapping_add(64);
        }
    }

    if apply {
        for (q, &b) in read.quality_scores_mut().iter_mut().zip(&bq) {
            *q = q.wrapping_sub(b.wrapping_sub(64));
        }

        read.append_string_tag(*b"ZQ", &bq);
    } else {
        read.append_string_tag(*b"BQ", &bq);
    }

    Ok(())
}

/// Whether the realignment put query base `i` on the reference base the CIGAR puts it on.
fn on_cigar(state: &[i32], i: usize, expected: i64) -> bool {
    (state[i] & 3) == 0 && i64::from(state[i] >> 2) == expected
}

/// BAQ capped by the base quality: `bq` holds the qualities and becomes min(quality, BAQ) in
/// matches, or 0 where the realignment disagrees with the CIGAR.
fn capped_baq(cigar: &[u32], pos: i64, xb: i64, state: &[i32], q: &[u8], bq: &mut [u8]) {
    let l_qseq = bq.len() as i64;
    let (mut x, mut y) = (pos, 0i64);

    for &op in cigar {
        let kind = op & 0xf;
        let mut l = i64::from(op >> 4);

        if l == 0 {
            continue;
        }

        match kind {
            0 | 7 | 8 => {
                if l > l_qseq - y {
                    l = l_qseq - y;
                }

                for i in y..y + l {
                    let i_us = i as usize;

                    bq[i_us] = if on_cigar(state, i_us, x - xb + (i - y)) {
                        bq[i_us].min(q[i_us])
                    } else {
                        0
                    };
                }

                x += l;
                y += l;
            }
            4 | 1 => {
                if l > l_qseq - y {
                    l = l_qseq - y;
                }

                y += l;
            }
            2 => x += l,
            _ => {}
        }
    }
}

/// Extended BAQ, which can exceed the base quality: in each run of matches, the BAQ of a base is
/// the smaller of the running maxima from either side.
fn extended_baq(cigar: &[u32], pos: i64, xb: i64, state: &[i32], q: &[u8], bq: &mut [u8]) {
    let l_qseq = bq.len() as i64;
    let mut left = vec![0u8; bq.len()];
    let mut rght = vec![0u8; bq.len()];
    let (mut x, mut y) = (pos, 0i64);
    let mut len = 0i64;

    for k in 0..cigar.len() {
        let op = cigar[k] & 0xf;
        let mut l = i64::from(cigar[k] >> 4);

        // Adjacent match ops count as one, so that 50M50M is the same as 100M.
        if matches!(op, 0 | 7 | 8) {
            if k + 1 < cigar.len() && matches!(cigar[k + 1] & 0xf, 0 | 7 | 8) {
                len += l;
                continue;
            }

            l += len;
            len = 0;
        }

        if l == 0 {
            continue;
        }

        match op {
            0 | 7 | 8 => {
                if l > l_qseq - y {
                    l = l_qseq - y;
                }

                // Only a CIGAR longer than the sequence leaves nothing here.
                if l <= 0 {
                    x += l;
                    y += l;
                    continue;
                }

                let (ys, ls) = (y as usize, l as usize);

                for i in ys..ys + ls {
                    bq[i] = if on_cigar(state, i, x - xb + (i as i64 - y)) {
                        q[i]
                    } else {
                        0
                    };
                }

                left[ys] = bq[ys];

                for i in ys + 1..ys + ls {
                    left[i] = bq[i].max(left[i - 1]);
                }

                rght[ys + ls - 1] = bq[ys + ls - 1];

                for i in (ys..ys + ls - 1).rev() {
                    rght[i] = bq[i].max(rght[i + 1]);
                }

                for i in ys..ys + ls {
                    bq[i] = left[i].min(rght[i]);
                }

                x += l;
                y += l;
            }
            4 | 1 => {
                if l > l_qseq - y {
                    l = l_qseq - y;
                }

                y += l;
            }
            2 => x += l,
            _ => {}
        }
    }
}

/// `sam_cap_mapq`: the capped mapping quality, or `None` to skip the read.
#[must_use]
pub fn cap_mapq(read: &Read, reference: &[u8], thres: i32) -> Option<i32> {
    let thres = if thres < 0 { 40 } else { thres };
    let ref_len = reference.len() as i64;
    let qual = read.quality_scores();

    let (mut mm, mut q, mut len, mut clip_q) = (0i32, 0i32, 0i32, 0i32);
    let (mut x, mut y) = (read.position(), 0usize);

    'cigar: for &op in read.cigar() {
        let (kind, l) = (op & 0xf, (op >> 4) as usize);

        match kind {
            0 | 7 | 8 => {
                for j in 0..l {
                    let pos = x + j as i64;

                    if pos >= ref_len || reference[pos as usize] == 0 {
                        break 'cigar;
                    }

                    let z = y + j;
                    let c1 = read.base(z);
                    let c2 = SEQ_NT16_TABLE[usize::from(reference[pos as usize])];

                    if c2 != 15 && c1 != 15 && qual[z] >= 13 {
                        len += 1;

                        if c1 != 0 && c1 != c2 && qual[z] >= 13 {
                            mm += 1;
                            q += i32::from(qual[z].min(33));
                        }
                    }
                }

                x += l as i64;
                y += l;
                len += l as i32;
            }
            2 => {
                for j in 0..l {
                    let pos = x + j as i64;

                    if pos >= ref_len || reference[pos as usize] == 0 {
                        break 'cigar;
                    }
                }

                x += l as i64;
            }
            4 => {
                clip_q += qual[y..y + l].iter().map(|&q| i32::from(q)).sum::<i32>();
                y += l;
            }
            5 => clip_q += 13 * l as i32,
            1 => y += l,
            3 => x += l as i64,
            _ => {}
        }
    }

    let mut t = 1.0f64;

    for i in 0..mm {
        t *= f64::from(len) / f64::from(i + 1);
    }

    let mut t = f64::from(q) - 4.343 * t.ln() + f64::from(clip_q) / 5.0;

    if t > f64::from(thres) {
        return None;
    }

    if t < 0.0 {
        t = 0.0;
    }

    let thres = f64::from(thres);
    let t = ((thres - t) / thres).sqrt() * thres;

    Some(c_int(t + 0.499))
}

/// A C `(int)` conversion as x86-64 performs it: values that are NaN or out of range become
/// `INT_MIN` (where a Rust cast would saturate). A certain BAQ base, whose phred value is
/// infinite, gets a quality of 0 this way, as in htslib.
fn c_int(value: f64) -> i32 {
    if value.is_nan() || !(-2_147_483_648.0..2_147_483_648.0).contains(&value) {
        i32::MIN
    } else {
        #[expect(clippy::cast_possible_truncation, reason = "the value is in range")]
        let n = value as i32;
        n
    }
}
