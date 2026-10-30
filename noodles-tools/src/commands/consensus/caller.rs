//! The consensus callers of `bam_consensus.c`: the Bayesian (Gap5-derived) caller and the simple
//! frequency-counting caller, with the per-read precomputation done as reads enter the pileup.
//!
//! The floating-point code follows the C operation order exactly, including the `fast_exp` and
//! `fast_log2` approximations, so calls and qualities match to the byte. C's double-to-int
//! conversions are done with [`c_int`], which (like x86's `cvttsd2si`) truncates and gives
//! `INT_MIN` for NaN or out-of-range values.

use std::sync::OnceLock;

use super::{
    pileup::{CHARD_CLIP, CSOFT_CLIP, Column, Entry},
    tables::{MQUAL_POW_1M, Q2P, QcalTable},
};

/// Consensus modes (`MODE_*`).
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub enum Mode {
    /// Frequency counting.
    Simple,
    /// samtools 1.16 compatible (no indel parameters).
    Bayes116,
    /// Gap5 parameters, favoring recall (the default).
    Recall,
    /// Parameters favoring precision.
    Precise,
    /// A blend of the precise and recall parameters.
    Mixed,
}

/// The options the callers use.
#[derive(Clone, Debug)]
pub struct CallerOptions {
    pub mode: Mode,
    pub use_qual: bool,
    pub min_qual: i32,
    pub adj_qual: bool,
    pub use_mqual: bool,
    pub scale_mqual: f64,
    pub nm_adjust: bool,
    pub nm_halo: i32,
    pub sc_cost: i32,
    pub low_mqual: i32,
    pub high_mqual: i32,
    pub min_depth: i32,
    pub call_fract: f64,
    pub het_fract: f64,
    pub cons_cutoff: i32,
    pub ambig: bool,
    pub default_qual: i32,
    pub homopoly_fix: f64,
}

/// `(int)x` as compiled for x86-64: truncation toward zero, and `INT_MIN` for NaN or values out
/// of range.
pub fn c_int(x: f64) -> i32 {
    if x.is_nan() || x >= 2_147_483_648.0 || x <= -2_147_483_649.0 {
        i32::MIN
    } else {
        x as i32
    }
}

/// `(int)(a + b * scale)` for ints `a` and `b`, without float arithmetic when `scale` is 1.
fn int_plus_scaled(a: i32, b: i32, scale: f64) -> i32 {
    if scale == 1.0 {
        // Exact in double precision too.
        a + b
    } else {
        c_int(f64::from(a) + f64::from(b) * scale)
    }
}

/// `(uint8_t)x` for a double, through a 32-bit truncation as x86-64 does it.
fn c_u8(x: f64) -> u8 {
    c_int(x) as u8
}

/// `e_tab`, `e_tab2`, and `e_log`, filled in once.
struct ExpTables {
    /// `exp(i)` for i in -500..=500, at `i + 500`.
    e_tab: Vec<f64>,
    /// `exp(i/10.)` for i in -500..=500, at `i + 500`.
    e_tab2: Vec<f64>,
}

fn exp_tables() -> &'static ExpTables {
    static TABLES: OnceLock<ExpTables> = OnceLock::new();

    TABLES.get_or_init(|| ExpTables {
        e_tab: (-500..=500).map(|i| f64::from(i).exp()).collect(),
        e_tab2: (-500..=500).map(|i| (f64::from(i) / 10.0).exp()).collect(),
    })
}

/// `fast_exp`: a table lookup. (C reads out of bounds for NaN; that reads 0 here.)
fn fast_exp(y: f64) -> f64 {
    let t = exp_tables();

    let (table, i) = if (-50.0..=50.0).contains(&y) {
        (&t.e_tab2, c_int(y * 10.0))
    } else {
        (&t.e_tab, c_int(y.clamp(-500.0, 500.0)))
    };

    usize::try_from(i64::from(i) + 500)
        .ok()
        .and_then(|i| table.get(i))
        .copied()
        .unwrap_or(0.0)
}

/// `fast_log2`: the exponent plus a cubic approximation of the log of the mantissa.
pub fn fast_log2(val: f64) -> f64 {
    let x = val.to_bits();
    let e = ((x >> 52) & 2047) as i32 - 1024;
    let m = f64::from_bits((x & !(2047u64 << 52)) + (1023u64 << 52));
    let val = ((-1.0 / 3.0) * m + 2.0) * m - 2.0 / 3.0;
    f64::from(e) + val
}

const TENLOG2OVERLOG10: f64 = 3.0103;

/// `ph_log`
fn ph_log(x: f64) -> f64 {
    -TENLOG2OVERLOG10 * fast_log2(x)
}

/// The fields of `cons_probs` in a flat array with the C struct's layout (16-byte aligned
/// arrays, zero padding), so that C's reads past the end of one table (qualities above 100)
/// read the next table as they do in C.
const PRIOR: usize = 0;
/// The highest quality the tables cover.
const MAX_QUAL: i32 = 100;
const LPRIOR15: usize = 26;
const PMM: usize = 42;
const PXX: usize = 144;
const PXM: usize = 246;
const POX: usize = 348;
const POM: usize = 450;
const POO: usize = 552;
const PUU: usize = 654;
const PUM: usize = 756;
const PMMD: usize = 858;
const POLY_MUL: usize = 959;
/// The flat table's length: the struct, then zeros for reads of up to 255 past any table.
const CONS_PROBS_PADDED: usize = PMMD + 256;

/// Precomputed matrices for the Bayesian caller (`cons_probs`).
#[derive(Clone, Debug)]
pub struct ConsProbs {
    /// Fixed-size, so that clamped table lookups need no bounds checks.
    v: Box<[f64; CONS_PROBS_PADDED]>,
    /// [`qual2`] for homopolymer lengths up to [`MAX_POLY`] and all qualities, at
    /// `poly * 256 + qual`.
    qual2: Box<[i32; (MAX_POLY + 1) * 256]>,
}

impl Default for ConsProbs {
    /// An uninitialized (zeroed) set, as the unused static is in C.
    fn default() -> Self {
        Self::from_values(vec![0.0; CONS_PROBS_PADDED])
    }
}

impl ConsProbs {
    /// `consensus_init`
    pub fn new(
        p_het: f64,
        p_indel: f64,
        het_scale: f64,
        poly_mul: f64,
        qcal: &QcalTable,
        mode: Mode,
    ) -> Self {
        let mut v = vec![0.0; CONS_PROBS_PADDED];
        let [smap, umap, omap] = qcal;

        v[POLY_MUL] = poly_mul;

        for i in 0..25 {
            v[PRIOR + i] = p_het / 6.0;
        }

        for i in [0, 6, 12, 18, 24] {
            v[PRIOR + i] = 1.0;
        }

        // Heterozygous deletions and insertions.
        for i in (4..24).step_by(5) {
            v[PRIOR + i] = p_indel / 6.0;
        }

        for i in 20..24 {
            v[PRIOR + i] = p_indel / 6.0;
        }

        for (j, i) in [0, 1, 2, 3, 4, 6, 7, 8, 9, 12, 13, 14, 18, 19, 24]
            .into_iter()
            .enumerate()
        {
            v[LPRIOR15 + j] = v[PRIOR + i].ln();
        }

        let phred = |q: i32| 1.0 - 10f64.powf(f64::from(-q) / 10.0);

        for i in 1..101 {
            let prob = phred(smap[i]);
            v[PMM + i] = prob.ln();
            v[PXX + i] = ((1.0 - prob) / 3.0).ln();
            v[PXM + i] = ((v[PMM + i].exp() + v[PXX + i].exp()) / 2.0).ln();
            v[PXM + i] += het_scale.ln();

            if mode == Mode::Bayes116 {
                v[PMMD + i] = v[PMM + i];
                v[POM + i] = v[PXM + i];
                v[PUM + i] = v[PXM + i];
                v[POX + i] = v[PXX + i];
                v[POO + i] = v[PXX + i];
                v[PUU + i] = v[PXX + i];
            } else {
                let prob = phred(omap[i]);
                v[POO + i] = ((1.0 - prob) / 3.0).ln();

                if v[POO + i] > v[PMM + i] - 0.5 {
                    v[POO + i] = v[PMM + i] - 0.5;
                }

                v[POX + i] = ((v[POO + i].exp() + v[PXX + i].exp()) / 2.0).ln();
                v[POM + i] = ((v[POO + i].exp() + v[PMM + i].exp()) / 2.0).ln();

                if v[POM + i] > v[PXM + i] + 0.5 {
                    v[POM + i] = v[PXM + i] + 0.5;
                }

                let prob = phred(umap[i]);
                v[PMMD + i] = prob.ln();
                v[PUU + i] = ((1.0 - prob) / 3.0).ln();

                if v[PUU + i] > v[PMM + i] - 0.5 {
                    v[PUU + i] = v[PMM + i] - 0.5;
                }

                v[PUM + i] = ((v[PUU + i].exp() + v[PMMD + i].exp()) / 2.0).ln();
            }
        }

        for table in [PMM, PXX, PXM, PMMD, POO, POX, POM, PUU, PUM] {
            v[table] = v[table + 1];
        }

        Self::from_values(v)
    }

    fn from_values(v: Vec<f64>) -> Self {
        let poly_mul = v[POLY_MUL];
        let qual2: Vec<i32> = (0..=MAX_POLY as i32)
            .flat_map(|poly| (0..=u8::MAX).map(move |q| qual2(q, poly, poly_mul)))
            .collect();

        Self {
            v: v.into_boxed_slice()
                .try_into()
                .expect("cons_probs has its fixed length"),
            qual2: qual2
                .into_boxed_slice()
                .try_into()
                .expect("one entry per poly and quality"),
        }
    }

    /// [`qual2`] with this set's homopolymer multiplier.
    fn qual2(&self, qual: u8, poly: i32) -> i32 {
        match usize::try_from(poly) {
            Ok(p) if p <= MAX_POLY => self.qual2[p * 256 + usize::from(qual)],
            _ => qual2(qual, poly, self.poly_mul()),
        }
    }

    /// `cp->table[q]`, with qualities clamped to the tables' range (samtools reads past a
    /// table, into the next one, for qualities above 100).
    fn at(&self, table: usize, q: i32) -> f64 {
        self.v[table + q.clamp(0, MAX_QUAL) as usize]
    }

    /// `cp->table[q]` for a byte quality, clamped as in [`Self::at`].
    fn at_u8(&self, table: usize, q: u8) -> f64 {
        self.v[table + usize::from(q.min(MAX_QUAL as u8))]
    }

    fn lprior15(&self, j: usize) -> f64 {
        self.v[LPRIOR15 + j]
    }

    fn poly_mul(&self) -> f64 {
        self.v[POLY_MUL]
    }
}

/// `q2p[q]`, with qualities clamped to 100 (samtools reads past the table above it).
fn q2p(q: u8) -> f64 {
    Q2P[usize::from(q.min(MAX_QUAL as u8))]
}

/// The quality adjusted by the mapping quality, `(uint8_t)ph_log(P+.75*M-P*M)`, for a base
/// quality `q` and the mapping quality table index `m`.
fn mq_adjusted(q: u8, m: usize) -> u8 {
    let pq = q2p(q);
    let m = MQUAL_POW_1M.get(m).copied().unwrap_or(0.0);
    c_u8(ph_log(pq + 0.75 * m - pq * m))
}

/// [`mq_adjusted`] for all base qualities and table indexes, at `m * 256 + q`.
fn mq_adjusted_table() -> &'static [u8; MQUAL_POW_1M.len() * 256] {
    static TABLE: OnceLock<Box<[u8; MQUAL_POW_1M.len() * 256]>> = OnceLock::new();

    TABLE.get_or_init(|| {
        let table: Vec<u8> = (0..MQUAL_POW_1M.len())
            .flat_map(|m| (0..=u8::MAX).map(move |q| mq_adjusted(q, m)))
            .collect();
        table
            .into_boxed_slice()
            .try_into()
            .expect("one entry per index and quality")
    })
}

/// The largest homopolymer length in [`ConsProbs`]' `qual2` table.
const MAX_POLY: usize = 100;

/// `MAX(1, qual-(poly-2)*cp->poly_mul)`, as an int.
fn qual2(qual: u8, poly: i32, poly_mul: f64) -> i32 {
    let q2 = f64::from(qual) - (f64::from(poly) - 2.0) * poly_mul;
    if 1.0 > q2 { 1 } else { c_int(q2) }
}

/// Scratch space for [`Caller::init_read`], reused across reads.
#[derive(Debug, Default)]
pub struct ReadScratch {
    /// The read's decoded bases.
    bases: Vec<u8>,
    /// Quality minima of the read's windows.
    window_min: Vec<u8>,
}

/// The result of the Bayesian caller (`consensus_t`).
#[derive(Clone, Copy, Debug, Default)]
pub struct Cons {
    /// The most likely call: A=0, C=1, G=2, T=3, *=4.
    pub call: i32,
    /// The most likely heterozygous call, as `"ACGT*"[het / 5]` and `"ACGT*"[het % 5]`.
    pub het_call: i32,
    pub het_logodd: i32,
    pub phred: i32,
    pub depth: i32,
}

/// Map the 15 possible combinations to 1-base or 2-base encodings.
const MAP_SING: [i32; 15] = [0, 5, 5, 5, 5, 1, 5, 5, 5, 2, 5, 5, 3, 5, 4];
const MAP_HET: [i32; 15] = [0, 1, 2, 3, 4, 6, 7, 8, 9, 12, 13, 14, 18, 19, 24];

/// From BAM nt16 codes (and 16 for `*`) to ACGT*N order.
const NT16_TO_ACGT: [usize; 32] = [
    5, 0, 1, 5, 2, 5, 5, 5, 3, 5, 5, 5, 5, 5, 5, 5, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
];

/// `nm_local` and `poly_len`: the local NM figure and the homopolymer length near `pos`, which
/// is a read offset plus 1.
fn local_nm(p: &Entry, pos: i32) -> (f64, i32) {
    const MASK: i32 = (1 << 24) - 1;

    let Some(nm) = p.local_nm() else {
        return (0.0, 0);
    };

    match usize::try_from(pos).ok().and_then(|i| nm.get(i)) {
        Some(&v) => (f64::from(v & MASK) / 10.0, v >> 24),
        None if pos < 0 => (f64::from(nm[0] & MASK), 0),
        None => (f64::from(nm[nm.len() - 1] & MASK), 0),
    }
}

/// `calculate_consensus_gap5`
fn consensus_gap5(column: &Column<'_>, opts: &CallerOptions, cp: &ConsProbs) -> Cons {
    let min_e_exp = f64::from(f64::MIN_EXP) * 2f64.ln() + 1.0;

    let mut s = [0.0f64; 15];
    let mut counts = [0i32; 6];

    let td = column.depth;
    let mut depth = 0;

    // The depth fudge of the mapping quality adjustment.
    let fudge = 1.0 + 2.0 * (0.5 - f64::from(td.min(30)) / 60.0);
    let low_mqual = f64::from(opts.low_mqual);
    let high_mqual = f64::from(opts.high_mqual);
    let mq_table = mq_adjusted_table();

    for p in column.reads() {
        if p.qual < opts.min_qual || p.ref_skip {
            continue;
        }

        let mut qual = p.qual as u8;

        if qual == 255 || (qual == 0 && p.qual_at(0) == 255) {
            qual = opts.default_qual as u8;
        }

        let base = NT16_TO_ACGT[(p.base4 & 31) as usize];
        let (nm, poly) = local_nm(p, p.seq_offset + 1);

        if opts.use_mqual {
            let mut mqual = f64::from(p.mapq);

            if opts.nm_adjust {
                mqual /= nm + 1.0;
                mqual *= fudge;
            }

            mqual *= opts.scale_mqual;

            if mqual < low_mqual {
                mqual = low_mqual;
            }

            if mqual > high_mqual {
                mqual = high_mqual;
            }

            // mqual is between the caps here unless it is NaN.
            let m = if mqual.is_nan() { -1 } else { mqual as i64 };

            qual = match usize::try_from(m) {
                Ok(m) if m < MQUAL_POW_1M.len() => mq_table[m * 256 + usize::from(qual)],
                _ => mq_adjusted(qual, usize::MAX),
            };
        }

        // Quality 0 would break the maths.
        if qual < 1 {
            qual = 1;
        }

        let qual2 = cp.qual2(qual, poly);

        let xx = cp.at_u8(PXX, qual);
        let mm_ = cp.at_u8(PMM, qual) - xx;
        let xm = cp.at_u8(PXM, qual) - xx;

        let (oo, om, ox, uu, um, mm) = match u8::try_from(qual2) {
            Ok(q) => (
                cp.at_u8(POO, q),
                cp.at_u8(POM, q),
                cp.at_u8(POX, q),
                cp.at_u8(PUU, q),
                cp.at_u8(PUM, q),
                cp.at_u8(PMMD, q),
            ),
            Err(_) => (
                cp.at(POO, qual2),
                cp.at(POM, qual2),
                cp.at(POX, qual2),
                cp.at(PUU, qual2),
                cp.at(PUM, qual2),
                cp.at(PMMD, qual2),
            ),
        };

        let (oo, om, ox) = (oo - xx, om - xx, ox - xx);
        let (uu, um, mm) = (uu - xx, um - xx, mm - xx);

        counts[base] += 1;

        match base {
            0 => {
                s[0] += mm_;
                s[1] += xm;
                s[2] += xm;
                s[3] += xm;
                s[4] += om;
                s[8] += ox;
                s[11] += ox;
                s[13] += ox;
                s[14] += oo;
            }
            1 => {
                s[1] += xm;
                s[5] += mm_;
                s[6] += xm;
                s[7] += xm;
                s[8] += om;
                s[4] += ox;
                s[11] += ox;
                s[13] += ox;
                s[14] += oo;
            }
            2 => {
                s[2] += xm;
                s[6] += xm;
                s[9] += mm_;
                s[10] += xm;
                s[11] += om;
                s[4] += ox;
                s[8] += ox;
                s[13] += ox;
                s[14] += oo;
            }
            3 => {
                s[3] += xm;
                s[7] += xm;
                s[10] += xm;
                s[12] += mm_;
                s[13] += om;
                s[4] += ox;
                s[8] += ox;
                s[11] += ox;
                s[14] += oo;
            }
            4 => {
                for j in [0, 1, 2, 3, 5, 6, 7, 9, 10, 12] {
                    s[j] += uu;
                }
                for j in [4, 8, 11, 13] {
                    s[j] += um;
                }
                s[14] += mm;
            }
            _ => {
                // N: equal weight to A, C, G, and T, but not a pad.
                for j in [0, 1, 2, 3, 5, 6, 7, 9, 10, 12] {
                    s[j] += mm_;
                }
                for j in [4, 8, 11, 13] {
                    s[j] += om;
                }
                s[14] += oo;
            }
        }

        depth += 1;
    }

    // Speculate on the consensus call.
    let mut shift = -f64::MAX;
    let mut max = -f64::MAX;
    let mut max_het = -f64::MAX;
    let mut call = 0;
    let mut het_call = 0;

    for (j, sj) in s.iter_mut().enumerate() {
        *sj += cp.lprior15(j);

        if shift < *sj {
            shift = *sj;
        }

        if !matches!(j, 0 | 5 | 9 | 12 | 14) {
            if max_het < *sj {
                max_het = *sj;
                het_call = j;
            }

            continue;
        }

        if max < *sj {
            max = *sj;
            call = j;
        }
    }

    // Shift and normalize.
    let mut norm = [0.0f64; 15];

    for sj in &mut s {
        *sj -= shift;
        let e = fast_exp(*sj);
        *sj = if *sj > min_e_exp {
            e
        } else {
            f64::MIN_POSITIVE
        };
    }

    let mut tot1 = 0.0;
    let mut tot2 = 0.0;

    for j in 0..15 {
        norm[j] += tot1;
        norm[14 - j] += tot2;
        tot1 += s[j];
        tot2 += s[14 - j];
    }

    if depth == 0 || depth == counts[5] {
        return Cons {
            call: 4,
            ..Cons::default()
        };
    }

    if norm[call] == 0.0 {
        norm[call] = f64::MIN_POSITIVE;
    }

    let ph = if s[call] == 1.0 && norm[call] < 0.01 {
        c_int(ph_log(norm[call]) + 0.5)
    } else {
        c_int(ph_log(1.0 - s[call] / (norm[call] + s[call])) + 0.5)
    };

    if norm[het_call] == 0.0 {
        norm[het_call] = f64::MIN_POSITIVE;
    }

    let het_ph =
        c_int(TENLOG2OVERLOG10 * (fast_log2(s[het_call]) - fast_log2(norm[het_call])) + 0.5);

    Cons {
        call: MAP_SING[call],
        het_call: MAP_HET[het_call],
        het_logodd: het_ph,
        phred: ph.max(0),
        depth,
    }
}

/// `calculate_consensus_gap5m`: the Bayesian caller, blending the precise and recall
/// parameters in mixed mode.
fn consensus_gap5m(
    column: &Column<'_>,
    opts: &CallerOptions,
    cp_r: &ConsProbs,
    cp_p: &ConsProbs,
) -> Cons {
    if opts.mode != Mode::Mixed {
        let cp = if opts.mode == Mode::Precise {
            cp_p
        } else {
            cp_r
        };
        return consensus_gap5(column, opts, cp);
    }

    let cons_p = consensus_gap5(column, opts, cp_p);
    let mut cons_r = consensus_gap5(column, opts, cp_r);

    let mut cons = cons_p;

    if cons_p.phred > 0 && cons_r.phred > 0 && cons_p.call == cons_r.call {
        // Both agree on a homozygous call.
        cons.phred += cons_r.phred.min(20);
    } else if cons_p.het_logodd >= 0 && cons_r.het_logodd >= 0 && cons_p.het_call == cons_r.het_call
    {
        // Both agree on a heterozygous call.
        cons.het_logodd += cons_r.het_logodd.min(20);
    } else if cons_p.het_logodd >= 0 {
        let q2 = cons_r.phred.max(cons_r.het_logodd);
        cons.het_logodd = 1.max(cons.het_logodd - q2 / 2);
    } else if cons_r.het_logodd >= 70 {
        let q1 = cons_p.phred;
        let q2 = cons_r.het_logodd;
        cons = cons_r;
        let a = f64::from((q2 - q1 * 2) / 2);
        let b = 1.0 + f64::from(q2) / (f64::from(q1) + 1.0);
        let hi = if a > b { a } else { b };
        cons.het_logodd = c_int(if 15.0 < hi { 15.0 } else { hi });
    } else if cons_r.het_logodd >= 0 {
        let q1 = cons_p.phred;
        let q2 = cons_r.het_logodd;
        cons = cons_r;
        let a = f64::from(q2) - 0.3 * f64::from(q1);
        let hi = if 1.0 > a { 1.0 } else { a };
        cons.het_logodd = c_int(hi + f64::from(5 * i32::from(cons_p.het_call == cons_r.het_call)));
        cons.phred = 0;
    } else {
        // Neither is heterozygous, but they differ in the call (very rare).
        cons_r.phred /= 2;

        if cons_r.phred > cons_p.phred {
            cons = cons_r;
        }

        cons.phred = cons.phred.max(10);
    }

    cons
}

/// `calculate_consensus_simple`: returns the call and its quality.
fn consensus_simple(column: &Column<'_>, opts: &CallerOptions) -> (u8, i32) {
    // nt16 to A, C, G, T compatibility, with weights on pure bases.
    const SEQI2A: [u64; 16] = [0, 8, 0, 4, 0, 4, 0, 2, 0, 4, 0, 2, 0, 2, 0, 1];
    const SEQI2C: [u64; 16] = [0, 0, 8, 4, 0, 0, 4, 2, 0, 0, 4, 2, 0, 0, 2, 1];
    const SEQI2G: [u64; 16] = [0, 0, 0, 0, 8, 4, 4, 1, 0, 0, 0, 0, 4, 2, 2, 1];
    const SEQI2T: [u64; 16] = [0, 0, 0, 0, 0, 0, 0, 0, 8, 4, 4, 2, 8, 2, 2, 1];
    const HET: &[u8; 32] = b"NACMGRSVTWYHKDBN*ac?g???t???????";

    let mut score = [0u64; 17];
    let mut tot_depth = 0;

    for p in column.reads() {
        let q = p.qual;

        if q < opts.min_qual {
            continue;
        }

        let w = if opts.use_qual { q as u64 } else { 1 };
        let b = p.base4;

        if let Ok(b) = usize::try_from(b)
            && b < 16
        {
            score[1] = score[1].wrapping_add(SEQI2A[b].wrapping_mul(w));
            score[2] = score[2].wrapping_add(SEQI2C[b].wrapping_mul(w));
            score[4] = score[4].wrapping_add(SEQI2G[b].wrapping_mul(w));
            score[8] = score[8].wrapping_add(SEQI2T[b].wrapping_mul(w));
        } else {
            score[16] = score[16].wrapping_add(8u64.wrapping_mul(w));
        }

        tot_depth += 1;
    }

    let tscore = (0..5).fold(0u64, |t, i| t.wrapping_add(score[1 << i]));

    let mut call1 = 15;
    let mut call2 = 15;
    let mut score1 = 0u64;
    let mut score2 = 0u64;

    for i in 0..5 {
        let c = 1 << i;

        if score1 < score[c] {
            score2 = score1;
            call2 = call1;
            score1 = score[c];
            call1 = c;
        } else if score2 < score[c] {
            score2 = score[c];
            call2 = c;
        }
    }

    let mut used_score = score1;
    let mut used_base = call1;

    if score2 as f64 >= opts.het_fract * score1 as f64 && opts.ambig {
        used_base |= call2;
        used_score = used_score.wrapping_add(score2);
    }

    if tot_depth < opts.min_depth || (used_score as f64) < opts.call_fract * tscore as f64 {
        used_base = if call1 == 16 { 16 } else { 0 };
    }

    let qual = if used_base != 0 {
        c_int(100.0 * used_score as f64 / tscore as f64)
    } else {
        0
    };

    (HET[used_base], qual)
}

/// The probability tables and options of the callers.
#[derive(Debug)]
pub struct Caller {
    pub opts: CallerOptions,
    recall: ConsProbs,
    precise: ConsProbs,
}

impl Caller {
    pub fn new(opts: CallerOptions, recall: ConsProbs, precise: ConsProbs) -> Self {
        Self {
            opts,
            recall,
            precise,
        }
    }

    /// `consensus_base`: the consensus base and quality of a column.
    pub fn call(&self, column: &Column<'_>) -> (u8, i32) {
        let opts = &self.opts;

        if opts.mode == Mode::Simple {
            return consensus_simple(column, opts);
        }

        let cons = consensus_gap5m(column, opts, &self.recall, &self.precise);

        let (mut cb, mut cq);

        if cons.depth < opts.min_depth && cons.call != 4 {
            cb = b'N';
            cq = 0;
        } else if cons.het_logodd > 0 && opts.ambig {
            const HET: &[u8; 25] = b"AMRWaMCSYcRSGKgWYKTtacgt*";
            cb = HET[cons.het_call as usize];
            cq = cons.het_logodd;
        } else {
            cb = b"ACGT*"[cons.call as usize];
            cq = cons.phred;
        }

        if cq < opts.cons_cutoff && cb != b'*' && cons.het_call % 5 != 4 && cons.het_call / 5 != 4 {
            // A het base/* keeps the base or *; otherwise N.
            cb = b'N';
            cq = 0;
        }

        (cb, cq)
    }

    /// `nm_init`: precomputes a read's local NM and homopolymer figures as it enters the
    /// pileup, and redistributes homopolymer qualities with `-p`. Returns whether to keep it.
    /// `bases` is scratch space for the read's decoded bases.
    pub fn init_read(&self, p: &mut Entry, scratch: &mut ReadScratch) -> bool {
        let opts = &self.opts;

        if !opts.use_mqual {
            return true;
        }

        let qlen = p.l_qseq;

        if qlen <= 0 {
            return false;
        }

        let n = qlen as usize;
        let ReadScratch { bases, window_min } = scratch;
        bases.clear();
        bases.extend((0..qlen).map(|i| p.seqi(i) as u8));
        let bases = bases.as_slice();

        let mut local_nm = p.take_nm_buffer(n);
        let poly_adj = if opts.homopoly_fix != 0.0 {
            opts.homopoly_fix
        } else {
            1.0
        };

        if opts.adj_qual {
            adjust_by_local_quality(p, bases, window_min, &mut local_nm, opts, poly_adj);
        }

        if opts.homopoly_fix != 0.0 {
            homopoly_qual_fix(p, bases);
        }

        // Homopolymer lengths in the top 8 bits.
        let mut i = 0;
        while i < qlen {
            let base = bases[i as usize];
            let mut j = i + 1;

            while j < qlen && bases[j as usize] == base {
                j += 1;
            }

            let poly = (j - i - 1).min(100);

            for nm in &mut local_nm[i as usize..j as usize] {
                *nm = (poly.max(*nm >> 24) << 24) | (*nm & ((1 << 24) - 1));
            }

            i = j;
        }

        adjust_by_md(p, &mut local_nm, opts);

        p.set_local_nm(local_nm);
        true
    }
}

/// The `adj_qual` part of `nm_init`: raises the local NM where a base's quality stands above
/// the local minimum.
fn adjust_by_local_quality(
    p: &Entry,
    bases: &[u8],
    window_min: &mut Vec<u8>,
    local_nm: &mut [i32],
    opts: &CallerOptions,
    poly_adj: f64,
) {
    const QHALO: i32 = 8;
    const QHALOP: i32 = 2;

    let qlen = p.l_qseq;
    let quals = p.quals();
    let qual = |i: i32| i32::from(quals[i as usize]);
    let bayes116 = opts.mode == Mode::Bayes116;

    // The minimum quality of each window of `2 * QHALO` bases, by its start.
    window_minima(quals, 2 * QHALO as usize, window_min);

    let mut qmin = qual(0);
    let mut qminp = qual(0);
    let seqi = |i: i32| bases[i as usize];
    let mut base = seqi(0);
    let mut polyl = 0;
    let mut polyr;

    // The minimum quality of the initial homopolymer.
    for i in 1..qlen {
        if seqi(i) != base {
            break;
        }

        if i < QHALOP && qminp > qual(i) {
            qminp = qual(i);
        }
    }

    // The minimum quality of the first bases.
    let mut i = 0;
    while i < qlen && i < QHALO {
        qmin = qmin.min(qual(i));
        i += 1;
    }

    while i < qlen - QHALO {
        if opts.homopoly_fix != 0.0 && seqi(i) != base {
            polyl = i;
            base = seqi(i);
            qminp = qual(i);

            let mut j = i + 1;
            while j < qlen {
                if seqi(j) != base {
                    break;
                }

                if i < QHALOP && qminp > qual(j) {
                    qminp = qual(j);
                }

                j += 1;
            }

            polyr = j - 1;
        } else {
            polyr = polyl;
        }

        let pl = polyr - polyl;

        let t = if bayes116 {
            (qual(i) + 5 * qmin) / 4
        } else {
            int_plus_scaled(qual(i) / 3, qminp - pl * 2, poly_adj)
        };

        local_nm[i as usize] += if t < qual(i) { qual(i) - t } else { 0 };

        qminp = qual(i);
        let (lo, hi) = (polyl.max(i - QHALOP), polyr.min(i + QHALOP));
        if lo <= hi {
            let window = &quals[lo as usize..=hi as usize];
            qminp = window.iter().fold(qminp, |m, &q| m.min(i32::from(q)));
        }

        if qmin > qual(i + QHALO) {
            qmin = qual(i + QHALO);
        } else if qmin <= qual(i - QHALO) {
            // The minimum of quals[i - QHALO + 1..=i + QHALO], capped at 99.
            qmin = 99.min(i32::from(window_min[(i - QHALO + 1) as usize]));
        }

        i += 1;
    }

    while i < qlen {
        let t = if bayes116 {
            (qual(i) + 5 * qmin) / 4
        } else {
            int_plus_scaled(qual(i) / 3, qminp, poly_adj)
        };

        local_nm[i as usize] += if t < qual(i) { qual(i) - t } else { 0 };
        i += 1;
    }
}

/// Writes the minimum of each `width`-long window of `values` (a power of two) to `dst`, indexed
/// by the window's start, by halving: windows of 2, then 4, and so on.
fn window_minima(values: &[u8], width: usize, dst: &mut Vec<u8>) {
    debug_assert!(width.is_power_of_two());
    dst.clear();

    if values.len() < width {
        return;
    }

    dst.extend_from_slice(values);
    let mut len = values.len();
    let mut span = 1;

    while span < width {
        len -= span;

        for i in 0..len {
            dst[i] = dst[i].min(dst[i + span]);
        }

        span *= 2;
    }

    dst.truncate(len);
}

/// The soft-clip and MD parts of `nm_init`: raises the local NM near soft-clips and
/// substitutions.
fn adjust_by_md(p: &Entry, local_nm: &mut [i32], opts: &CallerOptions) {
    let Ok(Some(md)) = crate::bam_raw::string_field(p.record(), *b"MD") else {
        return;
    };

    let qlen = p.l_qseq;
    let halo = opts.nm_halo;
    let mut add_range = |from: i32, to: i32, v: i32| {
        let from = usize::try_from(from).unwrap_or(0);
        let to = usize::try_from(to).unwrap_or(0).min(local_nm.len());

        if let Some(nms) = local_nm.get_mut(from..to) {
            for nm in nms {
                *nm = nm.wrapping_add(v);
            }
        }
    };
    let mut add = |i: i32, v: i32| add_range(i, i + 1, v);

    let ncig = p.n_cigar();
    let op = |i: u32| (p.cigar(i) & 0xf) as i32;

    if ncig > 0 {
        // The cost of being near a soft-clip.
        if op(0) == CSOFT_CLIP || (op(0) == CHARD_CLIP && ncig > 1 && op(1) == CSOFT_CLIP) {
            let mut i = 0;
            while i < halo && i < qlen {
                add(i, opts.sc_cost);
                i += 1;
            }
            while i < halo * 2 && i < qlen {
                add(i, opts.sc_cost >> 1);
                i += 1;
            }
        }

        let last = ncig - 1;
        if op(last) == CSOFT_CLIP
            || (op(last) == CHARD_CLIP && ncig > 1 && op(last - 1) == CSOFT_CLIP)
        {
            let mut i = qlen - 1;
            while i >= qlen - halo && i >= 0 {
                add(i, opts.sc_cost);
                i -= 1;
            }
            while i >= qlen - halo * 2 && i >= 0 {
                add(i, opts.sc_cost >> 1);
                i -= 1;
            }
        }
    }

    // Substitutions in the MD tag. (C never advances `pos` past them.)
    let mut pos: i32 = 0;
    let mut k = 0;

    while k < md.len() {
        let c = md[k];

        if c.is_ascii_digit() {
            let start = k;
            while k < md.len() && md[k].is_ascii_digit() {
                k += 1;
            }
            // strtol saturates at LONG_MAX; `pos += i` then truncates to int.
            let n = std::str::from_utf8(&md[start..k])
                .ok()
                .and_then(|s| s.parse::<i64>().ok())
                .unwrap_or(i64::MAX);
            pos = (i64::from(pos).wrapping_add(n)) as i32;
            continue;
        }

        if c == b'^' {
            k += 1;
            while k < md.len() && !md[k].is_ascii_digit() {
                k += 1;
            }
            continue;
        }

        // `i` runs on from one band to the next, as in the three C loops.
        let mut i = if pos - halo * 2 >= 0 {
            pos - halo * 2
        } else {
            0
        };

        for (bound, v) in [(pos - halo, 5), (pos + halo, 10), (pos + halo * 2, 5)] {
            let end = bound.min(qlen);

            if i < end {
                add_range(i, end, v);
                i = end;
            }
        }

        k += 1;
    }
}

/// `homopoly_qual_fix`: averages the error probabilities of the outermost pairs of bases in
/// each homopolymer, so the low qualities are spread to both ends.
fn homopoly_qual_fix(p: &mut Entry, bases: &[u8]) {
    static PH2ERR: OnceLock<Vec<f64>> = OnceLock::new();
    let ph2err =
        PH2ERR.get_or_init(|| (0..256).map(|i| 10f64.powf(f64::from(i) / -10.0)).collect());

    let qlen = p.l_qseq;
    let qual = p.quals_mut();

    let mut i = 0;
    while i < qlen {
        let s = i;
        let base = bases[i as usize];

        while i + 1 < qlen && bases[(i + 1) as usize] == base {
            i += 1;
        }

        let (mut j, mut k) = (s as usize, i as usize);
        while j < k {
            let e = ph2err[usize::from(qual[j])] + ph2err[usize::from(qual[k])];
            let q = c_u8(-fast_log2(e / 2.0) * 3.0104 + 0.49);
            qual[j] = q;
            qual[k] = q;
            j += 1;
            k -= 1;
        }

        i += 1;
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_window_minima() {
        let values: Vec<u8> = (0..40u32).map(|i| (i * 37 % 23) as u8).collect();
        let mut dst = Vec::new();
        window_minima(&values, 16, &mut dst);

        let expected: Vec<u8> = values
            .windows(16)
            .map(|w| *w.iter().min().unwrap())
            .collect();
        assert_eq!(dst, expected);

        window_minima(&values[..15], 16, &mut dst);
        assert!(dst.is_empty());
    }

    #[test]
    fn test_fast_log2() {
        // Exact at powers of two, close elsewhere.
        assert_eq!(fast_log2(1.0), 0.0);
        assert_eq!(fast_log2(8.0), 3.0);
        assert!((fast_log2(10.0) - 10f64.log2()).abs() < 0.01);
        assert!((fast_log2(0.3) - 0.3f64.log2()).abs() < 0.01);
    }

    #[test]
    fn test_fast_exp() {
        assert_eq!(fast_exp(0.0), 1.0);
        assert_eq!(fast_exp(-0.25), (-0.2f64).exp());
        assert_eq!(fast_exp(-60.5), (-60f64).exp());
        assert_eq!(fast_exp(-1000.0), (-500f64).exp());
    }

    #[test]
    fn test_c_int() {
        assert_eq!(c_int(2.9), 2);
        assert_eq!(c_int(-2.9), -2);
        assert_eq!(c_int(f64::NAN), i32::MIN);
        assert_eq!(c_int(1e10), i32::MIN);
        assert_eq!(c_u8(256.5), 0);
    }

    #[test]
    fn test_cons_probs_layout() {
        let qcal = super::super::tables::STATIC_QCAL[0];
        let cp = ConsProbs::new(1e-3, 2e-4, 1.0, 0.01, &qcal, Mode::Recall);
        // The flat table's last entry is 0, so log(1 - 10^0) = -inf at quality 100.
        assert_eq!(cp.at(PMM, 100), f64::NEG_INFINITY);
        // Qualities above the tables' range use quality 100.
        assert_eq!(cp.at(PMM, 102), cp.at(PMM, 100));
        assert_eq!(cp.at(PMMD, 101), cp.at(PMMD, 100));
        assert_eq!(cp.at_u8(PXX, 250), cp.at(PXX, 100));
    }
}
