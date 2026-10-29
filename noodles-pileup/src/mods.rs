//! Base modifications from MM and ML tags (htslib 1.24's `sam_mods.c`).
//!
//! [`BaseMods`] walks a read left to right, reporting the modifications called at each query
//! position, as `bam_parse_basemod`, `bam_mods_at_next_pos`, and `bam_mods_at_qpos` do.

use std::{error, fmt, io::Write};

use crate::Read;

/// The 4-bit code of the complement of each 4-bit base code (`seqi_rc`).
const SEQI_RC: [u8; 16] = [0, 8, 4, 12, 2, 10, 6, 14, 1, 9, 5, 13, 3, 11, 7, 15];
const SEQ_NT16_STR: &[u8; 16] = b"=ACMGRSVTWYHKDBN";
/// htslib's `MAX_BASE_MOD`.
const MAX_MODS: usize = 256;
const REVERSE: u16 = 0x10;

/// A modification called at a base (`hts_base_mod`).
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub struct BaseModification {
    /// The modification code (a character, such as `m`), or a negated ChEBI number.
    pub modified_base: i32,
    /// The canonical base (`A`, `C`, `G`, `T`, or `N`).
    pub canonical_base: u8,
    /// 0 for the `+` strand, 1 for `-`.
    pub strand: u8,
    /// The likelihood (0-255) from ML, or -1 if there is none.
    pub qual: i32,
}

/// An error parsing MM and ML tags.
#[derive(Clone, Debug, Eq, PartialEq)]
pub struct ParseError(&'static str);

impl error::Error for ParseError {}

impl fmt::Display for ParseError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.write_str(self.0)
    }
}

/// The state of one modification type (one code of an MM entry).
#[derive(Clone, Copy, Debug)]
struct ModState {
    ty: i32,
    canonical: u8,
    strand: u8,
    /// Canonical bases left until the next call (`MMcount`).
    count: i64,
    /// The position in the MM string of the next delta (`MM`).
    mm: usize,
    /// On the reverse strand, the end of the deltas not yet used (`MMend`).
    mm_end: usize,
    /// The index of the next likelihood in ML, if there is an ML tag.
    ml: Option<isize>,
    stride: usize,
}

/// The modification calls of a read, walked left to right (`hts_base_mod_state`).
#[derive(Clone, Debug, Default)]
pub struct BaseMods {
    mods: Vec<ModState>,
    seq_pos: usize,
    /// The MM value, NUL terminated.
    mm: Vec<u8>,
    ml: Vec<u8>,
}

impl BaseMods {
    /// `bam_parse_basemod`: parses a read's MM (or Mm) and ML (or Ml) tags. A read without MM
    /// has no modifications.
    ///
    /// # Errors
    ///
    /// Returns an error if the tags are malformed or do not fit the sequence.
    pub fn parse(read: &Read) -> Result<Self, ParseError> {
        let mut state = Self::default();

        let Some((ty, mm)) = read.tag(*b"MM").or_else(|| read.tag(*b"Mm")) else {
            return Ok(state);
        };

        if ty != b'Z' {
            return Err(ParseError("MM tag is not of type Z"));
        }

        let l_qseq = read.sequence_len();

        if let Some(n) = read.integer_tag(*b"MN")
            && n != l_qseq as i64
            && l_qseq > 0
        {
            return Err(ParseError(
                "MM/MN data length is incompatible with SEQ length",
            ));
        }

        let ml = match read.tag(*b"ML").or_else(|| read.tag(*b"Ml")) {
            Some((b'B', value)) if value.first() == Some(&b'C') && value.len() >= 5 => {
                let n = u32::from_le_bytes([value[1], value[2], value[3], value[4]]) as usize;
                Some(value[5..5 + n.min(value.len() - 5)].to_vec())
            }
            Some(_) => return Err(ParseError("ML tag is not of type B,C")),
            None => None,
        };

        state.mm = mm.to_vec();

        if state.mm.last() != Some(&0) {
            state.mm.push(0);
        }

        let reverse = read.flags() & REVERSE != 0;
        let mut freq = [0i64; 16];

        if reverse {
            for i in 0..l_qseq {
                freq[usize::from(read.base(i))] += 1;
            }

            freq[15] = l_qseq as i64;
        }

        let mm = &state.mm;
        let ml_len = ml.as_ref().map_or(0, Vec::len) as isize;
        let mut ml_pos: isize = 0;
        let mut cp = 0;

        while mm[cp] != 0 {
            while mm[cp] != 0 {
                // [ACGTNU][+-]([a-zA-Z]+|[0-9]+)[.?]?(,\d+)*;
                let mut btype = mm[cp];
                cp += 1;

                if !matches!(btype, b'A' | b'C' | b'G' | b'T' | b'U' | b'N') {
                    return Err(ParseError("malformed MM tag (base)"));
                }

                if btype == b'U' {
                    btype = b'T';
                }

                let btype = crate::baq::SEQ_NT16_TABLE[usize::from(btype)];

                let strand = match mm[cp] {
                    b'+' => 0,
                    b'-' => 1,
                    _ => return Err(ParseError("malformed MM tag (strand)")),
                };

                cp += 1;

                let mut ms = cp;
                let mut cp_end = None;
                let mut chebi = 0;

                if mm[cp].is_ascii_digit() {
                    let (n, end) = str2uint(mm, cp)
                        .filter(|&(_, end)| end != cp)
                        .ok_or(ParseError("malformed MM tag (invalid ChEBI code)"))?;
                    chebi = n as i32;
                    cp = end;
                    cp_end = Some(end);
                    ms = cp - 1;
                } else {
                    while mm[cp] != 0 && mm[cp].is_ascii_alphabetic() {
                        cp += 1;
                    }

                    if mm[cp] == 0 {
                        return Err(ParseError("malformed MM tag"));
                    }
                }

                let me = cp;

                match mm[cp] {
                    b'.' | b'?' => cp += 1,
                    b',' | b';' => {}
                    _ => return Err(ParseError("malformed MM tag")),
                }

                let stride = me - ms;
                let mut ndelta = 0isize;
                let delta: i64;

                if reverse {
                    // The deltas count from the right; the left-most one is the remainder.
                    let mut total = 0i64;

                    loop {
                        if mm[cp] == b',' {
                            cp += 1;
                        }

                        if mm[cp] == 0 || mm[cp] == b';' {
                            break;
                        }

                        let (n, end) = str2uint(mm, cp)
                            .filter(|&(_, end)| end != cp)
                            .ok_or(ParseError("Hit end of MM tag. Missing semicolon?"))?;
                        cp = end;
                        cp_end = Some(end);
                        total += n + 1;
                        ndelta += 1;
                    }

                    delta = freq[usize::from(SEQI_RC[usize::from(btype)])] - total;
                } else if mm[cp] == b',' {
                    let (n, end) = str2uint(mm, cp + 1)
                        .filter(|&(_, end)| end != cp + 1)
                        .ok_or(ParseError("Failed to parse integer from MM tag"))?;
                    delta = n;
                    cp_end = Some(end);
                } else {
                    delta = i64::from(i32::MAX);
                    cp_end = Some(cp);
                }

                for (n, m) in (ms..me).enumerate() {
                    if delta < 0 {
                        return Err(ParseError("MM tag refers to bases beyond sequence length"));
                    }

                    let ml_index = ml.as_ref().map(|_| {
                        if reverse {
                            ml_pos + n as isize + (ndelta - 1) * stride as isize
                        } else {
                            ml_pos + n as isize
                        }
                    });

                    // In htslib, a reverse entry without deltas has no end (`NULL`).
                    let end = cp_end.unwrap_or(0);

                    state.mods.push(ModState {
                        ty: if chebi != 0 { -chebi } else { i32::from(mm[m]) },
                        canonical: btype,
                        strand,
                        count: delta,
                        mm: if reverse { me + 1 } else { end },
                        mm_end: if reverse { end } else { 0 },
                        ml: ml_index,
                        stride,
                    });

                    if state.mods.len() >= MAX_MODS {
                        return Err(ParseError("Too many base modification types"));
                    }
                }

                // Skip the deltas.
                if ml.is_some() {
                    if reverse {
                        ml_pos += ndelta * stride as isize;
                    } else {
                        while mm[cp] != 0 && mm[cp] != b';' {
                            if mm[cp] == b',' {
                                ml_pos += stride as isize;
                            }

                            cp += 1;
                        }
                    }

                    if ml_pos > ml_len {
                        return Err(ParseError("Insufficient number of entries in ML tag"));
                    }
                } else if let Some(end) = cp_end
                    && reverse
                {
                    cp = end;
                } else {
                    while mm[cp] != 0 && mm[cp] != b';' {
                        cp += 1;
                    }
                }

                if mm[cp] == 0 {
                    return Err(ParseError("Hit end of MM tag. Missing semicolon?"));
                }

                cp += 1;
            }
        }

        if ml.is_some() && ml_pos != ml_len {
            return Err(ParseError("Too many entries in ML tag"));
        }

        state.ml = ml.unwrap_or_default();
        Ok(state)
    }

    /// Whether the read has no modification calls.
    pub fn is_empty(&self) -> bool {
        self.mods.is_empty()
    }

    /// `bam_mods_at_next_pos`: the modifications at the next query position, written to `out`.
    ///
    /// # Errors
    ///
    /// Returns an error if the MM tag is malformed.
    pub fn next_position(
        &mut self,
        read: &Read,
        out: &mut Vec<BaseModification>,
    ) -> Result<(), ParseError> {
        out.clear();

        if self.seq_pos >= read.sequence_len() {
            return Ok(());
        }

        let reverse = read.flags() & REVERSE != 0;
        let mut base = read.base(self.seq_pos);
        self.seq_pos += 1;

        if reverse {
            base = SEQI_RC[usize::from(base)];
        }

        let mut i = 0;

        while i < self.mods.len() {
            let state = &mut self.mods[i];

            if state.canonical != base && state.canonical != 15 {
                i += 1;
                continue;
            }

            let remaining = state.count;
            state.count -= 1;

            if remaining > 0 {
                i += 1;
                continue;
            }

            let mm_ptr = state.mm;
            out.push(self.call(i));

            self.step_ml(i, reverse);

            let mm = &self.mm;
            let state = &mut self.mods[i];

            if reverse {
                let mut cp = state.mm_end;

                if cp == 0 || cp - 1 < state.mm {
                    return Err(ParseError(
                        "Assert failed while processing base modification states",
                    ));
                }

                cp -= 1;

                while cp != state.mm && mm[cp] != b',' {
                    cp -= 1;
                }

                state.mm_end = cp;
                state.count = if cp != state.mm {
                    str2uint(mm, cp + 1)
                        .ok_or(ParseError("Error parsing unsigned integer from MM tag"))?
                        .0
                } else {
                    i64::from(i32::MAX)
                };
            } else if mm[state.mm] == b',' {
                let (n, end) = str2uint(mm, state.mm + 1)
                    .ok_or(ParseError("Error parsing unsigned integer from MM tag"))?;
                state.count = n;
                state.mm = end;
            } else {
                state.count = i64::from(i32::MAX);
            }

            let (count, mm_now) = (self.mods[i].count, self.mods[i].mm);

            // Other codes of the same MM entry are called at the same bases.
            let mut j = i + 1;

            while j < self.mods.len() && self.mods[j].mm == mm_ptr {
                out.push(self.call(j));
                self.mods[j].count = count;
                self.mods[j].mm = mm_now;
                self.step_ml(j, reverse);
                j += 1;
            }

            i = j;
        }

        Ok(())
    }

    /// `bam_mods_at_qpos`: the modifications at query position `qpos`, which must not be before
    /// the last one asked for. Positions in between are walked past.
    ///
    /// # Errors
    ///
    /// Returns an error if the MM tag is malformed.
    pub fn at_query_position(
        &mut self,
        read: &Read,
        qpos: usize,
        out: &mut Vec<BaseModification>,
    ) -> Result<(), ParseError> {
        out.clear();

        // htslib loops forever past the end of the sequence.
        let end = qpos.min(read.sequence_len().saturating_sub(1));

        while self.seq_pos <= end && self.seq_pos < read.sequence_len() {
            self.next_position(read, out)?;
        }

        Ok(())
    }

    fn call(&self, i: usize) -> BaseModification {
        let state = &self.mods[i];

        BaseModification {
            modified_base: state.ty,
            canonical_base: SEQ_NT16_STR[usize::from(state.canonical)],
            strand: state.strand,
            qual: state
                .ml
                .and_then(|k| usize::try_from(k).ok())
                .and_then(|k| self.ml.get(k))
                .map_or(-1, |&q| i32::from(q)),
        }
    }

    fn step_ml(&mut self, i: usize, reverse: bool) {
        let state = &mut self.mods[i];

        if let Some(k) = &mut state.ml {
            if reverse {
                *k -= state.stride as isize;
            } else {
                *k += state.stride as isize;
            }
        }
    }
}

/// Writes modification calls as samtools mpileup prints them: nothing if there are none, else
/// `[`, then for each call its strand (`+` or `-`), its code (or ChEBI number in parentheses),
/// and its likelihood if known, then `]`.
pub fn write_calls(dst: &mut Vec<u8>, calls: &[BaseModification]) {
    if calls.is_empty() {
        return;
    }

    dst.push(b'[');

    for call in calls {
        dst.push(if call.strand == 0 { b'+' } else { b'-' });

        if call.modified_base < 0 {
            write!(dst, "({})", -call.modified_base).expect("writing to a Vec cannot fail");
        } else {
            dst.push(call.modified_base as u8);
        }

        if call.qual >= 0 {
            write!(dst, "{}", call.qual).expect("writing to a Vec cannot fail");
        }
    }

    dst.push(b']');
}

/// htslib's `hts_str2uint(s, &end, 31, &failed)`: an optional `+` and decimal digits. Returns
/// the value and the end, or `None` if the value does not fit in 31 bits.
fn str2uint(s: &[u8], start: usize) -> Option<(i64, usize)> {
    const LIMIT: i64 = (1 << 31) - 1;

    let mut i = start + usize::from(s[start] == b'+');
    let mut n = 0i64;
    let mut failed = false;

    while s[i].is_ascii_digit() {
        n = n * 10 + i64::from(s[i] - b'0');
        failed |= n > LIMIT;
        n = n.min(LIMIT);
        i += 1;
    }

    (!failed).then_some((n, i))
}

#[cfg(test)]
mod tests {
    use super::*;

    fn read(seq: &[u8], reverse: bool, tags: &[u8]) -> Read {
        let codes: Vec<u8> = seq
            .iter()
            .map(|&b| crate::baq::SEQ_NT16_TABLE[usize::from(b)])
            .collect();

        let mut raw = Vec::new();
        raw.extend(0i32.to_le_bytes());
        raw.extend(0i32.to_le_bytes());
        raw.push(2);
        raw.push(60);
        raw.extend(0u16.to_le_bytes());
        raw.extend(1u16.to_le_bytes());
        raw.extend((if reverse { 0x10u16 } else { 0 }).to_le_bytes());
        raw.extend((seq.len() as u32).to_le_bytes());
        raw.extend((-1i32).to_le_bytes());
        raw.extend((-1i32).to_le_bytes());
        raw.extend(0i32.to_le_bytes());
        raw.extend(b"r\0");
        raw.extend(((seq.len() as u32) << 4).to_le_bytes());

        for pair in codes.chunks(2) {
            raw.push((pair[0] << 4) | pair.get(1).copied().unwrap_or(0));
        }

        raw.extend(std::iter::repeat_n(30, seq.len()));
        raw.extend_from_slice(tags);
        Read::new(raw).expect("valid record")
    }

    fn calls(read: &Read) -> Vec<(usize, Vec<(char, i32)>)> {
        let mut state = BaseMods::parse(read).expect("valid tags");
        let mut out = Vec::new();
        let mut calls = Vec::new();

        for qpos in 0..read.sequence_len() {
            state
                .at_query_position(read, qpos, &mut out)
                .expect("valid");

            if !out.is_empty() {
                calls.push((
                    qpos,
                    out.iter()
                        .map(|m| (m.modified_base as u8 as char, m.qual))
                        .collect(),
                ));
            }
        }

        calls
    }

    #[test]
    fn test_forward() {
        // C+m: skip one C, call the next; then call the next C.
        let mut tags = b"MMZC+m,1,0;\0".to_vec();
        tags.extend(b"MLBC");
        tags.extend(2u32.to_le_bytes());
        tags.extend([200, 50]);

        let r = read(b"ACGCCA", false, &tags);
        assert_eq!(calls(&r), [(3, vec![('m', 200)]), (4, vec![('m', 50)])]);
    }

    #[test]
    fn test_linked_codes() {
        let mut tags = b"MMZC+mh,0;\0".to_vec();
        tags.extend(b"MLBC");
        tags.extend(2u32.to_le_bytes());
        tags.extend([10, 20]);

        let r = read(b"ACG", false, &tags);
        assert_eq!(calls(&r), [(1, vec![('m', 10), ('h', 20)])]);
    }

    #[test]
    fn test_reverse() {
        // On the reverse strand, C+m counts Gs from the right: call the right-most G.
        let mut tags = b"MMZC+m,0;\0".to_vec();
        tags.extend(b"MLBC");
        tags.extend(1u32.to_le_bytes());
        tags.extend([99]);

        let r = read(b"AGTGA", true, &tags);
        assert_eq!(calls(&r), [(3, vec![('m', 99)])]);
    }

    #[test]
    fn test_parse_errors() {
        assert!(BaseMods::parse(&read(b"ACG", false, b"MMZX+m,0;\0")).is_err());
        assert!(BaseMods::parse(&read(b"ACG", false, b"MMZC+m,0\0")).is_err());
    }
}
