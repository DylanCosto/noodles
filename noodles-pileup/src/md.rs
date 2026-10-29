//! MD and NM tags from a reference sequence (samtools 1.24's `bam_fillmd1_core`).

use crate::{Read, baq::SEQ_NT16_TABLE};

const UNMAPPED: u16 = 0x4;

/// What [`fill_md`] changes besides computing the tags (samtools calmd's options).
#[derive(Clone, Copy, Debug)]
pub struct FillMdOptions {
    /// Update the MD and NM tags (`-N` turns this off).
    pub update_tags: bool,
    /// Change bases identical to the reference to `=` (`-e`).
    pub use_equal: bool,
    /// Mask the matching bases (as `N`, quality 0) of reads with at least this many mismatches;
    /// 0 disables (`-n`).
    pub max_nm: i32,
    /// Drop all tags but RG (`-d`).
    pub drop_tags: bool,
    /// Bin base qualities to `q / 10 * 10 + 7` (`-q`).
    pub bin_qualities: bool,
}

impl Default for FillMdOptions {
    fn default() -> Self {
        Self {
            update_tags: true,
            use_equal: false,
            max_nm: 0,
            drop_tags: false,
            bin_qualities: false,
        }
    }
}

/// Existing tags that [`fill_md`] replaced, as samtools reports them.
#[derive(Clone, Debug, Default, Eq, PartialEq)]
pub struct FillMdReport {
    /// The old and new NM values, if NM changed.
    pub nm: Option<(i64, i32)>,
    /// The old and new MD values, if MD changed (ignoring case).
    pub md: Option<(Vec<u8>, Vec<u8>)>,
}

/// Computes the MD string and edit distance (NM) of a read against `reference` (the whole
/// sequence of its reference sequence), and updates the tags and the read as `options` says.
/// Returns `None`, changing nothing, if the read has no sequence.
///
/// New tags are appended; changed ones are removed and appended, and unchanged ones are kept in
/// place.
pub fn fill_md(read: &mut Read, reference: &[u8], options: &FillMdOptions) -> Option<FillMdReport> {
    if read.l_seq == 0 {
        return None;
    }

    // The CIGAR is read while the sequence changes.
    let cigar = std::mem::take(&mut read.cigar);
    let report = fill_md_with(read, &cigar, reference, options);
    read.cigar = cigar;
    Some(report)
}

fn fill_md_with(
    read: &mut Read,
    cigar: &[u32],
    reference: &[u8],
    options: &FillMdOptions,
) -> FillMdReport {
    let ref_len = reference.len() as i64;
    let l_qseq = read.l_seq;
    let in_bounds = |rpos: i64| rpos < ref_len && reference[rpos as usize] != 0;

    // A match, as samtools counts it: the same unambiguous base, or a base already `=`.
    let is_match = |read: &Read, z: usize, r: i64| {
        let c1 = read.base(z);
        let c2 = SEQ_NT16_TABLE[usize::from(reference[r as usize])];
        (c1 == c2 && c1 != 15 && c2 != 15) || c1 == 0
    };

    let mut md: Vec<u8> = Vec::new();
    let mut matched = 0i32;
    let mut nm = 0i32;
    let (mut qpos, mut rpos) = (0usize, read.pos);

    'ops: for &op in cigar {
        let (kind, oplen) = (op & 0xf, (op >> 4) as usize);

        match kind {
            0 | 7 | 8 => {
                for j in 0..oplen {
                    let z = qpos + j;
                    let r = rpos + j as i64;

                    if !in_bounds(r) || z >= l_qseq {
                        break 'ops;
                    }

                    if is_match(read, z, r) {
                        if options.use_equal {
                            read.set_base(z, 0);
                        }

                        matched += 1;
                    } else {
                        push_int(&mut md, matched);
                        md.push(reference[r as usize].to_ascii_uppercase());
                        matched = 0;
                        nm += 1;
                    }
                }

                rpos += oplen as i64;
                qpos += oplen;
            }
            2 => {
                push_int(&mut md, matched);
                md.push(b'^');

                let mut j = 0;

                while j < oplen && in_bounds(rpos + j as i64) {
                    md.push(reference[(rpos + j as i64) as usize].to_ascii_uppercase());
                    j += 1;
                }

                matched = 0;
                rpos += j as i64;
                nm += j as i32;

                if j < oplen {
                    break 'ops;
                }
            }
            1 | 4 => {
                qpos += oplen;

                if kind == 1 {
                    nm += oplen as i32;
                }
            }
            3 => rpos += oplen as i64,
            _ => {}
        }
    }

    push_int(&mut md, matched);

    if options.max_nm > 0 && nm >= options.max_nm {
        let (mut qpos, mut rpos) = (0usize, read.pos);

        'mask: for &op in cigar {
            let (kind, oplen) = (op & 0xf, (op >> 4) as usize);

            match kind {
                0 | 7 | 8 => {
                    for j in 0..oplen {
                        let z = qpos + j;
                        let r = rpos + j as i64;

                        if !in_bounds(r) || z >= l_qseq {
                            break 'mask;
                        }

                        if is_match(read, z, r) {
                            read.set_base(z, 15);
                            read.quality_scores_mut()[z] = 0;
                        }
                    }

                    rpos += oplen as i64;
                    qpos += oplen;
                }
                2 | 3 => rpos += oplen as i64,
                1 | 4 => qpos += oplen,
                _ => {}
            }
        }
    }

    let mut report = FillMdReport::default();

    if options.update_tags && read.flag & UNMAPPED == 0 {
        match read.integer_tag(*b"NM") {
            None => read.append_tag(*b"NM", b'i', &nm.to_le_bytes()),
            Some(old) => {
                if i64::from(nm) != old {
                    report.nm = Some((old, nm));
                    read.remove_tag(*b"NM");
                    read.append_tag(*b"NM", b'i', &nm.to_le_bytes());
                }
            }
        }

        match read.tag(*b"MD") {
            None => read.append_string_tag(*b"MD", &md),
            Some((_, value)) => {
                let old = &value[..value.iter().position(|&b| b == 0).unwrap_or(value.len())];

                if !old.eq_ignore_ascii_case(&md) {
                    report.md = Some((old.to_vec(), md.clone()));
                    read.remove_tag(*b"MD");
                    read.append_string_tag(*b"MD", &md);
                }
            }
        }
    }

    if options.drop_tags {
        read.retain_only_tag(*b"RG");
    }

    if options.bin_qualities {
        for q in read.quality_scores_mut() {
            if *q >= 3 {
                *q = *q / 10 * 10 + 7;
            }
        }
    }

    report
}

fn push_int(dst: &mut Vec<u8>, n: i32) {
    dst.extend_from_slice(n.to_string().as_bytes());
}
