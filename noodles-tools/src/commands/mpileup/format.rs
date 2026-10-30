//! Formatting pileup columns.

use std::io::{self, Write};

use noodles_pileup as pileup;
use noodles_pileup::{Entry, Pileup, Read, baq::SEQ_NT16_TABLE, mods::BaseMods};

use super::{
    Mods, OutputOptions, PRINT_FLAG, PRINT_LAST, PRINT_MAPQ, PRINT_MAPQ_CHAR, PRINT_MODS,
    PRINT_PNEXT, PRINT_POS, PRINT_QNAME, PRINT_QPOS, PRINT_QPOS5, PRINT_RLEN, PRINT_RNAME,
    PRINT_RNEXT,
};

/// `print_empty_pileup`
pub(super) fn write_empty<W: Write>(
    out: &mut W,
    line: &mut Vec<u8>,
    options: &OutputOptions,
    name: &[u8],
    pos: i64,
    n: usize,
    reference: Option<&[u8]>,
) -> io::Result<()> {
    line.clear();
    line.extend_from_slice(name);
    line.push(b'\t');
    push_int(line, pos + 1);
    line.push(b'\t');
    line.push(match reference {
        Some(r) if pos < r.len() as i64 => r[pos as usize],
        _ => b'N',
    });

    for _ in 0..n {
        line.extend_from_slice(b"\t0\t*\t*");
        push_empty_columns(line, options);
    }

    line.push(b'\n');
    out.write_all(line)
}

pub(super) fn push_empty_columns(line: &mut Vec<u8>, options: &OutputOptions) {
    if options.flags == 0 && options.aux_tags.is_empty() {
        return;
    }

    let mut flag = PRINT_MAPQ_CHAR;

    while flag < PRINT_LAST {
        if flag != PRINT_MODS && options.flags & flag != 0 {
            line.extend_from_slice(b"\t*");
        }

        flag <<= 1;
    }

    for _ in &options.aux_tags {
        line.extend_from_slice(b"\t*");
    }
}

/// The base quality samtools filters and prints: that of the base at `qpos`, or 0 past the end.
pub(super) fn base_quality(read: &Read, entry: &Entry) -> i32 {
    if (entry.qpos as usize) < read.sequence_len() {
        i32::from(read.quality_scores()[entry.qpos as usize])
    } else {
        0
    }
}

/// The columns of one input at one position; `pileup` is `None` for an input without reads
/// there.
#[expect(clippy::too_many_arguments, reason = "the line's context")]
pub(super) fn write_sample(
    line: &mut Vec<u8>,
    pileup: Option<&mut Pileup>,
    options: &OutputOptions,
    names: &[Vec<u8>],
    pos: i64,
    reference: Option<&[u8]>,
    seq: &mut Vec<u8>,
    qual: &mut Vec<u8>,
    ins: &mut Vec<u8>,
) -> anyhow::Result<()> {
    seq.clear();
    qual.clear();

    let Some(pileup) = pileup else {
        line.extend_from_slice(b"\t0\t*\t*");
        push_empty_columns(line, options);
        return Ok(());
    };

    let mut cnt = 0i64;
    let n_entries = pileup.entries().len();

    if options.mods == Mods::None {
        for entry in &pileup.entries()[..n_entries] {
            let read = pileup.read(entry);
            let c = base_quality(read, entry);

            if c >= options.min_bq {
                push_base(seq, read, None, entry, pos, reference, options, ins)?;
                qual.push(if c + 33 < 126 { (c + 33) as u8 } else { 126 });
                cnt += 1;
            }
        }
    } else {
        // Modification states advance as bases are printed, so they are borrowed mutably.
        for k in 0..n_entries {
            let entry = pileup.entries()[k];
            let (read, mods) = pileup.read_and_base_mods(&entry);
            let c = base_quality(read, &entry);

            if c >= options.min_bq {
                push_base(seq, read, mods, &entry, pos, reference, options, ins)?;
                qual.push(if c + 33 < 126 { (c + 33) as u8 } else { 126 });
                cnt += 1;
            }
        }
    }

    line.push(b'\t');
    push_int(line, cnt);
    line.push(b'\t');

    if n_entries == 0 {
        line.extend_from_slice(b"*\t*");
        push_empty_columns(line, options);
        return Ok(());
    }

    if seq.is_empty() {
        line.push(b'*');
    } else {
        line.extend_from_slice(seq);
    }

    line.push(b'\t');

    if qual.is_empty() {
        line.push(b'*');
    } else {
        line.extend_from_slice(qual);
    }

    let mut flag = if options.flags == 0 {
        PRINT_LAST
    } else {
        PRINT_MAPQ_CHAR
    };

    while flag < PRINT_LAST {
        if flag != PRINT_MODS && options.flags & flag != 0 {
            line.push(b'\t');
            let mut n = 0;

            for entry in &pileup.entries()[..n_entries] {
                let read = pileup.read(entry);

                if base_quality(read, entry) < options.min_bq {
                    continue;
                }

                if n > 0 && flag != PRINT_MAPQ_CHAR {
                    line.push(b',');
                }

                n += 1;

                match flag {
                    PRINT_MAPQ_CHAR => {
                        line.push((read.mapping_quality() as u32 + 33).min(126) as u8)
                    }
                    PRINT_QPOS => push_int(line, entry.qpos + 1),
                    PRINT_QPOS5 => {
                        let pos5 = if read.is_reverse_complemented() {
                            read.sequence_len() as i64 - entry.qpos + i64::from(entry.is_del)
                        } else {
                            entry.qpos + 1
                        };
                        push_int(line, pos5);
                    }
                    PRINT_QNAME => line.extend_from_slice(read.name()),
                    PRINT_FLAG => push_int(line, i64::from(read.flags())),
                    PRINT_RNAME => {
                        if read.reference_sequence_id() >= 0 {
                            line.extend_from_slice(&names[read.reference_sequence_id() as usize]);
                        } else {
                            line.push(b'*');
                        }
                    }
                    PRINT_POS => push_int(line, read.position() + 1),
                    PRINT_MAPQ => push_int(line, i64::from(read.mapping_quality())),
                    PRINT_RNEXT => {
                        if read.mate_reference_sequence_id() >= 0 {
                            line.extend_from_slice(
                                &names[read.mate_reference_sequence_id() as usize],
                            );
                        } else {
                            line.push(b'*');
                        }
                    }
                    PRINT_PNEXT => push_int(line, read.mate_position() + 1),
                    PRINT_RLEN => push_int(line, read.sequence_len() as i64),
                    _ => {}
                }
            }

            if n == 0 {
                line.push(b'*');
            }
        }

        flag <<= 1;
    }

    for &tag in &options.aux_tags {
        line.push(b'\t');
        let mut n = 0;

        for entry in &pileup.entries()[..n_entries] {
            let read = pileup.read(entry);

            if base_quality(read, entry) < options.min_bq {
                continue;
            }

            if n > 0 {
                line.push(options.sep);
            }

            n += 1;
            push_tag_value(line, read, tag, options.empty);
        }

        if n == 0 {
            line.push(b'*');
        }
    }

    Ok(())
}

/// An --output-extra tag value, as samtools prints it.
pub(super) fn push_tag_value(line: &mut Vec<u8>, read: &Read, tag: [u8; 2], empty: u8) {
    let Some((ty, value)) = read.tag(tag) else {
        line.push(empty);
        return;
    };

    match ty {
        b'Z' | b'H' => line.extend_from_slice(&value[..value.len() - 1]),
        b'c' => push_int(line, i64::from(value[0].cast_signed())),
        b'C' => push_int(line, i64::from(value[0])),
        b's' => push_int(line, i64::from(i16::from_le_bytes([value[0], value[1]]))),
        b'S' => push_int(line, i64::from(u16::from_le_bytes([value[0], value[1]]))),
        b'i' => push_int(
            line,
            i64::from(i32::from_le_bytes([value[0], value[1], value[2], value[3]])),
        ),
        b'I' => push_int(
            line,
            i64::from(u32::from_le_bytes([value[0], value[1], value[2], value[3]])),
        ),
        b'f' => push_kputd(
            line,
            f64::from(f32::from_le_bytes([value[0], value[1], value[2], value[3]])),
        ),
        b'A' => line.push(value[0]),
        _ => line.push(b'*'),
    }
}

/// `pileup_seq`: one read's base column entry.
#[expect(clippy::too_many_arguments, reason = "the read's context")]
#[inline]
pub(super) fn push_base(
    seq: &mut Vec<u8>,
    read: &Read,
    mut mods: Option<&mut BaseMods>,
    entry: &Entry,
    pos: i64,
    reference: Option<&[u8]>,
    options: &OutputOptions,
    ins: &mut Vec<u8>,
) -> anyhow::Result<()> {
    pub(super) const LOWER: &[u8; 16] = b",acmgrsvtwyhkdbn";
    pub(super) const UPPER: &[u8; 16] = b".ACMGRSVTWYHKDBN";

    let reverse = read.is_reverse_complemented();

    if !options.no_ends && entry.is_head {
        seq.push(b'^');
        let mapq = read.mapping_quality();
        seq.push(if mapq > 93 { 126 } else { mapq + 33 });
    }

    if entry.is_del {
        seq.push(if entry.is_refskip {
            if reverse { b'<' } else { b'>' }
        } else if reverse && options.rev_del {
            b'#'
        } else {
            b'*'
        });
    } else {
        let mut c = if (entry.qpos as usize) < read.sequence_len() {
            read.base(entry.qpos as usize)
        } else {
            15
        };

        if let Some(r) = reference {
            let rb = if pos < r.len() as i64 {
                SEQ_NT16_TABLE[usize::from(r[pos as usize])]
            } else {
                15
            };

            if c == rb {
                c = 0;
            }
        }

        seq.push(if reverse { LOWER } else { UPPER }[usize::from(c)]);

        if let Some(state) = mods.as_deref_mut() {
            let mut calls = Vec::new();
            state.at_query_position(read, entry.qpos as usize, &mut calls)?;
            pileup::mods::write_calls(seq, &calls);
        }
    }

    let mut del_len = if entry.indel < 0 {
        (-entry.indel) as usize
    } else {
        0
    };

    if entry.indel > 0 {
        let mods = mods.filter(|_| options.mods == Mods::BasesAndInsertions);
        let (len, after) = pileup::insertion_with_mods(read, entry, mods, ins)?;
        del_len = after;

        if options.no_ins < 2 {
            seq.push(b'+');
            push_int(seq, len as i64);
        }

        if options.no_ins == 0 {
            // Modifications (in brackets) keep their case.
            let mut in_mod = false;
            let pad = if reverse && options.rev_del {
                b'#'
            } else {
                b'*'
            };

            seq.extend(ins.iter().map(|&b| {
                if b == b'[' {
                    in_mod = true;
                } else if b == b']' {
                    in_mod = false;
                }

                if in_mod {
                    b
                } else if !reverse {
                    b.to_ascii_uppercase()
                } else if b == b'*' {
                    pad
                } else {
                    b.to_ascii_lowercase()
                }
            }));
        }
    }

    if del_len > 0 {
        if options.no_del < 2 {
            push_int(seq, -(del_len as i64));
        }

        if options.no_del == 0 {
            for j in 1..=del_len as i64 {
                let c = match reference {
                    Some(r) if ((pos as i32 as i64) + j) < r.len() as i64 => r[(pos + j) as usize],
                    _ => b'N',
                };

                seq.push(if reverse {
                    c.to_ascii_lowercase()
                } else {
                    c.to_ascii_uppercase()
                });
            }
        }
    }

    if !options.no_ends && entry.is_tail {
        seq.push(b'$');
    }

    Ok(())
}

pub(super) fn push_int(dst: &mut Vec<u8>, n: i64) {
    let mut buf = [0; lexical_core::BUFFER_SIZE];
    dst.extend_from_slice(lexical_core::write(n, &mut buf));
}

/// htslib's `kputd`: up to six significant digits, trailing zeros dropped, and `%g` outside
/// [0.0001, 999999].
pub(super) fn push_kputd(dst: &mut Vec<u8>, d: f64) {
    if d == 0.0 {
        dst.extend_from_slice(if d.is_sign_negative() { b"-0" } else { b"0" });
        return;
    }

    let mut d = d;

    if d < 0.0 {
        dst.push(b'-');
        d = -d;
    }

    if !(0.0001..=999_999.0).contains(&d) {
        dst.extend_from_slice(crate::depth::format_g(d, 6).as_bytes());
        return;
    }

    // Six significant digits as an integer, and how many of them follow the decimal point.
    let (scale, decimals) = if d < 0.001 {
        (1e9, 9)
    } else if d < 0.01 {
        (1e8, 8)
    } else if d < 0.1 {
        (1e7, 7)
    } else if d < 1.0 {
        (1e6, 6)
    } else if d < 10.0 {
        (1e5, 5)
    } else if d < 100.0 {
        (1e4, 4)
    } else if d < 1000.0 {
        (1e3, 3)
    } else if d < 10_000.0 {
        (1e2, 2)
    } else if d < 100_000.0 {
        (1e1, 1)
    } else {
        (1.0, 0)
    };

    #[expect(
        clippy::cast_possible_truncation,
        clippy::cast_sign_loss,
        reason = "kputd stores rint() of the scaled value in a uint32_t"
    )]
    let i = (d * scale).round_ties_even() as u32;

    // The fixed-point text of i / scale with trailing zeros culled, as kputd's digit shuffling
    // produces.
    let digits = format!("{i:0width$}", width = decimals + 1);
    let (int_part, frac_part) = digits.split_at(digits.len() - decimals);
    let frac = frac_part.trim_end_matches('0');

    dst.extend_from_slice(int_part.as_bytes());

    if !frac.is_empty() {
        dst.push(b'.');
        // kputd keeps at most six significant digits after the leading zeros.
        dst.extend_from_slice(frac.as_bytes());
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_push_kputd() {
        let t = |d: f64| {
            let mut buf = Vec::new();
            push_kputd(&mut buf, d);
            String::from_utf8(buf).expect("ASCII")
        };

        assert_eq!(t(0.0), "0");
        assert_eq!(t(-0.0), "-0");
        assert_eq!(t(1.5), "1.5");
        assert_eq!(t(-2.25), "-2.25");
        assert_eq!(t(100.0), "100");
        assert_eq!(t(0.000_123_456_7), "0.000123457");
        assert_eq!(t(0.999_999_9), "1");
        assert_eq!(t(123_456.7), "123457");
        assert_eq!(t(1e-5), "1e-05");
        assert_eq!(t(1e7), "1e+07");
    }
}
