//! The keys that identify duplicates: one per single read, one per pair.

use anyhow::bail;
use noodles_bam as bam;
use noodles_sam as sam;
use noodles_sam::alignment::record::data::field::Value;

use super::{
    Key, MATE_CIGAR, Mode, O_FF, O_FR, O_RF, O_RR, R_LE, R_RI,
    record::{
        cigar_ops, flags_of, i32_field, parse_cigar_string, unclipped_end, unclipped_other_end,
        unclipped_other_start, unclipped_start,
    },
};
use crate::bam_raw::field;

/// `make_single_key`
pub(super) fn make_single_key(record: &bam::Record, read_group: i32, barcode: i32) -> Key {
    let ops = cigar_ops(record).unwrap_or_default();
    let raw = record.as_ref();
    let this_ref = i32_field(raw, field::REFERENCE_SEQUENCE_ID) + 1;

    let (this_coord, orientation) = if flags_of(raw).is_reverse_complemented() {
        (unclipped_end(record, &ops), O_RR)
    } else {
        (unclipped_start(record, &ops), O_FF)
    };

    Key {
        this_coord,
        this_ref,
        barcode,
        read_group,
        single: true,
        orientation,
        ..Key::default()
    }
}

/// `make_pair_key`
pub(super) fn make_pair_key(
    record: &bam::Record,
    mode: Mode,
    read_group: i32,
    barcode: i32,
) -> anyhow::Result<Key> {
    let ops = cigar_ops(record)?;
    let raw = record.as_ref();
    let flags = flags_of(raw);

    let this_ref = i32_field(raw, field::REFERENCE_SEQUENCE_ID) + 1;
    let other_ref = i32_field(raw, field::MATE_REFERENCE_SEQUENCE_ID) + 1;
    let mpos = i64::from(i32_field(raw, 24..28));
    let pos = i64::from(i32_field(raw, field::POSITION));

    let mut this_coord = unclipped_start(record, &ops);
    let this_end = unclipped_end(record, &ops);

    let data = record.data();
    let mate_cigar = match sam::alignment::record::Data::get(&data, &MATE_CIGAR).transpose()? {
        Some(Value::String(s)) => parse_cigar_string(s.as_ref()),
        Some(_) => bail!(
            "[markdup] error, MC tag wrong type. Please use the MC tag provided by samtools fixmate."
        ),
        None => bail!("[markdup] error, no MC tag. Please run samtools fixmate on file first."),
    };

    let other_end = unclipped_other_end(mpos, &mate_cigar);
    let mut other_coord = unclipped_other_start(mpos, &mate_cigar);

    let is_rev = flags.is_reverse_complemented();
    let is_mrev = flags.is_mate_reverse_complemented();
    let is_read1 = flags.is_first_segment();

    let leftmost;
    let orientation;

    match mode {
        Mode::T => {
            leftmost = if this_ref != other_ref {
                this_ref < other_ref
            } else if is_rev == is_mrev {
                if is_rev {
                    this_end <= other_end
                } else {
                    this_coord <= other_coord
                }
            } else if is_rev {
                this_end <= other_coord
            } else {
                this_coord <= other_end
            };

            if leftmost {
                if is_rev == is_mrev {
                    other_coord = other_end;
                    orientation = match (is_rev, is_read1) {
                        (false, true) | (true, false) => O_FF,
                        (false, false) | (true, true) => O_RR,
                    };
                } else if !is_rev {
                    orientation = O_FR;
                    other_coord = other_end;
                } else {
                    orientation = O_RF;
                    this_coord = this_end;
                }
            } else if is_rev == is_mrev {
                this_coord = this_end;
                orientation = match (is_rev, is_read1) {
                    (false, true) | (true, false) => O_RR,
                    (false, false) | (true, true) => O_FF,
                };
            } else if !is_rev {
                orientation = O_RF;
                other_coord = other_end;
            } else {
                orientation = O_FR;
                this_coord = this_end;
            }
        }
        Mode::S => {
            let difference = if this_ref != other_ref {
                i64::from(this_ref - other_ref)
            } else if is_rev == is_mrev {
                if is_rev {
                    this_end - other_end
                } else {
                    this_coord - other_coord
                }
            } else if is_rev {
                this_end - other_coord
            } else {
                this_coord - other_end
            };

            leftmost = match difference.cmp(&0) {
                std::cmp::Ordering::Less => true,
                std::cmp::Ordering::Greater => false,
                std::cmp::Ordering::Equal => {
                    if pos == mpos {
                        is_read1
                    } else {
                        pos < mpos
                    }
                }
            };

            orientation = match (leftmost, is_rev == is_mrev, is_rev) {
                (true, true, false) => O_FF,
                (true, true, true) => O_RR,
                (true, false, false) => O_FR,
                (true, false, true) => O_RF,
                (false, true, false) => O_RR,
                (false, true, true) => O_FF,
                (false, false, false) => O_RF,
                (false, false, true) => O_FR,
            };

            this_coord = if is_rev {
                unclipped_end(record, &ops)
            } else {
                unclipped_start(record, &ops)
            };

            other_coord = if is_mrev {
                unclipped_other_end(mpos, &mate_cigar)
            } else {
                unclipped_other_start(mpos, &mate_cigar)
            };
        }
    }

    Ok(Key {
        this_coord,
        other_coord,
        this_ref,
        other_ref,
        barcode,
        read_group,
        single: false,
        leftmost: if leftmost { R_LE } else { R_RI },
        orientation,
    })
}
