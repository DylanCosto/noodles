#![warn(missing_docs)]

//! **noodles-pileup** is an htslib-compatible pileup engine for encoded BAM records, with base
//! alignment quality (BAQ), mapping quality capping, and MD/NM calculation.
//!
//! These are ports of htslib 1.24's `bam_plp_*`, `sam_prob_realn`, and `sam_cap_mapq`, and
//! samtools 1.24's `bam_fillmd1`. They keep the same rules and floating-point operations, so the
//! reads and qualities at each position, adjusted qualities, and tags match samtools byte for
//! byte.
//!
//! A [`Read`] wraps an encoded BAM record (for example, from [`noodles_bam::Record`]). Push reads
//! in coordinate order into a [`Pileup`], and after each push, call [`Pileup::next_position`] until it
//! returns `false`; each `true` prepares the [`Entry`]s at one position.

pub mod baq;
mod hash;
pub mod md;
pub mod mods;
mod pileup;
mod read;

pub use self::{
    pileup::{Entry, Error, Overlaps, Pileup, insertion, insertion_with_mods},
    read::Read,
};
