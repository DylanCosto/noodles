//! Recording this tool's invocation as an `@PG` header line, as samtools does.

use std::env;

use anyhow::Context;
use noodles_sam::{
    self as sam,
    header::record::value::{
        Map,
        map::{Program, program::tag},
    },
};

const NAME: &str = "ntools";

/// Appends an `@PG` line for this invocation to every program chain in the header.
pub fn add(header: &mut sam::Header) -> anyhow::Result<()> {
    let command_line = env::args().collect::<Vec<_>>().join(" ");

    let map = Map::<Program>::builder()
        .insert(tag::NAME, NAME)
        .insert(tag::VERSION, env!("CARGO_PKG_VERSION"))
        .insert(tag::COMMAND_LINE, command_line)
        .build()
        .context("building @PG record")?;

    header
        .programs_mut()
        .add(NAME, map)
        .context("adding @PG record")
}
