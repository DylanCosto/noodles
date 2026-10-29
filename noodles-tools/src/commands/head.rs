//! `ntools head`: print the header and first records of a file as SAM (samtools head).

use std::{
    io::{self, Write},
    path::PathBuf,
};

use anyhow::Context;
use noodles_sam as sam;

use crate::{
    alignment::{AlignmentReader, ReadOptions},
    output::{AlignmentWriter, OutputFormat, WriteOptions},
};

#[derive(Debug, clap::Args)]
#[command(disable_help_flag = true)]
pub struct Args {
    /// Input file [default: stdin].
    #[arg(default_value = "-")]
    input: String,

    /// Print only the first N header lines [default: all].
    #[arg(short = 'h', long = "headers", value_parser = parse_count)]
    headers: Option<u64>,

    /// Print the first N records.
    #[arg(short = 'n', long = "records", default_value = "0", value_parser = parse_count)]
    records: u64,

    /// Reference FASTA (indexed) for CRAM input.
    #[arg(short = 'T', long = "reference")]
    reference: Option<PathBuf>,

    /// Additional threads for decompression.
    #[arg(short = '@', long = "threads", default_value_t = 0)]
    threads: usize,

    /// Print help.
    #[arg(long = "help", action = clap::ArgAction::Help)]
    help: Option<bool>,
}

pub fn run(args: Args) -> anyhow::Result<()> {
    let read_options = ReadOptions {
        threads: args.threads,
        reference: args.reference.clone(),
    };

    let (mut reader, _) = AlignmentReader::open(&args.input, &read_options)?;
    let header = reader
        .read_header()
        .with_context(|| format!("[head] failed to read the header from \"{}\"", args.input))?;

    let mut text = Vec::new();
    sam::io::Writer::new(&mut text).write_header(&header)?;

    if let Some(n) = args.headers {
        let end = text
            .iter()
            .enumerate()
            .filter(|&(_, &b)| b == b'\n')
            .nth(usize::try_from(n).unwrap_or(usize::MAX).wrapping_sub(1))
            .map_or(text.len(), |(i, _)| i + 1);

        text.truncate(if n == 0 { 0 } else { end });
    }

    {
        let mut stdout = io::stdout().lock();
        stdout.write_all(&text)?;
        stdout.flush()?;
    }

    if args.records == 0 {
        return Ok(());
    }

    let write_options = WriteOptions {
        format: OutputFormat::Sam,
        level: crate::output::parse_level(0)?,
        threads: 0,
        reference: args.reference.clone(),
    };

    let mut writer = AlignmentWriter::create(None, &write_options, &header)?;
    let mut remaining = args.records;

    reader.for_each_while(&header, |rec| {
        writer.write(&header, &rec)?;
        remaining -= 1;
        Ok(remaining > 0)
    })?;

    writer.finish(&header)?;

    Ok(())
}

/// C's `strtoull(s, NULL, 0)`: decimal, `0x` hexadecimal, or `0` octal digits up to the first
/// other character, 0 if there are none, and a leading `-` negating modulo 2^64.
fn parse_count(s: &str) -> Result<u64, String> {
    let s = s.trim_start();

    let (negative, s) = match s.as_bytes().first() {
        Some(b'-') => (true, &s[1..]),
        Some(b'+') => (false, &s[1..]),
        _ => (false, s),
    };

    let (radix, digits) = if let Some(hex) = s.strip_prefix("0x").or_else(|| s.strip_prefix("0X"))
        && hex.starts_with(|c: char| c.is_ascii_hexdigit())
    {
        (16, hex)
    } else if s.starts_with('0') {
        (8, s)
    } else {
        (10, s)
    };

    let n = digits
        .chars()
        .map_while(|c| c.to_digit(radix))
        .try_fold(0u64, |n, d| {
            n.checked_mul(u64::from(radix))?.checked_add(u64::from(d))
        })
        .unwrap_or(u64::MAX);

    Ok(if negative { n.wrapping_neg() } else { n })
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_parse_count() {
        assert_eq!(parse_count("12"), Ok(12));
        assert_eq!(parse_count("0x1f"), Ok(31));
        assert_eq!(parse_count("010"), Ok(8));
        assert_eq!(parse_count("7abc"), Ok(7));
        assert_eq!(parse_count("abc"), Ok(0));
        assert_eq!(parse_count("-1"), Ok(u64::MAX));
    }
}
