//! `ntools flags`: convert between numeric and textual flags (samtools flags).

use std::{
    io::{self, Write},
    process::ExitCode,
};

/// Flag bits, names, and descriptions, in bit order (`bam_flag2str`).
const FLAGS: [(u16, &str, &str); 12] = [
    (
        0x1,
        "PAIRED",
        "paired-end / multiple-segment sequencing technology",
    ),
    (
        0x2,
        "PROPER_PAIR",
        "each segment properly aligned according to aligner",
    ),
    (0x4, "UNMAP", "segment unmapped"),
    (0x8, "MUNMAP", "next segment in the template unmapped"),
    (0x10, "REVERSE", "SEQ is reverse complemented"),
    (
        0x20,
        "MREVERSE",
        "SEQ of next segment in template is rev.complemented",
    ),
    (0x40, "READ1", "the first segment in the template"),
    (0x80, "READ2", "the last segment in the template"),
    (0x100, "SECONDARY", "secondary alignment"),
    (
        0x200,
        "QCFAIL",
        "not passing quality controls or other filters",
    ),
    (0x400, "DUP", "PCR or optical duplicate"),
    (0x800, "SUPPLEMENTARY", "supplementary alignment"),
];

#[derive(Debug, clap::Args)]
pub struct Args {
    /// Flags, each an integer (decimal, 0x hexadecimal, or 0 octal) or a comma-separated list of
    /// names.
    #[arg(allow_hyphen_values = true)]
    flags: Vec<String>,
}

pub fn run(args: Args) -> anyhow::Result<ExitCode> {
    let mut out = io::stdout().lock();

    if args.flags.is_empty() {
        write_usage(&mut out)?;
        return Ok(ExitCode::SUCCESS);
    }

    for arg in &args.flags {
        let Some(mask) = parse(arg).filter(|&mask| mask >= 0) else {
            out.flush()?;
            eprintln!("[flags] Could not parse \"{arg}\"");
            write_usage(&mut io::stderr().lock())?;
            return Ok(ExitCode::FAILURE);
        };

        writeln!(out, "{mask:#x}\t{mask}\t{}", names(mask))?;
    }

    Ok(ExitCode::SUCCESS)
}

fn write_usage<W: Write>(out: &mut W) -> io::Result<()> {
    out.write_all(
        b"About: Convert between textual and numeric flag representation\n\
Usage: ntools flags FLAGS...\n\
\n\
Each FLAGS argument is either an INT (in decimal/hexadecimal/octal) representing\n\
a combination of the following numeric flag values, or a comma-separated string\n\
NAME,...,NAME representing a combination of the following flag names:\n\
\n",
    )?;

    for (bit, name, description) in FLAGS {
        writeln!(out, "{:>#6x} {bit:5}  {name:<15}{description}", bit)?;
    }

    Ok(())
}

/// `bam_str2flag`: an integer as C's `strtol(s, &end, 0)` reads it (truncated to `int`), or
/// names, matched without case; `None` if a name is unknown.
fn parse(s: &str) -> Option<i32> {
    if let Some(n) = strtol(s) {
        return Some(n as i32);
    }

    if s.is_empty() {
        return Some(0);
    }

    s.split(',').try_fold(0, |flags, name| {
        let (bit, _, _) = FLAGS
            .iter()
            .find(|(_, n, _)| n.eq_ignore_ascii_case(name))?;
        Some(flags | i32::from(*bit))
    })
}

/// C's `strtol(s, &end, 0)`, or `None` if it reads no digits. Out-of-range values saturate.
fn strtol(s: &str) -> Option<i64> {
    let t = s.trim_start_matches(|c: char| c.is_ascii_whitespace());

    let (negative, t) = match t.as_bytes().first() {
        Some(b'-') => (true, &t[1..]),
        Some(b'+') => (false, &t[1..]),
        _ => (false, t),
    };

    let (radix, digits) = if let Some(hex) = t.strip_prefix("0x").or_else(|| t.strip_prefix("0X"))
        && hex.starts_with(|c: char| c.is_ascii_hexdigit())
    {
        (16, hex)
    } else if t.starts_with('0') {
        (8, t)
    } else {
        (10, t)
    };

    let len = digits.chars().take_while(|c| c.is_digit(radix)).count();

    if len == 0 {
        return None;
    }

    let n = digits[..len].chars().fold(0i64, |n, c| {
        n.saturating_mul(i64::from(radix))
            .saturating_add(i64::from(c.to_digit(radix).expect("a digit")))
    });

    Some(if negative { n.saturating_neg() } else { n })
}

/// `bam_flag2str`: the names of the set flags, comma separated.
fn names(mask: i32) -> String {
    FLAGS
        .iter()
        .filter(|(bit, _, _)| mask & i32::from(*bit) != 0)
        .map(|(_, name, _)| *name)
        .collect::<Vec<_>>()
        .join(",")
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_parse() {
        assert_eq!(parse("99"), Some(99));
        assert_eq!(parse("0x10"), Some(16));
        assert_eq!(parse("010"), Some(8));
        assert_eq!(parse("12abc"), Some(12));
        assert_eq!(parse("paired,Read1"), Some(0x41));
        assert_eq!(parse(""), Some(0));
        assert_eq!(parse("PAIRED,"), None);
        assert_eq!(parse("PAIR"), None);
        assert_eq!(parse("-1"), Some(-1));
    }

    #[test]
    fn test_names() {
        assert_eq!(names(0x41), "PAIRED,READ1");
        assert_eq!(names(0), "");
    }
}
