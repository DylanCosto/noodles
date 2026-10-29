//! Parsing of samtools-style flag arguments.

use anyhow::{Context, bail};
use noodles_sam::alignment::record::Flags;

/// Flag names accepted by samtools, in bit order.
const NAMES: [(&str, Flags); 12] = [
    ("PAIRED", Flags::SEGMENTED),
    ("PROPER_PAIR", Flags::PROPERLY_SEGMENTED),
    ("UNMAP", Flags::UNMAPPED),
    ("MUNMAP", Flags::MATE_UNMAPPED),
    ("REVERSE", Flags::REVERSE_COMPLEMENTED),
    ("MREVERSE", Flags::MATE_REVERSE_COMPLEMENTED),
    ("READ1", Flags::FIRST_SEGMENT),
    ("READ2", Flags::LAST_SEGMENT),
    ("SECONDARY", Flags::SECONDARY),
    ("QCFAIL", Flags::QC_FAIL),
    ("DUP", Flags::DUPLICATE),
    ("SUPPLEMENTARY", Flags::SUPPLEMENTARY),
];

/// Parses a flag set given as a decimal, `0x` hexadecimal, or `0` octal integer, or as a
/// comma-separated list of names (`UNMAP,SECONDARY`).
pub fn parse(s: &str) -> anyhow::Result<Flags> {
    let s = s.trim();

    if s.starts_with(|c: char| c.is_ascii_digit()) {
        let n = if let Some(hex) = s.strip_prefix("0x").or_else(|| s.strip_prefix("0X")) {
            u16::from_str_radix(hex, 16)
        } else if s.len() > 1
            && let Some(oct) = s.strip_prefix('0')
        {
            u16::from_str_radix(oct, 8)
        } else {
            s.parse()
        }
        .with_context(|| format!("invalid flag value: {s}"))?;

        return Ok(Flags::from_bits_retain(n));
    }

    let mut flags = Flags::empty();

    for name in s.split(',') {
        let name = name.trim();

        match NAMES.iter().find(|(n, _)| n.eq_ignore_ascii_case(name)) {
            Some((_, flag)) => flags |= *flag,
            None => bail!("unknown flag name: {name}"),
        }
    }

    Ok(flags)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_parse() -> anyhow::Result<()> {
        assert_eq!(parse("4")?, Flags::UNMAPPED);
        assert_eq!(parse("0x904")?, Flags::from_bits_retain(0x904));
        assert_eq!(parse("010")?, Flags::MATE_UNMAPPED);
        assert_eq!(parse("0")?, Flags::empty());
        assert_eq!(
            parse("UNMAP,secondary")?,
            Flags::UNMAPPED | Flags::SECONDARY
        );
        assert!(parse("NOPE").is_err());
        assert!(parse("0xZZ").is_err());
        Ok(())
    }
}
