//! The statistics report (`-s`, `-f`, `--json`) and library size estimate.

use std::io::{self, Write};

use super::Stats;

/// The Picard library size estimate samtools reports.
pub(super) fn estimate_library_size(
    paired_reads: i64,
    paired_duplicate_reads: i64,
    optical: i64,
) -> i64 {
    let coverage = |x: f64, c: f64, n: f64| c / x - 1.0 + (-n / x).exp();

    let non_optical_pairs = (paired_reads - optical) / 2;
    let unique_pairs = (paired_reads - paired_duplicate_reads) / 2;
    let duplicate_pairs = (paired_duplicate_reads - optical) / 2;

    if non_optical_pairs == 0
        || duplicate_pairs == 0
        || unique_pairs == 0
        || non_optical_pairs <= duplicate_pairs
    {
        eprintln!(
            "[markdup] warning, unable to calculate estimated library size. Read pairs \
             {non_optical_pairs} should be greater than duplicate pairs {duplicate_pairs}, which \
             should both be non zero."
        );
        return 0;
    }

    let (u, n) = (unique_pairs as f64, non_optical_pairs as f64);
    let (mut m, mut big_m) = (1.0, 100.0);

    if coverage(m * u, u, n) < 0.0 {
        eprintln!("[markdup] warning, unable to calculate estimated library size.");
        return 0;
    }

    while coverage(big_m * u, u, n) > 0.0 {
        big_m *= 10.0;
    }

    for _ in 0..40 {
        let r = (m + big_m) / 2.0;
        let value = coverage(r * u, u, n);

        if value > 0.0 {
            m = r;
        } else if value < 0.0 {
            big_m = r;
        } else {
            break;
        }
    }

    #[expect(
        clippy::cast_possible_truncation,
        reason = "samtools truncates the estimate to an unsigned long the same way"
    )]
    let size = (u * (m + big_m) / 2.0) as i64;
    size
}

/// `write_stats`. The library size is estimated (possibly warning) before the title is written.
pub(super) fn write_stats<W: Write>(
    out: &mut W,
    title: &str,
    name: &str,
    stats: &Stats,
) -> io::Result<()> {
    let library_size = estimate_library_size(stats.pair, stats.duplicate, stats.optical);

    writeln!(out, "{title}{name}")?;

    for (label, value) in stat_fields(stats, library_size) {
        writeln!(out, "{label}: {value}")?;
    }

    Ok(())
}

/// The reported statistics, in samtools' order.
pub(super) fn stat_fields(stats: &Stats, library_size: i64) -> [(&'static str, i64); 15] {
    let primary_total = stats.single_duplicate + stats.duplicate;

    [
        ("READ", stats.reading),
        ("WRITTEN", stats.writing),
        ("EXCLUDED", stats.excluded),
        ("EXAMINED", stats.examined),
        ("PAIRED", stats.pair),
        ("SINGLE", stats.single),
        ("DUPLICATE PAIR", stats.duplicate),
        ("DUPLICATE SINGLE", stats.single_duplicate),
        ("DUPLICATE PAIR OPTICAL", stats.optical),
        ("DUPLICATE SINGLE OPTICAL", stats.single_optical),
        ("DUPLICATE NON PRIMARY", stats.non_primary_duplicate),
        ("DUPLICATE NON PRIMARY OPTICAL", stats.non_primary_optical),
        ("DUPLICATE PRIMARY TOTAL", primary_total),
        (
            "DUPLICATE TOTAL",
            primary_total + stats.non_primary_duplicate,
        ),
        ("ESTIMATED_LIBRARY_SIZE", library_size),
    ]
}

/// `write_json_stats`: one block of fields at `indent`, without a trailing newline.
pub(super) fn write_json_stats<W: Write>(
    out: &mut W,
    indent: &str,
    group_name: Option<&str>,
    stats: &Stats,
) -> io::Result<()> {
    let library_size = estimate_library_size(stats.pair, stats.duplicate, stats.optical);

    if let Some(name) = group_name {
        writeln!(out, "{indent}\"READ GROUP\": \"{name}\",")?;
    }

    let fields = stat_fields(stats, library_size);

    for (i, (label, value)) in fields.iter().enumerate() {
        let separator = if i + 1 < fields.len() { ",\n" } else { "" };
        write!(out, "{indent}\"{label}\": {value}{separator}")?;
    }

    Ok(())
}

/// The `--json` report. As in samtools, names are written without escaping.
pub(super) fn write_json_report<W: Write>(
    out: &mut W,
    command_line: &str,
    total: &Stats,
    use_read_groups: bool,
    group_stats: &[Stats],
    group_names: &[String],
) -> io::Result<()> {
    pub(super) const INDENT_1: &str = "    ";
    pub(super) const INDENT_2: &str = "        ";
    pub(super) const INDENT_3: &str = "            ";

    writeln!(out, "{{")?;
    writeln!(out, "{INDENT_1}\"COMMAND\": \"{command_line}\",")?;
    write_json_stats(out, INDENT_1, None, total)?;
    writeln!(out, "{}", if use_read_groups { "," } else { "" })?;

    if use_read_groups {
        writeln!(out, "{INDENT_1}\"READ GROUPS\": [")?;

        if group_stats[0].reading > 0 {
            writeln!(out, "{INDENT_2}{{")?;
            write_json_stats(out, INDENT_3, Some("ungrouped"), &group_stats[0])?;
            writeln!(out)?;
            writeln!(out, "{INDENT_2}}},")?;
        }

        for (i, (name, stats)) in group_names.iter().zip(&group_stats[1..]).enumerate() {
            writeln!(out, "{INDENT_2}{{")?;
            write_json_stats(out, INDENT_3, Some(name), stats)?;
            writeln!(out)?;
            let separator = if i + 1 < group_names.len() { "," } else { "" };
            writeln!(out, "{INDENT_2}}}{separator}")?;
        }

        writeln!(out, "{INDENT_1}]")?;
    }

    writeln!(out, "}}")
}
