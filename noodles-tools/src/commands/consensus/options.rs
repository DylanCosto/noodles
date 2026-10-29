//! Command-line options, parsed as samtools' `getopt_long` loop does.
//!
//! samtools applies options in command-line order and some of them interact (`-X` presets set
//! several parameters, `-f pileup` changes the default block size), so the arguments are parsed
//! here with GNU `getopt_long` rules: permuted operands, bundled short options, attached or
//! separate values, and unambiguous long option prefixes.

use std::{fs, path::PathBuf};

use anyhow::{Context, bail};

use super::{
    caller::{CallerOptions, Mode},
    pileup::Filters,
    tables::{QcalTable, STATIC_QCAL},
};

/// The output format.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub enum Format {
    Fasta,
    Fastq,
    Pileup,
}

/// The parsed options of `samtools consensus`.
#[derive(Debug)]
pub struct Options {
    pub caller: CallerOptions,
    pub fmt: Format,
    pub line_len: usize,
    /// `-a`, given once or twice.
    pub all_bases: u32,
    pub show_del: bool,
    pub show_ins: bool,
    pub mark_ins: bool,
    pub filters: Filters,
    pub ref_qual: i32,
    /// The block size for threads (`-Z`).
    pub span: i64,
    pub region: Option<String>,
    pub regions_file: Option<PathBuf>,
    pub threads: usize,
    pub reference: Option<PathBuf>,
    pub output: Option<PathBuf>,
    pub p_het: f64,
    pub p_indel: f64,
    pub het_scale: f64,
    pub homopoly_redux: f64,
    pub qcal: QcalTable,
    pub input: String,
}

/// What the command line asks for.
pub enum Parsed {
    Run(Box<Options>),
    /// Print the usage to stdout and succeed (no input file given).
    Usage,
}

pub const USAGE: &str = "\
Usage: ntools consensus [options] <in.bam>

Options:
  -r, --region REG      Limit query to REG. Requires an index
      --regions-file FILE
                        Output regions listed in the BED FILE
  -f, --format FMT      Output in format FASTA, FASTQ or PILEUP [FASTA]
  -l, --line-len INT    Wrap FASTA/Q at line length INT [70]
  -o, --output FILE     Output consensus to FILE
  -m, --mode STR        Switch consensus mode to \"simple\"/\"bayesian\" [bayesian]
  -a                    Output all bases (start/end of reference)
  --rf, --incl-flags STR|INT
                        Only include reads with any flag bit set [0]
  --ff, --excl-flags STR|INT
                        Exclude reads with any flag bit set
                        [UNMAP,SECONDARY,QCFAIL,DUP]
  --min-MQ INT          Exclude reads with mapping quality below INT [0]
  --min-BQ INT          Exclude reads with base quality below INT [0]
  --show-del yes/no     Whether to show deletion as \"*\" [no]
  --show-ins yes/no     Whether to show insertions [yes]
  --mark-ins            Add '+' before every inserted base/qual [off]
  -A, --ambig           Enable IUPAC ambiguity codes [off]
  -d, --min-depth INT   Minimum depth of INT [1]
  -Z, --block-size INT  Size of chromosome block (bp) when threading [100000]
      --ref-qual INT    QUAL to use for reference bases [0]

For simple consensus mode:
  -q, --(no-)use-qual   Use quality values in calculation [off]
  -c, --call-fract INT  At least INT portion of bases must agree [0.75]
  -H, --het-fract INT   Minimum fraction of 2nd-most to most common base [0.15]

For default \"Bayesian\" consensus mode:
  -C, --cutoff C        Consensus cutoff quality C [10]
      --(no-)adj-qual   Modify quality with local minima [on]
      --(no-)use-MQ     Use mapping quality in calculation [on]
      --(no-)adj-MQ     Modify mapping quality by local NM [on]
      --NM-halo INT     Size of window for NM count in --adj-MQ [50]
      --scale-MQ FLOAT  Scale mapping quality by FLOAT [1.00]
      --low-MQ  INT     Cap minimum mapping quality [1]
      --high-MQ INT     Cap maximum mapping quality [60]
      --P-het FLOAT     Probability of heterozygous site[1.0e-03]
      --P-indel FLOAT   Probability of indel sites[2.0e-04]
      --het-scale FLOAT Heterozygous SNP probability multiplier[1.0e+00]
  -p, --homopoly-fix    Spread low-qual bases to both ends of homopolymers
      --homopoly-score FLOAT
                        Qual fraction adjustment for -p option [0.5]
  -t, --qual-calibration FILE / :config (see man page)
                        Load quality calibration file

  -X, --config STR      Use pre-defined configuration set. STR from:
                        hiseq, hifi, r10.4_sup, r10.4_dup and ultima

Global options:
      --input-fmt-option OPT[=VAL]
               Specify a single input file format option in the form
               of OPTION or OPTION=VALUE
  -T, --reference FILE
               Reference sequence FASTA FILE [null]
  -@, --threads INT
               Number of additional decompression threads to use [0]
      --verbosity INT
               Set level of verbosity
";

/// Option keys, shared by the short and long forms.
const SHORT: &[(char, &str, bool)] = &[
    ('@', "threads", true),
    ('q', "use-qual", false),
    ('d', "min-depth", true),
    ('c', "call-fract", true),
    ('H', "het-fract", true),
    ('r', "region", true),
    ('5', "5", false),
    ('f', "format", true),
    ('C', "cutoff", true),
    ('a', "a", false),
    ('A', "ambig", false),
    ('l', "line-len", true),
    ('o', "output", true),
    ('m', "mode", true),
    ('p', "homopoly-fix", false),
    ('t', "qual-calibration", true),
    ('X', "config", true),
    ('T', "reference", true),
    ('Z', "block-size", true),
];

/// Long options: name, key, and whether it takes an argument. A key of `?` is rejected.
const LONG: &[(&str, &str, bool)] = &[
    ("input-fmt", "?", true),
    ("input-fmt-option", "input-fmt-option", true),
    ("output-fmt", "output-fmt", true),
    ("output-fmt-option", "?", true),
    ("reference", "reference", true),
    ("threads", "threads", true),
    ("write-index", "write-index", false),
    ("verbosity", "verbosity", true),
    ("use-qual", "use-qual", false),
    ("no-use-qual", "no-use-qual", false),
    ("adj-qual", "adj-qual", false),
    ("no-adj-qual", "no-adj-qual", false),
    ("use-MQ", "use-MQ", false),
    ("no-use-MQ", "no-use-MQ", false),
    ("adj-MQ", "adj-MQ", false),
    ("no-adj-MQ", "no-adj-MQ", false),
    ("NM-halo", "NM-halo", true),
    ("SC-cost", "SC-cost", true),
    ("scale-MQ", "scale-MQ", true),
    ("low-MQ", "low-MQ", true),
    ("high-MQ", "high-MQ", true),
    ("min-depth", "min-depth", true),
    ("call-fract", "call-fract", true),
    ("het-fract", "het-fract", true),
    ("region", "region", true),
    ("regions-file", "regions-file", true),
    ("format", "format", true),
    ("cutoff", "cutoff", true),
    ("ambig", "ambig", false),
    ("line-len", "line-len", true),
    ("default-qual", "default-qual", true),
    ("het-only", "het-only", false),
    ("show-del", "show-del", true),
    ("show-ins", "show-ins", true),
    ("mark-ins", "mark-ins", false),
    ("output", "output", true),
    ("incl-flags", "incl-flags", true),
    ("rf", "incl-flags", true),
    ("excl-flags", "excl-flags", true),
    ("ff", "excl-flags", true),
    ("min-MQ", "min-MQ", true),
    ("min-BQ", "min-BQ", true),
    ("P-het", "P-het", true),
    ("P-indel", "P-indel", true),
    ("het-scale", "het-scale", true),
    ("mode", "mode", true),
    ("homopoly-fix", "homopoly-fix", false),
    ("homopoly-score", "homopoly-score", true),
    ("homopoly-redux", "homopoly-redux", true),
    ("qual-calibration", "qual-calibration", true),
    ("config", "config", true),
    ("ref-qual", "ref-qual", true),
    ("block-size", "block-size", true),
];

/// GNU `getopt_long` over `args`: returns the options in order and the operands.
#[expect(clippy::type_complexity, reason = "a pair of plain lists")]
fn getopt(args: &[String]) -> anyhow::Result<(Vec<(&'static str, Option<String>)>, Vec<String>)> {
    let mut opts = Vec::new();
    let mut operands = Vec::new();
    let mut i = 0;

    while i < args.len() {
        let arg = &args[i];
        i += 1;

        if arg == "--" {
            operands.extend(args[i..].iter().cloned());
            break;
        }

        if let Some(long) = arg.strip_prefix("--") {
            let (name, value) = match long.split_once('=') {
                Some((name, value)) => (name, Some(value.to_owned())),
                None => (long, None),
            };

            let candidates: Vec<_> = match LONG.iter().find(|(n, _, _)| *n == name) {
                Some(exact) => vec![exact],
                None => LONG
                    .iter()
                    .filter(|(n, _, _)| n.starts_with(name))
                    .collect(),
            };

            let &(_, key, has_arg) = match candidates.as_slice() {
                [] => bail!("unrecognized option '--{name}'"),
                [first, rest @ ..]
                    if rest.iter().all(|(_, k, a)| *k == first.1 && *a == first.2) =>
                {
                    *first
                }
                _ => bail!("option '--{name}' is ambiguous"),
            };

            if key == "?" {
                bail!("unsupported option '--{name}'");
            }

            let value = if has_arg {
                match value {
                    Some(value) => Some(value),
                    None if i < args.len() => {
                        i += 1;
                        Some(args[i - 1].clone())
                    }
                    None => bail!("option '--{name}' requires an argument"),
                }
            } else if value.is_some() {
                bail!("option '--{name}' doesn't allow an argument");
            } else {
                None
            };

            opts.push((key, value));
        } else if let Some(cluster) = arg.strip_prefix('-').filter(|s| !s.is_empty()) {
            for (j, c) in cluster.char_indices() {
                let Some(&(_, key, has_arg)) = SHORT.iter().find(|(s, _, _)| *s == c) else {
                    bail!("invalid option -- '{c}'");
                };

                if has_arg {
                    let rest = &cluster[j + c.len_utf8()..];

                    let value = if !rest.is_empty() {
                        rest.to_owned()
                    } else if i < args.len() {
                        i += 1;
                        args[i - 1].clone()
                    } else {
                        bail!("option requires an argument -- '{c}'");
                    };

                    opts.push((key, Some(value)));
                    break;
                }

                opts.push((key, None));
            }
        } else {
            operands.push(arg.clone());
        }
    }

    Ok((opts, operands))
}

/// C's `atoi`: leading whitespace, an optional sign, and digits; 0 if there are none.
pub fn atoi(s: &str) -> i32 {
    let s = s.trim_start_matches(|c: char| c.is_ascii_whitespace());
    let (neg, digits) = match s.as_bytes().first() {
        Some(b'-') => (true, &s[1..]),
        Some(b'+') => (false, &s[1..]),
        _ => (false, s),
    };

    let n = digits
        .bytes()
        .take_while(u8::is_ascii_digit)
        .fold(0i64, |n, d| {
            n.saturating_mul(10).saturating_add(i64::from(d - b'0'))
        });

    (if neg { -n } else { n }) as i32
}

/// C's `atof`: the longest leading part that is a number; 0 if there is none.
pub fn atof(s: &str) -> f64 {
    let s = s.trim_start_matches(|c: char| c.is_ascii_whitespace());

    (1..=s.len())
        .rev()
        .filter(|&n| s.is_char_boundary(n))
        .find_map(|n| {
            let prefix = &s[..n];
            // Rust accepts "infinity" and "nan" like strtod, but not a trailing `e`.
            prefix.parse::<f64>().ok()
        })
        .unwrap_or(0.0)
}

/// `load_qcal`: a preset (`:hifi`, ...), `:flat`, or a calibration file of
/// `QUAL <q> <subst> <undercall> <overcall>` lines.
fn load_qcal(name: &str) -> anyhow::Result<QcalTable> {
    let preset = match name {
        ":hifi" => Some(1),
        ":hiseq" => Some(2),
        ":r10.4_sup" => Some(3),
        ":r10.4_dup" => Some(4),
        ":ultima" => Some(5),
        _ => None,
    };

    if let Some(i) = preset {
        return Ok(STATIC_QCAL[i]);
    }

    // The three maps are adjacent in C, so a quality above 100 writes into the next one.
    let mut flat = [0i32; 303];
    for (i, v) in flat.iter_mut().enumerate() {
        *v = (i % 101) as i32;
    }

    let to_table = |flat: &[i32; 303]| {
        let mut q: QcalTable = [[0; 101]; 3];
        for (k, map) in q.iter_mut().enumerate() {
            map.copy_from_slice(&flat[k * 101..(k + 1) * 101]);
        }
        q
    };

    if name == ":flat" {
        return Ok(to_table(&flat));
    }

    let text = fs::read(name)?;
    let mut max = 0i32;
    let mut last_qual = 0i32;

    let set = |flat: &mut [i32; 303], map: usize, q: i32, v: i32| {
        if let Some(slot) = usize::try_from(q)
            .ok()
            .and_then(|q| flat.get_mut(map * 101 + q))
        {
            *slot = v;
        }
    };
    let get = |flat: &[i32; 303], map: usize, q: i32| {
        usize::try_from(q)
            .ok()
            .and_then(|q| flat.get(map * 101 + q))
            .copied()
            .unwrap_or(0)
    };

    // Lines as `khgetline` reads them: no empty line after the last newline.
    let text = text.strip_suffix(b"\n").unwrap_or(&text);

    for line in text.split(|&c| c == b'\n') {
        let line = line.strip_suffix(b"\r").unwrap_or(line);

        if line.first() == Some(&b'#') {
            continue;
        }

        let Some([v, s, u, o]) = scan_qual_line(line) else {
            bail!("invalid quality calibration line");
        };

        while v > last_qual {
            for map in 0..3 {
                let prev = get(&flat, map, last_qual);
                set(&mut flat, map, last_qual + 1, prev);
            }
            last_qual += 1;
        }

        if (0..100).contains(&v) {
            set(&mut flat, 0, v, s);
            set(&mut flat, 1, v, u);
            set(&mut flat, 2, v, o);
        }

        if v < max {
            bail!("Qual calibration file is not in ascending order");
        }

        max = v;
    }

    for i in max + 1..101 {
        for map in 0..3 {
            let prev = get(&flat, map, max);
            set(&mut flat, map, i, prev);
        }
    }

    Ok(to_table(&flat))
}

/// `sscanf(line, "QUAL %d %d %d %d", ...)` succeeding with all four values.
fn scan_qual_line(line: &[u8]) -> Option<[i32; 4]> {
    let mut rest = std::str::from_utf8(line.strip_prefix(b"QUAL")?).ok()?;
    let mut values = [0; 4];

    for v in &mut values {
        let s = rest.trim_start_matches(|c: char| c.is_ascii_whitespace());
        let sign = usize::from(s.starts_with(['-', '+']));
        let digits = s[sign..].bytes().take_while(u8::is_ascii_digit).count();

        if digits == 0 {
            return None;
        }

        *v = atoi(&s[..sign + digits]);
        rest = &s[sign + digits..];
    }

    Some(values)
}

/// Parses the arguments after `consensus`.
pub fn parse(args: &[String]) -> anyhow::Result<Parsed> {
    let (list, operands) = getopt(args)?;

    let mut caller = CallerOptions {
        mode: Mode::Recall,
        use_qual: false,
        min_qual: 0,
        adj_qual: true,
        use_mqual: true,
        scale_mqual: 1.0,
        nm_adjust: true,
        nm_halo: 50,
        sc_cost: 60,
        low_mqual: 1,
        high_mqual: 60,
        min_depth: 1,
        call_fract: 0.75,
        het_fract: 0.5,
        cons_cutoff: 10,
        ambig: false,
        default_qual: 10,
        homopoly_fix: 0.0,
    };

    let mut o = Options {
        caller: caller.clone(),
        fmt: Format::Fasta,
        line_len: 70,
        all_bases: 0,
        show_del: false,
        show_ins: true,
        mark_ins: false,
        filters: Filters {
            incl_flags: 0,
            excl_flags: 0x4 | 0x100 | 0x200 | 0x400,
            min_mqual: 0,
        },
        ref_qual: 0,
        span: 500_000,
        region: None,
        regions_file: None,
        threads: 0,
        reference: None,
        output: None,
        p_het: 1e-3,
        p_indel: 2e-4,
        het_scale: 1.0,
        homopoly_redux: 0.01,
        qcal: STATIC_QCAL[0],
        input: String::new(),
    };

    for (key, value) in list {
        let v = value.as_deref().unwrap_or_default();
        let yes = || v.starts_with(['y', 'Y']);

        match key {
            "a" => o.all_bases += 1,
            "use-qual" => caller.use_qual = true,
            "no-use-qual" => caller.use_qual = false,
            "use-MQ" => caller.use_mqual = true,
            "no-use-MQ" => caller.use_mqual = false,
            "scale-MQ" => caller.scale_mqual = atof(v),
            "low-MQ" => caller.low_mqual = atoi(v),
            "high-MQ" => caller.high_mqual = atoi(v),
            "min-depth" => caller.min_depth = atoi(v),
            "call-fract" => caller.call_fract = atof(v),
            "het-fract" => caller.het_fract = atof(v),
            "region" => {
                if o.regions_file.is_some() {
                    bail!("option -r and --regions-file are incompatible");
                }
                o.region = Some(v.to_owned());
            }
            "regions-file" => {
                if o.region.is_some() {
                    bail!("option -r and --regions-file are incompatible");
                }
                o.regions_file = Some(PathBuf::from(v));
            }
            "cutoff" => caller.cons_cutoff = atoi(v),
            "ambig" => caller.ambig = true,
            "homopoly-fix" => caller.homopoly_fix = 0.5,
            "homopoly-score" => caller.homopoly_fix = atof(v),
            "homopoly-redux" => o.homopoly_redux = atof(v),
            "default-qual" => caller.default_qual = atoi(v),
            "het-only" | "5" | "output-fmt" | "input-fmt-option" | "write-index" | "verbosity" => {}
            "show-del" => o.show_del = yes(),
            "show-ins" => o.show_ins = yes(),
            "mark-ins" => o.mark_ins = true,
            "min-MQ" => o.filters.min_mqual = atoi(v),
            "min-BQ" => caller.min_qual = atoi(v),
            "P-het" => o.p_het = atof(v),
            "P-indel" => o.p_indel = atof(v),
            "het-scale" => o.het_scale = atof(v),
            "adj-qual" => caller.adj_qual = true,
            "no-adj-qual" => caller.adj_qual = false,
            "adj-MQ" => caller.nm_adjust = true,
            "no-adj-MQ" => caller.nm_adjust = false,
            "NM-halo" => caller.nm_halo = atoi(v),
            "SC-cost" => caller.sc_cost = atoi(v),
            "block-size" => o.span = i64::from(atoi(v).max(2)),
            "mode" => {
                caller.mode = match v.to_ascii_lowercase().as_str() {
                    "simple" => Mode::Simple,
                    "bayesian_m" => Mode::Mixed,
                    "bayesian_p" => Mode::Precise,
                    "bayesian_r" | "bayesian" => Mode::Recall,
                    "bayesian_116" => Mode::Bayes116,
                    _ => bail!("Unknown mode {v}"),
                };
            }
            "line-len" => {
                let n = atoi(v);
                o.line_len = if n <= 0 {
                    i32::MAX as usize
                } else {
                    n as usize
                };
            }
            "format" => {
                o.fmt = match v.to_ascii_lowercase().as_str() {
                    "fasta" => Format::Fasta,
                    "fastq" => Format::Fastq,
                    "pileup" => {
                        // Pileup uses more memory, so a smaller default block.
                        if o.span == 500_000 {
                            o.span = 100_000;
                        }
                        Format::Pileup
                    }
                    _ => bail!("Unknown format {v}"),
                };
            }
            "output" => o.output = Some(PathBuf::from(v)),
            "config" => apply_config(v, &mut caller, &mut o)?,
            "incl-flags" => {
                o.filters.incl_flags = crate::flags::parse(v)
                    .with_context(|| format!("could not parse --rf {v}"))?
                    .bits();
            }
            "excl-flags" => {
                o.filters.excl_flags = crate::flags::parse(v)
                    .with_context(|| format!("could not parse --ff {v}"))?
                    .bits();
            }
            "qual-calibration" => {
                o.qcal = load_qcal(v)
                    .with_context(|| format!("failed to load quality calibration '{v}'"))?;
            }
            "reference" => o.reference = Some(PathBuf::from(v)),
            "ref-qual" => o.ref_qual = atoi(v),
            "threads" => {
                let n: i64 = v.trim().parse().context("Invalid threads value.")?;
                o.threads = usize::try_from(n).unwrap_or(0);
            }
            _ => unreachable!("every option key is handled"),
        }
    }

    o.caller = caller;

    match operands.as_slice() {
        [] => Ok(Parsed::Usage),
        [input] => {
            o.input = input.clone();
            Ok(Parsed::Run(Box::new(o)))
        }
        _ => bail!("expected one input file"),
    }
}

/// `-X`: a configuration preset.
fn apply_config(name: &str, caller: &mut CallerOptions, o: &mut Options) -> anyhow::Result<()> {
    let name = name.to_ascii_lowercase();

    let qcal = match name.as_str() {
        "hifi" => 1,
        "hiseq" => 2,
        "r10.4_sup" => 3,
        "r10.4_dup" => 4,
        "ultima" => 5,
        _ => bail!("Unrecognised configuration name: \"{name}\""),
    };

    o.qcal = STATIC_QCAL[qcal];
    caller.mode = Mode::Recall;
    o.homopoly_redux = 0.01;

    match name.as_str() {
        "hiseq" => {}
        "ultima" => {
            caller.homopoly_fix = 0.3;
            o.het_scale = 0.37;
            caller.scale_mqual = 2.0;
            caller.low_mqual = 10;
        }
        _ => {
            caller.homopoly_fix = 0.3;
            caller.low_mqual = 5;
            caller.scale_mqual = 1.5;
            o.het_scale = 0.37;
        }
    }

    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    fn args(s: &str) -> Vec<String> {
        s.split_whitespace().map(String::from).collect()
    }

    #[test]
    fn test_getopt() -> anyhow::Result<()> {
        let (opts, operands) = getopt(&args("-aa in.bam -fpileup --show-d yes -@ 4 --rf=0x2"))?;

        assert_eq!(operands, ["in.bam"]);
        assert_eq!(
            opts,
            [
                ("a", None),
                ("a", None),
                ("format", Some("pileup".into())),
                ("show-del", Some("yes".into())),
                ("threads", Some("4".into())),
                ("incl-flags", Some("0x2".into())),
            ]
        );

        assert!(getopt(&args("--min 3")).is_err());
        assert!(getopt(&args("-h")).is_err());

        Ok(())
    }

    #[test]
    fn test_options_apply_in_order() -> anyhow::Result<()> {
        let Parsed::Run(o) = parse(&args("--low-MQ 3 -X hifi x.bam"))? else {
            bail!("expected options");
        };
        assert_eq!(o.caller.low_mqual, 5);

        let Parsed::Run(o) = parse(&args("-X hifi --low-MQ 3 -f pileup -Z 7 x.bam"))? else {
            bail!("expected options");
        };
        assert_eq!(o.caller.low_mqual, 3);
        assert_eq!(o.span, 7);

        Ok(())
    }

    #[test]
    fn test_atoi_and_atof() {
        assert_eq!(atoi(" 12x"), 12);
        assert_eq!(atoi("-3"), -3);
        assert_eq!(atoi("x"), 0);
        assert_eq!(atof("1e-3"), 1e-3);
        assert_eq!(atof("0.5x"), 0.5);
        assert_eq!(atof("2e"), 2.0);
        assert_eq!(atof("junk"), 0.0);
    }
}
