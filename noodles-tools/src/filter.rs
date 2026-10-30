//! samtools view record filters: flags, mapping quality, subsampling, names, tags, and regions.

use std::{collections::HashSet, fs, io, path::Path};

use anyhow::{Context, bail};
use bstr::BString;
use noodles_sam::{
    self as sam,
    alignment::record::{
        Flags,
        data::field::{Tag, Value},
    },
    header::record::value::map::read_group::tag::LIBRARY,
};

use crate::{
    alignment::Rec,
    regions::{RegionSet, record_interval},
};

/// Mapping quality samtools uses for "unavailable"; it passes every `-q` threshold.
const MISSING_MAPPING_QUALITY: u8 = 255;

/// A samtools view record filter. Every configured condition must hold.
#[derive(Debug, Default)]
pub struct Filter {
    pub require: Flags,
    pub exclude: Flags,
    pub exclude_all: Option<Flags>,
    pub min_mapping_quality: u8,
    /// Keep only records whose `RG` is one of these read groups (`-l`).
    pub read_groups: Option<HashSet<BString>>,
    pub subsample: Option<Subsample>,
    pub names: Option<HashSet<BString>>,
    pub tag: Option<TagFilter>,
    /// Keep only records overlapping these intervals (`-L`).
    pub regions: Option<RegionSet>,
}

impl Filter {
    pub fn matches(&self, record: &Rec<'_>, header: &sam::Header) -> io::Result<bool> {
        let flags = record.flags()?;

        if !flags.contains(self.require) || flags.intersects(self.exclude) {
            return Ok(false);
        }

        if let Some(all) = self.exclude_all
            && flags.contains(all)
        {
            return Ok(false);
        }

        if self.min_mapping_quality > 0 {
            let mapq = record.mapping_quality()?.unwrap_or(MISSING_MAPPING_QUALITY);

            if mapq < self.min_mapping_quality {
                return Ok(false);
            }
        }

        let record = record.as_dyn();

        if let Some(names) = &self.names {
            let name: &[u8] = record.name().map_or(b"*", |name| name.as_ref());
            if !names.contains(name) {
                return Ok(false);
            }
        }

        if let Some(read_groups) = &self.read_groups {
            let read_group = match record.data().get(&Tag::READ_GROUP).transpose()? {
                Some(Value::String(id)) => id,
                _ => return Ok(false),
            };

            if !read_groups.contains::<[u8]>(read_group.as_ref()) {
                return Ok(false);
            }
        }

        if let Some(tag) = &self.tag
            && !tag.matches(record)?
        {
            return Ok(false);
        }

        if let Some(regions) = &self.regions {
            match record_interval(record, header)? {
                Some((id, start, end)) if regions.first_overlap(id, start, end).is_some() => {}
                _ => return Ok(false),
            }
        }

        if let Some(subsample) = &self.subsample {
            let name: &[u8] = record.name().map_or(b"*", |name| name.as_ref());
            if !subsample.keeps(name) {
                return Ok(false);
            }
        }

        Ok(true)
    }
}

/// The IDs of the read groups in `library` (their `LB`), for `samtools view -l`.
pub fn library_read_groups(header: &sam::Header, library: &str) -> HashSet<BString> {
    header
        .read_groups()
        .iter()
        .filter(|(_, read_group)| {
            read_group
                .other_fields()
                .get(&LIBRARY)
                .is_some_and(|lb| lb == library.as_bytes())
        })
        .map(|(id, _)| id.clone())
        .collect()
}

/// Reads one read name per line (the first whitespace-delimited word), as `samtools view -N`.
pub fn read_names(path: &Path) -> anyhow::Result<HashSet<BString>> {
    let text = fs::read(path).with_context(|| format!("reading {}", path.display()))?;

    Ok(text
        .split(|&b| b == b'\n')
        .filter_map(|line| line.split(|b| b.is_ascii_whitespace()).next())
        .filter(|name| !name.is_empty())
        .map(BString::from)
        .collect())
}

/// `-d TAG[:VALUE]`: keep records with the tag, and optionally one of the given values.
#[derive(Debug)]
pub struct TagFilter {
    tag: Tag,
    /// `None` keeps any value.
    values: Option<HashSet<BString>>,
}

impl TagFilter {
    /// Builds a filter from `-d` arguments. As in samtools, every argument must name the same tag.
    pub fn from_args(args: &[String]) -> anyhow::Result<Option<Self>> {
        let mut filter: Option<Self> = None;

        for arg in args {
            let (raw_tag, value) = match arg.split_once(':') {
                Some((tag, value)) => (tag, Some(value)),
                None => (arg.as_str(), None),
            };

            let tag: [u8; 2] = raw_tag
                .as_bytes()
                .try_into()
                .with_context(|| format!("invalid tag: {raw_tag}"))?;
            let tag = Tag::from(tag);

            let filter = filter.get_or_insert_with(|| Self {
                tag,
                values: Some(HashSet::new()),
            });

            if filter.tag != tag {
                bail!("only one tag can be given with -d");
            }

            match (value, &mut filter.values) {
                (Some(value), Some(values)) => {
                    values.insert(BString::from(value));
                }
                _ => filter.values = None,
            }
        }

        Ok(filter)
    }

    fn matches(&self, record: &dyn sam::alignment::Record) -> io::Result<bool> {
        let Some(value) = record.data().get(&self.tag).transpose()? else {
            return Ok(false);
        };

        let Some(values) = &self.values else {
            return Ok(true);
        };

        // As in samtools, string-like values compare as text and integers by their decimal form.
        let text: BString = match value {
            Value::Character(c) => BString::from(vec![c]),
            Value::String(s) | Value::Hex(s) => BString::from(s.to_vec()),
            Value::Int8(n) => n.to_string().into(),
            Value::UInt8(n) => n.to_string().into(),
            Value::Int16(n) => n.to_string().into(),
            Value::UInt16(n) => n.to_string().into(),
            Value::Int32(n) => n.to_string().into(),
            Value::UInt32(n) => n.to_string().into(),
            Value::Float(_) | Value::Array(_) => return Ok(false),
        };

        Ok(values.contains(&text))
    }
}

/// `-s SEED.FRACTION`: keeps a pseudorandom fraction of templates, chosen by read name so that
/// mates are kept together, matching samtools for the same seed.
#[derive(Clone, Copy, Debug)]
pub struct Subsample {
    seed: u32,
    fraction: f64,
}

impl Subsample {
    /// Parses samtools' `INT.FRAC` form: the integer part seeds, the rest is the fraction kept.
    pub fn parse(s: &str) -> anyhow::Result<Self> {
        let (raw_seed, raw_fraction) = match s.find('.') {
            Some(i) => (&s[..i], &s[i..]),
            None => (s, ""),
        };

        let seed: i64 = if raw_seed.is_empty() {
            0
        } else {
            raw_seed
                .parse()
                .with_context(|| format!("invalid subsampling seed: {s}"))?
        };

        let fraction: f64 = if raw_fraction.is_empty() {
            1.0
        } else {
            raw_fraction
                .parse()
                .with_context(|| format!("invalid subsampling fraction: {s}"))?
        };

        // samtools turns small seeds (0, 1, 2, ...) into well-mixed ones with srand/rand.
        let seed = if seed == 0 {
            0
        } else {
            #[expect(
                clippy::cast_possible_truncation,
                reason = "srand takes an unsigned int, so samtools truncates the seed the same way"
            )]
            let seed = seed as u32;
            glibc_first_rand(seed)
        };

        Ok(Self { seed, fraction })
    }

    pub fn keeps(&self, name: &[u8]) -> bool {
        let k = wang_hash(x31_hash(name) ^ self.seed);
        f64::from(k & 0x00ff_ffff) / f64::from(0x0100_0000) < self.fraction
    }
}

/// khash's `__ac_X31_hash_string`.
fn x31_hash(s: &[u8]) -> u32 {
    let Some((&first, rest)) = s.split_first() else {
        return 0;
    };

    rest.iter().fold(u32::from(first), |h, &c| {
        (h << 5).wrapping_sub(h).wrapping_add(u32::from(c))
    })
}

/// khash's `__ac_Wang_hash`.
fn wang_hash(mut key: u32) -> u32 {
    key = key.wrapping_add(!(key << 15));
    key ^= key >> 10;
    key = key.wrapping_add(key << 3);
    key ^= key >> 6;
    key = key.wrapping_add(!(key << 11));
    key ^= key >> 16;
    key
}

/// The first value glibc's `rand()` returns after `srand(seed)` (the TYPE_3 additive generator).
fn glibc_first_rand(seed: u32) -> u32 {
    const MODULUS: i64 = 2_147_483_647;

    let mut r = [0i64; 344];
    // glibc treats a zero seed as one.
    r[0] = i64::from(if seed == 0 { 1 } else { seed });

    for i in 1..31 {
        // glibc computes 16807 * r[i - 1] % (2^31 - 1) with Schrage's method on a signed int,
        // whose result equals this, taking the seed as a signed 32-bit value.
        let prev = i64::from(u32::try_from(r[i - 1]).unwrap_or(0).cast_signed());
        let mut v = (16807 * prev) % MODULUS;
        if v < 0 {
            v += MODULUS;
        }
        r[i] = v;
    }

    for i in 31..34 {
        r[i] = r[i - 31];
    }

    for i in 34..344 {
        r[i] = i64::from((r[i - 31] as u32).wrapping_add(r[i - 3] as u32));
    }

    let next = (r[344 - 31] as u32).wrapping_add(r[344 - 3] as u32);
    next >> 1
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_glibc_first_rand() {
        // Known values of `srand(seed); rand()` on glibc.
        assert_eq!(glibc_first_rand(1), 1_804_289_383);
        assert_eq!(glibc_first_rand(42), 71_876_166);
    }

    #[test]
    fn test_subsample_parse() -> anyhow::Result<()> {
        let s = Subsample::parse("0.25")?;
        assert_eq!(s.seed, 0);
        assert!((s.fraction - 0.25).abs() < f64::EPSILON);

        let s = Subsample::parse("42.5")?;
        assert_eq!(s.seed, glibc_first_rand(42));
        assert!((s.fraction - 0.5).abs() < f64::EPSILON);

        Ok(())
    }
}
