//! samtools view record filters: flags, mapping quality, subsampling, names, read groups, tags,
//! regions, query lengths, and filter expressions.

use std::{cell::RefCell, collections::HashSet, fs, io, path::Path};

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
    alignment::{RawEncoder, Rec},
    bam_raw,
    expr::Expression,
    regions::{RegionSet, record_interval},
};

thread_local! {
    /// Encodes non-BAM records for the filters that read the BAM encoding.
    static ENCODER: RefCell<RawEncoder> = RefCell::new(RawEncoder::new());
}

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
    /// A filter expression (`-e`).
    pub expression: Option<Expression>,
    /// Exclude records whose query length, with hard clips, is below this (`-m`); 0 is off.
    pub min_query_length: i64,
    /// Keep or drop records by read group ID (`-r`, `-R`).
    pub read_group_ids: Option<ReadGroupFilter>,
    /// Whether records need an `RG` tag (`-n`).
    pub read_group_presence: ReadGroupPresence,
}

/// Whether records need an `RG` tag (`-n`).
#[derive(Clone, Copy, Debug, Default, Eq, PartialEq)]
pub enum ReadGroupPresence {
    #[default]
    Optional,
    Required,
}

impl Filter {
    /// Returns whether a record passes every filter.
    pub fn matches(&self, record: &Rec<'_>, header: &sam::Header) -> io::Result<bool> {
        Ok(self.matches_expression(record, header)? && self.matches_fields(record, header)?)
    }

    /// Returns whether a record passes the filter expression (`-e`), which samtools tests
    /// before the others (and before `-B` edits the record).
    pub fn matches_expression(&self, record: &Rec<'_>, header: &sam::Header) -> io::Result<bool> {
        let Some(expression) = &self.expression else {
            return Ok(true);
        };

        with_raw(record, header, |raw| {
            let record = bam_raw::Record::new(raw)?;
            expression
                .matches(record, header)
                .map_err(|e| io::Error::other(format!("{e:#}")))
        })
    }

    /// Returns whether a record passes every filter but the expression.
    pub fn matches_fields(&self, record: &Rec<'_>, header: &sam::Header) -> io::Result<bool> {
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

        if self.min_query_length > 0 && query_length(record)? < self.min_query_length {
            return Ok(false);
        }

        let record = record.as_dyn();

        if (self.read_group_ids.is_some()
            || self.read_group_presence == ReadGroupPresence::Required)
            && !self.matches_read_group(record)?
        {
            return Ok(false);
        }

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

    /// `-r`/`-R` and `-n`. As in samtools, records without an `RG` tag pass `-r`/`-R`; only `-n`
    /// excludes them.
    fn matches_read_group(&self, record: &dyn sam::alignment::Record) -> io::Result<bool> {
        let Some(value) = record.data().get(&Tag::READ_GROUP).transpose()? else {
            return Ok(self.read_group_presence == ReadGroupPresence::Optional);
        };

        let Some(read_groups) = &self.read_group_ids else {
            return Ok(true);
        };

        let listed = match value {
            Value::String(id) => read_groups.ids.contains::<[u8]>(id.as_ref()),
            _ => false,
        };

        Ok(listed == (read_groups.mode == ListMode::Include))
    }
}

/// Calls `f` with a record's BAM encoding.
fn with_raw<T>(
    record: &Rec<'_>,
    header: &sam::Header,
    f: impl FnOnce(&[u8]) -> io::Result<T>,
) -> io::Result<T> {
    if let Rec::Bam(record) = record {
        return f(record.as_ref());
    }

    ENCODER.with_borrow_mut(|encoder| f(record.encode_raw(header, encoder)?))
}

/// The query length samtools' `-m` compares: the bases of query-consuming CIGAR operations plus
/// hard clips.
fn query_length(record: &Rec<'_>) -> io::Result<i64> {
    use sam::alignment::record::cigar::op::Kind;

    const HARD_CLIP: u32 = 5;

    if let Rec::Bam(record) = record {
        let mut ops = Vec::new();
        bam_raw::Record::new(record.as_ref())?.read_cigar(&mut ops)?;

        return Ok(ops
            .iter()
            // `bam_cigar_type` bit 1: the operation consumes the query.
            .filter(|&&op| (0x3c1a7 >> ((op & 0xf) << 1)) & 1 != 0 || op & 0xf == HARD_CLIP)
            .map(|&op| i64::from(op >> 4))
            .sum());
    }

    let mut n = 0;

    for result in record.as_dyn().cigar().iter() {
        let op = result?;

        if op.kind().consumes_read() || op.kind() == Kind::HardClip {
            n += i64::try_from(op.len()).map_err(io::Error::other)?;
        }
    }

    Ok(n)
}

/// Whether a list of IDs selects or rejects what it names.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub enum ListMode {
    Include,
    /// The list file was given as `^FILE`.
    Exclude,
}

/// `-r RG` and `-R [^]FILE`: read group IDs to keep, or (with `^`) to drop.
#[derive(Debug)]
pub struct ReadGroupFilter {
    ids: HashSet<BString>,
    mode: ListMode,
}

impl ReadGroupFilter {
    /// Builds the filter from `-r` IDs and `-R` files. As in samtools, a file's IDs are its
    /// whitespace-separated words, and keep and drop lists cannot be mixed.
    pub fn from_args(ids: &[String], files: &[String]) -> anyhow::Result<Option<Self>> {
        if ids.is_empty() && files.is_empty() {
            return Ok(None);
        }

        let excludes = files.iter().filter(|f| f.starts_with('^')).count();

        if excludes > 0 && (excludes < files.len() || !ids.is_empty()) {
            bail!("cannot mix include and exclude read-group files in the same command line");
        }

        let mut set: HashSet<BString> = ids.iter().map(|id| BString::from(id.as_str())).collect();

        for file in files {
            let path = file.strip_prefix('^').unwrap_or(file);
            set.extend(read_words(Path::new(path))?);
        }

        let mode = if excludes > 0 {
            ListMode::Exclude
        } else {
            ListMode::Include
        };

        Ok(Some(Self { ids: set, mode }))
    }

    /// Removes the `@RG` lines of read groups whose records are dropped. samtools keeps only the
    /// listed read groups, even when the list names those to drop (`-R ^FILE`); here a drop list
    /// removes the lines it names instead.
    pub fn retain_header_lines(&self, header: &mut sam::Header) {
        header
            .read_groups_mut()
            .retain(|id, _| self.ids.contains(id) == (self.mode == ListMode::Include));
    }
}

/// The whitespace-separated words of a file, as samtools reads ID lists (`fscanf("%s")`).
fn read_words(path: &Path) -> anyhow::Result<HashSet<BString>> {
    let text = fs::read(path)
        .with_context(|| format!("failed to open \"{}\" for reading", path.display()))?;

    Ok(text
        .split(u8::is_ascii_whitespace)
        .filter(|word| !word.is_empty())
        .map(BString::from)
        .collect())
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

/// `-d TAG[:VALUE]` and `-D TAG:FILE`: keep records with the tag, and optionally one of the
/// given values.
#[derive(Debug)]
pub struct TagFilter {
    tag: Tag,
    /// `None` keeps any value.
    values: Option<HashSet<BString>>,
}

impl TagFilter {
    /// Builds a filter from `-d` and `-D` arguments. As in samtools, every argument must name
    /// the same tag, and values are required as soon as one is given (a `-D` file counts even
    /// if it is empty).
    pub fn from_args(tags: &[String], tag_files: &[String]) -> anyhow::Result<Option<Self>> {
        let mut filter: Option<Self> = None;

        for arg in tags {
            let bytes = arg.as_bytes();

            if bytes.len() < 2 || (bytes.len() > 2 && bytes[2] != b':') {
                bail!("Invalid \"tag:value\" option: \"{arg}\"");
            }

            let filter = Self::for_tag(&mut filter, bytes, arg)?;

            if bytes.len() > 3 {
                filter
                    .values
                    .get_or_insert_with(HashSet::new)
                    .insert(BString::from(&bytes[3..]));
            }
        }

        for arg in tag_files {
            let bytes = arg.as_bytes();

            // samtools also accepts `;`, for MinGW path translation.
            if bytes.len() < 4 || !matches!(bytes[2], b':' | b';') {
                bail!("Invalid \"tag:file\" option: \"{arg}\"");
            }

            let words = read_words(Path::new(&arg[3..]))?;
            Self::for_tag(&mut filter, bytes, arg)?
                .values
                .get_or_insert_with(HashSet::new)
                .extend(words);
        }

        Ok(filter)
    }

    /// The filter for the tag `arg` starts with, which must be the tag of any earlier argument.
    fn for_tag<'f>(
        filter: &'f mut Option<Self>,
        bytes: &[u8],
        arg: &str,
    ) -> anyhow::Result<&'f mut Self> {
        let tag = Tag::from([bytes[0], bytes[1]]);
        let filter = filter.get_or_insert(Self { tag, values: None });

        if filter.tag != tag {
            let previous: [u8; 2] = filter.tag.into();
            bail!(
                "Different tag \"{}\" was specified before: \"{arg}\"",
                String::from_utf8_lossy(&previous)
            );
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
