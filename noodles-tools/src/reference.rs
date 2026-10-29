//! Random access to an indexed (`.fai`) FASTA file.

use std::{collections::HashMap, ffi::OsString, fs::File, path::Path};

use anyhow::Context;
use noodles_fasta::{self as fasta, fai};

use crate::commands::faidx::{build_index, read_bases};

/// Builds and saves a FASTA file's `.fai` index if it has none, as HTSlib's `fai_load` does.
pub fn ensure_index(path: &Path) -> anyhow::Result<()> {
    let mut index = OsString::from(path);
    index.push(".fai");

    if !Path::new(&index).exists() {
        build_index(path).with_context(|| format!("Could not build faidx: {}", path.display()))?;
    }

    Ok(())
}

/// An indexed FASTA file, read through its index.
pub struct IndexedFasta {
    file: File,
    records: HashMap<Vec<u8>, fai::Record>,
}

impl IndexedFasta {
    pub fn open(path: &Path) -> anyhow::Result<Self> {
        ensure_index(path)?;
        let index = fasta::fs::read_associated_index(path)
            .with_context(|| format!("Could not load faidx: {}", path.display()))?;
        let records = index
            .as_ref()
            .iter()
            .map(|record| (record.name().to_vec(), record.clone()))
            .collect();
        let file = File::open(path).with_context(|| format!("opening {}", path.display()))?;
        Ok(Self { file, records })
    }

    /// `faidx_fetch_seq64` of the 0-based closed interval `[beg, end]`, clamped to the sequence
    /// as htslib does. Returns `None` if the sequence is not in the index.
    pub fn fetch(
        &self,
        name: &[u8],
        beg: i64,
        end: i64,
        dst: &mut Vec<u8>,
    ) -> anyhow::Result<bool> {
        let Some(record) = self.records.get(name) else {
            return Ok(false);
        };

        let length = i64::try_from(record.length())?;
        let mut beg = beg;
        let mut end = end;

        if end < beg {
            beg = end;
        }

        beg = beg.clamp(0, length);
        end = if end < 0 { 0 } else { end.min(length - 1) };

        dst.clear();

        if end >= beg {
            let mut raw = Vec::new();
            read_bases(
                &self.file,
                record,
                usize::try_from(beg + 1)?,
                usize::try_from(end + 1)?,
                &mut raw,
                dst,
            )?;
        }

        Ok(true)
    }
}
