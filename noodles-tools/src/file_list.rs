//! Lists of input file names (samtools' `read_file_list`, for `coverage -b` and `depth -f`).

use std::{fs, path::Path};

use anyhow::{Context, bail};

/// Reads file names, one per line. Trailing whitespace is removed and empty lines are skipped;
/// every name must be an existing file or a URL.
pub fn read(list: &Path) -> anyhow::Result<Vec<String>> {
    let text = fs::read_to_string(list).with_context(|| format!("{}", list.display()))?;

    let files: Vec<String> = text
        .lines()
        .map(str::trim_end)
        .filter(|line| !line.is_empty())
        .map(String::from)
        .collect();

    for file in &files {
        if !file.contains("://") && !Path::new(file).exists() {
            bail!(
                "The file list \"{}\" appears broken, could not locate: {file}",
                list.display()
            );
        }
    }

    if files.is_empty() {
        bail!("No files read from {}", list.display());
    }

    Ok(files)
}
