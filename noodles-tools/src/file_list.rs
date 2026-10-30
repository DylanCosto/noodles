//! Lists of input file names: samtools' `read_file_list` (`coverage -b`, `depth -f`) and
//! HTSlib's `hts_readlines` (`merge -b`, `cat -b`).

use std::{
    fs::{self, File},
    io::{BufRead, BufReader, Read},
    path::{Path, PathBuf},
};

use anyhow::{Context, bail};
use bstr::ByteSlice;
use flate2::read::MultiGzDecoder;

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

/// Reads a list of file names, one per line, skipping empty lines (`hts_readlines`). A gzipped
/// list is decompressed; a name that is not a file but starts with `:` is a comma-separated
/// list.
pub fn read_lines(src: &Path) -> anyhow::Result<Vec<PathBuf>> {
    let file = match File::open(src) {
        Ok(file) => file,
        Err(e) => {
            if let Some(list) = src.to_str().and_then(|s| s.strip_prefix(':')) {
                return Ok(list.split(',').map(PathBuf::from).collect());
            }

            return Err(e).with_context(|| format!("Invalid file list \"{}\"", src.display()));
        }
    };

    let mut reader = BufReader::new(file);
    let is_gzipped = reader.fill_buf()?.starts_with(&[0x1f, 0x8b]);

    let mut text = Vec::new();

    if is_gzipped {
        MultiGzDecoder::new(reader).read_to_end(&mut text)?;
    } else {
        reader.read_to_end(&mut text)?;
    }

    // Lines end at `\n`, with any `\r` before it removed.
    let paths = ByteSlice::lines(text.as_slice())
        .filter(|line| !line.is_empty())
        .map(|line| line.to_path().map(Path::to_path_buf))
        .collect::<Result<_, _>>()
        .with_context(|| format!("Invalid file list \"{}\"", src.display()))?;

    Ok(paths)
}
