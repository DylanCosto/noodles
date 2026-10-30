//! Piling up an indexed input one chromosome per thread, spilling finished chromosomes to temporary files until their turn to be written.

use std::{
    fs::File,
    io::{self, BufWriter, Write},
    path::Path,
    thread,
};

use noodles_sam as sam;

use super::{Args, PileupContext, check_reference_sequences, pileup_region};
use crate::alignment::{AlignmentReader, ReadOptions};

/// The output a chromosome keeps in memory before spilling the rest to a temporary file.
pub(super) const SPILL_THRESHOLD: usize = 64 << 20;

/// The inputs' headers and the reference sequence names, for [`pileup_by_chromosome`].
pub(super) type Chromosomes = (Vec<sam::Header>, Vec<Vec<u8>>);

/// Returns the headers and reference sequence names if the pileup can run one chromosome at a
/// time on several threads: every input is indexed, no region is given, and there are several
/// sequences, all with names that are valid regions.
///
/// Chromosomes are independent in the pileup (reads, overlapping mates, and the depth limit
/// never span two), and a whole-chromosome region gives exactly that chromosome's part of the
/// output, also with `-a`/`-aa`. So the output is the same as reading the inputs through.
pub(super) fn whole_chromosomes(
    args: &Args,
    inputs: &[String],
    indexes: Option<&[String]>,
    read_options: &ReadOptions,
) -> anyhow::Result<Option<Chromosomes>> {
    use crate::alignment::AlignmentIndex;

    if args.threads == 0 || args.region.is_some() {
        return Ok(None);
    }

    let mut headers = Vec::with_capacity(inputs.len());

    for (i, input) in inputs.iter().enumerate() {
        if input == "-" {
            return Ok(None);
        }

        let (mut reader, path) = AlignmentReader::open(input, read_options)?;
        let header = reader.read_header()?;

        let Some(path) = path else {
            return Ok(None);
        };

        let explicit = indexes.map(|indexes| Path::new(&indexes[i]));

        if AlignmentIndex::read(reader.format(), &path, explicit).is_err() {
            return Ok(None);
        }

        headers.push(header);
    }

    check_reference_sequences(&headers, inputs)?;

    let names: Vec<Vec<u8>> = headers[0]
        .reference_sequences()
        .keys()
        .map(|name| name.to_vec())
        .collect();

    if names.len() < 2
        || names
            .iter()
            .any(|name| name.contains(&b':') || name.is_empty())
    {
        return Ok(None);
    }

    Ok(Some((headers, names)))
}

/// Runs the pileup of each chromosome as a whole-chromosome region, `jobs` at a time, and writes
/// their output in order. Chromosomes ahead of the one being written keep their output in
/// memory, then in a temporary file.
pub(super) fn pileup_by_chromosome<W: Write>(
    context: &PileupContext<'_>,
    names: &[Vec<u8>],
    jobs: usize,
    out: &mut W,
) -> anyhow::Result<()> {
    let regions: Vec<String> = names
        .iter()
        .map(|name| String::from_utf8_lossy(name).into_owned())
        .collect();

    let jobs = jobs.max(1);

    // With fewer chromosomes than threads, each one realigns on several.
    let active = jobs.min(regions.len());
    let workers = if jobs / active > 1 { jobs / active } else { 0 };

    let read_options = ReadOptions {
        threads: 0,
        reference: context.read_options.reference.clone(),
    };

    let job_context = PileupContext {
        read_options: &read_options,
        workers,
        ..*context
    };
    let job_context = &job_context;

    thread::scope(|scope| -> anyhow::Result<()> {
        let mut pending = std::collections::VecDeque::new();
        let mut spawned = 0;

        for i in 0..regions.len() {
            while spawned < regions.len() && spawned < i + jobs {
                let region = &regions[spawned];

                pending.push_back(scope.spawn(move || -> anyhow::Result<SpillWriter> {
                    let mut writer = SpillWriter::default();
                    pileup_region(job_context, Some(region), None, &mut writer)?;
                    Ok(writer)
                }));

                spawned += 1;
            }

            let handle = pending.pop_front().expect("spawned ahead");
            let output = handle.join().expect("pileup thread panicked")?;
            output.copy_to(out)?;
        }

        Ok(())
    })
}

/// Output kept in memory up to [`SPILL_THRESHOLD`], then in a temporary file.
#[derive(Default)]
pub(super) struct SpillWriter {
    pub(super) memory: Vec<u8>,
    pub(super) file: Option<BufWriter<File>>,
}

impl SpillWriter {
    /// Writes everything to `out`, in order.
    pub(super) fn copy_to<W: Write>(self, out: &mut W) -> anyhow::Result<()> {
        use std::io::Seek as _;

        out.write_all(&self.memory)?;

        if let Some(file) = self.file {
            let mut file = file.into_inner().map_err(io::IntoInnerError::into_error)?;
            file.rewind()?;
            io::copy(&mut file, out)?;
        }

        Ok(())
    }
}

impl Write for SpillWriter {
    fn write(&mut self, buf: &[u8]) -> io::Result<usize> {
        if let Some(file) = &mut self.file {
            return file.write(buf);
        }

        if self.memory.len() + buf.len() > SPILL_THRESHOLD {
            let mut file = BufWriter::with_capacity(1 << 20, tempfile::tempfile()?);
            file.write_all(buf)?;
            self.file = Some(file);
        } else {
            self.memory.extend_from_slice(buf);
        }

        Ok(buf.len())
    }

    fn flush(&mut self) -> io::Result<()> {
        match &mut self.file {
            Some(file) => file.flush(),
            None => Ok(()),
        }
    }
}
