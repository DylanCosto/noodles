//! CRAM decoding across a worker pool, one container per task.
//!
//! A reader thread reads containers, workers decode them and fold each record into a
//! per-container value (SAM text, raw BAM, counts), and the calling thread receives the values in
//! container order. Containers are independent, so this is equivalent to decoding serially.

use std::{
    collections::BTreeMap,
    io::{self, Read},
    thread,
};

use anyhow::Context;
use noodles_cram::{self as cram, io::reader::MdNm};
use noodles_fasta as fasta;
use noodles_sam as sam;

/// Containers queued per worker; bounds memory while keeping workers busy.
const QUEUE_DEPTH_PER_WORKER: usize = 2;

/// Decodes every container on `workers` threads.
///
/// `fold` is called on a worker for each record with that container's accumulator; `merge` is
/// called on this thread with each container's accumulator, in file order. An error from `merge`
/// stops decoding.
pub fn decode<R, T, F, M>(
    reader: &mut cram::io::Reader<R>,
    header: &sam::Header,
    repository: &fasta::Repository,
    workers: usize,
    fold: F,
    mut merge: M,
) -> anyhow::Result<()>
where
    R: Read + Send,
    T: Default + Send,
    F: Fn(&mut T, &cram::Record<'_>) -> io::Result<()> + Sync,
    M: FnMut(T) -> anyhow::Result<()>,
{
    let workers = workers.max(1);
    let md_nm = reader.md_nm();

    thread::scope(|scope| {
        let (work_tx, work_rx) = crossbeam_channel::bounded::<(u64, cram::io::reader::Container)>(
            workers * QUEUE_DEPTH_PER_WORKER,
        );
        let (done_tx, done_rx) =
            crossbeam_channel::bounded::<(u64, io::Result<T>)>(workers * QUEUE_DEPTH_PER_WORKER);

        let read = scope.spawn(move || -> io::Result<()> {
            for id in 0.. {
                let mut container = cram::io::reader::Container::default();

                if reader.read_container(&mut container)? == 0 {
                    break;
                }

                if work_tx.send((id, container)).is_err() {
                    // The consumer stopped; it reports why.
                    break;
                }
            }

            Ok(())
        });

        for _ in 0..workers {
            let work_rx = work_rx.clone();
            let done_tx = done_tx.clone();
            let fold = &fold;

            scope.spawn(move || {
                for (id, container) in work_rx {
                    let result = decode_container(&container, header, repository, md_nm, fold);

                    if done_tx.send((id, result)).is_err() {
                        break;
                    }
                }
            });
        }

        // Only workers hold these now, so the channels close when the workers finish.
        drop(work_rx);
        drop(done_tx);

        let mut pending = BTreeMap::new();
        let mut next = 0;

        for (id, result) in &done_rx {
            pending.insert(id, result);

            while let Some(result) = pending.remove(&next) {
                merge(result.with_context(|| format!("decoding CRAM container {next}"))?)?;
                next += 1;
            }
        }

        read.join()
            .map_err(|_| anyhow::anyhow!("CRAM reader thread panicked"))?
            .context("reading CRAM containers")
    })
}

fn decode_container<T, F>(
    container: &cram::io::reader::Container,
    header: &sam::Header,
    repository: &fasta::Repository,
    md_nm: MdNm,
    fold: &F,
) -> io::Result<T>
where
    T: Default,
    F: Fn(&mut T, &cram::Record<'_>) -> io::Result<()>,
{
    let compression_header = container.compression_header()?;
    let mut acc = T::default();

    for result in container.slices() {
        let slice = result?;
        let (core_data_src, external_data_srcs) = slice.decode_blocks()?;

        let records = slice.records_with_md_nm(
            repository.clone(),
            header,
            &compression_header,
            &core_data_src,
            &external_data_srcs,
            md_nm,
        )?;

        for record in &records {
            fold(&mut acc, record)?;
        }
    }

    Ok(acc)
}
