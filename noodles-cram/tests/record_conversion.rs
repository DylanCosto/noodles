//! Reconstructing records in bulk matches the generic, per-base conversion.

use std::{error::Error, io};

use noodles_core::Position;
use noodles_cram as cram;
use noodles_fasta as fasta;
use noodles_sam::{
    self as sam,
    alignment::{
        RecordBuf,
        io::Write as _,
        record::{
            Flags,
            cigar::{Op, op::Kind},
            data::field::Tag,
        },
        record_buf::data::field::{Value, value::Array},
    },
};

const REFERENCE_LEN: usize = 1000;

#[test]
fn test_records_match_generic_conversion() -> Result<(), Box<dyn Error>> {
    let reference = reference_sequence(REFERENCE_LEN);
    let records = records(&reference);
    let header = header();
    let src = write(&header, &reference, &records)?;

    let mut reader = cram::io::reader::Builder::default()
        .set_reference_sequence_repository(repository(&reference))
        .build_from_reader(&src[..]);
    reader.read_header()?;
    let actual: Vec<RecordBuf> = reader.records(&header).collect::<io::Result<_>>()?;

    let mut expected = Vec::new();
    let mut sam_via_view = Vec::new();
    let mut sam_via_record = Vec::new();

    let mut reader = cram::io::Reader::new(&src[..]);
    reader.read_header()?;
    let mut container = cram::io::reader::Container::default();
    let mut sequence = Vec::new();

    while reader.read_container(&mut container)? != 0 {
        let compression_header = container.compression_header()?;

        for result in container.slices() {
            let slice = result?;
            let (core_data_src, external_data_srcs) = slice.decode_blocks()?;
            let records = slice.records(
                repository(&reference),
                &header,
                &compression_header,
                &core_data_src,
                &external_data_srcs,
            )?;

            for record in &records {
                expected.push(RecordBuf::try_from_alignment_record(&header, record)?);

                let mut writer = sam::io::Writer::new(&mut sam_via_record);
                writer.write_alignment_record(&header, record)?;

                let view = record.with_sequence(&mut sequence);
                let mut writer = sam::io::Writer::new(&mut sam_via_view);
                writer.write_alignment_record(&header, &view)?;
            }
        }
    }

    assert_eq!(actual.len(), records.len());
    assert_eq!(actual, expected);
    assert_eq!(sam_via_view, sam_via_record);

    Ok(())
}

/// Records exercising every kind of read feature: substitutions (including against lowercase
/// reference bases), insertions, deletions, skips, soft clips, missing sequences, and unmapped
/// reads.
fn records(reference: &[u8]) -> Vec<RecordBuf> {
    use Kind::{Deletion, Insertion, Match, Skip, SoftClip};

    let mut records = Vec::new();

    for (i, ops) in [
        &[(Match, 50)][..],
        &[
            (SoftClip, 3),
            (Match, 20),
            (Insertion, 2),
            (Match, 10),
            (Deletion, 4),
            (Match, 15),
        ],
        &[(Match, 10), (Skip, 30), (Match, 10), (SoftClip, 5)],
    ]
    .into_iter()
    .enumerate()
    {
        for start in [1, 100, 401, 900] {
            let mut sequence = sequence_for(reference, start, ops);

            // Mismatches against the reference become substitution features.
            for j in (2..sequence.len()).step_by(17) {
                sequence[j] = match sequence[j].to_ascii_uppercase() {
                    b'A' => b'C',
                    b'C' => b'G',
                    b'G' => b'T',
                    _ => b'A',
                };
            }

            records.push(
                RecordBuf::builder()
                    .set_name(format!("r{i}-{start}"))
                    .set_flags(Flags::empty())
                    .set_reference_sequence_id(0)
                    .set_alignment_start(Position::new(start).unwrap())
                    .set_mapping_quality(sam::alignment::record::MappingQuality::new(60).unwrap())
                    .set_cigar(ops.iter().map(|&(kind, len)| Op::new(kind, len)).collect())
                    .set_quality_scores(vec![30; sequence.len()].into())
                    .set_sequence(sequence.into())
                    .set_data(
                        [
                            (Tag::from(*b"XI"), Value::from(7)),
                            (
                                Tag::from(*b"XB"),
                                Value::Array(Array::UInt8(vec![1, 2, 255])),
                            ),
                        ]
                        .into_iter()
                        .collect(),
                    )
                    .build(),
            );
        }
    }

    records.sort_by_key(|record| record.alignment_start());

    records.push(
        RecordBuf::builder()
            .set_name("unmapped")
            .set_flags(Flags::UNMAPPED)
            .set_sequence(b"ACGTN".to_vec().into())
            .set_quality_scores(vec![20; 5].into())
            .build(),
    );

    records
}

fn sequence_for(reference: &[u8], start: usize, ops: &[(Kind, usize)]) -> Vec<u8> {
    let mut sequence = Vec::new();
    let mut i = start - 1;

    for &(kind, len) in ops {
        match kind {
            Kind::Match => {
                sequence.extend(&reference[i..i + len]);
                i += len;
            }
            Kind::Deletion | Kind::Skip => i += len,
            _ => sequence.extend(b"N".repeat(len)),
        }
    }

    sequence
}

/// A pseudorandom reference with a lowercase (soft-masked) stretch.
fn reference_sequence(len: usize) -> Vec<u8> {
    let mut state = 7u32;

    (0..len)
        .map(|i| {
            state = state.wrapping_mul(1_103_515_245).wrapping_add(12_345);
            let base = b"ACGT"[(state >> 16) as usize % 4];
            if (380..520).contains(&i) {
                base.to_ascii_lowercase()
            } else {
                base
            }
        })
        .collect()
}

fn header() -> sam::Header {
    format!("@SQ\tSN:sq0\tLN:{REFERENCE_LEN}\n")
        .parse()
        .expect("valid header")
}

fn repository(reference: &[u8]) -> fasta::Repository {
    let definition = fasta::record::Definition::new("sq0", None);
    let record = fasta::Record::new(definition, reference.to_vec().into());
    fasta::Repository::new(vec![record])
}

fn write(
    header: &sam::Header,
    reference: &[u8],
    records: &[RecordBuf],
) -> Result<Vec<u8>, Box<dyn Error>> {
    let mut writer = cram::io::writer::Builder::default()
        .set_reference_sequence_repository(repository(reference))
        .build_from_writer(Vec::new());

    writer.write_header(header)?;

    for record in records {
        writer.write_alignment_record(header, record)?;
    }

    writer.try_finish(header)?;

    Ok(writer.into_inner())
}
