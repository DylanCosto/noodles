//! Decoding fields that CRAM does not store for unmapped records.

use std::io;

use noodles_cram as cram;
use noodles_sam::{
    self as sam,
    alignment::{
        RecordBuf,
        io::Write as _,
        record::{Flags, MappingQuality},
        record_buf::{QualityScores, Sequence},
    },
};

// CRAM stores mapping qualities only for mapped records. As in HTSlib, unmapped records decode
// with a mapping quality of 0, so CRAM-to-SAM/BAM output matches samtools.
#[test]
fn test_unmapped_record_mapping_quality_is_zero() -> io::Result<()> {
    let header = sam::Header::default();

    let record = RecordBuf::builder()
        .set_name("r0")
        .set_flags(Flags::UNMAPPED)
        .set_sequence(Sequence::from(b"ACGT".to_vec()))
        .set_quality_scores(QualityScores::from(vec![30, 30, 30, 30]))
        .build();

    let mut writer = cram::io::Writer::new(Vec::new());
    writer.write_header(&header)?;
    writer.write_alignment_record(&header, &record)?;
    writer.finish(&header)?;
    let src = writer.get_ref().clone();

    let mut reader = cram::io::Reader::new(&src[..]);
    reader.read_header()?;
    let records: Vec<_> = reader.records(&header).collect::<io::Result<_>>()?;

    assert_eq!(records.len(), 1);
    assert_eq!(records[0].mapping_quality(), MappingQuality::new(0));

    Ok(())
}
