/// Whether decoding adds the MD and NM tags a CRAM writer left out.
///
/// HTSlib's CRAM writer drops MD and NM tags that match the alignment, and its reader rebuilds
/// them unless `CRAM_OPT_DECODE_MD` is set to 0 (as `samtools fastq`, `mpileup`, `depth`, and
/// `coverage` do, and as `--input-fmt-option decode_md=0` does).
#[derive(Clone, Copy, Debug, Default, Eq, PartialEq)]
pub enum MdNm {
    /// Add the MD and NM tags that mapped records with a reference lack (HTSlib's default).
    #[default]
    Regenerate,
    /// Keep the tags as stored.
    AsStored,
}
