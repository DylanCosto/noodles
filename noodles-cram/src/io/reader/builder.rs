use std::{
    fs::File,
    io::{self, Read},
    path::Path,
};

use noodles_fasta as fasta;

use super::{Fields, Reader};

/// A CRAM reader builder.
#[derive(Debug, Default)]
pub struct Builder {
    reference_sequence_repository: fasta::Repository,
    fields: Fields,
}

impl Builder {
    /// Sets the reference sequence repository.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_cram::io::reader::Builder;
    /// use noodles_fasta as fasta;
    ///
    /// let reference_sequence_repository = fasta::Repository::default();
    ///
    /// let builder = Builder::default()
    ///     .set_reference_sequence_repository(reference_sequence_repository);
    /// ```
    pub fn set_reference_sequence_repository(
        mut self,
        reference_sequence_repository: fasta::Repository,
    ) -> Self {
        self.reference_sequence_repository = reference_sequence_repository;
        self
    }

    /// Sets the record fields to decode.
    ///
    /// By default, every field is decoded.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_cram::{io::reader::Builder, io::reader::Fields};
    /// let builder = Builder::default().set_fields(Fields::default().set_quality_scores(false));
    /// ```
    pub fn set_fields(mut self, fields: Fields) -> Self {
        self.fields = fields;
        self
    }

    /// Builds a CRAM reader from a path.
    ///
    /// # Examples
    ///
    /// ```no_run
    /// use noodles_cram::io::reader::Builder;
    /// let reader = Builder::default().build_from_path("sample.cram")?;
    /// # Ok::<_, std::io::Error>(())
    /// ````
    pub fn build_from_path<P>(self, src: P) -> io::Result<Reader<File>>
    where
        P: AsRef<Path>,
    {
        File::open(src).map(|file| self.build_from_reader(file))
    }

    /// Builds a CRAM reader from a reader.
    ///
    /// # Examples
    ///
    /// ```
    /// # use std::io;
    /// use noodles_cram::io::reader::Builder;
    /// let reader = Builder::default().build_from_reader(io::empty());
    /// # Ok::<_, io::Error>(())
    /// ```
    pub fn build_from_reader<R>(self, reader: R) -> Reader<R>
    where
        R: Read,
    {
        Reader {
            inner: reader,
            reference_sequence_repository: self.reference_sequence_repository,
            fields: self.fields,
        }
    }
}
