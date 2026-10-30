# ntools

A samtools-compatible command-line toolkit built on this fork of [noodles]. It aims to produce
the same bytes as samtools 1.24 while being faster, especially with threads.

[noodles]: https://github.com/zaeleus/noodles

## Install

```sh
cargo install --path noodles-tools
```

This builds `ntools` with libdeflate (the `libdeflate` feature, on by default).

## Commands

| Command | samtools equivalent | Notes |
|---|---|---|
| `view` | `view` | SAM/BAM/CRAM conversion; `-c`, `-f/-F/-G`, `-q`, `-m`, `-e` (filter expressions), `-L`, `-M`, `-P`, `-s`, `-N`, `-r/-R/-n`, `-d/-D`, `-x`, `-U`, `-p`, `-B`, `-t`, `-X`, `-l` (library), `-h/-H`, `-b/-C/-u/-1`, `--output-fmt-option level=N`, `--write-index`, `-?`, `-S`, regions and `*` |
| `sort` | `sort` | coordinate, `-n` (natural), `-N` (bytewise); `-m`, `-T`, `-@`; stable, like samtools |
| `merge` | `merge -c -p` | BAM inputs; reference dictionaries are merged and IDs translated |
| `index` | `index` | BAI (default) or CSI (`-c`) for BAM; CRAI for CRAM |
| `faidx`, `fqidx` | `faidx`, `fqidx` | FASTA and FASTQ, plain or BGZF (`.fai` and `.gzi` as HTSlib builds them); regions as HTSlib parses them (`{name}` quoting, `1,000`, `2k`), `-r`, `-o` (BGZF for `.gz`), `-n`, `-c`, `-i`, `--mark-strand`, `--fai-idx/--gzi-idx`, `--write-index`, `-f` |
| `flagstat` | `flagstat` | default, `-O tsv`, `-O json` |
| `idxstats` | `idxstats` | from BAM index metadata; scans CRAM/SAM |
| `depth` | `depth` | `-a/-aa`, `-r`, `-b`, `-J`, `-s`, `-q/-Q`, `-g/-G`, `-l`, `-H`, several inputs |
| `coverage` | `coverage` | per reference or `-r` region; `-q/-Q`, `-l`, flag filters (no `-m` histogram) |
| `fastq`, `fasta` | `fastq`, `fasta` | `-0/-1/-2/-s/-o`, `-n/-N`, `-t/-T` (`*` for all tags), `-F/-f/-G/--rf`, `-d/-D` (tag filters), `-O` (OQ qualities), `-i` (CASAVA), `--i1/--i2`, `--index-format`, `--barcode-tag/--quality-tag`, `-v`, `-U/--UMI-tag`, `--no-sc/--sc-aux/--no-sc-bkp`, `-c`, `.gz` output |
| `quickcheck` | `quickcheck` | `-v`, `-q`, `-u` |
| `reheader` | `reheader` | BAM; copies compressed blocks |
| `cat` | `cat` | BAM; copies compressed blocks |
| `split` | `split` | by read group or `-d` tag; `-f`, `-p`, `-u`, `-M` |
| `collate` | `collate` | groups templates by name; `-f`, `-O/-o`, `-u`, `-l`, `-m`, `-T` |
| `fixmate` | `fixmate` | `-m`, `-c`, `-r`, `-p`, `-z`, `-M` (base modification tags); parallel with `-@` |
| `dict` | `dict` | `-a`, `-s`, `-u`, `-H`, `-o`, `-A` (`AN` aliases), `-l` (`AH` for alternate loci) |
| `bedcov` | `bedcov` | `-c`, `-d`, `-j`, `-Q`, `-g/-G`, `-H`; several inputs |
| `addreplacerg` | `addreplacerg` | `-r` (repeatable, `\t` escapes), `-R`, `-m overwrite_all\|orphan_only`, `-w` |
| `markdup` | `markdup` | `-r`, `-l`, `-m t\|s`, `-c`, `-S`, `-t`, `-d` (optical duplicates), `-T`, `-s/-f`, `--json`, `--no-multi-dup`, `--include-fails`, `--use-read-groups`, `--duplicate-count`, `--read-coords`/`--coords-order`, `--barcode-tag`/`-name`/`-rgx`, `--move-umi-to-tag`, `--umi-separator` |
| `stats` | `stats` | all sections; `-r`, `-t`, regions, `-f/-F/-d`, `-c`, `-i`, `-m`, `-q`, `-l`, `-I`, `-S/-P`, `-p`, `-g`, `-x`, `--GC-depth`, `--ref-stats`; per-cycle counts use `-@` threads |
| `mpileup` | `mpileup` | text pileup with BAQ (`-B`, `-E`), overlap removal (`-x`), `-f`, `-r`, `-l`, `-a/-aa`, `-q/-Q`, `-d`, `-A`, `-C`, `-G`, `-R`, `-6`, `--rf/--ff`, `-b`, `-X`, `-s`, `-O`, `--output-BP-5`, `--output-QNAME`, `--output-extra`, `--output-sep/-empty`, `--reverse-del`, `--no-output-ins/-del/-ends`, `-M` (base modifications from MM/ML, also in insertions unless `--no-output-ins-mods`); reads and BAQ use `-@` threads, and indexed inputs with several sequences pile up one chromosome per thread |
| `calmd` | `calmd` | `-e`, `-r` (BAQ), `-A`, `-E`, `-n`, `-C`, `-d`, `-q`, `-N`, `-Q`, `-b/-u`; records are processed on `-@` threads |
| `ampliconclip` | `ampliconclip` | soft or `--hard-clip`, `--both-ends`, `--strand`, `--tolerance`, `--original`, `--keep-tag`, `--clipped`, `--fail`, `--filter-len/--fail-len/--unmap-len`, `--no-excluded`, `--rejects-file`, `--primer-counts`, `-f` |
| `import` | `import` | FASTQ/FASTA (plain or gzipped) to unaligned SAM/BAM/CRAM; `-0/-1/-2/-s`, `--i1/--i2` (index reads to `BC`/`QT`), `-b`, `-i` (CASAVA), `-U` (UMIs), `-N` (SRA names), `-T` (tags from comments, such as `MM`/`ML`), `-r/-R`, `--order`, `-O`, `-u`, `-o` |
| `reset` | `reset` | `-x`/`--keep-tag` (including `^` lists), `--reject-PG`, `--no-RG`, `--dupflag`, `-O`, `-o` |
| `checksum` | `checksum` | `-F/-f/-b`, `-c`, `-t` (including `*` negation), `-O`, `-P/-C/-M`, `-z`, `-N`, `-q`, `-v`, `-a`, `-T`, `-B`, `-m` (merge reports), `-o` |
| `head` | `head` | `-h` (header lines), `-n` (records) |
| `flags` | `flags` | numbers and flag names, both ways |
| `samples` | `samples` | `-h`, `-i`, `-T`, `-f/-F` (matching references), `-X`, `-o`, paths from stdin |
| `consensus` | `consensus` | FASTA/FASTQ/pileup; `simple` and Bayesian modes (`bayesian`, `_r`, `_p`, `_m`, `_116`), `-X` presets and `-t` calibration files, `-a/-aa`, `--show-ins/--show-del/--mark-ins`, `-A`, `-r`, `--regions-file`, `-T`, every tuning option; `-@` on an indexed input computes `-Z` blocks in parallel, with the same output as without threads |
| `cram-size` | `cram-size` | sizes by content ID and expanded compression method (gzip level, rANS order and transforms, ...), with the data series using each block; `-v`, `-e` (encodings), `-o` |
| `reference` | `reference` | from `MD` tags (any input) or embedded CRAM references (`-e`); `-r` (indexed), `-q`, `-o` |
| `rmdup` | `rmdup` | paired-end, `-s`, `-S`; libraries from `@RG` `LB`; single-end mode holds only the reads that can still be replaced |

CRAM input without `--reference` uses the FASTA file named by the `@SQ` `UR` fields when it exists
locally with a `.fai` index, as HTSlib does.

`-@ N` adds worker threads throughout. Beyond BGZF compression and decompression, which is what
samtools mostly uses its threads for, ntools also formats SAM text, decodes CRAM containers,
sorts, and fixes mates in parallel.

A typical duplicate-marking pipeline is unchanged apart from the program name:

```sh
ntools collate -Ou -@ 8 in.bam | ntools fixmate -m -@ 8 - - | ntools sort -@ 8 - \
  | ntools markdup -@ 8 - out.bam
```

On the short-read corpus this pipeline takes 1.30 s, against 3.24 s for the same samtools
pipeline, and marks the same 26,749 duplicates in identical records.

## Compatibility

`scripts/compare-samtools.sh` checks ntools against samtools byte for byte on the benchmark
corpus (short paired-end, tagged, and nanopore long reads; BAM, SAM, and CRAM 3.0/3.1). All
checks pass.

Intentional differences:

- `@PG` lines name ntools. Use `--no-PG` (or `-P` for `reheader`) for identical headers.
- `view` on a CRAM `*` region returns every unplaced unmapped read; samtools 1.24 returns only
  some of them from the corpus CRAM. ntools matches samtools' output for the equivalent BAM.
- `view -e` evaluates HTSlib's filter expression language with HTSlib's semantics (undefined
  values, tag presence as truth, no short-circuiting), with these fixes: an undefined value
  `%` a number stays undefined (HTSlib turns `[XX] % 2` into 0 when `XX` is absent, so
  `[XX] % 2 == 0` selects every record without the tag), and `%` by a value between -1 and 1 is
  undefined (HTSlib crashes). A syntax error is reported before any record is read; samtools
  reports it at the first record, and not at all for an empty input. A regular expression
  computed from record fields is compiled for each record; HTSlib reuses the first pattern it
  saw.
- `view -e` regular expressions are POSIX extended expressions run by the Rust `regex` crate on
  bytes (matching for `=~`/`!~` only needs a yes or no, so leftmost-longest rules do not
  matter). Back-references (`\1`) are not supported. As in glibc, `\w`, `\W`, `\s`, `\S`,
  `\b`, `\B`, `\<`, `\>` keep their GNU meaning, a backslash is literal inside brackets, and a
  backslash before another letter or digit is that character (`\d` is `d`, not a digit).
- `view -R ^FILE` removes the `@RG` lines of the excluded read groups. samtools instead keeps
  only those lines and removes the lines of the read groups whose records it outputs.
- `view -P` applies `-p` and `-U` to the records it fetches (the reads in the regions and their
  mates); samtools ignores both options with `-P`.
- `view -B` on a record without base qualities leaves them all missing. HTSlib marks only the
  first one missing and keeps merged values in the others (SAM output shows `*` either way).
  noodles cannot parse `B` operations in SAM text, so `-B` needs BAM input.
- `view -T ref.fa` on SAM without `@SQ` lines uses `ref.fa.fai` only if it exists; HTSlib
  builds it.
- `merge` combines `@RG`/`@PG` lines with the same ID, like `samtools merge -c -p`, and keeps
  each record's tag order. Plain `samtools merge` renames clashing IDs with random suffixes and
  moves `RG` to the last tag.
- BAI/CSI bins are written in insertion order, not HTSlib's hash order. samtools reads them and
  returns identical query results.
- CRAM written by ntools keeps `MD`/`NM` tags as stored; HTSlib drops and regenerates them.
  Decoded records are the same apart from tag order. Reading CRAM, ntools regenerates `MD`/`NM`
  as HTSlib does: for mapped records with a reference, missing tags are added after the stored
  ones (a `cF` tag can turn this off, and is removed).
- Reading CRAM records that overhang the end of the reference, ntools fills the bases past the
  end with `N`, as HTSlib does when it can decode them; HTSlib fails on such slices when it
  must regenerate `MD`.
- `reheader` and `cat` recompress the block holding the end of the header with noodles' encoder,
  so files can differ in size by a few bytes; records are identical.
- `collate` groups each template's records together, in input order, as samtools does, but
  templates come out in a different (hash) order.
- Header lines are written grouped by type (`@HD`, `@SQ`, `@RG`, `@PG`, `@CO`). HTSlib keeps
  the input's line order and appends new lines at the end, so `addreplacerg` places a new `@RG`
  line before existing `@PG` lines where samtools places it after them.
- `ampliconclip` fixes samtools bugs: `OA` tags follow the SAM specification
  (`RNAME,POS,strand,CIGAR,MAPQ,NM;`), where samtools writes the read name as RNAME and drops
  the comma before MAPQ without an `NM` tag; clipping the 3' end of a read that ends in a
  deletion keeps its other CIGAR operations, which samtools drops (leaving, say, `1D` for 18
  bases, so the read is then unmapped); and clipped records get a recalculated BAM bin.
- `import` parses `-T` tags that follow a CASAVA barcode (`-i`), which samtools loses. It reads
  an empty input file as having no records; samtools fails on it ("truncated
  file"), as it does on the empty `-0` file `samtools fastq` often writes. An `-r` read group
  line without an ID, or with a malformed field, is an error; samtools writes nothing and exits
  successfully.
- `markdup`: with `-d`, duplicates whose names have no coordinates are not compared with each
  other; samtools compares them using uninitialized memory, so a chain mixing parseable and
  unparseable names can come out differently. `--read-coords` and `--barcode-rgx` take Rust
  `regex` syntax rather than POSIX extended regular expressions (the same for ordinary
  patterns). Without `-S`, `--move-umi-to-tag` moves the UMI of every read; samtools 1.24 moves
  only those of reads still buffered at the end of the input.
- `consensus` fixes samtools bugs instead of reproducing them. Pileup output with `-a` lists each
  column once; samtools repeats growing runs of depth-0 lines after positions whose consensus is
  a deletion (and prints such positions as depth 0). `-@` gives the serial output, and on an
  unindexed input writes pileup lines, which samtools drops. A read's qualities never come from
  other reads (samtools takes a deleted or padded position's quality from stale bytes of earlier
  records), and qualities above 100 use the tables' quality-100 values instead of reading past
  them. With `-@` and a `-T` FASTA that lacks a needed sequence, ntools stops with an error where
  samtools writes truncated or garbage sequence.
- `bedcov -d` adds its column without changing the depth sum, which counts deletions and
  reference skips unless `-j`, as documented; samtools also drops them from the sum when `-d`
  is given. The `-d` column counts aligned bases, as in samtools.
- `mpileup` stops with an error when an input's reference sequences differ from the first
  input's (positions are compared by reference index). samtools mixes up such inputs, or fails
  partway when one has more sequences.
- `reference` advances the reference position over `N` (reference skip) CIGAR operations;
  samtools does not, so it writes the bases after a splice at the wrong positions (76 wrong
  bases on the tagged corpus, against none for ntools). A region with no data is labelled with
  its range (`>chr:1-100`), as one with data is; samtools labels it as the whole sequence. An
  `MD` tag that runs past the CIGAR is reported as an error; samtools fails without a message.
- `fastq`/`fasta` group records by name. samtools compares each record with a name it can
  overwrite: after a copy of a read without qualities is replaced by one with them, the next
  template is merged into the group, losing reads or pairing reads of different templates.
  Index reads (`--i1/--i2`) end with the barcode tag: samtools reads on past its end (as with
  `i*i*` and a barcode without a separator) and writes the following tag bytes as an index
  read. Barcodes that are not strings give no index reads; samtools reads their binary values
  as text. As in samtools, reads without qualities are written with `B` qualities, and `-v`
  applies only to index reads.
- `fixmate -M`: a hard-clipped secondary or supplementary alignment whose `ML` tag is not a
  `B:C` array loses its base modification tags, as invalid ones are; samtools drops `ML` and
  every tag after it, and writes corrupt records when `ML` has fewer values than `MM` calls.

## Performance

`scripts/benchmark.sh` times both tools on the corpus (median of 5 warm runs).

Measured on an Intel Xeon Gold 5418Y VM (24 vCPUs) against samtools/HTSlib 1.24, both with
libdeflate. "short" is 1M synthetic 150 bp paired reads, "tagged" 250k reads with extra tags, and
"long" 9,784 nanopore reads with MM/ML modification tags. `-@ N` adds N threads to both tools.

| Operation | samtools (s) | ntools (s) | Speedup |
|---|---:|---:|---:|
| view BAM->SAM, short | 0.381 | 0.351 | 1.09x |
| view BAM->SAM, long | 2.135 | 1.706 | 1.25x |
| view BAM->SAM, long, -@4 | 0.568 | 0.350 | 1.62x |
| view BAM->SAM, long, -@8 | 0.300 | 0.242 | 1.24x |
| view CRAM->SAM, short | 0.886 | 0.969 | 0.91x |
| view CRAM->SAM, long | 2.898 | 2.482 | 1.17x |
| view CRAM->SAM, long, -@8 | 0.577 | 0.343 | 1.68x |
| view CRAM->BAM, long, -@8 | 2.199 | 1.445 | 1.52x |
| view BAM->BAM, long, -@4 | 3.448 | 2.353 | 1.47x |
| view -c -q 30, short | 0.144 | 0.149 | 0.97x |
| flagstat, long, -@4 | 0.330 | 0.340 | 0.97x |
| index, long, -@4 | 0.333 | 0.341 | 0.98x |
| sort, short, -@4 | 0.997 | 0.596 | 1.67x |
| sort -n, short, -@4 | 1.088 | 0.580 | 1.88x |
| collate -Ou, short, -@4 | 2.106 | 0.386 | 5.46x |
| fixmate -m, short | 3.427 | 3.349 | 1.02x |
| fixmate -m, short, -@4 | 0.882 | 0.450 | 1.96x |
| markdup, short | 3.817 | 2.458 | 1.55x |
| markdup, short, -@4 | 0.909 | 0.756 | 1.20x |
| fastq, short | 0.608 | 0.292 | 2.08x |
| depth, short | 0.282 | 0.254 | 1.11x |
| stats, short | 1.659 | 1.142 | 1.45x |
| stats, short, -@4 | 1.581 | 0.727 | 2.17x |
| stats, long, -@8 | 5.699 | 1.577 | 3.61x |
| mpileup, short | 2.899 | 2.529 | 1.15x |
| mpileup -f (BAQ), short, -@8 | 27.788 | 4.365 | 6.37x |
| calmd -b, short, -@4 | 0.687 | 0.643 | 1.07x |
| calmd -r -A -b (BAQ), short, -@8 | 26.262 | 4.617 | 5.69x |
| mpileup -M (base modifications), long | 16.429 | 10.986 | 1.50x |
| ampliconclip, short | 3.048 | 1.792 | 1.70x |
| ampliconclip --both-ends, short, -@4 | 0.809 | 0.395 | 2.05x |
| coverage, short | 1.145 | 0.413 | 2.77x |
| consensus, short | 6.383 | 5.562 | 1.15x |
| consensus, short, -@4 | 1.633 | 1.429 | 1.14x |
| consensus, long | 11.325 | 10.012 | 1.13x |
| consensus, long, -@4 | 3.467 | 3.022 | 1.15x |
| faidx 1024 regions | 0.018 | 0.020 | 0.90x |

Where ntools is slower, it is single-threaded (single-threaded `mpileup` without BAQ spends its
time in the pileup engine itself, as samtools does). Short-read CRAM decoding still trails HTSlib in
per-record tag decoding and allocation; its rANS decoding (AVX2) now keeps pace with HTSlib's.
Rows near 1.0x (flagstat, index, counting) are bound by BGZF decompression, which both tools do
with libdeflate. With threads, ntools parallelizes more of the work than samtools does, so the
gap reverses.

## Library changes in this fork

ntools relies on these changes to the noodles crates, which downstream tools get too:

- CRAM: correctness fixes, table-driven rANS decoding, direct tok3 name decoding, and
  reference-range retrieval (upstream PRs #440–#443).
- CRAM: unmapped records decode with a mapping quality of 0, as in HTSlib.
- CRAM: sequences are reconstructed in bulk (`Record::with_sequence`), and tag blocks are found
  without hashing.
- CRAM: `Record::encode_bam_record` encodes a decoded record as raw BAM directly, so CRAM input
  takes the same raw BAM paths as BAM input.
- CRAM: 32-state order-1 rANS Nx16 streams (quality scores) decode with AVX2 when the CPU has
  it, falling back to the scalar decoder otherwise.
- BAM: `Record` exposes its encoded bytes (`AsRef<[u8]>`, `TryFrom<Vec<u8>>`,
  `From<Record> for Vec<u8>`) for copying records without re-encoding.
- SAM: integer arrays are formatted through a buffer, which speeds up long-read tags (`MM`/`ML`).
- CRAM: decoding regenerates `MD`/`NM` tags that HTSlib's writer drops (as HTSlib's reader
  does), and fills bases past the end of the reference with `N` instead of panicking.
- CRAM: a missing reference sequence is an error rather than a panic, and
  `Reader::set_reference_sequence_repository` lets a reference be chosen after reading the
  header.
- CRAM: `io::reader::read_block`, `Block`, and `Container::src` expose a container's raw blocks
  (method as stored, content type and ID, sizes, compressed bytes) for block-level tools such
  as `cram-size`.
- BAM: a record whose CIGAR is held in a `CG` field (over 65,535 operations) no longer lists
  that field in its data, as in HTSlib; `Record::cigar` already returned its operations.
- New crate `noodles-pileup`: HTSlib-compatible pileup (`bam_plp`), BAQ, mapping quality capping,
  MD/NM calculation, and MM/ML base modifications over encoded BAM records, which `mpileup` and
  `calmd` use. Positions inside a CIGAR operation are resolved from a cached
  cursor without reading the CIGAR, the per-position scan reads a compact array rather than the
  buffered records, and overlapping mates' qualities are combined in lockstep along shared match
  operations.
