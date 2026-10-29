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
| `view` | `view` | SAM/BAM/CRAM conversion; `-c`, `-f/-F/-G`, `-q`, `-L`, `-M`, `-s`, `-N`, `-d`, `-x`, `-U`, `-X`, `-h/-H`, `-b/-C/-u/-1/-l`, `--write-index`, regions and `*` |
| `sort` | `sort` | coordinate, `-n` (natural), `-N` (bytewise); `-m`, `-T`, `-@`; stable, like samtools |
| `merge` | `merge -c -p` | BAM inputs; reference dictionaries are merged and IDs translated |
| `index` | `index` | BAI (default) or CSI (`-c`) for BAM; CRAI for CRAM |
| `faidx` | `faidx` | build `.fai`; extract regions (`-r`, `-n`, `-i`, `-o`) |
| `flagstat` | `flagstat` | default, `-O tsv`, `-O json` |
| `idxstats` | `idxstats` | from BAM index metadata; scans CRAM/SAM |
| `depth` | `depth` | `-a/-aa`, `-r`, `-b`, `-J`, `-q/-Q`, `-g/-G`, `-l`, `-H`, several inputs |
| `coverage` | `coverage` | per reference or `-r` region |
| `fastq`, `fasta` | `fastq`, `fasta` | `-0/-1/-2/-s/-o`, `-n/-N`, `-t/-T`, `-F/-f/-G`, `.gz` output |
| `quickcheck` | `quickcheck` | `-v`, `-q`, `-u` |
| `reheader` | `reheader` | BAM; copies compressed blocks |
| `cat` | `cat` | BAM; copies compressed blocks |

`-@ N` adds worker threads throughout. Beyond BGZF compression and decompression, which is what
samtools uses its threads for, ntools also formats SAM text and decodes CRAM containers in
parallel.

## Compatibility

`scripts/compare-samtools.sh` checks ntools against samtools byte for byte on the benchmark
corpus (short paired-end, tagged, and nanopore long reads; BAM, SAM, and CRAM 3.0/3.1). All 70
checks pass.

Intentional differences:

- `@PG` lines name ntools. Use `--no-PG` (or `-P` for `reheader`) for identical headers.
- `view` on a CRAM `*` region returns every unplaced unmapped read; samtools 1.24 returns only
  some of them from the corpus CRAM. ntools matches samtools' output for the equivalent BAM.
- `merge` combines `@RG`/`@PG` lines with the same ID, like `samtools merge -c -p`, and keeps
  each record's tag order. Plain `samtools merge` renames clashing IDs with random suffixes and
  moves `RG` to the last tag.
- BAI/CSI bins are written in insertion order, not HTSlib's hash order. samtools reads them and
  returns identical query results.
- CRAM written by ntools keeps `MD`/`NM` tags as stored; HTSlib drops and regenerates them.
  Decoded records are the same apart from tag order.
- `reheader` and `cat` recompress the block holding the end of the header with noodles' encoder,
  so files can differ in size by a few bytes; records are identical.

## Performance

`scripts/benchmark.sh` times both tools on the corpus (median of 5 warm runs).

BENCHMARK_TABLE

## Library changes in this fork

ntools relies on these changes to the noodles crates, which downstream tools get too:

- CRAM: correctness fixes, table-driven rANS decoding, direct tok3 name decoding, and
  reference-range retrieval (upstream PRs #440–#443).
- CRAM: unmapped records decode with a mapping quality of 0, as in HTSlib.
- CRAM: sequences are reconstructed in bulk (`Record::with_sequence`), and tag blocks are found
  without hashing.
- BAM: `Record` exposes its encoded bytes (`AsRef<[u8]>`, `TryFrom<Vec<u8>>`) for copying
  records without re-encoding.
- SAM: integer arrays are formatted through a buffer, which speeds up long-read tags (`MM`/`ML`).
