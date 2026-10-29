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
| `split` | `split` | by read group or `-d` tag; `-f`, `-p`, `-u`, `-M` |
| `collate` | `collate` | groups templates by name; `-f`, `-O/-o`, `-u`, `-l`, `-m`, `-T` |
| `fixmate` | `fixmate` | `-m`, `-c`, `-r`, `-p`, `-z`; parallel with `-@` |
| `dict` | `dict` | `-a`, `-s`, `-u`, `-H`, `-o` |
| `bedcov` | `bedcov` | `-c`, `-d`, `-j`, `-Q`, `-g/-G`, `-H`; several inputs |
| `addreplacerg` | `addreplacerg` | `-r` (repeatable, `\t` escapes), `-R`, `-m overwrite_all\|orphan_only`, `-w` |
| `markdup` | `markdup` | `-r`, `-l`, `-m t\|s`, `-c`, `-s/-f`, `--include-fails`, `--use-read-groups`, `--duplicate-count` |

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
- `merge` combines `@RG`/`@PG` lines with the same ID, like `samtools merge -c -p`, and keeps
  each record's tag order. Plain `samtools merge` renames clashing IDs with random suffixes and
  moves `RG` to the last tag.
- BAI/CSI bins are written in insertion order, not HTSlib's hash order. samtools reads them and
  returns identical query results.
- CRAM written by ntools keeps `MD`/`NM` tags as stored; HTSlib drops and regenerates them.
  Decoded records are the same apart from tag order.
- `reheader` and `cat` recompress the block holding the end of the header with noodles' encoder,
  so files can differ in size by a few bytes; records are identical.
- `collate` groups each template's records together, in input order, as samtools does, but
  templates come out in a different (hash) order.
- Header lines are written grouped by type (`@HD`, `@SQ`, `@RG`, `@PG`, `@CO`). HTSlib keeps
  the input's line order and appends new lines at the end, so `addreplacerg` places a new `@RG`
  line before existing `@PG` lines where samtools places it after them.
- `markdup` does not yet support `-S`, `-t`, `-d` (optical duplicates), or barcode options.

## Performance

`scripts/benchmark.sh` times both tools on the corpus (median of 5 warm runs).

Measured on an Intel Xeon Gold 5418Y VM (24 vCPUs) against samtools/HTSlib 1.24, both with
libdeflate. "short" is 1M synthetic 150 bp paired reads, "tagged" 250k reads with extra tags, and
"long" 9,784 nanopore reads with MM/ML modification tags. `-@ N` adds N threads to both tools.

| Operation | samtools (s) | ntools (s) | Speedup |
|---|---:|---:|---:|
| view BAM->SAM, short | 0.381 | 0.346 | 1.10x |
| view BAM->SAM, long | 2.129 | 1.699 | 1.25x |
| view BAM->SAM, long, -@4 | 0.568 | 0.347 | 1.64x |
| view BAM->SAM, long, -@8 | 0.300 | 0.234 | 1.28x |
| view CRAM->SAM, short | 0.880 | 1.513 | 0.58x |
| view CRAM->SAM, long | 2.891 | 3.361 | 0.86x |
| view CRAM->SAM, long, -@8 | 0.570 | 0.453 | 1.26x |
| view CRAM->BAM, long, -@8 | 2.132 | 1.441 | 1.48x |
| view BAM->BAM, long, -@4 | 3.372 | 2.357 | 1.43x |
| view -c -q 30, short | 0.144 | 0.150 | 0.96x |
| flagstat, long, -@4 | 0.330 | 0.340 | 0.97x |
| index, long, -@4 | 0.333 | 0.341 | 0.98x |
| sort, short, -@4 | 0.994 | 0.579 | 1.72x |
| sort -n, short, -@4 | 1.063 | 0.574 | 1.85x |
| collate -Ou, short, -@4 | 2.107 | 0.382 | 5.52x |
| fixmate -m, short | 3.409 | 3.340 | 1.02x |
| fixmate -m, short, -@4 | 0.891 | 0.457 | 1.95x |
| markdup, short | 3.812 | 2.408 | 1.58x |
| markdup, short, -@4 | 0.912 | 0.700 | 1.30x |
| fastq, short | 0.606 | 0.626 | 0.97x |
| depth, short | 0.283 | 0.319 | 0.89x |
| coverage, short | 1.134 | 0.403 | 2.81x |
| faidx 1024 regions | 0.018 | 0.020 | 0.90x |

Where ntools is slower, it is single-threaded: per-thread CRAM decoding still trails HTSlib,
whose rANS decoder uses SIMD. Rows near 1.0x (flagstat, index, counting) are bound by BGZF
decompression, which both tools do with libdeflate. With threads, ntools parallelizes more of
the work than samtools does, so the gap reverses.

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
