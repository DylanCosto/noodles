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
| `markdup` | `markdup` | `-r`, `-l`, `-m t\|s`, `-c`, `-S`, `-t`, `-d` (optical duplicates), `-T`, `-s/-f`, `--json`, `--no-multi-dup`, `--include-fails`, `--use-read-groups`, `--duplicate-count` |
| `stats` | `stats` | all sections; `-r`, `-t`, regions, `-f/-F/-d`, `-c`, `-i`, `-m`, `-q`, `-l`, `-I`, `-S/-P`, `-p`, `-g`, `-x`, `--GC-depth`, `--ref-stats`; per-cycle counts use `-@` threads |

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
- `markdup` does not yet support `--read-coords`/`--coords-order` (coordinates from a regex),
  barcode options, or `--move-umi-to-tag`. With `-d`, duplicates whose names have no
  coordinates are not compared with each other; samtools compares them using uninitialized
  memory, so a chain mixing parseable and unparseable names can come out differently.

## Performance

`scripts/benchmark.sh` times both tools on the corpus (median of 5 warm runs).

Measured on an Intel Xeon Gold 5418Y VM (24 vCPUs) against samtools/HTSlib 1.24, both with
libdeflate. "short" is 1M synthetic 150 bp paired reads, "tagged" 250k reads with extra tags, and
"long" 9,784 nanopore reads with MM/ML modification tags. `-@ N` adds N threads to both tools.

| Operation | samtools (s) | ntools (s) | Speedup |
|---|---:|---:|---:|
| view BAM->SAM, short | 0.381 | 0.352 | 1.08x |
| view BAM->SAM, long | 2.140 | 1.711 | 1.25x |
| view BAM->SAM, long, -@4 | 0.572 | 0.348 | 1.64x |
| view BAM->SAM, long, -@8 | 0.299 | 0.235 | 1.27x |
| view CRAM->SAM, short | 0.888 | 1.063 | 0.84x |
| view CRAM->SAM, long | 2.884 | 2.562 | 1.13x |
| view CRAM->SAM, long, -@8 | 0.568 | 0.354 | 1.60x |
| view CRAM->BAM, long, -@8 | 2.133 | 1.405 | 1.52x |
| view BAM->BAM, long, -@4 | 3.367 | 2.364 | 1.42x |
| view -c -q 30, short | 0.145 | 0.148 | 0.98x |
| flagstat, long, -@4 | 0.331 | 0.339 | 0.98x |
| index, long, -@4 | 0.333 | 0.340 | 0.98x |
| sort, short, -@4 | 0.997 | 0.586 | 1.70x |
| sort -n, short, -@4 | 1.070 | 0.579 | 1.85x |
| collate -Ou, short, -@4 | 2.087 | 0.370 | 5.64x |
| fixmate -m, short | 3.509 | 3.363 | 1.04x |
| fixmate -m, short, -@4 | 0.888 | 0.472 | 1.88x |
| markdup, short | 3.822 | 2.401 | 1.59x |
| markdup, short, -@4 | 0.911 | 0.725 | 1.26x |
| fastq, short | 0.607 | 0.649 | 0.94x |
| depth, short | 0.284 | 0.315 | 0.90x |
| stats, short | 1.658 | 1.085 | 1.53x |
| stats, short, -@4 | 1.584 | 0.661 | 2.40x |
| stats, long, -@8 | 5.674 | 1.484 | 3.82x |
| coverage, short | 1.143 | 0.446 | 2.56x |
| faidx 1024 regions | 0.018 | 0.022 | 0.82x |

Where ntools is slower, it is single-threaded. Short-read CRAM decoding still trails HTSlib in
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
