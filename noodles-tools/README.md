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
| `markdup` | `markdup` | `-r`, `-l`, `-m t\|s`, `-c`, `-S`, `-t`, `-d` (optical duplicates), `-T`, `-s/-f`, `--json`, `--no-multi-dup`, `--include-fails`, `--use-read-groups`, `--duplicate-count`, `--read-coords`/`--coords-order`, `--barcode-tag`/`-name`/`-rgx`, `--move-umi-to-tag`, `--umi-separator` |
| `stats` | `stats` | all sections; `-r`, `-t`, regions, `-f/-F/-d`, `-c`, `-i`, `-m`, `-q`, `-l`, `-I`, `-S/-P`, `-p`, `-g`, `-x`, `--GC-depth`, `--ref-stats`; per-cycle counts use `-@` threads |
| `mpileup` | `mpileup` | text pileup with BAQ (`-B`, `-E`), overlap removal (`-x`), `-f`, `-r`, `-l`, `-a/-aa`, `-q/-Q`, `-d`, `-A`, `-C`, `-G`, `-R`, `-6`, `--rf/--ff`, `-b`, `-X`, `-s`, `-O`, `--output-BP-5`, `--output-QNAME`, `--output-extra`, `--output-sep/-empty`, `--reverse-del`, `--no-output-ins/-del/-ends`, `-M` (base modifications from MM/ML, also in insertions unless `--no-output-ins-mods`); reads and BAQ use `-@` threads, and indexed inputs with several sequences pile up one chromosome per thread |
| `calmd` | `calmd` | `-e`, `-r` (BAQ), `-A`, `-E`, `-n`, `-C`, `-d`, `-q`, `-N`, `-Q`, `-b/-u`; records are processed on `-@` threads |
| `ampliconclip` | `ampliconclip` | soft or `--hard-clip`, `--both-ends`, `--strand`, `--tolerance`, `--original`, `--keep-tag`, `--clipped`, `--fail`, `--filter-len/--fail-len/--unmap-len`, `--no-excluded`, `--rejects-file`, `--primer-counts`, `-f` |
| `head` | `head` | `-h` (header lines), `-n` (records) |
| `flags` | `flags` | numbers and flag names, both ways |
| `samples` | `samples` | `-h`, `-i`, `-T`, `-f/-F` (matching references), `-X`, `-o`, paths from stdin |
| `consensus` | `consensus` | FASTA/FASTQ/pileup; `simple` and Bayesian modes (`bayesian`, `_r`, `_p`, `_m`, `_116`), `-X` presets and `-t` calibration files, `-a/-aa`, `--show-ins/--show-del/--mark-ins`, `-A`, `-r`, `--regions-file`, `-T`, every tuning option; `-@` on an indexed input computes `-Z` blocks in parallel, as samtools does |

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
- `ampliconclip` recalculates the BAM bin of clipped records; samtools keeps the original bin.
  It keeps samtools' quirks otherwise, such as the read name in the `OA` tag's reference field.
- `markdup`: with `-d`, duplicates whose names have no coordinates are not compared with each
  other; samtools compares them using uninitialized memory, so a chain mixing parseable and
  unparseable names can come out differently. `--read-coords` and `--barcode-rgx` take Rust
  `regex` syntax rather than POSIX extended regular expressions (the same for ordinary
  patterns). Without `-S`, `--move-umi-to-tag` moves the UMI of every read; samtools 1.24 moves
  only those of reads still buffered at the end of the input.
- `consensus` follows samtools where it reads past a record: a read padded or deleted after its
  last stored base (e.g. a mapped `SEQ *` read, kept by `-m simple` or `--no-use-MQ`) takes a
  quality from the bytes after the record, which in samtools can be uninitialized memory. ntools
  reads the stale bytes of earlier records there, as a reused `bam1_t` holds them, so those
  columns may differ. Base qualities above 100 read past samtools' probability tables; ntools
  reads the following tables within the same struct and 0 beyond it. With `-@` and a `-T` FASTA
  that lacks a needed sequence, samtools silently writes truncated or garbage sequence; ntools
  stops with an error. As in samtools, `-@` with pileup output on an unindexed input drops the
  pileup lines. CRAM files whose `MD`/`NM` tags HTSlib stripped (to regenerate on decoding) give
  a different consensus, because ntools does not regenerate the tags and `--adj-MQ` uses `MD`.

## Performance

`scripts/benchmark.sh` times both tools on the corpus (median of 5 warm runs).

Measured on an Intel Xeon Gold 5418Y VM (24 vCPUs) against samtools/HTSlib 1.24, both with
libdeflate. "short" is 1M synthetic 150 bp paired reads, "tagged" 250k reads with extra tags, and
"long" 9,784 nanopore reads with MM/ML modification tags. `-@ N` adds N threads to both tools.

| Operation | samtools (s) | ntools (s) | Speedup |
|---|---:|---:|---:|
| view BAM->SAM, short | 0.388 | 0.353 | 1.10x |
| view BAM->SAM, long | 2.129 | 1.742 | 1.22x |
| view BAM->SAM, long, -@4 | 0.569 | 0.353 | 1.61x |
| view BAM->SAM, long, -@8 | 0.299 | 0.246 | 1.22x |
| view CRAM->SAM, short | 0.883 | 0.970 | 0.91x |
| view CRAM->SAM, long | 2.901 | 2.504 | 1.16x |
| view CRAM->SAM, long, -@8 | 0.577 | 0.346 | 1.67x |
| view CRAM->BAM, long, -@8 | 2.167 | 1.456 | 1.49x |
| view BAM->BAM, long, -@4 | 3.409 | 2.373 | 1.44x |
| view -c -q 30, short | 0.145 | 0.150 | 0.97x |
| flagstat, long, -@4 | 0.330 | 0.342 | 0.96x |
| index, long, -@4 | 0.333 | 0.342 | 0.97x |
| sort, short, -@4 | 1.009 | 0.644 | 1.57x |
| sort -n, short, -@4 | 1.146 | 0.636 | 1.80x |
| collate -Ou, short, -@4 | 2.103 | 0.400 | 5.26x |
| fixmate -m, short | 3.496 | 3.352 | 1.04x |
| fixmate -m, short, -@4 | 0.886 | 0.445 | 1.99x |
| markdup, short | 3.820 | 2.384 | 1.60x |
| markdup, short, -@4 | 0.912 | 0.707 | 1.29x |
| fastq, short | 0.613 | 0.623 | 0.98x |
| depth, short | 0.285 | 0.317 | 0.90x |
| stats, short | 1.663 | 1.169 | 1.42x |
| stats, short, -@4 | 1.584 | 0.723 | 2.19x |
| stats, long, -@8 | 5.884 | 1.597 | 3.68x |
| mpileup, short | 2.906 | 2.989 | 0.97x |
| mpileup -f (BAQ), short, -@8 | 27.796 | 4.489 | 6.19x |
| coverage, short | 1.151 | 0.406 | 2.83x |
| consensus, short | 6.37 | 6.08 | 1.05x |
| consensus, short, -@4 | 1.61 | 1.54 | 1.05x |
| consensus, long | 11.51 | 10.89 | 1.06x |
| consensus, long, -@4 | 3.43 | 3.23 | 1.06x |
| faidx 1024 regions | 0.018 | 0.020 | 0.90x |

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
