# noodles-pileup

An htslib-compatible pileup engine for encoded BAM records, with base alignment quality (BAQ),
mapping quality capping, and MD/NM calculation.

These are ports of htslib 1.24's `bam_plp_*` (`sam.c`), `sam_prob_realn` and `sam_cap_mapq`
(`realn.c`, `probaln.c`), and samtools 1.24's `bam_fillmd1` (`bam_md.c`). They keep the same
rules and floating-point operations, so pileups, adjusted qualities, and tags match samtools and
htslib byte for byte. `ntools mpileup` and `ntools calmd` are built on them.

```rust,no_run
use noodles_pileup::{Overlaps, Pileup, Read};

# fn records() -> Vec<Vec<u8>> { Vec::new() }
let mut pileup = Pileup::new(8000, Overlaps::Adjust);

for raw in records() {
    pileup.push(Some(Read::new(raw)?))?;

    while pileup.next_position()? {
        let (tid, pos) = pileup.position();

        for entry in pileup.entries() {
            let read = pileup.read(entry);
            let _ = (tid, pos, read.name(), entry.qpos);
        }
    }
}

pileup.push(None)?;

while pileup.next_position()? {
    // ...
}
# Ok::<_, Box<dyn std::error::Error>>(())
```
