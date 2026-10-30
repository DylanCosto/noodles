"""Writes a SAM file with a read of 70,000 CIGAR operations, which BAM stores in a CG tag.

Usage: long-cigar-test-data.py > reads.sam

A second, ordinary read on the other strand shares the read group, so checksums and outputs
cover both kinds of record.
"""

import random

rng = random.Random(3)
OPS = 35_000

cigar = "1M1I" * OPS
seq = "".join(rng.choice("ACGT") for _ in range(2 * OPS))
qual = "".join(chr(33 + rng.randint(2, 40)) for _ in range(2 * OPS))

print("@HD\tVN:1.6\tSO:coordinate")
print("@SQ\tSN:c\tLN:100000")
print("@RG\tID:g")
print(f"long\t0\tc\t100\t60\t{cigar}\t*\t0\t0\t{seq}\t{qual}\tRG:Z:g\tXA:i:5")
print("short\t16\tc\t200\t60\t10M\t*\t0\t0\tACGTACGTAC\tIIIIIIIIII\tRG:Z:g")
