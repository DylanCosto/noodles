"""Writes a synthetic reference and SAM file for checking mpileup --output-mods against samtools.

Usage: basemod-test-data.py REFERENCE_FASTA [SEED] > reads.sam

Reads on both strands carry MM tags with single, linked (C+hm), ChEBI (C+76792), N (N+n), and
opposite-strand (G-m) modifications, with and without ML, and implicit (.) or explicit (?)
skips. Some have insertions (with modified inserted bases), deletions, or both.
"""

import random
import sys

rng = random.Random(int(sys.argv[2]) if len(sys.argv) > 2 else 7)
REF_LEN = 400
COMPLEMENT = {"A": "T", "C": "G", "G": "C", "T": "A", "N": "N"}

ref = "".join(rng.choice("ACGT") for _ in range(REF_LEN))

with open(sys.argv[1], "w") as fasta:
    fasta.write(">c1\n")
    fasta.write("".join(ref[i : i + 60] + "\n" for i in range(0, REF_LEN, 60)))


def modifications(seq, rev, specs):
    """MM and ML values calling random bases of each spec's canonical base."""
    original = "".join(COMPLEMENT[c] for c in reversed(seq)) if rev else seq
    mm, ml = [], []

    for base, strand, codes, skip in specs:
        n = sum(1 for c in original if c == base or base == "N")
        deltas, called = [], 0

        while True:
            d = rng.randint(0, 3)
            if called + d + 1 > n:
                break
            deltas.append(d)
            called += d + 1

        stride = 1 if codes.isdigit() else len(codes)
        mm.append(f"{base}{strand}{codes}{skip}" + "".join(f",{d}" for d in deltas) + ";")
        ml += [rng.randint(0, 255) for _ in range(len(deltas) * stride)]

    return "".join(mm), ml


SPECS = [
    [("C", "+", "m", ".")],
    [("C", "+", "hm", "?")],
    [("C", "+", "76792", "."), ("A", "+", "a", "")],
    [("N", "+", "n", "?")],
    [("G", "-", "m", "."), ("C", "+", "m", "")],
    [("T", "+", "g", "."), ("C", "+", "mh", ".")],
]

reads = []

for k in range(300):
    pos = rng.randint(1, 300)
    length = rng.randint(30, 80)
    seq = [b if rng.random() >= 0.05 else rng.choice("ACGT") for b in ref[pos - 1 : pos - 1 + length]]
    at = rng.randint(5, length - 15)
    kind = rng.random()

    if kind < 0.4:
        ins = [rng.choice("ACGT") for _ in range(rng.randint(1, 5))]
        seq = seq[:at] + ins + seq[at:]
        cigar = f"{at}M{len(ins)}I{length - at}M"
    elif kind < 0.6:
        d = rng.randint(1, 4)
        seq = seq[:at] + seq[at + d :]
        cigar = f"{at}M{d}D{length - at - d}M"
    elif kind < 0.7:
        seq = seq[:at] + [rng.choice("ACGT") for _ in range(2)] + seq[at + 2 :]
        cigar = f"{at}M2I2D{length - at - 2}M"
    else:
        cigar = f"{length}M"

    seq = "".join(seq)
    rev = rng.random() < 0.5
    mm, ml = modifications(seq, rev, rng.choice(SPECS))
    tags = [f"MM:Z:{mm}"]

    if ml and rng.random() < 0.85:
        tags.append("ML:B:C," + ",".join(map(str, ml)))

    qual = "".join(chr(33 + rng.randint(5, 40)) for _ in seq)
    fields = [f"r{k}", 16 if rev else 0, "c1", pos, rng.randint(0, 60), cigar, "*", 0, 0, seq, qual]
    reads.append((pos, "\t".join(map(str, fields + tags))))

sys.stdout.write(f"@HD\tVN:1.6\tSO:coordinate\n@SQ\tSN:c1\tLN:{REF_LEN}\n")

for _, line in sorted(reads, key=lambda r: r[0]):
    sys.stdout.write(line + "\n")
