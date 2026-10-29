"""Writes a synthetic SAM file to stdout for checking markdup -S/-t/-d against samtools.

Usage: markdup-test-data.py [SEED] [--barcodes]

Reads have Illumina-style names (three formats) with optical and library duplicates,
supplementary alignments (SA), XA tags, mate-unmapped pairs, secondary alignments, single reads,
QC fails, and a few names whose coordinates cannot be parsed. Sort it by name, run
`samtools fixmate -m`, and sort it by coordinate before marking duplicates.

With --barcodes, names with UMIs get one of three, and templates get a BC:Z tag of one of two
barcodes (some none, and some a BC:i tag of the wrong type), for --barcode-tag/-name/-rgx.
"""

import random
import sys
import zlib

args = [a for a in sys.argv[1:] if a != "--barcodes"]
barcodes = "--barcodes" in sys.argv[1:]
rng = random.Random(int(args[0]) if args else 7)
REF_LEN = 200_000
READ_LEN = 100

out = sys.stdout
out.write("@HD\tVN:1.6\tSO:unsorted\n")
out.write(f"@SQ\tSN:chr1\tLN:{REF_LEN}\n@SQ\tSN:chr2\tLN:{REF_LEN}\n")
out.write("@RG\tID:a\tSM:s\n@RG\tID:b\tSM:s\n")

lines = []
serial = 0


def qual():
    base = rng.randint(10, 40)
    return "".join(chr(33 + max(2, min(41, base + rng.randint(-8, 8)))) for _ in range(READ_LEN))


def seq():
    return "".join(rng.choice("ACGT") for _ in range(READ_LEN))


def name(tile, x, y, style):
    global serial
    serial += 1
    if style == 0:
        return f"M{serial % 3}:7:FC1:1:{tile}:{x}:{y}"
    if style == 1:
        return f"HS:{tile % 3}:{tile}:{x}:{y}"  # HiSeq style, 4 colons
    if style == 2:
        return f"M:7:FC1:1:{tile}:{x}:{y}:ACGTACGT"  # with UMI, 7 colons
    return f"odd_name_{serial}"


def cigar():
    r = rng.random()
    if r < 0.15:
        clip = rng.randint(1, 10)
        return f"{clip}S{READ_LEN - clip}M"
    if r < 0.25:
        clip = rng.randint(1, 10)
        return f"{READ_LEN - clip}M{clip}S"
    return f"{READ_LEN}M"


def pair(qname, ref, pos1, pos2, rev1, rg, kind):
    """kind: normal, supp, munmap, xa, qcfail"""
    f1 = 1 | 2 | 64 | (16 if rev1 else 0) | (32 if not rev1 else 0)
    f2 = 1 | 2 | 128 | (0 if rev1 else 16) | (16 * 0) | (32 if rev1 else 0)
    if not rev1:
        f2 |= 16
    if kind == "qcfail":
        f1 |= 512
        f2 |= 512
    c1, c2 = cigar(), cigar()
    tags1 = [f"RG:Z:{rg}"]
    tags2 = [f"RG:Z:{rg}"]
    if kind == "munmap":
        f1 = 1 | 8 | 64 | (16 if rev1 else 0)
        f2 = 1 | 4 | 128
        lines.append([qname, f1, ref, pos1, 60, c1, "=", pos1, 0, seq(), qual()] + tags1)
        lines.append([qname, f2, ref, pos1, 0, "*", "=", pos1, 0, seq(), qual()] + tags2)
        return
    if kind == "supp":
        spos = rng.randint(1, REF_LEN - READ_LEN)
        sa = f"SA:Z:chr2,{spos},+,50M50S,30,0;"
        tags1.append(sa)
        lines.append([qname, 2048 | (f1 & ~(2)), "chr2", spos, 30, "50M50H", ref, pos2, 0,
                      seq()[:50], qual()[:50], f"SA:Z:{ref},{pos1},{'-' if rev1 else '+'},{c1},60,0;",
                      f"RG:Z:{rg}"])
    if kind == "xa":
        tags2.append(f"XA:Z:chr2,+{rng.randint(1, REF_LEN)},100M,1;")
    if kind == "secondary":
        lines.append([qname, 256 | f1, "chr2", rng.randint(1, REF_LEN - READ_LEN), 0, "100M", ref, pos2,
                      0, "*", "*", f"RG:Z:{rg}"])
    tlen = pos2 + READ_LEN - pos1
    lines.append([qname, f1, ref, pos1, 60, c1, "=", pos2, tlen, seq(), qual()] + tags1)
    lines.append([qname, f2, ref, pos2, 60, c2, "=", pos1, -tlen, seq(), qual()] + tags2)


for fragment in range(4000):
    ref = rng.choice(["chr1", "chr2"])
    pos1 = rng.randint(1, REF_LEN - 1000)
    pos2 = pos1 + rng.randint(0, 400)
    rev1 = rng.random() < 0.5
    rg = rng.choice("ab")
    copies = 1 + (rng.random() < 0.35) * rng.randint(1, 5)
    tile = rng.choice([1101, 1102, 2101])
    x, y = rng.randint(1000, 30000), rng.randint(1000, 30000)
    style = rng.choices([0, 1, 2, 3], weights=[80, 8, 8, 4])[0]
    single = rng.random() < 0.08

    for copy in range(copies):
        if copy == 0 or rng.random() < 0.5:
            # optical neighbour (within 100 px) or far away
            dx, dy = (rng.randint(-120, 120), rng.randint(-120, 120)) if copy else (0, 0)
            ctile = tile
        else:
            dx, dy = rng.randint(-20000, 20000), rng.randint(-20000, 20000)
            ctile = rng.choice([1101, 1102, 2101])
        qname = name(ctile, max(1, x + dx), max(1, y + dy), style)
        if single:
            flag = 16 if rev1 else 0
            if rng.random() < 0.1:
                flag |= 512
            tags = [f"RG:Z:{rg}"]
            if rng.random() < 0.2:
                spos = rng.randint(1, REF_LEN - READ_LEN)
                tags.append(f"SA:Z:chr2,{spos},+,40M60S,20,0;")
                lines.append([qname, 2048, "chr2", spos, 20, "40M60H", "*", 0, 0, seq()[:40], qual()[:40],
                              f"SA:Z:{ref},{pos1},+,{READ_LEN}M,60,0;", f"RG:Z:{rg}"])
            lines.append([qname, flag, ref, pos1, 60, cigar(), "*", 0, 0, seq(), qual()] + tags)
        else:
            kind = rng.choices(["normal", "supp", "munmap", "xa", "qcfail", "secondary"],
                               weights=[70, 10, 6, 6, 4, 4])[0]
            pair(qname, ref, pos1, pos2, rev1, rg, kind)

for line in lines:
    if barcodes:
        h = zlib.crc32(line[0].encode())
        if line[0].endswith(":ACGTACGT"):
            line[0] = line[0][:-8] + ["AAAA", "CCGG", "TTAG"][h % 3]
        if h % 17 == 0:
            line.append("BC:i:5")
        elif h % 5 != 0:
            line.append("BC:Z:" + ["GATT", "CAGA"][h % 2])
    out.write("\t".join(str(f) for f in line) + "\n")
