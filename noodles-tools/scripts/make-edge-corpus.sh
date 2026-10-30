#!/usr/bin/env bash
# Builds a small, deterministic corpus of edge-case alignment files for compare-edge.sh.
#
# Usage: make-edge-corpus.sh OUTDIR [SAMTOOLS]
#
# Each case is written as NAME.sam (hand-made, in the order written), then converted with samtools
# to NAME.bam (coordinate-sorted, with .bai and .csi), NAME.name.bam (name-sorted),
# NAME.3.0.cram, NAME.3.1.cram, NAME.2.1.cram, and NAME.embed.cram (3.0 with an embedded
# reference), each with a .crai. Cases that samtools cannot sort or convert keep only their SAM
# file; the list of converted cases is written to OUTDIR/cases.
#
# The reference is ref.fa (chr1 200 bp, chr2 120 bp, chrL 40,000 bp) and manysq.fa (400 tiny
# sequences) for the manysq case, both with .fai indexes. edge.fa and edge.fq hold FASTA and
# FASTQ edge cases for faidx, fqidx, dict, and import.
set -euo pipefail

out=$1
st=${2:-samtools}

mkdir -p "$out"
cd "$out"

# Writes a pseudo-random DNA sequence of length $2 seeded with $1 (a 32-bit LCG, so the output
# does not depend on the awk implementation's random number generator).
dna() {
    awk -v seed="$1" -v n="$2" 'BEGIN {
        x = seed; b = "ACGT"; s = ""
        for (i = 0; i < n; i++) {
            x = (x * 69069 + 1) % 4294967296
            s = s substr(b, int(x / 1073741824) + 1, 1)
        }
        print s
    }'
}

# Writes a FASTA record wrapped at $3 columns.
fasta() {
    printf '>%s\n' "$1"
    printf '%s\n' "$2" | fold -w "$3"
}

# Repeats string $1 $2 times.
rep() {
    local s=$1 n=$2 r=""
    while [ "$n" -gt 0 ]; do
        if [ $((n & 1)) -eq 1 ]; then r=$r$s; fi
        s=$s$s
        n=$((n >> 1))
    done
    printf '%s' "$r"
}

chr1=$(dna 1 200)
chr2=$(dna 2 120)
chrL=$(dna 3 40000)

{
    fasta chr1 "$chr1" 60
    fasta chr2 "$chr2" 60
    fasta chrL "$chrL" 60
} >ref.fa
"$st" faidx ref.fa

# Prints $3 bases of reference sequence $1 from 1-based position $2.
sub() {
    local s
    case $1 in
        chr1) s=$chr1 ;;
        chr2) s=$chr2 ;;
        chrL) s=$chrL ;;
    esac
    printf '%s' "${s:$(($2 - 1)):$3}"
}

# Prints $1 copies of quality character $2 (default I).
qual() {
    rep "${2:-I}" "$1"
}

# Prints a tab-separated SAM line from its arguments.
rec() {
    local IFS=$'\t'
    printf '%s\n' "$*"
}

hd=$'@HD\tVN:1.6\tSO:unsorted'
sq=$'@SQ\tSN:chr1\tLN:200\n@SQ\tSN:chr2\tLN:120\n@SQ\tSN:chrL\tLN:40000'

# empty_none: no header and no records.
: >empty_none.sam

# empty_nosq: a header without @SQ lines and no records.
printf '%s\n' "$hd" >empty_nosq.sam

# empty_sq: @SQ lines but no records.
printf '%s\n%s\n' "$hd" "$sq" >empty_sq.sam

# unmapped: only unmapped reads, placed and unplaced, single and paired.
{
    printf '%s\n%s\n' "$hd" "$sq"
    rec u_unplaced 4 '*' 0 0 '*' '*' 0 0 ACGTACGTAC "$(qual 10)"
    rec u_placed 4 chr1 50 0 '*' '*' 0 0 "$(sub chr1 50 12)" "$(qual 12 '#')"
    rec u_placed_end 4 chr2 120 0 '*' '*' 0 0 ACGTA IIIII
    rec u_pair 77 '*' 0 0 '*' '*' 0 0 AAAAACCCCC "$(qual 10)"
    rec u_pair 141 '*' 0 0 '*' '*' 0 0 GGGGGTTTTT "$(qual 10)"
    rec u_ppair 77 chr1 30 0 '*' = 30 0 "$(sub chr1 30 8)" "$(qual 8)"
    rec u_ppair 141 chr1 30 0 '*' = 30 0 "$(sub chr1 40 8)" "$(qual 8)"
    rec u_rev 20 chr1 60 0 '*' '*' 0 0 ACGTN IIIII
} >unmapped.sam

# noseq: SEQ and/or QUAL missing.
{
    printf '%s\n%s\n' "$hd" "$sq"
    rec s_both 0 chr1 10 60 10M '*' 0 0 '*' '*'
    rec s_qual 0 chr1 20 60 10M '*' 0 0 "$(sub chr1 20 10)" '*'
    rec s_qual_rev 16 chr1 25 60 4M1I5M '*' 0 0 "$(sub chr1 25 4)A$(sub chr1 29 5)" '*'
    rec s_unmapped 4 '*' 0 0 '*' '*' 0 0 '*' '*'
    rec s_unmapped_seq 4 '*' 0 0 '*' '*' 0 0 ACGT '*'
    rec s_secondary 256 chr1 15 0 10M '*' 0 0 '*' '*' NM:i:0
    rec s_lowq 0 chr2 5 60 6M '*' 0 0 "$(sub chr2 5 6)" '!!!!!!'
    rec s_highq 0 chr2 7 60 6M '*' 0 0 "$(sub chr2 7 6)" '~~~~~~'
} >noseq.sam

# secsup: secondary and supplementary alignments.
{
    printf '%s\n%s\n' "$hd" "$sq"
    rec chim 0 chr1 10 60 10M10S '*' 0 0 "$(sub chr1 10 10)$(sub chr2 50 10)" "$(qual 20)" SA:Z:chr2,50,+,10H10M,50,0\;
    rec chim 2048 chr2 50 50 10H10M '*' 0 0 "$(sub chr2 50 10)" "$(qual 10)" SA:Z:chr1,10,+,10M10S,60,0\;
    rec chim 256 chr1 100 0 10M10S '*' 0 0 '*' '*'
    rec secsup 2304 chr2 70 0 5M '*' 0 0 '*' '*'
    rec sec_unmapped 260 '*' 0 0 '*' '*' 0 0 ACGT IIII
    rec sup_rev 2064 chr1 150 10 3S7M '*' 0 0 "ACG$(sub chr1 150 7)" "$(qual 10)"
    rec pair_sup 2145 chr1 120 40 5M5H chr2 30 0 "$(sub chr1 120 5)" IIIII
    rec pair_sup 99 chr1 120 60 5M5S = 140 30 "$(sub chr1 120 5)GGGGG" "$(qual 10)"
    rec pair_sup 147 chr1 140 60 10M = 120 -30 "$(sub chr1 140 10)" "$(qual 10)"
} >secsup.sam

# longcigar: 65535, 65536, and 70,000 CIGAR operations (BAM stores the last two in a CG tag).
{
    printf '%s\n%s\n@RG\tID:g\tSM:x\n' "$hd" "$sq"
    n=32767
    rec lc_65535 0 chrL 100 60 "$(rep 1M1I $n)1M" '*' 0 0 "$(dna 11 $((2 * n + 1)))" "$(qual $((2 * n + 1)))" XA:i:1
    n=32768
    rec lc_65536 16 chrL 200 60 "$(rep 1M1I $n)" '*' 0 0 "$(dna 12 $((2 * n)))" "$(qual $((2 * n)) 5)"
    n=35000
    rec lc_70000 0 chrL 300 60 "$(rep 1M1I $n)" '*' 0 0 "$(dna 13 $((2 * n)))" '*' RG:Z:g
    rec lc_short 0 chrL 1000 60 10M '*' 0 0 "$(sub chrL 1000 10)" "$(qual 10)"
} >longcigar.sam

# lengths: one-base, sequence-less, and very long reads.
{
    printf '%s\n%s\n' "$hd" "$sq"
    rec l_one 0 chr1 5 60 1M '*' 0 0 "$(sub chr1 5 1)" I
    rec l_zero 4 '*' 0 0 '*' '*' 0 0 '*' '*'
    rec l_zero_mapped 0 chr1 7 60 5M '*' 0 0 '*' '*'
    rec l_long 0 chrL 5 60 30000M '*' 0 0 "$(sub chrL 5 30000)" "$(qual 30000)"
    rec l_long_unmapped 4 '*' 0 0 '*' '*' 0 0 "$(dna 21 100000)" "$(qual 100000 '+')"
} >lengths.sam

# ends: reads at position 1, at the reference end, and overhanging it.
{
    printf '%s\n%s\n' "$hd" "$sq"
    rec e_first 0 chr1 1 60 10M '*' 0 0 "$(sub chr1 1 10)" "$(qual 10)"
    rec e_first_clip 0 chr1 1 60 5S10M '*' 0 0 "AAAAA$(sub chr1 1 10)" "$(qual 15)"
    rec e_first_ins 0 chr1 1 60 2I8M '*' 0 0 "TT$(sub chr1 1 8)" "$(qual 10)"
    rec e_last 0 chr1 191 60 10M '*' 0 0 "$(sub chr1 191 10)" "$(qual 10)"
    rec e_last_one 16 chr2 120 60 1M '*' 0 0 "$(sub chr2 120 1)" I
    rec e_over 0 chr1 195 60 10M '*' 0 0 "$(sub chr1 195 6)ACGT" "$(qual 10)"
    rec e_over_del 0 chr2 115 60 3M10D2M '*' 0 0 "$(sub chr2 115 3)AC" IIIII
    rec e_end_ins 0 chr1 196 60 5M3I '*' 0 0 "$(sub chr1 196 5)GGG" "$(qual 8)"
    rec e_chrL_end 0 chrL 39991 60 10M '*' 0 0 "$(sub chrL 39991 10)" "$(qual 10)"
} >ends.sam

# beyond: a read placed entirely past the end of its reference.
{
    printf '%s\n%s\n' "$hd" "$sq"
    rec b_in 0 chr1 100 60 10M '*' 0 0 "$(sub chr1 100 10)" "$(qual 10)"
    rec b_past 0 chr1 250 60 10M '*' 0 0 ACGTACGTAC "$(qual 10)"
} >beyond.sam

# cigarops: =, X, N, P, H, S, I, D in various places.
{
    printf '%s\n%s\n' "$hd" "$sq"
    s=$(sub chr1 10 10)
    rec c_eqx 0 chr1 10 60 5=1X4= '*' 0 0 "${s:0:5}$([ "${s:5:1}" = A ] && echo C || echo A)${s:6:4}" "$(qual 10)" NM:i:1 MD:Z:5"${s:5:1}"4
    rec c_skip 0 chr1 30 60 3M20N5M '*' 0 0 "$(sub chr1 30 3)$(sub chr1 53 5)" "$(qual 8)"
    rec c_pad 0 chr1 60 60 2M1P1I3M '*' 0 0 "$(sub chr1 60 2)A$(sub chr1 62 3)" "$(qual 6)"
    rec c_hard 16 chr1 70 60 3H5M2H '*' 0 0 "$(sub chr1 70 5)" "$(qual 5)"
    rec c_soft_del 0 chr1 80 60 2S3M1D3M2S '*' 0 0 "GG$(sub chr1 80 3)$(sub chr1 84 3)CC" "$(qual 10)"
    rec c_lead_del 0 chr1 90 60 2D5M '*' 0 0 "$(sub chr1 92 5)" "$(qual 5)"
    rec c_trail_skip 0 chr1 100 60 5M5N '*' 0 0 "$(sub chr1 100 5)" "$(qual 5)"
    rec c_del_ins 0 chr1 110 60 3M1D1I3M '*' 0 0 "$(sub chr1 110 3)T$(sub chr1 114 3)" "$(qual 7)"
    rec c_only_ins 0 chr1 120 60 4I '*' 0 0 ACGT IIII
    rec c_xx 0 chr1 130 60 3X '*' 0 0 NNN III
    rec c_sh 0 chr1 140 60 2H2S4M2S2H '*' 0 0 "AA$(sub chr1 140 4)TT" "$(qual 8)"
} >cigarops.sam

# mates: mates on different references, mates absent, overlapping pairs.
{
    printf '%s\n%s\n' "$hd" "$sq"
    rec m_xref 97 chr1 10 60 10M chr2 20 0 "$(sub chr1 10 10)" "$(qual 10)" MC:Z:10M
    rec m_xref 145 chr2 20 60 10M chr1 10 0 "$(sub chr2 20 10)" "$(qual 10)" MC:Z:10M
    rec m_alone 65 chr1 40 60 10M chr1 300 0 "$(sub chr1 40 10)" "$(qual 10)"
    rec m_overlap 99 chr1 50 60 20M = 60 30 "$(sub chr1 50 20)" "$(qual 20)"
    rec m_overlap 147 chr1 60 60 20M = 50 -30 "$(sub chr1 60 20)" "$(qual 20 '5')"
    rec m_munmap 73 chr1 80 60 10M = 80 0 "$(sub chr1 80 10)" "$(qual 10)"
    rec m_munmap 133 chr1 80 0 '*' = 80 0 ACGTACGTAC "$(qual 10)"
    rec m_badtlen 99 chr1 100 60 10M = 150 999 "$(sub chr1 100 10)" "$(qual 10)"
    rec m_badtlen 147 chr1 150 60 10M = 100 5 "$(sub chr1 150 10)" "$(qual 10)"
    rec m_three 65 chr2 5 60 5M chr2 50 0 "$(sub chr2 5 5)" IIIII
    rec m_three 129 chr2 50 60 5M chr2 5 0 "$(sub chr2 50 5)" IIIII
    rec m_three 129 chr2 60 60 5M chr2 5 0 "$(sub chr2 60 5)" IIIII
    rec m_same 163 chr1 170 60 10M = 170 -10 "$(sub chr1 170 10)" "$(qual 10)"
    rec m_same 83 chr1 170 60 10M = 170 10 "$(sub chr1 170 10)" "$(qual 10)"
    rec m_samestart 163 chr1 180 60 5M = 180 -10 "$(sub chr1 180 5)" IIIII
    rec m_samestart 83 chr1 180 60 10M = 180 10 "$(sub chr1 180 10)" "$(qual 10)"
    rec m_wrongpos 99 chr2 70 60 10M = 75 20 "$(sub chr2 70 10)" "$(qual 10)"
    rec m_wrongpos 147 chr2 80 60 10M = 70 -20 "$(sub chr2 80 10)" "$(qual 10)"
} >mates.sam

# headers: unusual @HD/@SQ/@RG/@PG/@CO lines.
{
    printf '@HD\tVN:1.6\tSO:unsorted\tGO:query\tzz:x\n'
    printf '@SQ\tSN:chr1\tLN:200\tAS:asm:1\tAN:one,1\tSP:human\n'
    printf '@SQ\tSN:chr2\tLN:120\tM5:%s\n' "$(printf '%s' "$chr2" | md5sum | cut -c1-32)"
    printf '@SQ\tSN:chrL\tLN:40000\tDS:with spaces and : colons\n'
    printf '@RG\tID:rg1\tSM:s1\tLB:lib1\tDS:a:b:c\tPL:ILLUMINA\n'
    printf '@RG\tID:rg2\tSM:s 2\tLB:lib1\tPU:x:y\n'
    printf '@RG\tID:rg:3\tSM:s1\n'
    printf '@PG\tID:p1\tPN:aligner\tVN:1:2\tCL:aligner -x --y=a:b\n'
    printf '@PG\tID:p2\tPN:sorter\tPP:p1\n'
    printf '@PG\tID:p3\tPN:marker\tPP:p2\tCL:a:b\n'
    printf '@PG\tID:p2.1\tPN:other\tPP:p1\n'
    printf '@CO\tfree\ttext:with\ttabs\n'
    printf '@CO\t\n'
    printf '@CO\tsecond: comment\n'
    rec h_rg1 0 chr1 10 60 10M '*' 0 0 "$(sub chr1 10 10)" "$(qual 10)" RG:Z:rg1 PG:Z:p3
    rec h_rg2 16 chr1 20 60 10M '*' 0 0 "$(sub chr1 20 10)" "$(qual 10)" RG:Z:rg2
    rec h_rg3 0 chr2 30 60 10M '*' 0 0 "$(sub chr2 30 10)" "$(qual 10)" RG:Z:rg:3 PG:Z:p2.1
    rec h_rgmissing 0 chr2 40 60 10M '*' 0 0 "$(sub chr2 40 10)" "$(qual 10)" RG:Z:nosuch
    rec h_norg 4 '*' 0 0 '*' '*' 0 0 ACGT IIII
} >headers.sam

# badheader: a header field without a colon (an error for samtools).
{
    printf '%s\n%s\n' "$hd" "$sq"
    printf '@PG\tID:p1\tPN:aligner\tCL:aligner -x\t--y\n'
    rec bh 0 chr1 10 60 10M '*' 0 0 "$(sub chr1 10 10)" "$(qual 10)"
} >badheader.sam

# badpp: a @PG line whose PP names a missing program (HTSlib warns and keeps it).
{
    printf '%s\n%s\n' "$hd" "$sq"
    printf '@PG\tID:p1\tPN:one\tPP:nosuch\n@PG\tID:p2\tPN:two\tPP:p1\n'
    rec bp 0 chr1 10 60 10M '*' 0 0 "$(sub chr1 10 10)" "$(qual 10)"
} >badpp.sam

# duprg: duplicate @RG and @PG IDs (samtools may reject these).
{
    printf '%s\n%s\n' "$hd" "$sq"
    printf '@RG\tID:a\tSM:x\n@RG\tID:a\tSM:y\n'
    printf '@PG\tID:q\tPN:one\n@PG\tID:q\tPN:two\n@PG\tID:r\tPN:three\tPP:nosuch\n'
    rec d_a 0 chr1 10 60 10M '*' 0 0 "$(sub chr1 10 10)" "$(qual 10)" RG:Z:a
} >duprg.sam

# aux: every tag type, each B subtype, empty arrays, and integer limits.
{
    printf '%s\n%s\n' "$hd" "$sq"
    s=$(sub chr1 10 10)
    q=$(qual 10)
    rec a_scalar 0 chr1 10 60 10M '*' 0 0 "$s" "$q" XA:A:x Xc:i:-128 XC:i:255 Xs:i:-32768 XS:i:65535 \
        Xi:i:-2147483648 XI:i:4294967295 Xf:f:3.14159 Xz:f:-0 Xe:f:1e-30 XE:f:1e+30 XZ:Z:hello\ world XH:H:1AE301
    rec a_empty 0 chr1 11 60 10M '*' 0 0 "$s" "$q" XZ:Z: XH:H: Bc:B:c BC:B:C Bs:B:s BS:B:S Bi:B:i BI:B:I Bf:B:f
    rec a_arrays 0 chr1 12 60 10M '*' 0 0 "$s" "$q" Bc:B:c,-128,0,127 BC:B:C,0,255 Bs:B:s,-32768,32767 \
        BS:B:S,0,65535 Bi:B:i,-2147483648,2147483647 BI:B:I,0,4294967295 Bf:B:f,0.5,-1.25,1e10
    rec a_ints 0 chr1 13 60 10M '*' 0 0 "$s" "$q" X0:i:0 X1:i:127 X2:i:128 X3:i:-129 X4:i:32768 X5:i:-32769 X6:i:65536
    rec a_special 0 chr1 14 60 10M '*' 0 0 "$s" "$q" Xn:f:nan Xp:f:inf Xm:f:-inf XA:A:~ XB:Z:'tab-free:colons;semi'
    rec a_md 0 chr1 15 60 10M '*' 0 0 "$s" "$q" MD:Z:10 NM:i:0 MC:Z:10M AS:i:10 XS:i:5
    rec a_unmapped 4 '*' 0 0 '*' '*' 0 0 ACGT IIII BZ:B:c,1 ZZ:Z:x
} >aux.sam

# names: one-character and 254-character names, and pair suffixes.
{
    printf '%s\n%s\n' "$hd" "$sq"
    long=$(rep abcdefghij 25)abcd
    rec a 0 chr1 10 60 10M '*' 0 0 "$(sub chr1 10 10)" "$(qual 10)"
    rec "$long" 0 chr1 20 60 10M '*' 0 0 "$(sub chr1 20 10)" "$(qual 10)"
    rec r/1 0 chr1 30 60 10M '*' 0 0 "$(sub chr1 30 10)" "$(qual 10)"
    rec r/2 16 chr1 40 60 10M '*' 0 0 "$(sub chr1 40 10)" "$(qual 10)"
    rec 'x:1:2#0/1' 4 '*' 0 0 '*' '*' 0 0 ACGT IIII
    rec '!~' 4 '*' 0 0 '*' '*' 0 0 ACGT IIII
    rec "$long" 4 '*' 0 0 '*' '*' 0 0 ACGT IIII
} >names.sam

# flags: odd flag combinations.
{
    printf '%s\n%s\n' "$hd" "$sq"
    s=$(sub chr1 10 10)
    q=$(qual 10)
    rec f_paired_only 1 chr1 10 60 10M '*' 0 0 "$s" "$q"
    rec f_qcfail 512 chr1 10 60 10M '*' 0 0 "$s" "$q"
    rec f_dup 1024 chr1 10 60 10M '*' 0 0 "$s" "$q"
    rec f_dup_qc 1536 chr1 10 60 10M '*' 0 0 "$s" "$q"
    rec f_mate_unmapped_unpaired 8 chr1 10 60 10M '*' 0 0 "$s" "$q"
    rec f_both_ends 193 chr1 10 60 10M = 10 0 "$s" "$q"
    rec f_neither_end 33 chr1 10 60 10M = 10 0 "$s" "$q"
    rec f_all 4095 chr1 10 60 10M = 10 0 "$s" "$q"
    rec f_unmapped_all 3583 chr1 10 0 '*' = 10 0 "$s" "$q"
    rec f_proper_single 3 chr1 10 60 10M '*' 0 0 "$s" "$q"
    rec f_rev_unmapped 52 '*' 0 0 '*' '*' 0 0 "$s" "$q"
} >flags.sam

# mapq: MAPQ 255 (unavailable), 254, and 0.
{
    printf '%s\n%s\n' "$hd" "$sq"
    s=$(sub chr1 10 10)
    q=$(qual 10)
    rec q255 0 chr1 10 255 10M '*' 0 0 "$s" "$q"
    rec q254 0 chr1 10 254 10M '*' 0 0 "$s" "$q"
    rec q0 0 chr1 10 0 10M '*' 0 0 "$s" "$q"
    rec q255u 4 chr1 10 255 '*' '*' 0 0 "$s" "$q"
    rec q255s 256 chr1 10 255 10M '*' 0 0 '*' '*'
} >mapq.sam

# manysq: 400 tiny reference sequences.
{
    for i in $(seq 0 399); do
        fasta "s$i" "$(dna $((100 + i)) $((i % 37 + 1)))" 60
    done
} >manysq.fa
"$st" faidx manysq.fa
{
    printf '%s\n' "$hd"
    awk '{ printf "@SQ\tSN:%s\tLN:%s\n", $1, $2 }' manysq.fa.fai
    awk '$2 >= 3 && NR % 7 == 1 { print $1, $2 }' manysq.fa.fai | while read -r name len; do
        s=$("$st" faidx manysq.fa "$name" | tail -n +2 | tr -d '\n')
        rec "m_$name" 0 "$name" 1 60 "${len}M" '*' 0 0 "$s" "$(qual "$len")"
        rec "m_${name}_2" 16 "$name" 2 60 "$((len - 1))M" '*' 0 0 "${s:1}" "$(qual $((len - 1)))"
    done
    rec m_last 0 s399 1 60 5M '*' 0 0 "$("$st" faidx manysq.fa s399:1-5 | tail -n +2)" IIIII
    rec m_unmapped 4 '*' 0 0 '*' '*' 0 0 ACGT IIII
} >manysq.sam

# FASTA edge cases: an empty sequence, lowercase and IUPAC bases, short final lines, no final
# newline.
{
    printf '>empty\n'
    printf '>lower desc with spaces\nacgtnACGTN\nacg\n'
    printf '>iupac\nRYKMSWBDHVN\n'
    printf '>one\nA\n'
    printf '>nofinal\nACGTACGT'
} >edge.fa

# FASTQ edge cases: empty reads, '@' and '+' starting quality lines, comments, tags.
{
    printf '@r0\n\n+\n\n'
    printf '@r1 comment here\nACGT\n+\n@@@@\n'
    printf '@r2\tBC:Z:ACGT\tRG:Z:x\nACGTN\n+r2\n+++++\n'
    printf '@r3/1\nA\n+\nI\n'
    printf '@r4\nNNNNNNNNNN\n+\n!!!!!!!!!!\n'
} >edge.fq

cases=""
for sam in *.sam; do
    name=${sam%.sam}
    ref=ref.fa
    [ "$name" = manysq ] && ref=manysq.fa
    if ! "$st" sort --no-PG -o "$name.bam" "$sam" 2>"$name.convert.err"; then
        rm -f "$name.bam"
        continue
    fi
    "$st" index "$name.bam"
    "$st" index -c "$name.bam"
    "$st" sort --no-PG -n -o "$name.name.bam" "$sam"
    ok=1
    for v in 3.0 3.1 2.1; do
        if "$st" view --no-PG -T "$ref" -O "cram,version=$v" -o "$name.$v.cram" "$name.bam" 2>>"$name.convert.err" &&
            "$st" index "$name.$v.cram" 2>>"$name.convert.err"; then
            :
        else
            ok=0
        fi
    done
    if "$st" view --no-PG -T "$ref" -O cram,embed_ref=1 -o "$name.embed.cram" "$name.bam" 2>>"$name.convert.err" &&
        "$st" index "$name.embed.cram" 2>>"$name.convert.err"; then
        :
    else
        ok=0
    fi
    [ "$ok" = 1 ] || echo "warning: some CRAM conversions failed for $name" >&2
    cases="$cases $name"
done
printf '%s\n' $cases >cases
