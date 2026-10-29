#!/usr/bin/env bash
# Compares ntools output with samtools, byte for byte, on a benchmark corpus.
#
# Usage: compare-samtools.sh NTOOLS SAMTOOLS CORPUS_DIR [WORK_DIR]
#
# CORPUS_DIR holds short/, tagged/, and long/ inputs (input.bam with .bai/.csi, input.sam,
# input.{3.0,3.1,3.1.multislice,3.1.noref}.cram with .crai, queries.*.regions/.bed, and
# reference.path), plus synthetic-reference.fa and fragmented-reference.fa with .fai indexes,
# as produced by the noodles benchmark preparation scripts.
#
# Prints one OK or DIFF line per check and exits nonzero if any check differs.
set -uo pipefail

nt=$1
st=$2
corpus=$3
work=${4:-$(mktemp -d)}

mkdir -p "$work"
cd "$work" || exit 2

failures=0

check() {
    local name=$1 expected=$2 actual=$3

    if cmp -s <(eval "$expected" 2>/dev/null) <(eval "$actual" 2>/dev/null); then
        echo "OK   $name"
    else
        echo "DIFF $name"
        failures=$((failures + 1))
    fi
}

ref=$corpus/synthetic-reference.fa
long_ref=$(cat "$corpus/long/reference.path")
short=$corpus/short/input.bam
tagged=$corpus/tagged/input.bam
long=$corpus/long/input.bam

head -40 "$corpus/short/queries.100.bed" > q.bed
regions=$(head -20 "$corpus/short/queries.100.regions" | tr '\n' ' ')
"$st" sort -n --no-PG -o short.name.bam "$short"

# view: format conversion
for set in short tagged long; do
    check "view $set BAM->SAM" "$st view --no-PG -h $corpus/$set/input.bam" "$nt view --no-PG -h $corpus/$set/input.bam"
    check "view $set BAM->SAM -@4" "$st view --no-PG -h $corpus/$set/input.bam" "$nt view --no-PG -h -@ 4 $corpus/$set/input.bam"
done

for set in short tagged; do
    for c in 3.0 3.1 3.1.multislice 3.1.noref; do
        cram=$corpus/$set/input.$c.cram
        check "view $set CRAM $c" "$st view --no-PG -T $ref $cram" "$nt view --no-PG -T $ref $cram"
        check "view $set CRAM $c -@4" "$st view --no-PG -T $ref $cram" "$nt view --no-PG -@ 4 -T $ref $cram"
    done
done

check "view long CRAM -@8" "$st view --no-PG -T $long_ref $corpus/long/input.3.1.cram" "$nt view --no-PG -@ 8 -T $long_ref $corpus/long/input.3.1.cram"

# view: filters, counts, and regions
for f in "-c" "-c -f 2" "-c -F 0x904" "-c -q 30" "-c -G 12" "-c -F UNMAP,DUP -q 1"; do
    check "view $f" "$st view $f $short" "$nt view $f $short"
done

check "view regions" "$st view --no-PG $short $regions" "$nt view --no-PG $short $regions"
check "view CRAM regions" "$st view --no-PG -T $ref $corpus/short/input.3.1.cram $regions" "$nt view --no-PG -T $ref $corpus/short/input.3.1.cram $regions"
check "view unmapped region" "$st view --no-PG $short '*'" "$nt view --no-PG $short '*'"
check "view -L" "$st view --no-PG -L q.bed $short" "$nt view --no-PG -L q.bed $short"
check "view -M" "$st view --no-PG -M $short chrS:1000-5000 chrS:3000-9000" "$nt view --no-PG -M $short chrS:1000-5000 chrS:3000-9000"

for s in 0.1 1.25 42.3; do
    check "view -s $s" "$st view --no-PG -s $s $short" "$nt view --no-PG -s $s $short"
done

check "view -d NH:2" "$st view --no-PG -d NH:2 $tagged" "$nt view --no-PG -d NH:2 $tagged"
check "view -x" "$st view --no-PG -x NM -x MD,AS $short" "$nt view --no-PG -x NM -x MD,AS $short"
check "view BAM->BAM" "$st view --no-PG -h $short" "$nt view --no-PG -b -@ 4 $short | $st view --no-PG -h -"

# statistics
for set in short tagged long; do
    for fmt in default tsv json; do
        check "flagstat $set $fmt" "$st flagstat -O $fmt $corpus/$set/input.bam" "$nt flagstat -O $fmt -@ 4 $corpus/$set/input.bam"
    done
    check "idxstats $set" "$st idxstats $corpus/$set/input.bam" "$nt idxstats $corpus/$set/input.bam"
done

check "depth" "$st depth $short" "$nt depth $short"
check "depth -a" "$st depth -a $short" "$nt depth -a $short"
check "depth -J -q 30 -Q 20" "$st depth -J -q 30 -Q 20 $short" "$nt depth -J -q 30 -Q 20 $short"
check "depth long" "$st depth $long" "$nt depth $long"
check "depth two inputs" "$st depth -H $short $tagged" "$nt depth -H $short $tagged"
check "coverage" "$st coverage $short" "$nt coverage $short"
check "coverage long" "$st coverage $long" "$nt coverage $long"

# sort and merge
check "sort" "$st sort --no-PG -o - short.name.bam | $st view --no-PG -h -" "$nt sort --no-PG -m 1M -@ 3 -o - short.name.bam | $st view --no-PG -h -"
check "sort -n" "$st sort -n --no-PG -o - $short | $st view --no-PG -h -" "$nt sort -n --no-PG -o - $short | $st view --no-PG -h -"
check "sort -N" "$st sort -N --no-PG -o - $short | $st view --no-PG -h -" "$nt sort -N --no-PG -m 1M -o - $short | $st view --no-PG -h -"

"$st" view --no-PG -b -f 64 -o p1.bam "$short"
"$st" view --no-PG -b -F 64 -o p2.bam "$short"
check "merge" "$st merge -c -p --no-PG -f -o - p1.bam p2.bam | $st view --no-PG - | cut -f1-11" "$nt merge --no-PG -f -o /dev/stdout p1.bam p2.bam | $st view --no-PG - | cut -f1-11"

# duplicate marking pipeline
check "fixmate" "$st fixmate --no-PG short.name.bam - | $st view --no-PG -h -" "$nt fixmate --no-PG short.name.bam /dev/stdout | $st view --no-PG -h -"
check "fixmate -m -c -r -@4" "$st fixmate --no-PG -m -c -r short.name.bam - | $st view --no-PG -h -" "$nt fixmate --no-PG -m -c -r -@ 4 short.name.bam /dev/stdout | $st view --no-PG -h -"
"$st" fixmate -m --no-PG short.name.bam - | "$st" sort --no-PG -o fixed.bam -
check "markdup" "$st markdup --no-PG fixed.bam - | $st view --no-PG -h -" "$nt markdup --no-PG fixed.bam /dev/stdout | $st view --no-PG -h -"
check "markdup -c -m s" "$st markdup --no-PG -c -m s fixed.bam - | $st view --no-PG -h -" "$nt markdup --no-PG -c -m s fixed.bam /dev/stdout | $st view --no-PG -h -"
check "collate grouping" "$st view $short | cut -f1 | sort -u | wc -l" "$nt collate --no-PG -O $short | $st view - | awk -F'\\t' '\$1 != p { print \$1; p = \$1 }' | sort | uniq | wc -l"

# reads
check "fastq" "$st fastq short.name.bam" "$nt fastq short.name.bam"
check "fastq -n -t" "$st fastq -n -t $tagged" "$nt fastq -n -t $tagged"
check "fasta" "$st fasta short.name.bam" "$nt fasta short.name.bam"

# indexes and FASTA
cp "$corpus/fragmented-reference.fa" frag.fa
"$st" faidx frag.fa && mv frag.fa.fai expected.fai
"$nt" faidx frag.fa
check "faidx index" "cat expected.fai" "cat frag.fa.fai"
check "faidx regions" "$st faidx -r $corpus/fasta-regions.1024.txt frag.fa" "$nt faidx -r $corpus/fasta-regions.1024.txt frag.fa"
"$nt" index -o ntools.bai "$short"
check "index BAI queries" "$st view --no-PG -X $short $corpus/short/input.bam.bai $regions" "$st view --no-PG -X $short ntools.bai $regions"

# split
mkdir -p split-s split-n
"$st" view -h "$short" chrS:1-200000 | awk 'BEGIN{OFS="\t"} /^@RG/ {print "@RG\tID:A\tSM:a"; print "@RG\tID:B\tSM:b"; next} /^@/ {print; next} {sub(/\tRG:Z:synthetic/, ""); n++; print $0 (n%3==0 ? "" : (n%3==1 ? "\tRG:Z:A" : "\tRG:Z:B"))}' | "$st" view -b -o multi.bam -
(cd split-s && "$st" split --no-PG -u u.bam ../multi.bam)
(cd split-n && "$nt" split --no-PG -u u.bam ../multi.bam)
for f in multi_0.bam multi_1.bam u.bam; do
    check "split $f" "$st view --no-PG -h split-s/$f" "$st view --no-PG -h split-n/$f"
done

# dict, bedcov, addreplacerg
check "dict" "$st dict -a A -s S -u U frag.fa" "$nt dict -a A -s S -u U frag.fa"
check "bedcov -c -d 5" "$st bedcov -c -d 5 q.bed $short" "$nt bedcov -c -d 5 q.bed $short"
check "bedcov long -j" "$st bedcov -j $(head -30 "$corpus/long/queries.100.bed" > lq.bed; echo lq.bed) $long" "$nt bedcov -j lq.bed $long"
check "addreplacerg records" "$st addreplacerg --no-PG -r 'ID:z\tSM:z' multi.bam | grep -v '^@'" "$nt addreplacerg --no-PG -r 'ID:z\tSM:z' multi.bam | grep -v '^@'"

echo "$failures difference(s)"
[ "$failures" -eq 0 ]
