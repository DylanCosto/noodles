#!/usr/bin/env bash
# Times ntools against samtools on a benchmark corpus and prints a Markdown table.
#
# Usage: benchmark.sh NTOOLS SAMTOOLS CORPUS_DIR [RUNS] [WORK_DIR]
#
# Each command runs RUNS times (default 5) with warm caches; the table reports the median wall
# time. Output goes to /dev/null or WORK_DIR, so disk writes are included for file outputs.
set -euo pipefail

nt=$1
st=$2
corpus=$3
runs=${4:-5}
work=${5:-$(mktemp -d)}

mkdir -p "$work"

ref=$corpus/synthetic-reference.fa
long_ref=$(cat "$corpus/long/reference.path")
short=$corpus/short/input.bam
tagged=$corpus/tagged/input.bam
long=$corpus/long/input.bam
long_cram=$corpus/long/input.3.1.cram
short_cram=$corpus/short/input.3.1.cram

"$st" sort -n --no-PG -o "$work/short.name.bam" "$short"
"$st" fixmate -m --no-PG "$work/short.name.bam" - | "$st" sort --no-PG -o "$work/fixed.bam" -

median() {
    local cmd=$1 times=()

    for _ in $(seq "$runs"); do
        local start end
        start=$(date +%s.%N)
        bash -c "$cmd" >/dev/null 2>&1
        end=$(date +%s.%N)
        times+=("$(echo "$end - $start" | bc)")
    done

    printf '%s\n' "${times[@]}" | sort -g | awk '{a[NR]=$1} END {printf "%.3f", a[int((NR+1)/2)]}'
}

row() {
    local name=$1 args_st=$2 args_nt=$3 s n
    s=$(median "$st $args_st")
    n=$(median "$nt $args_nt")
    printf '| %s | %s | %s | %.2fx |\n' "$name" "$s" "$n" "$(echo "$s / $n" | bc -l)"
}

echo "| Operation | samtools (s) | ntools (s) | Speedup |"
echo "|---|---:|---:|---:|"

row "view BAM->SAM, short" "view --no-PG $short" "view --no-PG $short"
row "view BAM->SAM, long" "view --no-PG $long" "view --no-PG $long"
row "view BAM->SAM, long, -@4" "view --no-PG -@ 4 $long" "view --no-PG -@ 4 $long"
row "view BAM->SAM, long, -@8" "view --no-PG -@ 8 $long" "view --no-PG -@ 8 $long"
row "view CRAM->SAM, short" "view --no-PG -T $ref $short_cram" "view --no-PG -T $ref $short_cram"
row "view CRAM->SAM, long" "view --no-PG -T $long_ref $long_cram" "view --no-PG -T $long_ref $long_cram"
row "view CRAM->SAM, long, -@8" "view --no-PG -@ 8 -T $long_ref $long_cram" "view --no-PG -@ 8 -T $long_ref $long_cram"
row "view CRAM->BAM, long, -@8" "view --no-PG -@ 8 -b -o $work/x.bam -T $long_ref $long_cram" "view --no-PG -@ 8 -b -o $work/x.bam -T $long_ref $long_cram"
row "view BAM->BAM, long, -@4" "view --no-PG -@ 4 -b -o $work/x.bam $long" "view --no-PG -@ 4 -b -o $work/x.bam $long"
row "view -c -q 30, short" "view -c -q 30 $short" "view -c -q 30 $short"
row "flagstat, long, -@4" "flagstat -@ 4 $long" "flagstat -@ 4 $long"
row "index, long, -@4" "index -@ 4 -o $work/x.bai $long" "index -@ 4 -o $work/x.bai $long"
row "sort, short, -@4" "sort --no-PG -@ 4 -o $work/x.bam $work/short.name.bam" "sort --no-PG -@ 4 -o $work/x.bam $work/short.name.bam"
row "sort -n, short, -@4" "sort -n --no-PG -@ 4 -o $work/x.bam $short" "sort -n --no-PG -@ 4 -o $work/x.bam $short"
row "collate -Ou, short, -@4" "collate -Ou -@4 $short" "collate -Ou -@ 4 $short"
row "fixmate -m, short" "fixmate -m $work/short.name.bam $work/x.bam" "fixmate -m $work/short.name.bam $work/x.bam"
row "fixmate -m, short, -@4" "fixmate -m -@4 $work/short.name.bam $work/x.bam" "fixmate -m -@ 4 $work/short.name.bam $work/x.bam"
row "markdup, short" "markdup $work/fixed.bam $work/x.bam" "markdup $work/fixed.bam $work/x.bam"
row "markdup, short, -@4" "markdup -@4 $work/fixed.bam $work/x.bam" "markdup -@ 4 $work/fixed.bam $work/x.bam"
row "fastq, short" "fastq $short" "fastq $short"
row "depth, short" "depth $short" "depth $short"
row "stats, short" "stats $short" "stats $short"
row "stats, short, -@4" "stats -@ 4 $short" "stats -@ 4 $short"
row "stats, long, -@8" "stats -@ 8 $long" "stats -@ 8 $long"
row "mpileup, short" "mpileup $short" "mpileup $short"
row "mpileup -f (BAQ), short, -@8" "mpileup -f $ref $short" "mpileup -@ 8 -f $ref $short"
row "coverage, short" "coverage $short" "coverage $short"
row "faidx 1024 regions" "faidx -r $corpus/fasta-regions.1024.txt $corpus/fragmented-reference.fa" "faidx -r $corpus/fasta-regions.1024.txt $corpus/fragmented-reference.fa"
