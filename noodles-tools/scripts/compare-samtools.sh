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
# Prints one OK or DIFF line per check, then the slowest checks, and exits nonzero if any check
# differs. Per-check wall times are written to WORK_DIR/timings.tsv.
#
# samtools outputs are cached as checksums in $COMPARE_CACHE (default
# ${XDG_CACHE_HOME:-~/.cache}/ntools-compare-samtools), keyed on the samtools build, the corpus files, the command,
# and the sizes of the work files it names, so later runs only run ntools. Set COMPARE_CACHE=
# to run samtools every time.
set -uo pipefail

nt=$1
st=$2
corpus=$3
work=${4:-$(mktemp -d)}
scripts=$(cd "$(dirname "$0")" && pwd)

mkdir -p "$work"
cd "$work" || exit 2

failures=0

cache=${COMPARE_CACHE-${XDG_CACHE_HOME:-$HOME/.cache}/ntools-compare-samtools}
if [ -n "$cache" ]; then
    mkdir -p "$cache"
    cache_id=$({
        readlink -f "$(command -v "$st")"
        "$st" --version | head -2
        find "$corpus" -type f -printf '%P %s %T@\n' | sort
    } | sha256sum | cut -c1-64)
fi

timings=$work/timings.tsv
: >"$timings"

# Prints the cache key of a samtools command: the command and the sizes of the work files it
# names. Files named st.* or *.st.* are samtools outputs, absent when a cached result is used.
cache_key() {
    local - word
    set -f
    {
        echo "$cache_id"
        printf '%s\n' "$1"
        for word in $1; do
            case $word in
                st.* | *.st.*) ;;
                *) [ -f "$word" ] && printf '%s %s\n' "$word" "$(stat -c %s "$word")" ;;
            esac
        done
    } | sha256sum | cut -c1-64
}

check() {
    local name=$1 expected=$2 actual=$3
    local start=$EPOCHREALTIME key="" source=run expected_sum actual_sum pid

    # Expectations that run ntools are never cached.
    if [ -n "$cache" ] && [[ $expected != *"$nt"* ]]; then
        key=$(cache_key "$expected")
    fi

    if [ -n "$key" ] && [ -f "$cache/$key" ]; then
        source=cache
        expected_sum=$(<"$cache/$key")
        actual_sum=$(eval "$actual" 2>/dev/null | sha256sum)
    else
        (eval "$expected" 2>/dev/null | sha256sum >.expected.sum) &
        pid=$!
        actual_sum=$(eval "$actual" 2>/dev/null | sha256sum)
        wait "$pid"
        expected_sum=$(<.expected.sum)
        [ -n "$key" ] && echo "$expected_sum" >"$cache/$key"
    fi

    printf '%s\t%s\t%s\n' "$(awk -v a="$start" -v b="$EPOCHREALTIME" 'BEGIN { printf "%.3f", b - a }')" "$source" "$name" >>"$timings"

    if [ "$expected_sum" = "$actual_sum" ]; then
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

# Read groups a and c are in library x; d has no library; one read has an integer RG.
printf '@HD\tVN:1.6\n@SQ\tSN:c\tLN:100\n@RG\tID:a\tLB:x\n@RG\tID:b\tLB:y\n@RG\tID:c\tLB:x\n@RG\tID:d\n' > lib.sam
for rg in RG:Z:a RG:Z:b RG:Z:c RG:Z:d RG:Z:zz RG:i:5 XX:i:0; do
    printf 'r%s\t0\tc\t1\t60\t4M\t*\t0\t0\tACGT\tIIII\t%s\n' "${rg//:/}" "$rg" >> lib.sam
done
for opts in "-l x" "-l y" "-l q" "-c -l x"; do
    check "view $opts" "$st view --no-PG $opts lib.sam" "$nt view --no-PG $opts lib.sam"
done
check "view -l -U" "$st view --no-PG -l x -U /dev/stdout -o /dev/null lib.sam" "$nt view --no-PG -l x -U /dev/stdout -o /dev/null lib.sam"
check "view --output-fmt-option level=1" "$st view --no-PG -b --output-fmt-option level=1 $short | $st view --no-PG -" "$nt view --no-PG -b --output-fmt-option level=1 $short | $st view --no-PG -"
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
# markdup -S/-t/-d on synthetic Illumina-style names (the corpus names have no coordinates)
if command -v python3 >/dev/null; then
    python3 "$scripts/markdup-test-data.py" 7 > md.sam
    "$st" sort -n -o md.name.bam md.sam
    "$st" fixmate -m md.name.bam - | "$st" sort -o md.bam -

    for opts in "-S -t -d 100" "-t -d 100 -m s" "-S -d 2500 -r --duplicate-count" "-S -t -d 100 --use-read-groups --include-fails" "-t --no-multi-dup -d 100"; do
        check "markdup $opts" \
            "$st markdup --no-PG -f st.md.stats $opts md.bam - | $st view --no-PG -h -; tail -n +2 st.md.stats" \
            "$nt markdup --no-PG -f nt.md.stats $opts md.bam - | $st view --no-PG -h -; tail -n +2 nt.md.stats"
    done

    "$st" markdup --no-PG -S -t -d 100 md.bam md.marked.bam 2>/dev/null
    check "markdup -c -S -t -d 50 (remarking)" "$st markdup --no-PG -c -S -t -d 50 md.marked.bam - | $st view --no-PG -h -" "$nt markdup --no-PG -c -S -t -d 50 md.marked.bam - | $st view --no-PG -h -"
    check "markdup --json" "$st markdup --no-PG --json -S -d 100 --use-read-groups md.bam md.st.bam 2>&1 | grep '^ ' | grep -v COMMAND" "$nt markdup --no-PG --json -S -d 100 --use-read-groups md.bam md.nt.bam 2>&1 | grep '^ ' | grep -v COMMAND"

    # Barcodes from tags and names, UMIs moved to RX, and coordinates from a regex
    python3 "$scripts/markdup-test-data.py" 11 --barcodes > mdbc.sam
    "$st" sort -n -o mdbc.name.bam mdbc.sam
    "$st" fixmate -m mdbc.name.bam - | "$st" sort -o mdbc.bam -
    "$st" view -h mdbc.bam | grep -v odd_name | "$st" view -b -o mdbc.named.bam -

    # Options are eval'ed, so regexes are quoted inside them.
    for opts in "--barcode-tag BC -S -t -d 100" "--barcode-name -t -d 100 --use-read-groups" "--barcode-name --move-umi-to-tag -S -t" "--barcode-rgx ':([ACGT]+)\$' --move-umi-to-tag -S"; do
        check "markdup $opts" \
            "$st markdup --no-PG -f st.md.stats $opts mdbc.bam - | $st view --no-PG -h -; tail -n +2 st.md.stats" \
            "$nt markdup --no-PG -f nt.md.stats $opts mdbc.bam - | $st view --no-PG -h -; tail -n +2 nt.md.stats"
    done

    # Names without coordinates are removed: samtools fails on kept reads the regex does not match.
    for opts in "--read-coords '([0-9]+):([0-9]+):([0-9]+)(:[ACGT]+)?\$' -d 100 -t" "--read-coords '^([^:]+:[^:]+):.*:([0-9]+):([0-9]+)(:[ACGT]+)?\$' --coords-order tyx -d 2500 -t -S --barcode-name"; do
        check "markdup $opts" \
            "$st markdup --no-PG -f st.md.stats $opts mdbc.named.bam - | $st view --no-PG -h -; tail -n +2 st.md.stats" \
            "$nt markdup --no-PG -f nt.md.stats $opts mdbc.named.bam - | $st view --no-PG -h -; tail -n +2 nt.md.stats"
    done
fi
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
check "faidx -i -n 25 --mark-strand sign -c" "$st faidx -i -n 25 --mark-strand sign -c frag.fa nope $(head -1 frag.fa.fai | cut -f1):5-300" "$nt faidx -i -n 25 --mark-strand sign -c frag.fa nope $(head -1 frag.fa.fai | cut -f1):5-300"
bgzip=$(dirname "$st")/bgzip
if [ -x "$bgzip" ]; then
    "$bgzip" -c frag.fa > frag.fa.gz
    "$st" faidx frag.fa.gz && mv frag.fa.gz.fai expected.gz.fai && mv frag.fa.gz.gzi expected.gzi
    "$nt" faidx frag.fa.gz
    check "faidx index (bgzip)" "cat expected.gz.fai expected.gzi" "cat frag.fa.gz.fai frag.fa.gz.gzi"
    check "faidx regions (bgzip)" "$st faidx -r $corpus/fasta-regions.1024.txt frag.fa.gz" "$nt faidx -r $corpus/fasta-regions.1024.txt frag.fa.gz"
fi

# fqidx, on FASTQ written by samtools fastq (single-line) and a multi-line FASTQ
"$st" fastq "$tagged" 2>/dev/null | head -40000 > fq.fq
awk 'NR % 4 == 1 { print; next } NR % 4 == 3 { print "+"; next } { while (length($0) > 60) { print substr($0, 1, 60); $0 = substr($0, 61) } print }' fq.fq > fqml.fq
for f in fq.fq fqml.fq; do
    "$st" fqidx $f && mv $f.fai expected.fq.fai
    "$nt" fqidx $f
    check "fqidx index $f" "cat expected.fq.fai" "cat $f.fai"
    names=$(awk 'NR % 997 == 1 { printf "%s ", $1 } NR % 1999 == 1 { printf "%s:5-90 ", $1 }' $f.fai)
    check "fqidx regions $f -i" "$st fqidx -i $f $names" "$nt fqidx -i $f $names"
done
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

# stats, without the version and command lines
strip_stats="grep -v '^# This file was produced\|^# The command line'"
awk 'BEGIN { OFS = "\t" } { print $1, $2 + 1, $3 }' "$corpus/short/queries.100.bed" | sort -k1,1 -k2,2n > targets.txt
check "stats" "$st stats $short | $strip_stats" "$nt stats $short | $strip_stats"
check "stats -r -@4" "$st stats -r $ref $short | $strip_stats" "$nt stats -@ 4 -r $ref $short | $strip_stats"
check "stats long -r -@8" "$st stats -r $long_ref $long | $strip_stats" "$nt stats -@ 8 -r $long_ref $long | $strip_stats"
check "stats tagged" "$st stats $tagged | $strip_stats" "$nt stats $tagged | $strip_stats"
check "stats CRAM" "$st stats --reference $ref $corpus/short/input.3.1.cram | $strip_stats" "$nt stats --reference $ref $corpus/short/input.3.1.cram | $strip_stats"
check "stats unsorted" "$st stats short.name.bam | $strip_stats" "$nt stats short.name.bam | $strip_stats"
check "stats -t -g 3 -p" "$st stats -t targets.txt -g 3 -p $short | $strip_stats" "$nt stats -t targets.txt -g 3 -p $short | $strip_stats"
check "stats regions --ref-stats" "$st stats --ref-stats -r $ref $short $regions | $strip_stats" "$nt stats --ref-stats -r $ref $short $regions | $strip_stats"
check "stats -d -q 20 -i 500 -c 5,200,5 -x" "$st stats -d -q 20 -i 500 -c 5,200,5 -x $short | $strip_stats" "$nt stats -d -q 20 -i 500 -c 5,200,5 -x $short | $strip_stats"

# Barcode tags with separators, length mismatches, and invalid characters
"$st" view -h "$short" | head -40000 | awk 'BEGIN { OFS = "\t"; srand(3); b = "ACGTN" } /^@/ { print; next } {
    bc = ""; q = ""
    for (i = 0; i < 8; i++) { bc = bc substr(b, int(rand() * 5) + 1, 1); q = q sprintf("%c", 35 + int(rand() * 40)) }
    r = rand()
    if (r < 0.5) { $0 = $0 "\tBC:Z:" bc "-" substr(bc, 1, 6) "\tQT:Z:" q "!" substr(q, 1, 6) }
    else if (r < 0.6) { $0 = $0 "\tBC:Z:" substr(bc, 1, 5) "\tQT:Z:" substr(q, 1, 4) }
    else if (r < 0.65) { $0 = $0 "\tBC:Z:" bc "+X" }
    if (rand() < 0.4) { $0 = $0 "\tRX:Z:" bc "\tQX:Z:" q }
    if (rand() < 0.2) { $0 = $0 "\tCR:Z:" substr(bc, 1, 4) }
    print }' > barcodes.sam
check "stats barcodes" "$st stats barcodes.sam 2>/dev/null | $strip_stats" "$nt stats barcodes.sam 2>/dev/null | $strip_stats"

if [ -f md.marked.bam ]; then
    check "stats -I a (marked)" "$st stats -I a md.marked.bam | $strip_stats" "$nt stats -I a md.marked.bam | $strip_stats"
fi

# mpileup (BAQ with -f), on part of the short reads, long reads, and synthetic overlapping pairs
"$st" view -b -o mp.short.bam "$short" chrS:1-300000 && "$st" index mp.short.bam
long_name=$("$st" idxstats "$long" | head -1 | cut -f1)
"$st" view -b -o mp.long.bam "$long" "$long_name:1-300000" && "$st" index mp.long.bam
check "mpileup" "$st mpileup mp.short.bam 2>/dev/null" "$nt mpileup mp.short.bam 2>/dev/null"
check "mpileup -f (BAQ) -@4" "$st mpileup -f $ref mp.short.bam 2>/dev/null" "$nt mpileup -@ 4 -f $ref mp.short.bam 2>/dev/null"
# Long-read BAQ takes ~2 minutes single-threaded; single-threaded BAQ is covered on short reads.
check "mpileup long -f -@8" "$st mpileup -f $long_ref mp.long.bam 2>/dev/null" "$nt mpileup -@ 8 -f $long_ref mp.long.bam 2>/dev/null"
check "mpileup CRAM -f" "$st mpileup -f $ref -r chrS:1-300000 $corpus/short/input.3.1.cram 2>/dev/null" "$nt mpileup -f $ref -r chrS:1-300000 $corpus/short/input.3.1.cram 2>/dev/null"
check "mpileup -aa -r -l -s -O --output-extra" "$st mpileup -aa -r chrS:1000-60000 -l targets.txt -s -O --output-BP-5 --output-extra QNAME,FLAG,RG,NM mp.short.bam 2>/dev/null" "$nt mpileup -aa -r chrS:1000-60000 -l targets.txt -s -O --output-BP-5 --output-extra QNAME,FLAG,RG,NM mp.short.bam 2>/dev/null"
check "mpileup long ins/del options -@8" "$st mpileup --reverse-del --no-output-ins --no-output-del --no-output-ends -f $long_ref mp.long.bam 2>/dev/null" "$nt mpileup -@ 8 --reverse-del --no-output-ins --no-output-del --no-output-ends -f $long_ref mp.long.bam 2>/dev/null"

if [ -f md.bam ]; then
    check "mpileup overlapping pairs -A -Q 0 -d 3" "$st mpileup -A -Q 0 -d 3 md.bam 2>/dev/null" "$nt mpileup -A -Q 0 -d 3 md.bam 2>/dev/null"
    check "mpileup -C 50 -E -x two files" "$st mpileup -C 50 -E -x -f $ref md.bam md.bam 2>/dev/null" "$nt mpileup -C 50 -E -x -f $ref md.bam md.bam 2>/dev/null"
    # Indexed, with several sequences: one chromosome per thread
    "$st" index md.bam
    check "mpileup by chromosome -aa -@4" "$st mpileup -aa -A md.bam 2>/dev/null" "$nt mpileup -@ 4 -aa -A md.bam 2>/dev/null"
fi

# mpileup base modifications (MM/ML), on long reads and synthetic edge cases
check "mpileup long -M" "$st mpileup -M mp.long.bam 2>/dev/null" "$nt mpileup -M mp.long.bam 2>/dev/null"

if command -v python3 >/dev/null; then
    python3 "$scripts/basemod-test-data.py" mods.fa 7 > mods.sam
    "$st" view -b -o mods.bam mods.sam

    for opts in "-M" "-M -f mods.fa -B -aa" "-M --no-output-ins-mods --reverse-del" "-M --no-output-ins -Q 0"; do
        check "mpileup $opts (synthetic mods)" "$st mpileup $opts mods.bam 2>/dev/null" "$nt mpileup $opts mods.bam 2>/dev/null"
    done
fi

# ampliconclip, with primer pairs every 350 bp (and a few overlapping ones)
awk 'BEGIN { print "# primers"; for (s = 0; s < 300000; s += 350) { printf "chrS\t%d\t%d\tamp%d_LEFT\t60\t+\n", s, s + 22, s; printf "chrS\t%d\t%d\tamp%d_RIGHT\t60\t-\n", s + 328, s + 350, s; if (s % 2450 == 0) printf "chrS\t%d\t%d\talt%d\t1\t+\n", s + 5, s + 30, s } }' > primers.bed
for opts in "" "--hard-clip --original" "--both-ends --strand --keep-tag" "--both-ends --hard-clip --fail --fail-len 145 --filter-len 140" "--clipped --no-excluded --tolerance 0 --unmap-len 130"; do
    # OA tags differ: samtools writes the read name where the reference name belongs.
    check "ampliconclip $opts" \
        "$st ampliconclip --no-PG -b primers.bed $opts -O sam -f st.amp.stats mp.short.bam | sed 's/\tOA:Z:[^\t]*//'; tail -n +2 st.amp.stats" \
        "$nt ampliconclip --no-PG -b primers.bed $opts -O sam -f nt.amp.stats mp.short.bam | sed 's/\tOA:Z:[^\t]*//'; tail -n +2 nt.amp.stats"
done
check "ampliconclip --original (SAM-spec OA)" "echo 'OA:Z:chrS,9,+,150M,34,0;'" "$nt ampliconclip --no-PG -b primers.bed --original -O sam mp.short.bam | grep -o 'OA:Z:[^[:space:]]*' | head -1"
check "ampliconclip --primer-counts" \
    "$st ampliconclip --no-PG -b primers.bed --strand --primer-counts /dev/stdout -o /dev/null mp.short.bam" \
    "$nt ampliconclip --no-PG -b primers.bed --strand --primer-counts /dev/stdout -o /dev/null mp.short.bam"

# CRAMs written with HTSlib's defaults drop MD and NM tags that match; both tools regenerate them
"$st" view -C -T "$ref" -o mp.short.cram mp.short.bam
"$st" view -C -T "$ref" -x MD -x NM -o mp.short.nomd.cram mp.short.bam
"$st" view -C -T "$long_ref" -o mp.long.cram mp.long.bam
check "view CRAM with dropped MD/NM" "$st view -T $ref mp.short.cram" "$nt view -T $ref mp.short.cram"
check "view CRAM without MD/NM" "$st view -T $ref mp.short.nomd.cram" "$nt view -T $ref mp.short.nomd.cram"
check "view long CRAM with dropped MD/NM" "$st view -T $long_ref mp.long.cram" "$nt view -T $long_ref mp.long.cram"
check "consensus CRAM with dropped MD/NM" "$st consensus mp.short.cram" "$nt consensus mp.short.cram"

# import: FASTQ written by samtools fastq (with MM/ML tags, paired files, interleaved names)
"$st" fastq -T MM,ML mp.long.bam > imp.long.fq 2>/dev/null
"$st" collate -u -O mp.short.bam 2>/dev/null | "$st" fastq -1 imp.R1.fq -2 imp.R2.fq -s imp.S.fq -0 /dev/null - 2>/dev/null
"$st" collate -u -O mp.short.bam 2>/dev/null | "$st" fastq -N - > imp.inter.fq 2>/dev/null
printf '@M:7:FC:1:1101:10:20:ACGT 1:Y:0:GGTT\tXA:i:5\nACGTN\n+\nIII#I\n@M:7:FC:1:1101:11:21:A+G 2:N:18:NN\nAC\n+\nII\n' > imp.casava.fq
check "import -T * (long, MM/ML)" "$st import --no-PG -T '*' imp.long.fq" "$nt import --no-PG -T '*' imp.long.fq"
check "import -1 -2 -s -r --order" "$st import --no-PG -1 imp.R1.fq -2 imp.R2.fq -s imp.S.fq -r ID:a -r SM:b --order XO:5" "$nt import --no-PG -1 imp.R1.fq -2 imp.R2.fq -s imp.S.fq -r ID:a -r SM:b --order XO:5"
check "import interleaved -O bam" "$st import --no-PG -O bam imp.inter.fq | $st view -h -" "$nt import --no-PG -O bam imp.inter.fq | $st view -h -"
check "import --i1 -1 -2 -b" "$st import --no-PG --i1 imp.R2.fq -1 imp.R1.fq -2 imp.R2.fq -b" "$nt import --no-PG --i1 imp.R2.fq -1 imp.R1.fq -2 imp.R2.fq -b"
# samtools loses tags after a CASAVA barcode; ntools keeps them.
check "import -i -U -T *" "$st import --no-PG -i -U -T '*' imp.casava.fq | sed 's/\tXA:i:5//'" "$nt import --no-PG -i -U -T '*' imp.casava.fq | sed 's/\tXA:i:5//'"
check "import -i -T *: tags after a barcode" "echo 1" "$nt import --no-PG -i -T '*' imp.casava.fq | grep -c 'XA:i:5'"

# reset and checksum
check "reset" "$st reset --no-PG $tagged" "$nt reset --no-PG $tagged"
check "reset CRAM" "$st reset --no-PG $corpus/short/input.3.1.cram" "$nt reset --no-PG $corpus/short/input.3.1.cram"
check "checksum" "$st checksum $tagged" "$nt checksum $tagged"
check "checksum -a long" "$st checksum -a $long" "$nt checksum -a $long"

if [ -f md.bam ]; then
    check "reset --dupflag -x ^RG,ms --reject-PG" "$st reset --no-PG --dupflag -x '^RG,ms' --reject-PG samtools.1 md.bam" "$nt reset --no-PG --dupflag -x '^RG,ms' --reject-PG samtools.1 md.bam"
    check "checksum -q -v -T -P -C -M -OO" "$st checksum -q -v -T -P -C -M -OO md.bam" "$nt checksum -q -v -T -P -C -M -OO md.bam"
    "$st" split -f 'chk_%!.bam' md.bam
    for f in chk_*.bam; do "$st" checksum -q -o "$f.chk" "$f"; done
    check "checksum -m (merge)" "$st checksum -q -m chk_*.bam.chk" "$nt checksum -q -m chk_*.bam.chk"
fi

# head, flags, samples
check "head" "$st head $tagged" "$nt head $tagged"
check "head -n 5 -h 2" "$st head -n 5 -h 2 $long" "$nt head -n 5 -h 2 $long"
check "head CRAM (reference from @SQ UR)" "$st head -h 0 -n 1000 $corpus/short/input.3.1.cram" "$nt head -h 0 -n 1000 $corpus/short/input.3.1.cram"
check "flags" "$st flags 99 0x10 010 PAIRED,read1 12abc 70000 ''" "$nt flags 99 0x10 010 PAIRED,read1 12abc 70000 ''"
# Many read groups, so samples come out in khash order after several resizes
awk 'BEGIN { print "@HD\tVN:1.6"; print "@SQ\tSN:c\tLN:100"; for (i = 0; i < 600; i++) printf "@RG\tID:r%d\tSM:s%d\tLB:l%d\n", i, (i * 7919) % 451, i % 3 }' > rg.sam
"$st" view -b -o rg.bam rg.sam
check "samples" "$st samples rg.bam $short $tagged" "$nt samples rg.bam $short $tagged"
check "samples -h -i -T LB -f" "$st samples -h -i -T LB -f $ref rg.bam $short $corpus/short/input.3.1.cram" "$nt samples -h -i -T LB -f $ref rg.bam $short $corpus/short/input.3.1.cram"

# consensus: every mode and -X preset, output formats, -a/-aa, insertions and deletions, regions,
# and the threaded (-@) block layout, which samtools computes differently from the serial path
# samtools' pileup output with -a repeats depth-0 lines after positions whose consensus is a
# deletion (and gives such positions a depth-0 line); ntools fixes this, so with -a only covered
# columns are compared with samtools, and -a output is checked to be the same with threads.
covered="awk -F'\t' '\$4 != 0'"
compare_consensus() {
    local name=$1 opts=$2 input=$3
    case "$opts" in
        *pileup*-a*|*-a*pileup*)
            check "$name" "$st consensus $opts $input 2>/dev/null | $covered" "$nt consensus $opts $input 2>/dev/null | $covered" ;;
        *)
            check "$name" "$st consensus $opts $input 2>/dev/null" "$nt consensus $opts $input 2>/dev/null" ;;
    esac
}

for set in short long; do
    for opts in "" "-f fastq" "-f pileup" "-a -f fastq --show-del yes --show-ins no" "-A --mark-ins" \
        "-m simple" "-m simple -q -A -d 3 -c 0.6 -H 0.3" "-m bayesian_m -A" "-m bayesian_p" "-m bayesian_116" \
        "-X hiseq" "-X hifi -A" "-X r10.4_sup" "-X r10.4_dup -f fastq" "-X ultima" \
        "-p --no-adj-MQ --min-MQ 10 --min-BQ 5" "--no-use-MQ -C 20" "-t :hifi --P-het 1e-2 --het-scale 0.5" \
        "-@4" "-@4 -f pileup -a -Z 20000" "-@4 -f fastq --mark-ins -Z 777"; do
        compare_consensus "consensus $set $opts" "$opts" "mp.$set.bam"
    done
done
check "consensus -f pileup -aa: threads give the serial output" "$nt consensus -f pileup -aa mp.long.bam" "$nt consensus -@4 -Z 777 -f pileup -aa mp.long.bam"
check "consensus -f pileup -aa: one line per column" "echo 0" "$nt consensus -f pileup -aa mp.long.bam | cut -f1-3 | sort | uniq -d | wc -l"
check "consensus long -r -a -T" "$st consensus -r $long_name:20000-60000 -a -f fastq -T $long_ref --ref-qual 5 mp.long.bam" "$nt consensus -r $long_name:20000-60000 -a -f fastq -T $long_ref --ref-qual 5 mp.long.bam"
compare_consensus "consensus -r -a -f pileup" "-r chrS:1000-60000 -a -f pileup" mp.short.bam
compare_consensus "consensus -aa -f pileup -T" "-aa -f pileup -T $ref" mp.short.bam
printf "chrS\t5000\t9000\nchrS\t100\t2000\nchrS\t1500\t2500\n" > cons.bed
check "consensus --regions-file" "$st consensus --regions-file cons.bed -a mp.short.bam" "$nt consensus --regions-file cons.bed -a mp.short.bam"
check "consensus --regions-file -@2" "$st consensus -@2 -Z 700 --regions-file cons.bed -f pileup mp.short.bam" "$nt consensus -@2 -Z 700 --regions-file cons.bed -f pileup mp.short.bam"
check "consensus CRAM -r" "$st consensus -r chrS:1-300000 -f fastq $corpus/short/input.3.1.cram" "$nt consensus -r chrS:1-300000 -f fastq $corpus/short/input.3.1.cram"

# consensus edge cases: insertions at read starts, pads, deletions, reference skips, clips,
# heterozygous sites, N bases, low depth, and references with no reads
{
    printf "@HD\tVN:1.6\tSO:coordinate\n@SQ\tSN:a\tLN:60\n@SQ\tSN:b\tLN:30\n@SQ\tSN:c\tLN:50\n"
    printf "r1\t0\ta\t3\t60\t2I10M2D5M\t*\t0\t0\tGGACGTACGTACACGTA\tIIIIIIIIIIIIIIIII\tMD:Z:10^TT5\n"
    printf "r2\t16\ta\t3\t30\t3S6M1P2I4M3N6M\t*\t0\t0\tTTTACGTACTTGTACACGTAC\t#####55555??????IIIII\n"
    printf "r3\t0\ta\t5\t60\t8M1D9M2H\t*\t0\t0\tGTACGTTCACGTACNNA\t;;;;;;;;;;;;;;;;;\tMD:Z:6A1^G9\n"
    printf "r4\t1024\ta\t6\t60\t10M\t*\t0\t0\tTACGTACGTA\tIIIIIIIIII\n"
    printf "r5\t0\ta\t7\t0\t1S12M4I\t*\t0\t0\tAACGTACGTACGTCCCC\t*\n"
    printf "r6\t16\ta\t9\t60\t5=1X4=2D3M\t*\t0\t0\tGTACGCACGTGTA\t+++++++++++++\n"
    printf "r7\t0\tc\t10\t60\t15M\t*\t0\t0\tACGTTGCAACGTTGC\tIIIII#####IIIII\n"
    printf "r8\t0\tc\t12\t60\t4M1I8M\t*\t0\t0\tGTTGTCAACGTTG\tIIIIIIIIIIIII\n"
} > cons.sam
"$st" view -b -o cons.bam cons.sam && "$st" index cons.bam
for opts in "" "-f pileup" "-aa -f fastq --show-del yes" "-aa -f pileup" "-A --mark-ins -d 1 -C 0" "-m simple -A" "-m simple -q --show-del yes -f pileup" "--no-use-MQ -A" "-@2 -aa -Z 7 -f fastq"; do
    compare_consensus "consensus synthetic $opts" "$opts" cons.bam
done

# calmd
check "calmd" "$st calmd --no-PG mp.short.bam $ref 2>/dev/null" "$nt calmd --no-PG mp.short.bam $ref 2>/dev/null"
check "calmd -e -r -A -E -@4" "$st calmd --no-PG -e -r -A -E mp.short.bam $ref 2>/dev/null" "$nt calmd --no-PG -@ 4 -e -r -A -E mp.short.bam $ref 2>/dev/null"
check "calmd -r -n 1 -d -q -C 50" "$st calmd --no-PG -r -n 1 -d -q -C 50 mp.short.bam $ref 2>/dev/null" "$nt calmd --no-PG -r -n 1 -d -q -C 50 mp.short.bam $ref 2>/dev/null"
check "calmd long -e -r -A -@8" "$st calmd --no-PG -e -r -A mp.long.bam $long_ref 2>/dev/null" "$nt calmd --no-PG -@ 8 -e -r -A mp.long.bam $long_ref 2>/dev/null"
check "calmd CRAM -u" "$st calmd --no-PG -u $corpus/short/input.3.1.cram $ref 2>/dev/null | $st view -h -" "$nt calmd --no-PG -u $corpus/short/input.3.1.cram $ref 2>/dev/null | $st view -h -"

# cram-size
for set in short tagged long; do
    for c in 3.0 3.1 3.1.multislice 3.1.noref; do
        cram=$corpus/$set/input.$c.cram
        check "cram-size $set $c" "$st cram-size $cram" "$nt cram-size $cram"
        check "cram-size $set $c -v -e" "$st cram-size -v -e $cram" "$nt cram-size -v -e $cram"
    done
done
"$st" view -T "$ref" -O cram,embed_ref=1,version=3.1,archive -o embed.cram "$tagged" 2>/dev/null && "$st" index embed.cram
check "cram-size embedded reference, archive codecs" "$st cram-size -v -e embed.cram" "$nt cram-size -v -e embed.cram"

# reference
check "reference short" "$st reference -q $short" "$nt reference -q $short"
check "reference long" "$st reference -q $long" "$nt reference -q $long"
check "reference CRAM" "$st reference -q $corpus/short/input.3.1.cram" "$nt reference -q $corpus/short/input.3.1.cram"
check "reference -r" "$st reference -q -r chrS:1000000-1003000 $short" "$nt reference -q -r chrS:1000000-1003000 $short"
check "reference -e" "$st reference -q -e embed.cram" "$nt reference -q -e embed.cram"
check "reference -e -r" "$st reference -q -e -r chrS:1000000-1003000 embed.cram" "$nt reference -q -e -r chrS:1000000-1003000 embed.cram"
# samtools misplaces bases after N CIGAR operations; ntools must agree with the true reference.
check "reference tagged (spliced reads) matches the reference" "echo 0" \
    "cmp -l <($nt reference -q $tagged | grep -v '>' | tr -d '\\n') <($st faidx $ref chrS | grep -v '>' | tr -d '\\n') | awk '\$2 != 116' | wc -l"

# rmdup, on md.bam's duplicates as one library and with its read groups as two
"$st" view --no-PG -h md.bam | sed 's/^@RG\tID:a\t.*/&\tLB:la/; s/^@RG\tID:b\t.*/&\tLB:lb/' | "$st" view -b -o md.lb.bam -
for input in md.bam md.lb.bam; do
    for opts in "" "-s" "-S"; do
        check "rmdup $opts $input" "$st rmdup $opts $input - 2>/dev/null | $st view --no-PG -h -" "$nt rmdup $opts $input - 2>/dev/null | $st view --no-PG -h -"
    done
done
check "rmdup -s long" "$st rmdup -s $long - 2>/dev/null | $st view --no-PG -h -" "$nt rmdup -s $long - 2>/dev/null | $st view --no-PG -h -"

# Records whose CIGAR (over 65,535 operations) is held in a CG tag
python3 "$scripts/long-cigar-test-data.py" > cg.sam
"$st" view -b -o cg.bam cg.sam && "$st" index cg.bam
for opts in "-C" "-t '*' -a" "-P -M -C -B"; do
    check "checksum $opts (CG CIGAR)" "$st checksum $opts cg.bam" "$nt checksum $opts cg.bam"
done
check "view SAM (CG CIGAR)" "$st view --no-PG cg.bam" "$nt view --no-PG cg.bam"
check "sort -O sam (CG CIGAR)" "$st sort --no-PG -O sam cg.bam" "$nt sort --no-PG -O sam cg.bam"

# bedcov: bad BED lines are reported and skipped; -d counts aligned bases (as samtools; ntools
# keeps deletions and skips in the depth sum unless -j)
printf "chrS 1000  5000 extra\nchrX\t1\t5\nchrS\t9000\t8000\nchrS\nbrowser x\n#c\nchrS\t20000x\t21000y\nchrS\t-5\t100\n" > bad.bed
check "bedcov bad BED lines" "$st bedcov bad.bed $short 2>&1; echo \$?" "$nt bedcov bad.bed $short 2>&1; echo \$?"
printf "chrS\t0\t300000\nchrS\t5000\t9000\n" > dcov.bed
check "bedcov -d (threshold column)" "$st bedcov -d 5 dcov.bed $tagged | cut -f5" "$nt bedcov -d 5 dcov.bed $tagged | cut -f5"
check "bedcov -j -d" "$st bedcov -j -d 5 dcov.bed $tagged" "$nt bedcov -j -d 5 dcov.bed $tagged"

# depth: -s (overlapping mates), -l (soft clips excluded), deletions' qualities with -J -q, and
# spliced reads' skips
for opts in "-s" "-s -J -q 30 -Q 20" "-l 100" "-J -q 30 -Q 20"; do
    for input in $tagged $long; do
        check "depth $opts $(basename "$(dirname "$input")")" "$st depth $opts $input" "$nt depth $opts $input"
    done
done

# coverage -q is the minimum mapping quality, -Q the minimum base quality
check "coverage -q 20 -Q 10" "$st coverage -q 20 -Q 10 $short" "$nt coverage -q 20 -Q 10 $short"
check "coverage -l 120 -q 10 long" "$st coverage -l 120 -q 10 $long" "$nt coverage -l 120 -q 10 $long"

# mpileup refuses inputs whose reference sequences differ from the first's
check "mpileup mismatched references" "echo 1" "$nt mpileup mp.short.bam md.bam > /dev/null 2>&1; echo \$?"

echo "Slowest checks (seconds, samtools run or cached):"
sort -rn "$timings" | head -15
awk -F'\t' '{ t += $1; if ($2 == "cache") c++ } END { printf "Total %.0f s over %d checks, %d from cache\n", t, NR, c }' "$timings"
echo "$failures difference(s)"
[ "$failures" -eq 0 ]
