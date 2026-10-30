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

# view: filter expressions, read groups, tag files, -m, -p, -P, -B, -t
for e in 'mapq >= 30 && flag.proper_pair && !flag.reverse' 'rname == "chrS" && pos > 100000 && endpos < 200000' \
    'qname =~ "1$" || cigar =~ "S"' 'sclen > 5 || hclen > 0 || ncigar > 1' 'tlen < -300 || mpos - pos > 200' \
    'avg(qual) > 30 && min(qual) < 10' 'seq =~ "^A+C" || length(seq) != qlen' '[NM] > 2 || ![MD]' \
    'default([XS], -1) < 0 && exists([AS])' '(flag & 0x50) == 0x50 && ~flag & 16' 'sqrt(mapq) > 5 && pow(2, 3) == 8' \
    '[XX] || mapq > 50' '-tlen > 0' 'rlen % 7 == 3' '"a" < "b" && library == ""'; do
    check "view -e '$e'" "$st view --no-PG -e '$e' $short" "$nt view --no-PG -e '$e' $short"
done
check "view -c -e" "$st view -c -e 'mapq > 30 && [NM] > 1' $short" "$nt view -c -e 'mapq > 30 && [NM] > 1' $short"
check "view -e tagged" "$st view --no-PG -e '[NH] > 1 && [HI] == 1' $tagged" "$nt view --no-PG -e '[NH] > 1 && [HI] == 1' $tagged"
check "view -e long" "$st view --no-PG -e 'qlen > 10000 && sclen > 100' $long" "$nt view --no-PG -e 'qlen > 10000 && sclen > 100' $long"
check "view -e CRAM -@4" "$st view --no-PG -e 'mapq > 30 && [NM] > 1' -T $ref $corpus/short/input.3.1.cram" "$nt view --no-PG -@ 4 -e 'mapq > 30 && [NM] > 1' -T $ref $corpus/short/input.3.1.cram"
check "view -c -e CRAM -@4" "$st view -c -e 'mapq > 30' -T $ref $corpus/short/input.3.1.cram" "$nt view -c -@ 4 -e 'mapq > 30' -T $ref $corpus/short/input.3.1.cram"
check "view -e -U" "$st view --no-PG -e 'mapq > 30' -U /dev/stdout -o /dev/null $short" "$nt view --no-PG -e 'mapq > 30' -U /dev/stdout -o /dev/null $short"
for m in 1 150 151; do
    check "view -c -m $m" "$st view -c -m $m $short" "$nt view -c -m $m $short"
done
check "view -m long" "$st view --no-PG -m 5000 $long" "$nt view --no-PG -m 5000 $long"
printf 'a\n' > rg.a.txt
printf 'a c\n zz\n' > rg.ac.txt
: > rg.empty.txt
for opts in "-r a" "-r a -r c -r zz" "-R rg.a.txt" "-R rg.ac.txt" "-R rg.empty.txt" "-n" "-n -r a" "-c -r a" "-l x -r a"; do
    check "view -h $opts" "$st view --no-PG -h $opts lib.sam" "$nt view --no-PG -h $opts lib.sam"
done
# samtools keeps the @RG lines of the dropped read groups with -R ^FILE; compare records only.
check "view -R ^FILE" "$st view --no-PG -R ^rg.ac.txt lib.sam" "$nt view --no-PG -R ^rg.ac.txt lib.sam"
printf '2\n3\n' > nh.txt
for opts in "-D NH:nh.txt" "-d NH:1 -D NH:nh.txt" "-D NH:rg.empty.txt" "-d NH -d NH:2" "-c -D NH:nh.txt"; do
    check "view $opts" "$st view --no-PG $opts $tagged" "$nt view --no-PG $opts $tagged"
done
check "view -p" "$st view --no-PG -p -q 30 $short" "$nt view --no-PG -p -q 30 $short"
check "view -p -x region" "$st view --no-PG -p -f 2 -x NM $short chrS:100000-200000" "$nt view --no-PG -p -f 2 -x NM $short chrS:100000-200000"
check "view -p -b" "$st view --no-PG -p -e 'mapq < 20' $short" "$nt view --no-PG -b -p -e 'mapq < 20' $short | $st view --no-PG -"
check "view -c -p" "$st view -c -p -q 30 $short" "$nt view -c -p -q 30 $short"
head -5 "$corpus/short/queries.100.bed" > pairs.bed
for opts in "chrS:100000-110000" "chrS:100000-110000 chrS:500000-505000" "-L pairs.bed" "-L pairs.bed chrS:1-3000000" "-q 30 chrS:100000-110000" "-f 64 chrS:100000-110000"; do
    check "view -P $opts" "$st view --no-PG -P $short $opts" "$nt view --no-PG -P $short $opts"
done
check "view -P long" "$st view --no-PG -P $long $(head -1 "$corpus/long/queries.100.regions")" "$nt view --no-PG -P $long $(head -1 "$corpus/long/queries.100.regions")"
printf '@SQ\tSN:c\tLN:100\nr\t0\tc\t1\t60\t5M2B5M\t*\t0\t0\tACGTAGTCCA\t(((((5((((\nr2\t0\tc\t1\t60\t2B5M\t*\t0\t0\tACGTA\t*\nr3\t0\tc\t1\t60\t3M1B2M\t*\t0\t0\tACGTA\t*\nr4\t0\tc\t1\t60\t3M1I2M\t*\t0\t0\tACGTAC\t*\n' > back.sam
"$st" view -b -o back.bam back.sam
for opts in "-B" "-B -m 5" "-B -F 4" "-B -e 'qlen > 5' -U /dev/stdout -o /dev/null"; do
    check "view $opts" "$st view --no-PG $opts back.bam" "$nt view --no-PG $opts back.bam"
done
"$st" view "$short" chrS:1-20000 > nosq.sam
# ntools writes the @SQ lines from -t first; HTSlib appends them after the others.
check "view -t" "$st view --no-PG -h -t $ref.fai nosq.sam | sort" "$nt view --no-PG -h -t $ref.fai nosq.sam | sort"
check "view -T (fai for SAM)" "$st view --no-PG -b -T $ref nosq.sam | $st view --no-PG -h -" "$nt view --no-PG -b -T $ref nosq.sam | $st view --no-PG -h -"
check "view -S" "$st view -S -c $short" "$nt view -S -c $short"

# statistics
for set in short tagged long; do
    for fmt in default tsv json; do
        check "flagstat $set $fmt" "$st flagstat -O $fmt $corpus/$set/input.bam" "$nt flagstat -O $fmt -@ 4 $corpus/$set/input.bam"
    done
    check "idxstats $set" "$st idxstats $corpus/$set/input.bam" "$nt idxstats $corpus/$set/input.bam"
done

# idxstats -X: the index file follows the input
ln -sf "$short" idxstats.bam
for index in "$corpus/short/input.bam.bai" "$corpus/short/input.bam.csi"; do
    check "idxstats -X ${index##*.}" "$st idxstats -X idxstats.bam $index" "$nt idxstats -X idxstats.bam $index"
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

# sort options: -t (integer, float, string, and character tags), -n tie-breaks, -H/-R (minimiser
# options, accepted and ignored)
check "sort -n long (secondary/supplementary order)" "$st sort -n --no-PG -o - $long | $st view --no-PG -h -" "$nt sort -n --no-PG -o - $long | $st view --no-PG -h -"
check "sort -t CB" "$st sort -t CB --no-PG -o - $tagged | $st view --no-PG -h -" "$nt sort -t CB -m 1M --no-PG -o - $tagged | $st view --no-PG -h -"
check "sort -t NH -n" "$st sort -t NH -n --no-PG -o - $short | $st view --no-PG -h -" "$nt sort -t NH -n -@ 3 --no-PG -o - $short | $st view --no-PG -h -"
check "sort -t AS -N" "$st sort -t AS -N --no-PG -o - $short | $st view --no-PG -h -" "$nt sort -t AS -N -m 1M --no-PG -o - $short | $st view --no-PG -h -"
check "sort -t qs (float)" "$st sort -t qs --no-PG -o - $long | $st view --no-PG -h -" "$nt sort -t qs --no-PG -o - $long | $st view --no-PG -h -"
{
    printf '@HD\tVN:1.6\tSO:unsorted\tGO:query\n@SQ\tSN:c1\tLN:1000\n'
    for t in 'XX:i:5' 'XX:f:4.5' 'XX:Z:abc' 'XX:A:b' 'XX:i:-3' '' 'XX:f:5' 'XX:Z:ab' 'XX:A:a' 'XX:B:c,1' 'XX:H:1F' 'XX:i:4'; do
        for pos in 30 10; do
            printf 'r%s\t0\tc1\t%s\t60\t4M\t*\t0\t0\tACGT\tIIII%s\n' "$pos" "$pos" "${t:+	$t}"
        done
    done
} > sort.tags.sam
check "sort -t (mixed types)" "$st sort -t XX --no-PG -O sam sort.tags.sam" "$nt sort -t XX --no-PG -O sam sort.tags.sam"
check "sort -t -n (mixed types)" "$st sort -t XX -n --no-PG -O sam sort.tags.sam" "$nt sort -t XX -n --no-PG -O sam sort.tags.sam"
check "sort -H -R" "$st sort -H -R --no-PG -o - $short | $st view --no-PG -h -" "$nt sort -H -R --no-PG -o - $short | $st view --no-PG -h -"

# sort --template-coordinate: fixmate -m adds the MC tags it needs (short reads whose mates are
# missing from the corpus, which samtools rejects, are dropped); hand-made edge cases (soft and
# hard clips, unmapped reads and mates, unpaired reads on both strands, other references,
# secondary and supplementary records, same-position ties, libraries, CB, MI with /A suffixes,
# MC of *); precedence with -n/-N/-t; tag errors. --write-index (coordinate order only), with
# FILE##idx##INDEX names.
"$st" fixmate -m --no-PG short.name.bam - | "$st" view -b --no-PG -e 'exists([MC]) || !flag.paired || flag.munmap' -o sort.tc.short.bam -
"$st" sort -n --no-PG "$long" | "$st" fixmate -m --no-PG - sort.tc.long.bam
check "sort --template-coordinate" "$st sort --template-coordinate --no-PG -o - sort.tc.short.bam | $st view --no-PG -h -" "$nt sort --template-coordinate --no-PG -m 1M -@ 3 -o - sort.tc.short.bam | $st view --no-PG -h -"
check "sort --template-coordinate long" "$st sort --template-coordinate --no-PG -o - sort.tc.long.bam | $st view --no-PG -h -" "$nt sort --template-coordinate --no-PG -o - sort.tc.long.bam | $st view --no-PG -h -"
{
    s=ACGTACGTAC q=IIIIIIIIII
    printf '@HD\tVN:1.6\tSO:queryname\n@SQ\tSN:c1\tLN:10000\n@SQ\tSN:c2\tLN:10000\n@RG\tID:g1\tSM:s\tLB:libB\n@RG\tID:g2\tSM:s\tLB:libA\n@RG\tID:g3\tSM:s\n'
    printf "a\t99\tc1\t100\t60\t2S8M\t=\t300\t210\t$s\t$q\tRG:Z:g1\tMC:Z:10M2S\n"
    printf "a\t147\tc1\t300\t60\t10M2S\t=\t100\t-210\t${s}AC\t${q}II\tRG:Z:g1\tMC:Z:2S8M\n"
    printf "b\t99\tc1\t100\t60\t2S8M\t=\t300\t210\t$s\t$q\tRG:Z:g2\tMC:Z:10M2S\n"
    printf "b\t147\tc1\t300\t60\t10M2S\t=\t100\t-210\t${s}AC\t${q}II\tRG:Z:g2\tMC:Z:2S8M\n"
    printf "c\t73\tc1\t50\t60\t10M\t=\t50\t0\t$s\t$q\n"
    printf "c\t133\tc1\t50\t0\t*\t=\t50\t0\t$s\t$q\tMC:Z:10M\n"
    printf "d\t77\t*\t0\t0\t*\t*\t0\t0\t$s\t$q\n"
    printf "d\t141\t*\t0\t0\t*\t*\t0\t0\t$s\t$q\n"
    printf "e\t0\tc2\t20\t60\t10M\t*\t0\t0\t$s\t$q\n"
    printf "e2\t16\tc2\t11\t60\t10M\t*\t0\t0\t$s\t$q\n"
    printf "f\t97\tc2\t500\t60\t3S7M\tc1\t40\t0\t$s\t$q\tMC:Z:5H5M\n"
    printf "f\t145\tc1\t40\t60\t5H5M\tc2\t500\t0\tACGTA\tIIIII\tMC:Z:3S7M\n"
    printf "f\t353\tc1\t800\t0\t10M\tc1\t40\t0\t$s\t$q\tMC:Z:5H5M\n"
    printf "f\t2129\tc1\t900\t0\t4S6M\tc1\t40\t0\t$s\t$q\tMC:Z:5H5M\n"
    printf "g\t65\tc1\t100\t60\t10M\t=\t100\t0\t$s\t$q\tMC:Z:10M\tCB:Z:AAA\tMI:Z:1/A\n"
    printf "g\t129\tc1\t100\t60\t10M\t=\t100\t0\t$s\t$q\tMC:Z:10M\tCB:Z:AAA\tMI:Z:1/B\n"
    printf "h\t65\tc1\t100\t60\t10M\t=\t100\t0\t$s\t$q\tMC:Z:10M\tCB:Z:AAA\tMI:Z:1\n"
    printf "h\t129\tc1\t100\t60\t10M\t=\t100\t0\t$s\t$q\tMC:Z:10M\tCB:H:AA\tMI:Z:10\n"
    printf "i\t115\tc1\t100\t60\t10M\t=\t100\t0\t$s\t$q\tMC:Z:10M\tRG:Z:g3\n"
    printf "i\t179\tc1\t100\t60\t10M\t=\t100\t0\t$s\t$q\tMC:Z:10M\tRG:Z:gX\n"
    printf "j\t97\tc1\t200\t60\t3M2I3M2D2M\t=\t250\t0\t$s\t$q\tMC:Z:2H3S4M1D2N3M2S4H\n"
    printf "j\t145\tc1\t250\t60\t3S5M2S\t=\t200\t0\t$s\t$q\tMC:Z:*\n"
    printf "k\t69\tc1\t250\t0\t*\t=\t250\t0\t$s\t$q\tMC:Z:3S5M2S\n"
    printf "k\t137\tc1\t250\t60\t3S5M2S\t=\t250\t0\t$s\t$q\n"
    printf "l\t355\tc2\t70\t0\t10M\t=\t90\t0\t$s\t$q\tMC:Z:10M\n"
    printf "l\t2147\tc2\t70\t0\t10M\t=\t90\t0\t$s\t$q\tMC:Z:10M\n"
    printf "l\t99\tc2\t70\t0\t10M\t=\t90\t0\t$s\t$q\tMC:Z:10M\n"
    printf "l\t147\tc2\t90\t0\t10M\t=\t70\t0\t$s\t$q\tMC:Z:10M\n"
} > sort.tc.sam
check "sort --template-coordinate (edge cases)" "$st sort --template-coordinate --no-PG -O sam sort.tc.sam" "$nt sort --template-coordinate --no-PG -O sam sort.tc.sam"
check "sort --template-coordinate (edge cases, spilled)" "$st sort --template-coordinate --no-PG -O sam sort.tc.sam" "$nt sort --template-coordinate -m 1K --no-PG -O sam sort.tc.sam"
for opts in "-n --template-coordinate" "--template-coordinate -n" "--template-coordinate -N" "-t XX --template-coordinate" "--template-coordinate -n -t XX" "-n -N"; do
    check "sort $opts (precedence)" "$st sort $opts --no-PG -O sam sort.tc.sam" "$nt sort $opts --no-PG -O sam sort.tc.sam"
done
for tag in MC CB MI none; do
    case $tag in
        MC) fields='\tMC:i:5' ;;
        none) fields='' ;;
        *) fields="\\tMC:Z:10M\\t$tag:i:5" ;;
    esac
    printf "@SQ\tSN:c1\tLN:1000\nx\t65\tc1\t100\t60\t10M\t=\t100\t0\tACGTACGTAC\tIIIIIIIIII$fields\n" > "sort.tc.bad.$tag.sam"
    check "sort --template-coordinate bad $tag (exit status)" "$st sort --template-coordinate --no-PG -O sam sort.tc.bad.$tag.sam; echo \$?" "$nt sort --template-coordinate --no-PG -O sam sort.tc.bad.$tag.sam 2>/dev/null; echo \$?"
done
sort_region=$(head -1 "$corpus/short/queries.100.regions")
check "sort --write-index (##idx## BAI)" "$st sort --write-index --no-PG -o 'st.sort.wi.bam##idx##st.sort.wi.bai' short.name.bam && head -c 4 st.sort.wi.bai | od -c && $st view -X st.sort.wi.bam st.sort.wi.bai $sort_region | md5sum" "$nt sort --write-index --no-PG -o 'nt.sort.wi.bam##idx##nt.sort.wi.bai' short.name.bam && head -c 4 nt.sort.wi.bai | od -c && $st view -X nt.sort.wi.bam nt.sort.wi.bai $sort_region | md5sum"
check "sort --write-index (CSI)" "$st sort --write-index --no-PG -o st.sort.wi2.bam short.name.bam && head -c 4 st.sort.wi2.bam.csi | od -c && $st idxstats st.sort.wi2.bam && $st view st.sort.wi2.bam $sort_region | md5sum" "$nt sort --write-index --no-PG -o nt.sort.wi2.bam short.name.bam && head -c 4 nt.sort.wi2.bam.csi | od -c && $st idxstats nt.sort.wi2.bam && $st view nt.sort.wi2.bam $sort_region | md5sum"
check "sort --template-coordinate --write-index (ignored)" "$st sort --template-coordinate --write-index --no-PG -o st.sort.wi3.bam sort.tc.sam; ls st.sort.wi3.bam*" "$nt sort --template-coordinate --write-index --no-PG -o nt.sort.wi3.bam sort.tc.sam 2>/dev/null; ls nt.sort.wi3.bam* | sed s/^nt/st/"

"$st" view --no-PG -b -f 64 -o p1.bam "$short"
"$st" view --no-PG -b -F 64 -o p2.bam "$short"
check "merge" "$st merge -c -p --no-PG -f -o - p1.bam p2.bam | $st view --no-PG - | cut -f1-11" "$nt merge --no-PG -f -o /dev/stdout p1.bam p2.bam | $st view --no-PG - | cut -f1-11"

# merge options. The corpus inputs share their @RG and @PG IDs, which plain merge renames (with
# -s for reproducible suffixes) and rewrites in records, moving RG and PG to the end.
"$st" index p1.bam
"$st" index p2.bam
cp p1.bam.bai merge.1.bai
cp p2.bam.bai merge.2.bai
printf 'p1.bam\n\np2.bam\n' > merge.list
merge_region=$(head -1 "$corpus/short/queries.100.regions")
check "merge (renamed IDs)" "$st merge -s 7 --no-PG -o - p1.bam p2.bam | $st view --no-PG -h -" "$nt merge -s 7 --no-PG -o - p1.bam p2.bam | $st view --no-PG -h -"
check "merge -c -p" "$st merge -c -p --no-PG -o - p1.bam p2.bam | $st view --no-PG -h -" "$nt merge -c -p --no-PG -o - p1.bam p2.bam | $st view --no-PG -h -"
check "merge -r -1 -s" "$st merge -r -1 -s 3 --no-PG -o - p1.bam p2.bam | $st view --no-PG -h -" "$nt merge -r -1 -s 3 --no-PG -o - p1.bam p2.bam | $st view --no-PG -h -"
check "merge -b" "$st merge -b merge.list -c -p --no-PG -o - | $st view --no-PG -h -" "$nt merge -b merge.list -c -p --no-PG -o - | $st view --no-PG -h -"
check "merge -R" "$st merge -R $merge_region -c -p --no-PG -o - p1.bam p2.bam | $st view --no-PG -h -" "$nt merge -R $merge_region -c -p --no-PG -o - p1.bam p2.bam | $st view --no-PG -h -"
check "merge -X -R" "$st merge -X -R $merge_region -c -p --no-PG -o - p1.bam p2.bam merge.1.bai merge.2.bai | $st view --no-PG -h -" "$nt merge -X -R $merge_region -c -p --no-PG -o - p1.bam p2.bam merge.1.bai merge.2.bai | $st view --no-PG -h -"
check "merge -L" "$st merge -L q.bed -c -p --no-PG -o - p1.bam p2.bam | $st view --no-PG -h -" "$nt merge -L q.bed -c -p --no-PG -o - p1.bam p2.bam | $st view --no-PG -h -"
check "merge -h" "$st merge -h $tagged -s 5 --no-PG -o - p1.bam p2.bam | $st view --no-PG -h -" "$nt merge -h $tagged -s 5 --no-PG -o - p1.bam p2.bam | $st view --no-PG -h -"
"$st" view --no-PG -b -f 64 -o merge.n1.bam short.name.bam
"$st" view --no-PG -b -F 64 -o merge.n2.bam short.name.bam
check "merge -n" "$st merge -n -s 1 --no-PG -o - merge.n1.bam merge.n2.bam | $st view --no-PG -h -" "$nt merge -n -s 1 --no-PG -o - merge.n1.bam merge.n2.bam | $st view --no-PG -h -"
"$st" sort -t CB --no-PG -o merge.cb.bam "$tagged"
"$st" view --no-PG -b -f 64 -o merge.t1.bam merge.cb.bam
"$st" view --no-PG -b -F 64 -o merge.t2.bam merge.cb.bam
check "merge -t CB" "$st merge -t CB -c -p --no-PG -o - merge.t1.bam merge.t2.bam | $st view --no-PG -h -" "$nt merge -t CB -c -p --no-PG -o - merge.t1.bam merge.t2.bam | $st view --no-PG -h -"
"$st" view --no-PG -C -T "$ref" -o merge.2.cram p2.bam
"$st" index merge.2.cram
check "merge CRAM input -R" "$st merge -R $merge_region -c -p --reference $ref --no-PG -o - p1.bam merge.2.cram | $st view --no-PG -h -" "$nt merge -R $merge_region -c -p --reference $ref --no-PG -o - p1.bam merge.2.cram | $st view --no-PG -h -"
# Colliding @RG/@PG IDs with PG/PP links, unknown RG/PG tags (dropped), @CO lines (all kept),
# and a reference order that loses coordinate sorting
printf '@HD\tVN:1.6\tSO:coordinate\n@SQ\tSN:c1\tLN:1000\n@SQ\tSN:c2\tLN:1000\n@RG\tID:g1\tSM:s1\tPG:p1\n@RG\tID:g2\tSM:s2\n@PG\tID:p1\tPN:bwa\n@PG\tID:p2\tPN:x\tPP:p1\n@CO\tone\n' > merge.a.sam
printf 'a1\t0\tc1\t10\t60\t5M\t*\t0\t0\tACGTA\tIIIII\tRG:Z:g1\tXX:i:1\tPG:Z:p2\na2\t0\tc2\t10\t60\t5M\t*\t0\t0\tACGTA\tIIIII\tRG:Z:g9\tXX:i:2\na3\t4\t*\t0\t0\t*\t*\t0\t0\tACGTA\tIIIII\tRG:Z:g2\n' >> merge.a.sam
printf '@HD\tVN:1.6\tSO:coordinate\tSS:coordinate:x\n@SQ\tSN:c2\tLN:1000\n@SQ\tSN:c3\tLN:500\n@RG\tID:g1\tSM:s9\tPG:p1\n@PG\tID:p1\tPN:bowtie\n@CO\tone\n@CO\ttwo\n' > merge.b.sam
printf 'b1\t0\tc2\t5\t60\t5M\t*\t0\t0\tACGTA\tIIIII\tXY:Z:q\tRG:Z:g1\tPG:Z:p1\tXX:i:1\nb2\t0\tc3\t10\t60\t5M\t*\t0\t0\tACGTA\tIIIII\tPG:Z:zz\tRG:Z:g1\n' >> merge.b.sam
for opts in "-s 42" "-c -s 42" "-p -s 42" "-r -s 42" "-h merge.b.sam -s 42" "-r -s 42 merge.a.sam"; do
    check "merge $opts (header edge cases)" "$st merge $opts --no-PG -o - merge.a.sam merge.b.sam | $st view --no-PG -h -" "$nt merge $opts --no-PG -o - merge.a.sam merge.b.sam | $st view --no-PG -h -"
done

# merge --template-coordinate: corpus halves, the sort edge cases, colliding read group IDs whose
# libraries differ (so -c and -r change the order; samtools keys records by their translated
# read group before -r replaces it, and accepts H-typed RG values), precedence with -n/-N/-t, and
# a missing MC. --write-index: CSI, BAI or CSI named with FILE##idx##INDEX, a positional output,
# and the errors for SAM, stdout, and output that is not coordinate-sorted.
"$st" sort --template-coordinate --no-PG -o merge.tc.bam sort.tc.short.bam
"$st" view --no-PG -b -f 64 -o merge.tc1.bam merge.tc.bam
"$st" view --no-PG -b -F 64 -o merge.tc2.bam merge.tc.bam
check "merge --template-coordinate" "$st merge --template-coordinate -s 1 --no-PG -o - merge.tc1.bam merge.tc2.bam | $st view --no-PG -h -" "$nt merge --template-coordinate -s 1 --no-PG -o - merge.tc1.bam merge.tc2.bam | $st view --no-PG -h -"
check "merge --template-coordinate -c -p -r" "$st merge --template-coordinate -c -p -r --no-PG -o - merge.tc1.bam merge.tc2.bam | $st view --no-PG -h -" "$nt merge --template-coordinate -c -p -r --no-PG -o - merge.tc1.bam merge.tc2.bam | $st view --no-PG -h -"
"$st" sort --template-coordinate --no-PG -O sam -o merge.tc.edge.sam sort.tc.sam
check "merge --template-coordinate (edge cases)" "$st merge --template-coordinate --no-PG -o - merge.tc.edge.sam merge.tc.edge.sam | $st view --no-PG -h -" "$nt merge --template-coordinate --no-PG -o - merge.tc.edge.sam merge.tc.edge.sam | $st view --no-PG -h -"
{
    s=ACGTACGTAC q=IIIIIIIIII
    printf '@HD\tVN:1.6\tSO:unsorted\tGO:query\tSS:unsorted:template-coordinate\n@SQ\tSN:c1\tLN:1000\n@RG\tID:g1\tSM:s\tLB:libZ\n'
    printf "a1\t99\tc1\t100\t60\t10M\t=\t200\t110\t$s\t$q\tMC:Z:10M\tRG:Z:g1\n"
    printf "a1\t147\tc1\t200\t60\t10M\t=\t100\t-110\t$s\t$q\tMC:Z:10M\tRG:Z:g1\n"
    printf "a2\t0\tc1\t300\t60\t10M\t*\t0\t0\t$s\t$q\tRG:Z:g9\n"
} > merge.tc.a.sam
{
    s=ACGTACGTAC q=IIIIIIIIII
    printf '@HD\tVN:1.6\tSO:unsorted\tGO:query\tSS:unsorted:template-coordinate\n@SQ\tSN:c1\tLN:1000\n@RG\tID:g1\tSM:s\tLB:libA\n@RG\tID:AB\tSM:s\tLB:libC\n'
    printf "b1\t99\tc1\t100\t60\t10M\t=\t200\t110\t$s\t$q\tMC:Z:10M\tRG:Z:g1\n"
    printf "b1\t147\tc1\t200\t60\t10M\t=\t100\t-110\t$s\t$q\tMC:Z:10M\tRG:Z:g1\n"
    printf "b0\t0\tc1\t300\t60\t10M\t*\t0\t0\t$s\t$q\tRG:H:AB\n"
} > merge.tc.b.sam
for opts in "-s 3" "-c -s 3" "-r -s 3" "-c -r" "-h merge.tc.b.sam -s 3"; do
    check "merge --template-coordinate $opts (libraries)" "$st merge --template-coordinate $opts --no-PG -o - merge.tc.a.sam merge.tc.b.sam | $st view --no-PG -h -" "$nt merge --template-coordinate $opts --no-PG -o - merge.tc.a.sam merge.tc.b.sam | $st view --no-PG -h -"
done
for opts in "-n --template-coordinate" "--template-coordinate -n" "--template-coordinate -N" "-t XX --template-coordinate"; do
    check "merge $opts (precedence)" "$st merge $opts --no-PG -o - merge.tc.edge.sam merge.tc.edge.sam | $st view --no-PG -h -" "$nt merge $opts --no-PG -o - merge.tc.edge.sam merge.tc.edge.sam | $st view --no-PG -h -"
done
check "merge --template-coordinate missing MC (exit status)" "$st merge --template-coordinate --no-PG -o - sort.tc.bad.none.sam > /dev/null; echo \$?" "$nt merge --template-coordinate --no-PG -o - sort.tc.bad.none.sam 2>/dev/null > /dev/null; echo \$?"
for idx in '' merge.wi1.bai merge.wi2.csi merge.wi3.idx; do
    if [ -z "$idx" ]; then
        st_out=st.merge.wi0.bam nt_out=nt.merge.wi0.bam st_idx=st.merge.wi0.bam.csi nt_idx=nt.merge.wi0.bam.csi
    else
        st_out="st.$idx.bam##idx##st.$idx" nt_out="nt.$idx.bam##idx##nt.$idx" st_idx=st.$idx nt_idx=nt.$idx
    fi
    check "merge --write-index (${idx:-CSI})" "$st merge --write-index -c -p --no-PG -o '$st_out' p1.bam p2.bam && head -c 4 $st_idx | od -c && $st view -X ${st_out%%#*} $st_idx $merge_region | md5sum" "$nt merge --write-index -c -p --no-PG -o '$nt_out' p1.bam p2.bam && head -c 4 $nt_idx | od -c && $st view -X ${nt_out%%#*} $nt_idx $merge_region | md5sum"
done
check "merge --write-index (positional output)" "$st merge --write-index -c -p --no-PG 'st.merge.wi4.bam##idx##st.merge.wi4.bai' p1.bam p2.bam && $st idxstats st.merge.wi4.bam && $st view -X st.merge.wi4.bam st.merge.wi4.bai $merge_region | md5sum" "$nt merge --write-index -c -p --no-PG 'nt.merge.wi4.bam##idx##nt.merge.wi4.bai' p1.bam p2.bam && $st idxstats nt.merge.wi4.bam && $st view -X nt.merge.wi4.bam nt.merge.wi4.bai $merge_region | md5sum"
check "merge --write-index SAM (exit status)" "$st merge --write-index --no-PG -o st.merge.wi5.sam p1.bam p2.bam; echo \$?" "$nt merge --write-index --no-PG -o nt.merge.wi5.sam p1.bam p2.bam 2>/dev/null; echo \$?"
check "merge --write-index stdout (exit status)" "$st merge --write-index --no-PG -o - p1.bam p2.bam > /dev/null; echo \$?" "$nt merge --write-index --no-PG -o - p1.bam p2.bam 2>/dev/null > /dev/null; echo \$?"
check "merge -n --write-index (exit status)" "$st merge -n --write-index --no-PG -o st.merge.wi6.bam merge.n1.bam merge.n2.bam; echo \$?" "$nt merge -n --write-index --no-PG -o nt.merge.wi6.bam merge.n1.bam merge.n2.bam 2>/dev/null; echo \$?"

# duplicate marking pipeline
check "fixmate" "$st fixmate --no-PG short.name.bam - | $st view --no-PG -h -" "$nt fixmate --no-PG short.name.bam /dev/stdout | $st view --no-PG -h -"
check "fixmate -m -c -r -@4" "$st fixmate --no-PG -m -c -r short.name.bam - | $st view --no-PG -h -" "$nt fixmate --no-PG -m -c -r -@ 4 short.name.bam /dev/stdout | $st view --no-PG -h -"

# fixmate -M: base modification tags on long reads, and synthetic primaries, hard-clipped
# supplementary and secondary alignments (both strands), draft Mm/Ml tags, and invalid tags
fm_long_name=$("$st" idxstats "$long" | head -1 | cut -f1)
"$st" view -u "$long" "$fm_long_name:1-300000" | "$st" sort -n --no-PG -o fm.long.bam -
{
    printf '@HD\tVN:1.6\tSO:queryname\n@SQ\tSN:c1\tLN:1000\n'
    s=ACGTCCGACGTCAACGCCGA q=IIIIIIIIIIIIIIIIIIII mm='MM:Z:C+m,0,2,1;\tML:B:C,10,20,30'
    printf "t1\t0\tc1\t100\t60\t20M\t*\t0\t0\t$s\t$q\t$mm\tXX:i:1\n"
    printf "t1\t2048\tc1\t300\t60\t5H10M5H\t*\t0\t0\tCGACGTCAAC\tIIIIIIIIII\t$mm\tXX:i:2\n"
    printf "t1\t2064\tc1\t500\t60\t4H13M3H\t*\t0\t0\tCGTTGACGTCGGA\tIIIIIIIIIIIII\t$mm\tXX:i:3\n"
    printf "t1\t256\tc1\t600\t60\t20M\t*\t0\t0\t$s\t$q\t$mm\tMN:i:20\n"
    printf "t1\t256\tc1\t700\t60\t3H10M\t*\t0\t0\tTCCGACGTCA\tIIIIIIIIII\t$mm\n"
    printf "t2\t0\tc1\t100\t60\t20M\t*\t0\t0\t$s\t$q\tMm:Z:C+m,0,2,1;\tMl:B:C,10,20,30\tXX:i:1\n"
    printf "t2\t2048\tc1\t300\t60\t5H10M5H\t*\t0\t0\tCGACGTCAAC\tIIIIIIIIII\tMm:Z:C+m,0,2,1;\tMl:B:C,10,20,30\tMn:i:20\n"
    printf "t3\t0\tc1\t100\t60\t5H15M\t*\t0\t0\tCGACGTCAACGCCGA\tIIIIIIIIIIIIIII\t$mm\n"
    printf "t4\t0\tc1\t100\t60\t5H15M\t*\t0\t0\tCGACGTCAACGCCGA\tIIIIIIIIIIIIIII\tMM:Z:C+m,0,1;\tML:B:C,10,20\tMN:i:15\n"
    printf "t5\t0\tc1\t100\t60\t20M\t*\t0\t0\t$s\t$q\tMM:Z:C+m,0,2,9;\tML:B:C,10,20,30\n"
    printf "t6\t0\tc1\t100\t60\t20M\t*\t0\t0\t$s\t$q\tML:B:C,10,20,30\tMN:i:20\tXX:i:6\n"
    printf "t7\t65\tc1\t100\t60\t20M\t=\t200\t0\t$s\t$q\tMM:Z:C+m?,0,2,1;A+a,1;\tML:B:C,10,20,30,40\n"
    printf "t7\t129\tc1\t200\t60\t20M\t=\t100\t0\t$s\t$q\t$mm\n"
    printf "t7\t2113\tc1\t400\t60\t10H10M\t=\t200\t0\tTCAACGCCGA\tIIIIIIIIII\tMM:Z:C+m?,0,2,1;A+a,1;\tML:B:C,10,20,30,40\n"
    printf "t7\t2177\tc1\t400\t60\t10M10H\t=\t100\t0\tACGTCCGACG\tIIIIIIIIII\t$mm\n"
    printf "t8\t256\tc1\t100\t60\t20M\t*\t0\t0\t$s\t$q\t$mm\n"
    printf "t8\t256\tc1\t100\t60\t20M\t*\t0\t0\t$s\t$q\tMM:Z:C+m,0,2,9;\tML:B:C,10,20,30\n"
    printf "t9\t16\tc1\t100\t60\t20M\t*\t0\t0\tTCGGCGTTGACGTCGGACGT\t$q\t$mm\n"
    printf "t9\t2064\tc1\t300\t60\t5H10M5H\t*\t0\t0\tGTTGACGTCG\tIIIIIIIIII\t$mm\n"
    printf "t10\t0\tc1\t100\t60\t20M\t*\t0\t0\t$s\t*\tMM:Z:C+m.,0,2,1;\tML:B:C,10,20,30\n"
    printf 't11\t4\t*\t0\t0\t*\t*\t0\t0\t*\t*\tMM:Z:C+m;\tMN:i:0\n'
} > fm-mods.sam
for opts in "-M" "-M -r -m -c -@2"; do
    check "fixmate $opts (synthetic mods)" "$st fixmate --no-PG $opts fm-mods.sam - | $st view --no-PG -h -" "$nt fixmate --no-PG $opts fm-mods.sam /dev/stdout | $st view --no-PG -h -"
done
check "fixmate -M (long)" "$st fixmate --no-PG -M fm.long.bam - | $st view --no-PG -h -" "$nt fixmate --no-PG -M fm.long.bam /dev/stdout | $st view --no-PG -h -"
check "fixmate -M -@4 (long)" "$st fixmate --no-PG -M fm.long.bam - | $st view --no-PG -h -" "$nt fixmate --no-PG -M -@ 4 fm.long.bam /dev/stdout | $st view --no-PG -h -"
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
check "collate -n -r (accepted, ignored)" "$st view $short | cut -f1 | sort -u | wc -l" "$nt collate --no-PG -n 8 -r 100 -O $short | $st view - | awk -F'\\t' '\$1 != p { print \$1; p = \$1 }' | sort | uniq | wc -l"
check "collate header (GO, SS kept)" "$st collate --no-PG -O short.name.bam | $st view --no-PG -H -" "$nt collate --no-PG -O short.name.bam | $st view --no-PG -H -"

# reheader: -c commands and CRAM (copied, and edited in place with -i)
"$st" view -H "$short" | sed 's/SM:synthetic/SM:renamed/' > reheader.sam
check "reheader -c (BAM)" "$st reheader -P -c 'sed s/synthetic/xyz/g' $short | $st view --no-PG -h -" "$nt reheader -P -c 'sed s/synthetic/xyz/g' $short | $st view --no-PG -h -"
for v in 3.0 3.1; do
    cram=$corpus/short/input.$v.cram
    check "reheader CRAM $v" "$st reheader -P reheader.sam $cram | $st view --no-PG -h -T $ref -" "$nt reheader -P reheader.sam $cram | $st view --no-PG -h -T $ref -"
    check "reheader -c CRAM $v" "$st reheader -P -c 'grep -v ^@PG' $cram | $st view --no-PG -h -T $ref -" "$nt reheader -P -c 'grep -v ^@PG' $cram | $st view --no-PG -h -T $ref -"
    check "reheader -i CRAM $v" "cp $cram st.reheader.cram && $st reheader -P -i reheader.sam st.reheader.cram && $st view --no-PG -h -T $ref st.reheader.cram" "cp $cram nt.reheader.cram && $nt reheader -P -i reheader.sam nt.reheader.cram && $st view --no-PG -h -T $ref nt.reheader.cram"
done
"$st" view -C -T "$ref" -O cram,version=2.1 -o reheader.21.cram p1.bam
check "reheader -i CRAM 2.1" "cp reheader.21.cram st.reheader.21.cram && $st reheader -P -i reheader.sam st.reheader.21.cram && $st view --no-PG -h -T $ref st.reheader.21.cram" "cp reheader.21.cram nt.reheader.21.cram && $nt reheader -P -i reheader.sam nt.reheader.21.cram && $st view --no-PG -h -T $ref nt.reheader.21.cram"

# cat: @RG lines of later inputs are added (HTSlib appends them after @PG lines, so headers are
# compared sorted), -b lists, and CRAM: concatenation, read group renumbering, and the -q, -r,
# -p, and -f container selections
"$st" addreplacerg -r 'ID:other\tSM:x' -o cat.rg.bam p2.bam
"$st" view -C -T "$ref" -o cat.p1.cram p1.bam
"$st" view -C -T "$ref" -o cat.rg.cram cat.rg.bam
printf 'p1.bam\ncat.rg.bam\n' > cat.list
check "cat (@RG merge)" "$st cat --no-PG p1.bam cat.rg.bam | $st view --no-PG -h - | sort" "$nt cat --no-PG p1.bam cat.rg.bam | $st view --no-PG -h - | sort"
check "cat -b -h" "$st cat --no-PG -b cat.list -h cat.rg.bam p1.bam | $st view --no-PG -h - | sort" "$nt cat --no-PG -b cat.list -h cat.rg.bam p1.bam | $st view --no-PG -h - | sort"
check "cat CRAM (@RG renumbered)" "$st cat --no-PG cat.p1.cram cat.rg.cram | $st view --no-PG -h -T $ref - | sort" "$nt cat --no-PG cat.p1.cram cat.rg.cram | $st view --no-PG -h -T $ref - | sort"
for v in 3.0 3.1; do
    cram=$corpus/short/input.$v.cram
    check "cat CRAM $v" "$st cat --no-PG $cram $cram | $st view --no-PG -h -T $ref -" "$nt cat --no-PG $cram $cram | $st view --no-PG -h -T $ref -"
    check "cat -q CRAM $v" "$st cat -q $cram" "$nt cat -q $cram"
    for r in chrS:1-100000 chrS:500000-600000 '*'; do
        check "cat -q -r $r CRAM $v" "$st cat -q -r '$r' $cram" "$nt cat -q -r '$r' $cram"
        check "cat -r $r CRAM $v" "$st cat --no-PG -r '$r' $cram | $st view --no-PG -h -T $ref -" "$nt cat --no-PG -r '$r' $cram | $st view --no-PG -h -T $ref -"
        check "cat -f -r $r CRAM $v" "$st cat --no-PG -f -r '$r' $cram | $st view --no-PG -h -T $ref -" "$nt cat --no-PG -f -r '$r' $cram | $st view --no-PG -h -T $ref -"
        check "cat -r $r -p 2/3 CRAM $v" "$st cat --no-PG -r '$r' -p 2/3 $cram | $st view --no-PG -h -T $ref -" "$nt cat --no-PG -r '$r' -p 2/3 $cram | $st view --no-PG -h -T $ref -"
    done
    check "cat -r #:3-5 CRAM $v" "$st cat --no-PG -r '#:3-5' $cram | $st view --no-PG -h -T $ref -" "$nt cat --no-PG -r '#:3-5' $cram | $st view --no-PG -h -T $ref -"
done
long_cram=$corpus/long/input.3.0.cram
for r in $(sed -n '1p;2p' "$corpus/long/queries.100.regions"); do
    check "cat -r $r CRAM long" "$st cat --no-PG -r $r $long_cram | $st view --no-PG -h -T $long_ref -" "$nt cat --no-PG -r $r --reference $long_ref $long_cram | $st view --no-PG -h -T $long_ref -"
done

# reads
check "fastq" "$st fastq short.name.bam" "$nt fastq short.name.bam"
check "fastq -n -t" "$st fastq -n -t $tagged" "$nt fastq -n -t $tagged"
check "fasta" "$st fasta short.name.bam" "$nt fasta short.name.bam"

# fastq/fasta options: tag filters, OQ, soft-clip removal, routing, CASAVA and index reads, UMIs.
# samtools reads index reads past the end of the barcode tag ("i*i*" after a barcode without a
# separator), so only formats that stay within it are compared.
{
    printf '@HD\tVN:1.6\tSO:queryname\n@SQ\tSN:c1\tLN:1000\n'
    printf 'a1\t99\tc1\t1\t60\t2S6M2S\t=\t20\t30\tACGTACGTAC\tABCDEFGHIJ\tRG:Z:g1\tXI:i:5\tXF:f:0.5\tOQ:Z:##########\ts0:Z:old\tBC:Z:AAAA-CCCC\tQT:Z:12345678.\n'
    printf 'a1\t147\tc1\t20\t60\t3S5M1I1S\t=\t1\t-30\tTTTTGGGGCA\tJIHGFEDCBA\tRG:Z:g2\tXI:i:5\tOQ:Z:!!!\tBC:Z:AAAA-CCCC\n'
    printf 'a2\t77\t*\t0\t0\t*\t*\t0\t0\tACGTA\t*\tRG:Z:g1\tXA:A:x\tXH:H:1AFF\n'
    printf 'a2\t141\t*\t0\t0\t*\t*\t0\t0\tACGTA\t*\tRG:Z:g3\tXI:i:-7\n'
    printf 'a3\t0\tc1\t5\t60\t10S\t*\t0\t0\tACGTACGTAC\tABCDEFGHIJ\tRG:Z:g1\tXI:i:200\n'
    printf 'a4\t16\tc1\t5\t60\t4M2D4M2S\t*\t0\t0\tACGTACGTAC\tABCDEFGHIJ\tRG:Z:g2\tXI:i:5\tXF:f:0.5\n'
    printf 'a5\t1\t*\t0\t0\t*\t*\t0\t0\t*\t*\tRG:Z:g1\n'
    printf 'a6\t129\tc1\t5\t60\t4M\t*\t0\t0\tACGT\tABCD\tXI:i:-300\tXF:f:1e-07\n'
    printf 'a7\t65\tc1\t5\t60\t4M\t*\t0\t0\tACGT\tABCD\tXI:i:5\tXF:f:3.25\tBC:Z:GGTT\tQT:Z:9999\n'
    printf 'a7\t129\tc1\t5\t60\t4M\t*\t0\t0\tACGT\tABCD\tXF:f:0.5\tBC:Z:TT\n'
    printf 'a8\t1105\tc1\t5\t60\t4M\t*\t0\t0\tACGT\tABCD\tXB:B:f,1.5,0.25\tXC:B:S,1,65535\n'
    printf 'p3\t589\t*\t0\t0\t*\t*\t0\t0\tACGT\tIIII\tBC:Z:acgtNN\tQT:Z:######\n'
    printf 'p4\t4\t*\t0\t0\t*\t*\t0\t0\tACGT\tIIII\tBC:Z:1CG\n'
    printf 'p6#ACG\t4\t*\t0\t0\t*\t*\t0\t0\tACGT\tIIII\tOX:Z:TT\tRX:Z:aa-cc\n'
    printf 'p7:x\t4\t*\t0\t0\t*\t*\t0\t0\tACGT\tIIII\tRX:Z:GG\n'
} > fq-opts.sam
printf 'g1\ng3\r\n' > fq-rg.txt
for opts in "-T '*'" "-t -T XI,XF" "-d RG:g1 -d RG:g2" "-d XF:0.500000 -T XF" "-d XA" "-D RG:fq-rg.txt" \
    "-O -T '*'" "--no-sc -T '*'" "--no-sc --no-sc-bkp --sc-aux xx -O" "--no-sc --sc-aux xx -T xx" \
    "--rf 0x10 -F 0" "--incl-flags 0x80 -G 0xc0" "-U" "-U --UMI-tag RX -N" \
    "-i --index-format 'i4n1i*'" "-i --barcode-tag RG --if 'n1i*'"; do
    check "fastq $opts (synthetic)" "$st fastq $opts fq-opts.sam" "$nt fastq $opts fq-opts.sam"
done
check "fasta -T * --no-sc (synthetic)" "$st fasta -T '*' --no-sc fq-opts.sam" "$nt fasta -T '*' --no-sc fq-opts.sam"
check "fastq -1 -2 -0 -s (synthetic)" \
    "$st fastq -1 st.fq1 -2 st.fq2 -0 st.fq0 -s st.fqs fq-opts.sam; cat st.fq1 st.fq2 st.fq0 st.fqs" \
    "$nt fastq -1 nt.fq1 -2 nt.fq2 -0 nt.fq0 -s nt.fqs fq-opts.sam; cat nt.fq1 nt.fq2 nt.fq0 nt.fqs"
check "fastq --i1 --i2 -i -T BC,QT (synthetic)" \
    "$st fastq -i --index-format 'i4n1i*' --i1 st.i1.fq --i2 st.i2.fq -T BC,QT fq-opts.sam; cat st.i1.fq st.i2.fq" \
    "$nt fastq -i --index-format 'i4n1i*' --i1 nt.i1.fq --i2 nt.i2.fq -T BC,QT fq-opts.sam; cat nt.i1.fq nt.i2.fq"
check "fastq --i1 -1 shared, --quality-tag (synthetic)" \
    "$st fastq --IF 'i*' --quality-tag OQ --I1 st.ix.fq -1 st.ix.fq -v 30 fq-opts.sam; cat st.ix.fq" \
    "$nt fastq --IF 'i*' --quality-tag OQ --I1 nt.ix.fq -1 nt.ix.fq -v 30 fq-opts.sam; cat nt.ix.fq"
check "fastq -T * (long)" "$st fastq -T '*' $long" "$nt fastq -T '*' $long"
check "fastq --no-sc -O -t (tagged)" "$st fastq --no-sc -O -t $tagged" "$nt fastq --no-sc -O -t $tagged"
check "fastq -i -d RG -U (short)" "$st fastq -t -i --index-format 'i*' -d RG:synthetic -U short.name.bam" "$nt fastq -t -i --index-format 'i*' -d RG:synthetic -U short.name.bam"

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

# index: -m (CSI with a minimum shift; 0 means BAI), -M (several inputs), and the legacy
# `index in.bam out.index`, compared through idxstats and region queries with each index
ln -sf "$short" index.a.bam
ln -sf "$tagged" index.b.bam
for opts in "-m 12" "-m 17" "-c" "-m 0" "-m 12 -b" "-b -m 13"; do
    check "index $opts" "$st index $opts -o st.index.idx index.a.bam && head -c 3 st.index.idx | od -An -c && $st idxstats -X index.a.bam st.index.idx && $st view --no-PG -X index.a.bam st.index.idx $regions" "$nt index $opts -o nt.index.idx index.a.bam && head -c 3 nt.index.idx | od -An -c && $st idxstats -X index.a.bam nt.index.idx && $st view --no-PG -X index.a.bam nt.index.idx $regions"
done
cp "$short" st.index.a.bam && cp "$tagged" st.index.b.bam
cp "$short" nt.index.a.bam && cp "$tagged" nt.index.b.bam
check "index -M" "$st index -M -m 13 st.index.a.bam st.index.b.bam && $st idxstats -X st.index.b.bam st.index.b.bam.csi && $st view --no-PG -X st.index.a.bam st.index.a.bam.csi $regions" "$nt index -M -m 13 nt.index.a.bam nt.index.b.bam && $st idxstats -X nt.index.b.bam nt.index.b.bam.csi && $st view --no-PG -X nt.index.a.bam nt.index.a.bam.csi $regions"
check "index without -M" "$st index index.a.bam index.b.bam; echo \$?" "$nt index index.a.bam index.b.bam 2>/dev/null; echo \$?"
check "index in.bam out.index" "$st index index.a.bam st.index.legacy && $st idxstats -X index.a.bam st.index.legacy" "$nt index index.a.bam nt.index.legacy && $st idxstats -X index.a.bam nt.index.legacy"

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

# dict -A (chr aliases) and -l (alternate loci, plain or gzipped)
printf '>chr1 desc\nACGT\n>chrM\nAC\n>MT\nGG\n>X\nTT\n>chr_alt1\nAA\n>M\nCC\n' > alias.fa
printf '@SQ\tSN:chr1\tLN:4\nchr_alt1\t0\tchr1\t1\t60\t2M\t*\t0\t0\tAA\t*\nX\n\nMT\t4\r\n' > alias.alt
gzip -c alias.alt > alias.alt.gz
check "dict -A" "$st dict -A -u U alias.fa" "$nt dict -A -u U alias.fa"
check "dict -l" "$st dict --alt alias.alt -H -u U alias.fa" "$nt dict --alt alias.alt -H -u U alias.fa"
check "dict -l (gzip) --alias" "$st dict -l alias.alt.gz --alias -u U alias.fa" "$nt dict -l alias.alt.gz --alias -u U alias.fa"
check "dict -A frag.fa" "$st dict -A -u U frag.fa" "$nt dict -A -u U frag.fa"
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

# bedcov -X: index files follow the inputs
ln -sf "$short" bedcov.bam
check "bedcov -X" "$st bedcov -H -c -X q.bed bedcov.bam $tagged $corpus/short/input.bam.csi $corpus/tagged/input.bam.bai" "$nt bedcov -H -c -X q.bed bedcov.bam $tagged $corpus/short/input.bam.csi $corpus/tagged/input.bam.bai"

# depth: -s (overlapping mates), -l (soft clips excluded), deletions' qualities with -J -q, and
# spliced reads' skips
for opts in "-s" "-s -J -q 30 -Q 20" "-l 100" "-J -q 30 -Q 20"; do
    for input in $tagged $long; do
        check "depth $opts $(basename "$(dirname "$input")")" "$st depth $opts $input" "$nt depth $opts $input"
    done
done

# depth: input lists (-f), index files after the inputs (-X, used with -r), --incl-flags, and
# the ignored -d/-m. With -X, samtools labels the -H columns with the index names; ntools uses
# the inputs'.
printf "%s\n%s\n" "$short" "$tagged" > depth.list
ln -sf "$short" depth.bam
check "depth -f" "$st depth -H -f depth.list" "$nt depth -H -f depth.list"
check "depth -f -r" "$st depth -a -f depth.list -r chrS:100000-110000" "$nt depth -a -f depth.list -r chrS:100000-110000"
check "depth -X" "$st depth -X -r chrS:100000-200000 depth.bam $corpus/short/input.bam.csi" "$nt depth -X -r chrS:100000-200000 depth.bam $corpus/short/input.bam.csi"
check "depth -X two inputs" "$st depth -a -X -r chrS:5000-60000 depth.bam $tagged $corpus/short/input.bam.bai $corpus/tagged/input.bam.bai" "$nt depth -a -X -r chrS:5000-60000 depth.bam $tagged $corpus/short/input.bam.bai $corpus/tagged/input.bam.bai"
check "depth --incl-flags" "$st depth --incl-flags READ1,REVERSE $tagged" "$nt depth --incl-flags READ1,REVERSE $tagged"
check "depth -d -m (ignored)" "$st depth -d 5 -m 10 $short" "$nt depth -d 5 -m 10 $short"

# coverage -q is the minimum mapping quality, -Q the minimum base quality
check "coverage -q 20 -Q 10" "$st coverage -q 20 -Q 10 $short" "$nt coverage -q 20 -Q 10 $short"
check "coverage -l 120 -q 10 long" "$st coverage -l 120 -q 10 $long" "$nt coverage -l 120 -q 10 $long"

# coverage: histograms (-m, -D, -A, -w; the bin count follows the terminal unless -w or -o
# is given), the pileup's depth limit (-d), --min-depth, input lists (-b), several inputs, and
# the long-only options. samtools reads CRAM base qualities as 255 unless -Q needs them, so that
# column is left out for CRAM.
for opts in "-m" "-m -o /dev/stdout" "-D -w 70" "-A -w 25" "-A -D" "-m -r chrS:100000-300000" \
    "-w 3 -r chrS:100-101" "-d 5" "-d 2 -D -w 30" "-d 0" "-d -1" "--min-depth 20" \
    "--min-depth 15 -m -w 30" "--incl-flags READ2 --excl-flags DUP,SECONDARY" "--rf 0x40 --ff 0" \
    "-H -r chrS:1000-9000" "-r chrS:7999000-8000500" "--min-mq 30 --min-bq 20"; do
    check "coverage $opts" "$st coverage $opts $short" "$nt coverage $opts $short"
done
for opts in "-d 3" "-m -d 4 -w 90"; do
    check "coverage $opts tagged" "$st coverage $opts $tagged" "$nt coverage $opts $tagged"
done
printf "%s\n\n%s  \n" "$short" "$tagged" > coverage.list
check "coverage -b" "$st coverage -b coverage.list" "$nt coverage -b coverage.list"
check "coverage -b -D" "$st coverage -D -w 50 -b coverage.list" "$nt coverage -D -w 50 -b coverage.list"
check "coverage two inputs -m" "$st coverage -m -w 40 $short $tagged" "$nt coverage -m -w 40 $short $tagged"
check "coverage CRAM --reference" "$st coverage -d 4 --reference $ref $corpus/short/input.3.1.cram | cut -f1-7,9" "$nt coverage -d 4 --reference $ref $corpus/short/input.3.1.cram | cut -f1-7,9"
# reads that start together beyond the limit, empty reads, deletions, and skips
printf "@HD\tVN:1.6\tSO:coordinate\n@SQ\tSN:a\tLN:40\n@SQ\tSN:b\tLN:30\n" > coverage.cap.sam
for r in "5 5S AAAAA" "5 3M AAA" "5 3M AAA" "5 2M1D1M AAA" "5 4I AAAA" "6 1M A" "6 6M AAAAAA" \
    "6 2M3N2M AAAA" "8 3M AAA" "8 3M AAA" "8 1M2D1M AA" "8 3M AAA" "9 2S2M AAAA"; do
    read -r pos cigar seq <<<"$r"
    printf "r\t0\ta\t%s\t30\t%s\t*\t0\t0\t%s\t%s\n" "$pos" "$cigar" "$seq" "$(echo "$seq" | tr A I)" >> coverage.cap.sam
done
printf "s\t0\tb\t3\t30\t2M\t*\t0\t0\tAA\tII\n" >> coverage.cap.sam
for d in 1 2 3; do
    check "coverage -d $d (reads left out)" "$st coverage -d $d coverage.cap.sam; $st coverage -d $d -D -w 40 coverage.cap.sam" "$nt coverage -d $d coverage.cap.sam; $nt coverage -d $d -D -w 40 coverage.cap.sam"
done

# mpileup refuses inputs whose reference sequences differ from the first's
check "mpileup mismatched references" "echo 1" "$nt mpileup mp.short.bam md.bam > /dev/null 2>&1; echo \$?"

echo "Slowest checks (seconds, samtools run or cached):"
sort -rn "$timings" | head -15
awk -F'\t' '{ t += $1; if ($2 == "cache") c++ } END { printf "Total %.0f s over %d checks, %d from cache\n", t, NR, c }' "$timings"
echo "$failures difference(s)"
[ "$failures" -eq 0 ]
