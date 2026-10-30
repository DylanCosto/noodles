#!/usr/bin/env bash
# Compares ntools with samtools on the edge-case corpus written by make-edge-corpus.sh.
#
# Usage: compare-edge.sh NTOOLS SAMTOOLS EDGE_DIR [WORK_DIR]
#
# Every check runs a samtools command and an ntools command and compares their standard output
# (binary outputs are decoded to SAM by samtools first), their exit status, and whether each
# wrote anything to standard error. A check prints OK, or DIFF with the parts that differ
# (out, rc=ST/NT, err=ST/NT with e for empty and x for not). An ntools panic is reported as
# PANIC whatever samtools did. A difference that is known and intended (see `known`) is reported
# as KNOWN with its reason instead. The outputs of differing checks are kept in WORK_DIR/diffs/N/
# (cmd, st.out, nt.out, st.err, nt.err). Exits nonzero if any check differs (other than KNOWN).
#
# Checks named "@1 vs @4" compare ntools with itself, with one and four threads.
set -uo pipefail

nt=$1
st=$2
edge=$(cd "$3" && pwd)
work=${4:-$(mktemp -d)}

mkdir -p "$work"
work=$(cd "$work" && pwd)
cd "$work" || exit 2
rm -rf diffs
mkdir -p diffs refpath

# Keep HTSlib from looking up CRAM references by MD5 on the network.
export REF_PATH=$work/refpath REF_CACHE=$work/refpath

# BAQ on the 70,000-operation CIGAR asks for a 58 GB matrix, which HTSlib's calloc gets lazily
# from the kernel and then fills; cap memory so such allocations fail instead.
ulimit -v 6000000

failures=0
known_count=0
checks=0

# Prints why a difference in check $1 is expected, or nothing. These are the deliberate or
# harmless differences of the triage (samtools bugs that ntools fixes, inputs one tool accepts
# and the other rejects, and messages); panics are never expected.
known() {
    case $1 in
        *" 2.1.cram"* | *"2.1 "*) echo "ntools reads only CRAM 3; it rejects CRAM 2.1 with an error" ;;
        *duprg*) echo "duplicate @RG/@PG IDs: samtools warns and keeps both lines, ntools rejects the header" ;;
        *empty_none*) echo "a 0-byte SAM file: ntools reads it as an empty header, samtools fails" ;;
        "idxstats "*" sam") echo "unsorted SAM: samtools idxstats requires position order" ;;
        *" ends 3.0.cram" | *" ends 3.1.cram" | *" ends embed.cram" | *" ends 3.1.cram ("*)
            echo "HTSlib fails on CRAM slices with reads past the reference end (README)" ;;
        "view regions ends 3.1.cram" | "bedcov ends 3.1.cram" | "coverage -m ends 3.1.cram")
            echo "HTSlib fails on CRAM slices with reads past the reference end (README)" ;;
        "coverage"*".cram" | "coverage"*" noseq bam" | "coverage"*" longcigar bam")
            echo "coverage meanbaseq leaves out missing qualities and decodes CRAM qualities (README)" ;;
        "reference cigarops 3.1.cram") echo "reference advances over N operations (README)" ;;
        *" flags 3.0.cram" | *" flags 3.1.cram" | *" flags embed.cram" | "view regions flags 3.1.cram")
            echo "HTSlib's CRAM decoder sets PAIRED on a detached read with only MREVERSE (52 -> 53)" ;;
        "reheader "*" cram") echo "samtools warns that it cannot load the reference (message only)" ;;
        "addreplacerg"*" headers"* | "addreplacerg"*" badpp"*)
            echo "a new @RG line goes before @PG lines (README)" ;;
        "view (with @PG) badpp sam") echo "samtools warns about a PP link to a missing @PG (message only)" ;;
        "cat badpp headers") echo "header lines are grouped by type (README)" ;;
        "ampliconclip --hard-clip --both-ends mapq "* | "ampliconclip --hard-clip --both-ends noseq "*)
            echo "samtools ampliconclip crashes (SIGSEGV) hard-clipping reads without a sequence" ;;
        "fixmate"*" ends") echo "samtools soft clips deletion bases past the reference end (3M3D9S for 5 bases)" ;;
        "calmd -e -r -A beyond "*) echo "samtools BAQ reads uninitialized memory for a read past the reference end" ;;
        "cram-size -v -e empty_"*) echo "samtools cram-size fails on a CRAM file without data containers" ;;
        "view -C "*" "*) ;;
        "view -C "*)
            echo "CRAM written by ntools has no @SQ UR field, and decoding adds NM before MD (md_nm.rs)" ;;
        "import edge.fa") echo "samtools import fails on an empty FASTA record or a missing final newline" ;;
    esac
}

check() {
    local name=$1 expected=$2 actual=$3
    local st_rc nt_rc st_err nt_err what="" panic=""
    checks=$((checks + 1))
    (eval "$expected") >st.out 2>st.err </dev/null &
    local pid=$!
    (eval "$actual") >nt.out 2>nt.err </dev/null
    nt_rc=$?
    wait "$pid"
    st_rc=$?
    st_err=e
    nt_err=e
    [ -s st.err ] && st_err=x
    [ -s nt.err ] && nt_err=x
    cmp -s st.out nt.out || what="$what out"
    [ "$st_rc" = "$nt_rc" ] || what="$what rc=$st_rc/$nt_rc"
    [ "$st_err" = "$nt_err" ] || what="$what err=$st_err/$nt_err"
    grep -q 'panicked' nt.err && panic=PANIC
    if [ -z "$what" ] && [ -z "$panic" ]; then
        echo "OK   $name"
        return
    fi

    local d=diffs/$checks reason
    mkdir -p "$d"
    printf '%s\nST: %s\nNT: %s\n' "$name" "$expected" "$actual" >"$d/cmd"
    mv st.out st.err nt.out nt.err "$d/"
    reason=$(known "$name")

    if [ -z "$panic" ] && [ -n "$reason" ]; then
        echo "KNOWN $name ($checks:$what): $reason"
        known_count=$((known_count + 1))
    else
        echo "DIFF $name ($checks:${panic:+ $panic}$what)"
        failures=$((failures + 1))
    fi
}

# Decodes alignment output to SAM with samtools.
sam() {
    "$st" view --no-PG -h --reference "$ref" "$@"
}

cases=$(<"$edge/cases")
sams=$(cd "$edge" && ls ./*.sam | sed 's|^\./||; s|\.sam$||')
formats="sam bam 3.0.cram 3.1.cram 2.1.cram embed.cram"
ref=$edge/ref.fa

# Sets $in to the path of case $1 in format $2 and $ref to its reference.
input() {
    if [ "$2" = sam ]; then in=$edge/$1.sam; else in=$edge/$1.$2; fi
    ref=$edge/ref.fa
    [ "$1" = manysq ] && ref=$edge/manysq.fa
    return 0
}

stats_strip="grep -v '^# This file was produced\|^# The command line'"

# view: decoding every format, SAM inputs samtools cannot convert, and writing BAM and CRAM.
for c in $sams; do
    input "$c" sam
    check "view $c sam" "$st view --no-PG -h $in" "$nt view --no-PG -h $in"
    check "view -b $c sam" "$st view --no-PG -b -o st.bam $in && sam st.bam" "$nt view --no-PG -b -o nt.bam $in && sam nt.bam"
    # With a @PG line of its own, linked to the leaves of the input's @PG chains.
    check "view (with @PG) $c sam" "$st view -h $in | grep -v '^@PG.ID:samtools'" "$nt view -h $in | grep -v '^@PG.ID:ntools'"
done
for c in $cases; do
    for f in $formats; do
        input "$c" "$f"
        check "view $c $f" "$st view --no-PG -h --reference $ref $in" "$nt view --no-PG -h --reference $ref $in"
    done
    input "$c" bam
    check "view -@1 vs -@4 $c bam" "$nt view --no-PG -h -@ 1 $in" "$nt view --no-PG -h -@ 4 $in"
    check "view -c $c" "$st view -c $in" "$nt view -c $in"
    check "view -c -f/-F $c" "$st view -c -f 1 -F 0x904 -q 1 $in" "$nt view -c -f 1 -F 0x904 -q 1 $in"
    check "view regions $c" "$st view --no-PG $in chr1 chr2:100-200 chrL:39990-40000 s399 '*'" "$nt view --no-PG $in chr1 chr2:100-200 chrL:39990-40000 s399 '*'"
    input "$c" 3.1.cram
    check "view regions $c 3.1.cram" "$st view --no-PG --reference $ref $in chr1 chr2:100-200 chrL:39990-40000 s399 '*'" "$nt view --no-PG --reference $ref $in chr1 chr2:100-200 chrL:39990-40000 s399 '*'"
    check "view -@1 vs -@4 $c 3.1.cram" "$nt view --no-PG -h -@ 1 --reference $ref $in" "$nt view --no-PG -h -@ 4 --reference $ref $in"
    input "$c" bam
    check "view -C $c" "$st view --no-PG -C -T $ref -o st.cram $in && sam st.cram" "$nt view --no-PG -C -T $ref -o nt.cram $in && sam nt.cram"
    # Both CRAM files decoded by ntools, so only the writers differ.
    check "view -C round trip (fields 1-11) $c" "$st view --no-PG -C -T $ref -o st.cram $in && $nt view --no-PG -T $ref st.cram | cut -f1-11" "$nt view --no-PG -C -T $ref -o nt.cram $in && $nt view --no-PG -T $ref nt.cram | cut -f1-11"
    check "view -C -@1 vs -@4 $c" "$nt view --no-PG -C -@ 1 -T $ref $in | sam -" "$nt view --no-PG -C -@ 4 -T $ref $in | sam -"
    check "view -b -@1 vs -@4 $c" "$nt view --no-PG -b -@ 1 $in | sam -" "$nt view --no-PG -b -@ 4 $in | sam -"
    check "view -e $c" "$st view --no-PG -e 'length(seq) > 5 && mapq < 255 || [XA]' $in" "$nt view --no-PG -e 'length(seq) > 5 && mapq < 255 || [XA]' $in"
done

# flagstat, idxstats, quickcheck, samples, head, checksum: every format.
for c in $cases; do
    for f in sam bam 3.1.cram 2.1.cram; do
        input "$c" "$f"
        check "flagstat $c $f" "$st flagstat $in" "$nt flagstat $in"
        check "idxstats $c $f" "$st idxstats $in" "$nt idxstats $in"
        check "quickcheck $c $f" "$st quickcheck -v $in" "$nt quickcheck -v $in"
        check "samples $c $f" "$st samples -h -i $in" "$nt samples -h -i $in"
        check "head $c $f" "$st head --reference $ref -n 3 $in" "$nt head --reference $ref -n 3 $in"
        check "checksum $c $f" "$st checksum --reference $ref -a $in" "$nt checksum --reference $ref -a $in"
    done
    input "$c" bam
    check "flagstat -O json $c" "$st flagstat -O json $in" "$nt flagstat -O json $in"
    check "flagstat -O tsv $c" "$st flagstat -O tsv $in" "$nt flagstat -O tsv $in"
    check "flagstat -@1 vs -@4 $c" "$nt flagstat -@ 1 $in" "$nt flagstat -@ 4 $in"
    check "checksum -@1 vs -@4 $c" "$nt checksum -a -@ 1 $in" "$nt checksum -a -@ 4 $in"
done
for c in $sams; do
    input "$c" sam
    check "flagstat $c raw sam" "$st flagstat $in" "$nt flagstat $in"
    check "quickcheck $c raw sam" "$st quickcheck -v $in" "$nt quickcheck -v $in"
done
for c in $cases; do
    input "$c" bam
    head -c 100 "$in" >trunc.$c.bam
    check "quickcheck truncated $c" "$st quickcheck -v trunc.$c.bam" "$nt quickcheck -v trunc.$c.bam"
    check "view truncated $c" "$st view --no-PG trunc.$c.bam" "$nt view --no-PG trunc.$c.bam"
done

# index
for c in $cases; do
    input "$c" bam
    cp "$in" idx.$c.bam
    check "index $c" "$st index -o st.bai idx.$c.bam && $st idxstats -X idx.$c.bam st.bai && $st view --no-PG -X idx.$c.bam st.bai chr1 chr2:100-200 chrL:39990-40000 s399 '*'" \
        "$nt index -o nt.bai idx.$c.bam && $st idxstats -X idx.$c.bam nt.bai && $st view --no-PG -X idx.$c.bam nt.bai chr1 chr2:100-200 chrL:39990-40000 s399 '*'"
    check "index -c $c" "$st index -c -o st.csi idx.$c.bam && $st idxstats -X idx.$c.bam st.csi && $st view --no-PG -X idx.$c.bam st.csi chr1 chr2:100-200 '*'" \
        "$nt index -c -o nt.csi idx.$c.bam && $st idxstats -X idx.$c.bam nt.csi && $st view --no-PG -X idx.$c.bam nt.csi chr1 chr2:100-200 '*'"
    check "index -@1 vs -@4 $c" "$nt index -@ 1 -o nt1.bai idx.$c.bam && cat nt1.bai" "$nt index -@ 4 -o nt4.bai idx.$c.bam && cat nt4.bai"
    for v in 3.0 3.1 2.1; do
        input "$c" $v.cram
        cp "$in" idx.$c.$v.cram
        check "index $c $v.cram" "$st index -o st.crai idx.$c.$v.cram && gzip -dc st.crai" "$nt index -o nt.crai idx.$c.$v.cram && gzip -dc nt.crai"
    done
done
check "index unsorted" "$st index $edge/mates.name.bam -o st.bai" "$nt index $edge/mates.name.bam -o nt.bai"

# sort and merge
for c in $cases; do
    for f in sam bam 3.1.cram; do
        input "$c" "$f"
        check "sort $c $f" "$st sort --no-PG --reference $ref -O sam $in" "$nt sort --no-PG --reference $ref -O sam $in"
        check "sort -n $c $f" "$st sort --no-PG --reference $ref -n -O sam $in" "$nt sort --no-PG --reference $ref -n -O sam $in"
    done
    input "$c" bam
    check "sort -N $c" "$st sort --no-PG -N -O sam $in" "$nt sort --no-PG -N -O sam $in"
    check "sort -t RG $c" "$st sort --no-PG -t RG -O sam $in" "$nt sort --no-PG -t RG -O sam $in"
    check "sort -@1 vs -@4 $c" "$nt sort --no-PG -@ 1 -O sam $in" "$nt sort --no-PG -@ 4 -O sam $in"
    check "merge $c $c" "$st merge --no-PG -s 1 -O sam -o - $in $in" "$nt merge --no-PG -s 1 -O sam -o - $in $in"
    check "merge $c mates" "$st merge --no-PG -s 1 -O sam -o - $in $edge/mates.bam" "$nt merge --no-PG -s 1 -O sam -o - $in $edge/mates.bam"
    check "merge -n $c" "$st merge --no-PG -n -s 1 -O sam -o - $edge/$c.name.bam $edge/mates.name.bam" "$nt merge --no-PG -n -s 1 -O sam -o - $edge/$c.name.bam $edge/mates.name.bam"
    check "merge -R $c" "$st merge --no-PG -s 1 -R chr1:1-150 -O sam -o - $in $edge/mates.bam" "$nt merge --no-PG -s 1 -R chr1:1-150 -O sam -o - $in $edge/mates.bam"
done
check "merge empty_sq manysq" "$st merge --no-PG -O sam -o - $edge/empty_sq.bam $edge/manysq.bam" "$nt merge --no-PG -O sam -o - $edge/empty_sq.bam $edge/manysq.bam"

# fastq and fasta
for c in $cases; do
    for f in sam bam 3.1.cram; do
        input "$c" "$f"
        check "fastq $c $f" "$st fastq --reference $ref $in" "$nt fastq --reference $ref $in"
        check "fasta $c $f" "$st fasta --reference $ref $in" "$nt fasta --reference $ref $in"
    done
    input "$c" name.bam
    check "fastq -1 -2 -s -0 $c" "$st fastq -N -T '*' -1 st.1.fq -2 st.2.fq -s st.s.fq -0 st.0.fq $in; cat st.1.fq st.2.fq st.s.fq st.0.fq" "$nt fastq -N -T '*' -1 nt.1.fq -2 nt.2.fq -s nt.s.fq -0 nt.0.fq $in; cat nt.1.fq nt.2.fq nt.s.fq nt.0.fq"
    check "fastq -@1 vs -@4 $c" "$nt fastq -@ 1 -T '*' $in" "$nt fastq -@ 4 -T '*' $in"
done

# reheader, cat
for c in $cases; do
    input "$c" bam
    "$st" view -H "$edge/headers.bam" >newhdr.sam
    check "reheader $c" "$st reheader -P newhdr.sam $in | sam -" "$nt reheader -P newhdr.sam $in | sam -"
    check "cat $c" "$st cat --no-PG $in $in | sam -" "$nt cat --no-PG $in $in | sam -"
    check "cat $c headers" "$st cat --no-PG $in $edge/headers.bam | sam -" "$nt cat --no-PG $in $edge/headers.bam | sam -"
    input "$c" 3.1.cram
    check "reheader $c cram" "$st reheader -P newhdr.sam $in | sam -" "$nt reheader -P newhdr.sam $in | sam -"
    check "cat $c cram" "$st cat --no-PG $in $in | sam -" "$nt cat --no-PG $in $in | sam -"
done

# depth, coverage, bedcov, mpileup, stats, consensus, reference
printf 'chr1\t0\t200\nchr2\t100\t120\nchrL\t39980\t40000\ns5\t0\t3\nchr1\t190\t250\n' >edge.bed
for c in $cases; do
    for f in bam 3.1.cram; do
        input "$c" "$f"
        check "depth $c $f" "$st depth --reference $ref $in" "$nt depth --reference $ref $in"
        check "depth -a -s $c $f" "$st depth --reference $ref -a -s $in" "$nt depth --reference $ref -a -s $in"
        check "depth -aa -J -Q 0 -q 0 $c $f" "$st depth --reference $ref -aa -J -g 0x900 $in" "$nt depth --reference $ref -aa -J -g 0x900 $in"
        check "coverage $c $f" "$st coverage --reference $ref $in" "$nt coverage --reference $ref $in"
        check "coverage -m $c $f" "$st coverage --reference $ref -m -w 40 $in" "$nt coverage --reference $ref -m -w 40 $in"
        check "bedcov $c $f" "$st bedcov --reference $ref -c -j edge.bed $in" "$nt bedcov --reference $ref -c -j edge.bed $in"
        check "mpileup $c $f" "$st mpileup --reference $ref -f $ref $in" "$nt mpileup --reference $ref -f $ref $in"
        check "mpileup -B -a -Q 0 --ff 0 -x $c $f" "$st mpileup --reference $ref -B -a -Q 0 --ff 0 -x --output-QNAME --output-extra FLAG,RG,MD $in" "$nt mpileup --reference $ref -B -a -Q 0 --ff 0 -x --output-QNAME --output-extra FLAG,RG,MD $in"
        check "stats $c $f" "$st stats --reference $ref $in | $stats_strip" "$nt stats --reference $ref $in | $stats_strip"
        check "stats -r $c $f" "$st stats -r $ref --reference $ref $in | $stats_strip" "$nt stats -r $ref --reference $ref $in | $stats_strip"
        check "consensus $c $f" "$st consensus --reference $ref -a $in" "$nt consensus --reference $ref -a $in"
        check "consensus -f fastq -m simple $c $f" "$st consensus --reference $ref -f fastq -m simple --show-ins yes --show-del yes $in" "$nt consensus --reference $ref -f fastq -m simple --show-ins yes --show-del yes $in"
        check "reference $c $f" "$st reference -q $in" "$nt reference -q $in"
    done
    input "$c" bam
    check "depth -@1 vs -@4 $c" "$nt depth -a -@ 1 $in" "$nt depth -a -@ 4 $in"
    check "coverage -@1 vs -@4 $c" "$nt coverage -@ 1 $in" "$nt coverage -@ 4 $in"
    check "mpileup -@1 vs -@4 $c" "$nt mpileup -f $ref -@ 1 $in" "$nt mpileup -f $ref -@ 4 $in"
    check "stats -@1 vs -@4 $c" "$nt stats -@ 1 $in | $stats_strip" "$nt stats -@ 4 $in | $stats_strip"
    check "consensus -@1 vs -@4 $c" "$nt consensus -a -@ 1 $in" "$nt consensus -a -@ 4 $in"
    check "bedcov -@1 vs -@4 $c" "$nt bedcov -@ 1 edge.bed $in" "$nt bedcov -@ 4 edge.bed $in"
    input "$c" embed.cram
    check "reference -e $c" "$st reference -q -e $in" "$nt reference -q -e $in"
done

# calmd, fixmate, markdup, rmdup, collate, split, addreplacerg, ampliconclip, reset
printf 'chr1\t0\t12\tp1_LEFT\t1\t+\nchr1\t180\t200\tp1_RIGHT\t1\t-\nchr2\t0\t10\tp2_LEFT\t2\t+\nchrL\t39980\t40000\tp3_RIGHT\t3\t-\n' >primers.bed
for c in $cases; do
    for f in sam bam 3.1.cram; do
        input "$c" "$f"
        check "calmd $c $f" "$st calmd --no-PG --reference $ref $in $ref" "$nt calmd --no-PG --reference $ref $in $ref"
        check "calmd -e -r -A $c $f" "$st calmd --no-PG --reference $ref -e -r -A $in $ref" "$nt calmd --no-PG --reference $ref -e -r -A $in $ref"
        check "addreplacerg $c $f" "$st addreplacerg --no-PG --reference $ref -r ID:new -r SM:x -O sam $in" "$nt addreplacerg --no-PG --reference $ref -r ID:new -r SM:x -O sam $in"
        check "reset $c $f" "$st reset --no-PG --reference $ref -O sam $in" "$nt reset --no-PG --reference $ref -O sam $in"
        check "ampliconclip $c $f" "$st ampliconclip --no-PG --reference $ref -b primers.bed -O sam $in" "$nt ampliconclip --no-PG --reference $ref -b primers.bed -O sam $in"
        check "ampliconclip --hard-clip --both-ends $c $f" "$st ampliconclip --no-PG --reference $ref --hard-clip --both-ends --strand -b primers.bed -O sam $in" "$nt ampliconclip --no-PG --reference $ref --hard-clip --both-ends --strand -b primers.bed -O sam $in"
    done
    input "$c" bam
    check "calmd -@1 vs -@4 $c" "$nt calmd --no-PG -@ 1 $in $ref" "$nt calmd --no-PG -@ 4 $in $ref"
    check "addreplacerg -m orphan_only $c" "$st addreplacerg --no-PG -m orphan_only -r ID:o -r SM:s -O sam -o - $in" "$nt addreplacerg --no-PG -m orphan_only -r ID:o -r SM:s -O sam -o - $in"
    check "ampliconclip -@1 vs -@4 $c" "$nt ampliconclip --no-PG -@ 1 -b primers.bed -O sam $in" "$nt ampliconclip --no-PG -@ 4 -b primers.bed -O sam $in"
    check "reset -@1 vs -@4 $c" "$nt reset --no-PG -@ 1 -O sam $in" "$nt reset --no-PG -@ 4 -O sam $in"
    check "rmdup $c" "$st rmdup $in st.rmdup.bam && sam st.rmdup.bam" "$nt rmdup $in nt.rmdup.bam && sam nt.rmdup.bam"
    check "rmdup -s $c" "$st rmdup -s $in st.rmdup.bam && sam st.rmdup.bam" "$nt rmdup -s $in nt.rmdup.bam && sam nt.rmdup.bam"
    input "$c" name.bam
    check "fixmate $c" "$st fixmate --no-PG -O sam $in -" "$nt fixmate --no-PG -O sam $in -"
    check "fixmate -m -c -r -p $c" "$st fixmate --no-PG -m -c -r -p -O sam $in -" "$nt fixmate --no-PG -m -c -r -p -O sam $in -"
    check "fixmate -@1 vs -@4 $c" "$nt fixmate --no-PG -m -@ 1 -O sam $in -" "$nt fixmate --no-PG -m -@ 4 -O sam $in -"
    "$st" fixmate -m "$in" - 2>/dev/null | "$st" sort -o md.$c.bam - 2>/dev/null
    check "markdup $c" "$st markdup --no-PG -s -O sam md.$c.bam -" "$nt markdup --no-PG -s -O sam md.$c.bam -"
    check "markdup -r -S -d 100 --duplicate-count $c" "$st markdup --no-PG -r -S -d 100 --duplicate-count -O sam md.$c.bam -" "$nt markdup --no-PG -r -S -d 100 --duplicate-count -O sam md.$c.bam -"
    check "markdup -@1 vs -@4 $c" "$nt markdup --no-PG -@ 1 -O sam md.$c.bam -" "$nt markdup --no-PG -@ 4 -O sam md.$c.bam -"
    check "collate $c" "$st collate --no-PG -O $in | sam - | sort" "$nt collate --no-PG -O $in | sam - | sort"
    check "collate -@1 vs -@4 $c" "$nt collate --no-PG -@ 1 -O $in | sam - | sort" "$nt collate --no-PG -@ 4 -O $in | sam - | sort"
    rm -rf split-st split-nt
    mkdir split-st split-nt
    check "split $c" "$st split --no-PG -u split-st/u.bam -f 'split-st/%!.%.' $in && for f in \$(ls split-st); do echo \$f; sam split-st/\$f; done" \
        "$nt split --no-PG -u split-nt/u.bam -f 'split-nt/%!.%.' $in && for f in \$(ls split-nt); do echo \$f; sam split-nt/\$f; done"
    check "split -d RG $c" "rm -rf split-st/*; $st split --no-PG -d RG -f 'split-st/%!.%.' $in && for f in \$(ls split-st); do echo \$f; sam split-st/\$f; done" \
        "rm -rf split-nt/*; $nt split --no-PG -d RG -f 'split-nt/%!.%.' $in && for f in \$(ls split-nt); do echo \$f; sam split-nt/\$f; done"
done

# cram-size
for c in $cases; do
    for f in 3.0.cram 3.1.cram 2.1.cram embed.cram; do
        input "$c" "$f"
        check "cram-size -v -e $c $f" "$st cram-size -v -e $in" "$nt cram-size -v -e $in"
    done
done

# faidx, fqidx, dict, import, flags
for fa in edge.fa ref.fa manysq.fa; do
    cp "$edge/$fa" st.$fa
    cp "$edge/$fa" nt.$fa
    check "faidx $fa" "$st faidx st.$fa && cat st.$fa.fai" "$nt faidx nt.$fa && cat nt.$fa.fai"
    check "dict $fa" "$st dict -u U $edge/$fa" "$nt dict -u U $edge/$fa"
done
check "faidx edge.fa regions" "$st faidx st.edge.fa empty lower one:1-1 nofinal:5-20 iupac:3 lower:0-2 one:2-5" "$nt faidx nt.edge.fa empty lower one:1-1 nofinal:5-20 iupac:3 lower:0-2 one:2-5"
check "faidx edge.fa missing" "$st faidx st.edge.fa nosuch" "$nt faidx nt.edge.fa nosuch"
check "faidx ref.fa past end" "$st faidx st.ref.fa chr1:195-300 chr2:121-130 chrL:39999" "$nt faidx nt.ref.fa chr1:195-300 chr2:121-130 chrL:39999"
check "faidx -i -n 7" "$st faidx -i -n 7 st.edge.fa lower iupac" "$nt faidx -i -n 7 nt.edge.fa lower iupac"
cp "$edge/edge.fq" st.edge.fq
cp "$edge/edge.fq" nt.edge.fq
check "fqidx edge.fq" "$st fqidx st.edge.fq && cat st.edge.fq.fai" "$nt fqidx nt.edge.fq && cat nt.edge.fq.fai"
check "fqidx edge.fq regions" "$st fqidx st.edge.fq r0 r1 r2:2-3 r3/1 r4:10-20" "$nt fqidx nt.edge.fq r0 r1 r2:2-3 r3/1 r4:10-20"
check "import edge.fq" "$st import --no-PG -0 $edge/edge.fq" "$nt import --no-PG -0 $edge/edge.fq"
check "import -T '*' edge.fq" "$st import --no-PG -T '*' -0 $edge/edge.fq" "$nt import --no-PG -T '*' -0 $edge/edge.fq"
check "import -i -N -1 -2" "$st import --no-PG -N -1 $edge/edge.fq -2 $edge/edge.fq" "$nt import --no-PG -N -1 $edge/edge.fq -2 $edge/edge.fq"
check "import edge.fa" "$st import --no-PG -0 $edge/edge.fa" "$nt import --no-PG -0 $edge/edge.fa"
for fl in 0 4 4095 4096 65535 -1 0x900 0o20 PAIRED,UNMAP 'PAIRED,nosuch' '' SECONDARY,SUPPLEMENTARY,DUP,QCFAIL; do
    check "flags '$fl'" "$st flags '$fl'" "$nt flags '$fl'"
done

echo "$checks checks, $failures differ, $known_count known differences"
[ "$failures" -eq 0 ]
