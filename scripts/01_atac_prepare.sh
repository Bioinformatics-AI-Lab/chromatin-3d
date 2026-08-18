#!/usr/bin/env bash
# ENCODE filtered alignments -> Tn5-shifted cut sites -> peaks.
#
# We start from ENCODE's filtered, deduplicated BAM rather than from FASTQ.
# Read alignment is the most resource-hungry step in this project and the least
# assay-specific one -- ENCODE's pipeline is the reference standard for ATAC-seq
# and rerunning bowtie2 would demonstrate nothing about chromatin. Every step
# that is actually specific to ATAC-seq is below.
set -euo pipefail
THREADS=${THREADS:-8}
BAM=$1                                  # ENCODE filtered alignments (hg38)
SAMPLE=${2:-GM12878}
OUT=../results/atac; mkdir -p "$OUT"

samtools index -@ "$THREADS" "$BAM" 2>/dev/null || true

# Mitochondrial fraction is an ATAC QC metric in its own right (chrM routinely
# consumes 20-80% of an ATAC library and is pure background). Record it before
# discarding chrM, even though ENCODE has usually removed it already.
samtools idxstats "$BAM" > "$OUT/${SAMPLE}.idxstats.txt"
awk '$1=="chrM"{m=$3} {t+=$3} END{printf "chrM fraction: %.4f\n", (t?m/t:0)}' \
  "$OUT/${SAMPLE}.idxstats.txt"

# Keep properly paired, primary, non-QC-fail, MAPQ>=30 on the main assembly.
# ENCODE has already filtered and deduplicated, so this is a guard rather than
# a transformation -- it should remove very little.
MAIN=$(cut -f1 ../refs/hg38.chrom.sizes | grep -E '^chr([0-9]+|X)$' | tr '\n' ' ')
samtools view -@ "$THREADS" -b -f 2 -F 1804 -q 30 "$BAM" $MAIN \
| samtools sort -@ "$THREADS" -o "$OUT/${SAMPLE}.dedup.bam" -
samtools index -@ "$THREADS" "$OUT/${SAMPLE}.dedup.bam"

# Tn5 correction: the transposase inserts as a dimer and leaves a 9 bp stagger,
# so the true insertion point is +4 on the plus strand and -5 on the minus
# strand. Skipping this offsets every peak summit and footprint by ~4-5 bp.
# Emit 1 bp cut sites -- peaks are called on insertions, not on read pileup.
bedtools bamtobed -i "$OUT/${SAMPLE}.dedup.bam" \
| awk 'BEGIN{OFS="\t"}
       $6=="+" {s=$2+4; print $1, s, s+1, ".", ".", $6}
       $6=="-" {s=$3-5; if (s>1) print $1, s-1, s, ".", ".", $6}' \
| sort -k1,1 -k2,2n > "$OUT/${SAMPLE}.cutsites.bed"
bgzip -f -c "$OUT/${SAMPLE}.cutsites.bed" > "$OUT/${SAMPLE}.cutsites.bed.gz"
tabix -f -p bed "$OUT/${SAMPLE}.cutsites.bed.gz"
echo "cut sites: $(wc -l < "$OUT/${SAMPLE}.cutsites.bed")"

# --nomodel because MACS2's fragment-length model is built for ChIP-seq and is
# meaningless for transposition. --shift/--extsize centre a 150 bp window on
# each cut site. --keep-dup all because duplicates were removed upstream.
macs2 callpeak -t "$OUT/${SAMPLE}.cutsites.bed" -f BED -n "$SAMPLE" \
  --outdir "$OUT" --nomodel --shift -75 --extsize 150 \
  --keep-dup all -q 0.01 -g hs --call-summits

echo "peaks: $(wc -l < "$OUT/${SAMPLE}_peaks.narrowPeak")"

# The BAM is only needed for the fragment-size histogram in 02_atac_qc.py.
# After that it can be deleted -- see RUNBOOK cleanup.
