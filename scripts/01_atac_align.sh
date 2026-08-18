#!/usr/bin/env bash
# FASTQ -> filtered alignments -> Tn5-shifted cut sites -> peaks.
set -euo pipefail
THREADS=${THREADS:-8}
IDX=../refs/bowtie2/GRCh38
OUT=../results/atac; mkdir -p "$OUT"
R1=$1; R2=$2; SAMPLE=${3:-GM12878}

# -X 2000 is not optional. bowtie2 defaults to a 500 bp maximum insert, which
# silently discards the di- and tri-nucleosomal fragments the QC depends on.
bowtie2 --very-sensitive -X 2000 -p "$THREADS" -x "$IDX" -1 "$R1" -2 "$R2" \
  2> "$OUT/${SAMPLE}.bowtie2.log" \
| samtools sort -@ "$THREADS" -o "$OUT/${SAMPLE}.raw.bam" -
samtools index "$OUT/${SAMPLE}.raw.bam"

# Mitochondrial fraction is an ATAC QC metric in its own right -- record it
# before discarding chrM.
samtools idxstats "$OUT/${SAMPLE}.raw.bam" > "$OUT/${SAMPLE}.idxstats.txt"

# Properly paired, uniquely mapped, primary, not QC-fail (-F 1804), MAPQ>=30,
# then drop chrM and unplaced contigs.
samtools view -@ "$THREADS" -b -f 2 -F 1804 -q 30 "$OUT/${SAMPLE}.raw.bam" \
  $(cut -f1 ../refs/hg38.chrom.sizes | grep -E '^chr([0-9]+|X)$' | tr '\n' ' ') \
| samtools sort -@ "$THREADS" -o "$OUT/${SAMPLE}.filt.bam" -
samtools index "$OUT/${SAMPLE}.filt.bam"

sambamba markdup -r -t "$THREADS" "$OUT/${SAMPLE}.filt.bam" "$OUT/${SAMPLE}.dedup.bam" \
  2> "$OUT/${SAMPLE}.markdup.log"
samtools index "$OUT/${SAMPLE}.dedup.bam"

# Tn5 correction: the transposase inserts as a dimer leaving a 9 bp stagger, so
# the true insertion point is +4 on the plus strand and -5 on the minus strand.
# Emit 1 bp cut sites; peaks are called on insertions, not on read pileup.
bedtools bamtobed -i "$OUT/${SAMPLE}.dedup.bam" \
| awk 'BEGIN{OFS="\t"}
       $6=="+" {s=$2+4; print $1, s, s+1, ".", ".", $6}
       $6=="-" {s=$3-5; if (s>0) print $1, s-1, s, ".", ".", $6}' \
| sort -k1,1 -k2,2n > "$OUT/${SAMPLE}.cutsites.bed"
bgzip -f -@ "$THREADS" -c "$OUT/${SAMPLE}.cutsites.bed" > "$OUT/${SAMPLE}.cutsites.bed.gz"
tabix -f -p bed "$OUT/${SAMPLE}.cutsites.bed.gz"

# --nomodel because MACS2's fragment-length model is built for ChIP-seq and is
# meaningless for transposition. --shift/-extsize centre a 150 bp window on each
# cut site. --keep-dup all because duplicates were already removed properly.
macs2 callpeak -t "$OUT/${SAMPLE}.cutsites.bed" -f BED -n "$SAMPLE" \
  --outdir "$OUT" --nomodel --shift -75 --extsize 150 \
  --keep-dup all -q 0.01 -g hs --call-summits

echo "Peaks: $(wc -l < "$OUT/${SAMPLE}_peaks.narrowPeak")"
