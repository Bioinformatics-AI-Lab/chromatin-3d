#!/usr/bin/env bash
# Fetch inputs. We take ENCODE's filtered, deduplicated alignments rather than
# FASTQ -- see 01_atac_prepare.sh for why. This saves ~55 GB and ~4 hours and
# costs nothing that is specific to ATAC-seq.
set -euo pipefail
mkdir -p ../data/atac ../data/hic ../data/rna ../refs
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

ENCODE_EXP=ENCSR637XSC
MCOOL_URL="https://4dn-open-data-public.s3.amazonaws.com/fourfront-webprod/wfoutput/d6abea45-b0bb-4154-9854-1d3075b98097/4DNFIXP4QG5B.mcool"

# --- ATAC alignments -------------------------------------------------------
# Resolve the BAM from the ENCODE API rather than hardcoding an ENCFF id, so
# this keeps working when the experiment is re-released.
echo "== BAMs available for ${ENCODE_EXP} =="
python3 "$HERE/encode_files.py" "$ENCODE_EXP" --file-format bam --assembly GRCh38 \
  | tee ../data/atac/bam_manifest.tsv

# Pinned by accession for the same reason as the peaks below. This experiment
# has three isogenic replicates and therefore three filtered BAMs; we use one.
# A single replicate is sufficient for QC and for the integration analysis,
# which does not compare conditions.
BAM_ACC=ENCFF962FMH
curl -L --retry 3 -C - -o ../data/atac/${BAM_ACC}.bam \
  "https://www.encodeproject.org/files/${BAM_ACC}/@@download/${BAM_ACC}.bam"
echo "ATAC alignments: data/atac/${BAM_ACC}.bam"

# ENCODE's own peak calls, to cross-check ours against.
#
# Pinned by accession, not by output_type. ENCODE stores several files per
# output type (per replicate pair, and per pipeline version), so selecting by
# type alone is non-deterministic -- a later run can silently pick a different
# file. ENCFF748UZH is the experiment-level pseudoreplicated peak set for
# GM12878 and is the one commonly cited for this cell line.
PEAK_ACC=ENCFF748UZH
curl -L --retry 3 -o ../data/atac/${PEAK_ACC}.bed.gz \
  "https://www.encodeproject.org/files/${PEAK_ACC}/@@download/${PEAK_ACC}.bed.gz"

# --- Hi-C contact matrix ---------------------------------------------------
curl -L --retry 3 -C - -o ../data/hic/GM12878.mcool "$MCOOL_URL"
cooler ls ../data/hic/GM12878.mcool

# --- References ------------------------------------------------------------
cd ../refs
[ -f hg38.chrom.sizes ] || curl -fL --retry 3 -O https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.chrom.sizes
[ -f gencode.v44.annotation.gtf.gz ] || \
  curl -fL --retry 3 -O https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_44/gencode.v44.annotation.gtf.gz
# The FASTA is needed only for GC-phasing the compartment eigenvector (Stage 2).
if [ ! -f hg38.fa ]; then
  curl -L --retry 3 -C - -O https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz
  gunzip hg38.fa.gz
fi
[ -f hg38.fa.fai ] || samtools faidx hg38.fa
cd ..

# Index the BAM -- without it every samtools range query falls back to a full
# scan of a 7 GB file.
[ -f data/atac/${BAM_ACC}.bam.bai ] || samtools index -@ "${THREADS:-8}" data/atac/${BAM_ACC}.bam

# Verify everything landed. curl -f above makes HTTP errors fatal, but a
# truncated or empty file still needs catching before hours of downstream work.
for f in data/atac/${BAM_ACC}.bam data/hic/GM12878.mcool \
         refs/hg38.chrom.sizes refs/gencode.v44.annotation.gtf.gz refs/hg38.fa; do
  [ -s "../$f" ] || [ -s "$f" ] || { echo "MISSING or empty: $f"; exit 1; }
done
gunzip -t refs/gencode.v44.annotation.gtf.gz
echo "Downloads complete and verified."
