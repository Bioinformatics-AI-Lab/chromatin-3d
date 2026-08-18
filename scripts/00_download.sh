#!/usr/bin/env bash
# Fetch inputs. We take ENCODE's filtered, deduplicated alignments rather than
# FASTQ -- see 01_atac_prepare.sh for why. This saves ~55 GB and ~4 hours and
# costs nothing that is specific to ATAC-seq.
set -euo pipefail
mkdir -p ../data/{atac,hic,rna} ../refs

ENCODE_EXP=ENCSR637XSC
MCOOL_URL="https://4dn-open-data-public.s3.amazonaws.com/fourfront-webprod/wfoutput/d6abea45-b0bb-4154-9854-1d3075b98097/4DNFIXP4QG5B.mcool"

# --- ATAC alignments -------------------------------------------------------
# Resolve the BAM URL from the ENCODE API rather than hardcoding an ENCFF id,
# so this keeps working when the experiment is re-released. We want the
# "alignments" output_type on GRCh38 (ENCODE's filtered, deduplicated BAM).
curl -sL "https://www.encodeproject.org/experiments/${ENCODE_EXP}/?format=json" \
| python3 - <<'PY' > ../data/atac/bam_manifest.txt
import json, sys
exp = json.load(sys.stdin)
for f in exp["files"]:
    if (f.get("file_format") == "bam"
            and f.get("status") == "released"
            and f.get("assembly") == "GRCh38"):
        print(f["accession"], f.get("output_type"), f.get("file_size", 0),
              "https://www.encodeproject.org" + f["href"], sep="\t")
PY
echo "Available BAMs:"; cat ../data/atac/bam_manifest.txt
BAM_URL=$(grep -m1 -P '\talignments\t' ../data/atac/bam_manifest.txt | cut -f4)
[ -n "$BAM_URL" ] || { echo "No 'alignments' BAM found -- inspect bam_manifest.txt"; exit 1; }
curl -L -o ../data/atac/${ENCODE_EXP}.bam "$BAM_URL"

# ENCODE's own peak calls, for comparison against ours.
curl -sL -o ../data/atac/ENCFF748UZH.bed.gz \
  "https://www.encodeproject.org/files/ENCFF748UZH/@@download/ENCFF748UZH.bed.gz"

# --- Hi-C contact matrix ---------------------------------------------------
curl -L -o ../data/hic/GM12878.mcool "$MCOOL_URL"
cooler ls ../data/hic/GM12878.mcool

# --- References (no genome FASTA or bowtie2 index needed now) --------------
cd ../refs
[ -f hg38.chrom.sizes ] || curl -sLO https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.chrom.sizes
[ -f gencode.v44.annotation.gtf.gz ] || \
  curl -sLO https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_44/gencode.v44.annotation.gtf.gz
# Compartment phasing needs GC content, so the FASTA is still required for
# Stage 2. It is the one large reference left (~3 GB).
[ -f hg38.fa ] || { curl -sLO https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz; gunzip hg38.fa.gz; }
samtools faidx hg38.fa
echo "Downloads complete."
