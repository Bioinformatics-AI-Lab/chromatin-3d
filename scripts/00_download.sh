#!/usr/bin/env bash
# Fetch inputs. ENCODE serves files directly; 4DN requires login on the portal
# but mirrors processed files on a public S3 bucket, which is what we use.
set -euo pipefail
mkdir -p ../data/{atac,hic,rna} ../refs

ENCODE_EXP=ENCSR637XSC
MCOOL_URL="https://4dn-open-data-public.s3.amazonaws.com/fourfront-webprod/wfoutput/d6abea45-b0bb-4154-9854-1d3075b98097/4DNFIXP4QG5B.mcool"

# --- ATAC FASTQ ------------------------------------------------------------
# Resolve FASTQ URLs from the ENCODE API rather than hardcoding ENCFF ids,
# so this keeps working when the experiment is re-released.
curl -sL "https://www.encodeproject.org/experiments/${ENCODE_EXP}/?format=json" \
| python3 - <<'PY' > ../data/atac/fastq_urls.txt
import json, sys
exp = json.load(sys.stdin)
for f in exp["files"]:
    if f.get("file_format") == "fastq" and f.get("status") == "released":
        print("https://www.encodeproject.org" + f["href"], f.get("paired_end"), f["accession"])
PY
echo "FASTQ manifest written to data/atac/fastq_urls.txt"
awk '{print $1}' ../data/atac/fastq_urls.txt | xargs -P 4 -n 1 -I{} curl -sL -O --output-dir ../data/atac {}

# --- Hi-C contact matrix ---------------------------------------------------
curl -L -o ../data/hic/GM12878.mcool "$MCOOL_URL"
cooler ls ../data/hic/GM12878.mcool   # sanity check: prints available resolutions

# --- References ------------------------------------------------------------
cd ../refs
[ -f hg38.chrom.sizes ] || curl -sLO https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.chrom.sizes
[ -f gencode.v44.annotation.gtf.gz ] || \
  curl -sLO https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_44/gencode.v44.annotation.gtf.gz
echo "Build the bowtie2 index separately (bowtie2-build is slow); see README."
