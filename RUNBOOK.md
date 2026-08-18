# RUNBOOK

Concrete order of operations, with realistic time and disk estimates.

## 0. Where to run

| | Cores | RAM | Disk | Wall time |
|---|---|---|---|---|
| Stage 1 (ATAC from FASTQ) | 8–16 | 32 GB | ~150 GB | 3–6 h |
| Stage 2 (Hi-C features) | 4 | 16 GB | ~30 GB | 30–60 min |
| Stage 3 (integration) | 2 | 16 GB | — | 20–40 min |

Stages 2–3 run on a laptop. Stage 1 wants a workstation or a cloud VM
(`c6i.4xlarge` or similar, a few dollars for the run). Use compute you own —
a personal portfolio project doesn't belong on an employer's cluster.

## 1. Environment

```bash
git clone <your-repo-url> chromatin-3d && cd chromatin-3d
mamba env create -f environment.yml     # or: conda env create -f environment.yml
conda activate chromatin-3d
macs2 --version && cooler --version && pairtools --version
```

## 2. References (~25 GB, do this first — it's the long pole)

```bash
mkdir -p refs && cd refs

curl -LO https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.chrom.sizes
curl -LO https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz
gunzip hg38.fa.gz && samtools faidx hg38.fa
curl -LO https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_44/gencode.v44.annotation.gtf.gz

# Prebuilt bowtie2 index. Building it yourself takes ~3 h and ~8 GB RAM;
# there is no reason to. Confirm the current link on the Bowtie 2 homepage.
mkdir -p bowtie2 && cd bowtie2
curl -LO https://genome-idx.s3.amazonaws.com/bt/GRCh38_noalt_as.zip
unzip GRCh38_noalt_as.zip && cd ../..
```

Then point `config/config.yaml` at what you actually downloaded:

```yaml
genome:
  fasta: refs/hg38.fa
  bowtie2_index: refs/bowtie2/GRCh38_noalt_as/GRCh38_noalt_as
  chromsizes: refs/hg38.chrom.sizes
  gtf: refs/gencode.v44.annotation.gtf.gz
```

The no-alt analysis set uses UCSC-style `chr` names, which matches
`hg38.chrom.sizes`, the GENCODE GTF, and the 4DN matrix. Mixing an Ensembl-named
index (`1`, `2`, …) into this stack is the most common way to get a silently
empty TSS profile.

## 3. Download data

```bash
cd scripts && bash 00_download.sh && cd ..
```

Sanity check before committing to Stage 1:

```bash
cooler ls data/hic/GM12878.mcool      # should list several resolutions
ls -lh data/atac/*.fastq.gz
```

## 4. Stage 1 — ATAC

```bash
cd scripts
THREADS=8 bash 01_atac_align.sh \
  ../data/atac/<R1>.fastq.gz ../data/atac/<R2>.fastq.gz GM12878
cd ..
python scripts/02_atac_qc.py --config config/config.yaml --sample GM12878
cat results/atac_qc.json
```

**Stop and read the QC before going further.** Targets: TSS enrichment ≥ 5
(ideal ≥ 7), FRiP ≥ 0.2, and a fragment-size plot with a sub-100 bp peak plus a
visible ~200 bp shoulder. If TSS enrichment comes back near 1.0, the cause is
almost always chromosome-name mismatch between the GTF and the BAM, not a bad
library — check `samtools idxstats` output against the GTF's first column.

## 5. Stage 2 — Hi-C

```bash
python scripts/04_hic_features.py --mcool data/hic/GM12878.mcool
cat results/hic_qc.json
```

If it errors on a missing `weight` column, the matrix is unbalanced:
`cooler balance -p 8 data/hic/GM12878.mcool::/resolutions/10000`.

Optional read-level demo (download a `.pairs` file from the same 4DN experiment
set first):

```bash
cd scripts && THREADS=8 bash 03_hic_pairs_demo.sh ../data/hic/<accession>.pairs.gz && cd ..
```

## 6. Stage 3 — Integration

Set `rna.encode_quantification` in `config.yaml` to a downloaded GM12878 polyA
RNA-seq gene-quantification TSV first — without it the "expressed promoter"
framing isn't supported and the script says so.

```bash
python scripts/05_integrate.py --mcool data/hic/GM12878.mcool
cat results/integration.json
```

Expected shape of the result: median `log2_vs_null` clearly above 0, and
`median_within_boundary` > `median_across_boundary`. If enrichment is ~0, check
that `expected_cis.parquet` was written by Stage 2 and that `min_peak_tss_distance`
is at least two bins.

## 7. Commit

```bash
git add -A
git commit -m "ATAC + Hi-C processing, QC, and peak-to-promoter contact analysis"
git push
```

`.gitignore` excludes raw data, BAMs, and matrices but **keeps**
`results/*.json` and `results/figures/`. Commit those — a reviewer should see
your TSS enrichment plot and your integration statistics without downloading
40 GB or running anything.
