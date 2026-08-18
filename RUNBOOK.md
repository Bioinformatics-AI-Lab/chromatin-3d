# RUNBOOK

Everything runs locally. No cloud, no cost.

## Requirements

- ~40 GB free disk during the run; ~5 MB kept afterward
- A few GB of RAM (cooler reads matrix ranges out of HDF5)
- macOS (Apple Silicon or Intel) or Linux
- No admin rights needed — conda installs into your home directory

Total wall time: about 2 hours, most of it downloading.

## 1. Environment

On Apple Silicon, create the environment under x86 emulation. Several bioconda
packages here (`macs2` in particular) have no arm64 build, and the emulation
penalty is irrelevant now that read alignment is out of scope.

```
CONDA_SUBDIR=osx-64 mamba env create -f environment.yml
conda activate chromatin-3d
conda config --env --set subdir osx-64
```

On Intel Mac or Linux, drop the `CONDA_SUBDIR` prefix and the last line.

Verify:

```
samtools --version | head -1
macs2 --version
cooler --version
python -c "import cooltools, bioframe, pysam; print('ok')"
```

## 2. Download (~35 GB, the long pole)

```
cd scripts && bash 00_download.sh && cd ..
```

This fetches ENCODE's filtered ATAC alignments, the 4DN Hi-C matrix, hg38
chrom.sizes, the GENCODE GTF, and the genome FASTA (needed only for GC-phasing
the compartment eigenvector).

Sanity check before proceeding:

```
cooler ls data/hic/GM12878.mcool
samtools idxstats data/atac/ENCSR637XSC.bam | head -3
zcat refs/gencode.v44.annotation.gtf.gz | grep -v '^#' | cut -f1 | sort -u | head -3
```

The BAM and the GTF must agree on chromosome naming — both `chr1`-style. A
mismatch here is the single most common cause of a flat TSS profile.

## 3. Stage 1 — ATAC

```
cd scripts
THREADS=8 bash 01_atac_prepare.sh ../data/atac/ENCSR637XSC.bam GM12878
cd ..
python scripts/02_atac_qc.py --sample GM12878
cat results/atac_qc.json
```

Read the QC before continuing. Targets: TSS enrichment >= 5 (ideal >= 7),
FRiP >= 0.2, and a fragment-size plot with a sub-100 bp peak plus a visible
~200 bp shoulder.

If TSS enrichment comes back near 1.0, it is almost always chromosome naming,
not a bad library.

Optional cross-check against ENCODE's own peak calls:

```
bedtools jaccard -a <(sort -k1,1 -k2,2n results/atac/GM12878_peaks.narrowPeak) \
                 -b <(zcat data/atac/ENCFF748UZH.bed.gz | sort -k1,1 -k2,2n)
```

## 4. Stage 2 — Hi-C

```
python scripts/04_hic_features.py --mcool data/hic/GM12878.mcool
cat results/hic_qc.json
```

If it errors on a missing `weight` column the matrix is unbalanced:

```
cooler balance -p 8 data/hic/GM12878.mcool::/resolutions/10000
```

## 5. Stage 3 — Integration

Set `rna.encode_quantification` in `config/config.yaml` to a downloaded GM12878
polyA RNA-seq gene-quantification TSV first. Without it the script falls back to
all protein-coding TSS and says so.

```
python scripts/05_integrate.py --mcool data/hic/GM12878.mcool
cat results/integration.json
```

Expected shape: median `log2_vs_null` clearly above 0, and
`median_within_boundary` > `median_across_boundary`.

## 6. Commit results

```
git add results/*.json results/figures
git commit -m "Add ATAC QC, Hi-C features, and peak-to-promoter integration results"
git push
```

Then paste the actual numbers into the README under a short "Results" heading.
That is what makes the repo evidence rather than scaffolding.

## 7. Reclaim disk

```
rm -rf data refs
rm -f results/atac/*.bam results/atac/*.bai results/atac/*.bed results/atac/*.bed.gz*
rm -f results/hic/*.parquet
conda clean -a -y
du -sh .
```

Should come back under 5 MB. Everything removed is regenerable from the scripts,
which is the point of keeping the pipeline in version control.
