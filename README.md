# chromatin-3d

End-to-end processing, QC, and integration of chromatin accessibility (ATAC-seq) and 3D
genome architecture (Hi-C) data, using public ENCODE and 4D Nucleome datasets in GM12878.

This is a self-directed learning project. It is **not** wet-lab work and **not** production
infrastructure — it exists to demonstrate how to take two assay types from raw
data (or from published contact matrices, where raw depth is prohibitive) through
assay-specific QC to a biologically meaningful integrated result.

---

## Why this project

I work primarily in transcriptomics and methylation at cohort scale. ATAC-seq and Hi-C sit
adjacent to that work and answer questions expression alone cannot: *which* regulatory
elements are open, *which* transcription factors are acting, and *which* elements are
physically close enough to a promoter to plausibly regulate it.

The interesting part is not running the tools — it is that the statistical structure of
these assays is genuinely different from RNA-seq, and the analysis has to change to match:

| | scRNA-seq | ATAC-seq | Hi-C |
|---|---|---|---|
| Unit of observation | transcript (UMI) | Tn5 insertion site | *pair* of loci |
| Features | given (~20k genes) | inferred (peaks called from data) | combinatorial (bin × bin) |
| Sparsity origin | undersampling | ploidy-bounded (0/1/2 per locus) | combinatorially bounded |
| Dominant nuisance | depth, mito%, cell cycle | Tn5 sequence bias, CNV | **genomic distance decay** |
| Normalization | size factors, log/NB | TF-IDF, quantile | matrix balancing (ICE) |

Every stage below is an attempt to handle one of those differences correctly rather than
reaching for an RNA-seq idiom that does not transfer.

---

## Results

All numbers below are from the committed JSON in `results/`. Figures are in
`results/figures/`.

### ATAC-seq QC (GM12878, ENCODE `ENCFF962FMH`)

| Metric | Value | Threshold |
|---|---|---|
| TSS enrichment | **11.88** | ENCODE hg38: ≥5 acceptable, ≥7 ideal |
| FRiP | **0.439** | >0.2 typical for bulk |
| Nucleosome signal (mono/NFR) | 0.792 | — |
| Peaks called | 236,397 | ENCODE's own set: ~278,000 |
| Tn5-shifted cut sites | 163M | — |

The fragment-size distribution shows the expected nucleosome-free peak below
100 bp with ~200 bp periodicity, so nucleosome structure survived transposition.
Peak count is within 15% of ENCODE's published set for this cell line, called
independently with different parameters.

### Hi-C structure (4DN `4DNFIXP4QG5B`, 10 kb)

| Metric | Value |
|---|---|
| Bins / non-zero pixels | 308,837 / 1.57 × 10⁹ |
| A-compartment fraction | 0.43 |
| Insulation boundaries (100 / 200 / 500 kb window) | 7,465 / 5,203 / 3,024 |

Boundary counts are reported at three window sizes rather than one. Boundary
calls are strongly scale-dependent and concordance between callers is poor;
presenting a single window as ground truth would claim more than the data
supports.

### Integration: do accessible peaks contact expressed promoters?

2,281,273 peak–TSS pairs within 1 Mb, each compared against 20 distance-matched
null pairs drawn at the same bin separation on the same chromosome.

| | median log2(O/E vs null) | n |
|---|---|---|
| **All pairs** | **+0.147** | 2,277,901 |
| Within an insulation boundary | **+0.539** | 589,802 |
| Crossing an insulation boundary | **−0.021** | 1,688,099 |

Accessible regions contact promoters roughly 45% above distance-matched
background when no insulation boundary intervenes, and at essentially background
when one does (Mann–Whitney p < 1e-300). The insulation boundaries and the ATAC
peaks were called independently in this pipeline, so the agreement between them
is not built in by construction.

The distance-matched null is what makes this interpretable. Contact frequency
falls steeply with genomic separation and ATAC peaks cluster near promoters, so
peak–TSS pairs are systematically closer than random pairs. Comparing raw
contact frequency would recover the distance-decay curve and mislabel it
regulation.

**A null result, reported as such.** Peak signal strength does not predict
contact enrichment: Spearman rho = −0.024. The p-value is 1e-293, which is a
consequence of n ≈ 2.3 million rather than evidence of an effect. Accessibility
magnitude and 3D proximity appear to be close to independent here.

### Limitations

- **One cell line, no biological replicates.** No differential analysis is
  attempted. Comparative Hi-C has no settled framework equivalent to DESeq2, and
  presenting feature-set overlap as a statistical test would be misleading.
- **RNA and ATAC come from different ENCODE eras.** The quantification
  (`ENCFF783VBA`, ENCODE2-era total RNA-seq) retained only 9,383 of ~20,000
  protein-coding TSS at TPM ≥ 1 — lower than expected for a lymphoblastoid line,
  and most likely a protocol-era mismatch rather than biology. A matched
  ENCODE4 polyA quantification would be the better input.
- **Zero-contact pairs.** 3,372 pairs (0.15%) had zero observed contact and are
  excluded from log-ratio statistics; the Wilcoxon test runs on raw O/E and
  includes them. At this fraction the exclusion does not move the estimates, but
  dropping them silently would bias enrichment upward.
- **One Hi-C replicate, one ATAC replicate.** The ATAC experiment has three
  isogenic replicates; this analysis uses one.
- **Public data only.** No wet-lab work, and no chromatin experiment designed at
  the bench.

---

## What's here

```
chromatin-3d/
├── config/config.yaml          # accessions, genome, resolutions, thresholds
├── scripts/
│   ├── 00_download.sh          # fetch ENCODE + 4DN inputs
│   ├── 01_atac_prepare.sh      # ENCODE BAM → Tn5-shifted cut sites → peaks
│   ├── 02_atac_qc.py           # TSS enrichment, FRiP, fragment-size periodicity
│   ├── 03_hic_pairs_demo.sh    # pairtools parse/sort/dedup on subsampled pairs
│   ├── 04_hic_features.py      # ICE balancing, distance decay, compartments, insulation
│   └── 05_integrate.py         # peak→TSS contact enrichment vs distance-matched null
├── notebooks/01_report.ipynb   # figures and narrative
└── results/                    # outputs (gitignored except summary JSON + figures)
```

---

## Stages

Each stage is independently runnable and each has an explicit claim it supports. If you
only have a day, do Stage 1 and Stage 3.

### Stage 1 — ATAC-seq (`01_atac_prepare.sh`, `02_atac_qc.py`)

Input: ENCODE `ENCSR637XSC` (bulk ATAC-seq, GM12878, ENCODE4, Snyder lab), starting from
the **filtered, deduplicated alignments** rather than FASTQ.

That choice is deliberate. Read alignment is the most resource-hungry step here and the
least assay-specific one — ENCODE's pipeline is the reference standard for ATAC-seq, and
rerunning `bowtie2` would demonstrate nothing about chromatin while costing ~55 GB and
several hours. Everything below *is* specific to ATAC-seq.

1. Guard filter: MAPQ ≥ 30, properly paired, main assembly only, drop `chrM`
   (mitochondrial reads routinely consume 20–80% of an ATAC library and are pure
   background). Mitochondrial fraction is recorded as a QC metric before removal.
2. **Tn5 shift**: `+4` on the plus strand, `−5` on the minus strand. The informative
   coordinate is the 5′ end of the read, offset for the 9 bp staggered cut the Tn5 dimer
   leaves. Skipping this shifts every footprint and peak summit by ~4–5 bp.
3. Peak calling on cut sites: `macs2 callpeak --nomodel --shift -75 --extsize 150
   --keep-dup all`. `--nomodel` because MACS2's fragment-length model is built for
   ChIP-seq and is meaningless for transposition; `--keep-dup all` because duplicates were
   already removed properly upstream.

QC produced (`02_atac_qc.py`):

- **TSS enrichment score** — mean insertion density in the TSS ±100 bp core divided by the
  mean over the ±1900–2000 bp flanks. ENCODE's threshold for hg38 is ≥ 5 (ideal ≥ 7).
- **FRiP** — fraction of filtered reads in peaks; > 0.2 is a reasonable bulk target.
- **Fragment-size distribution** — should show a nucleosome-free peak below ~100 bp and
  clear ~200 bp periodicity. A featureless smear means the library is over-transposed or
  the nuclei ruptured, and no downstream analysis will rescue it.

These three replace the RNA-seq QC panel entirely. There is no mito% / gene-count analogue.

### Stage 2 — Hi-C (`03_hic_pairs_demo.sh`, `04_hic_features.py`)

Hi-C resolution scales as √N reads, so a kilobase-resolution map is a billions-of-reads
proposition. This stage therefore splits:

- **Read-level competence** (`03_hic_pairs_demo.sh`) on a subsampled `.pairs` file:
  `pairtools parse` → `sort` → `dedup` → `select`, showing pair classification (valid pair
  vs. dangling end vs. self-circle) and deduplication on the position+strand 4-tuple, then
  `cooler cload pairs` → `cooler zoomify` to build a multi-resolution `.mcool`.
- **Feature calling** (`04_hic_features.py`) on the published deep matrix, 4DN
  `4DNFIXP4QG5B.mcool` (in situ Hi-C, MboI, GM12878).

Features computed:

- **ICE balancing** (`cooler balance`) — iterative row/column rescaling to a doubly
  stochastic matrix, on the assumption of equal per-bin visibility. This is a matrix
  problem, not a per-sample size factor.
- **Distance-decay curve and expected** (`cooltools.expected_cis`) — P(contact) vs genomic
  separation, spanning several orders of magnitude. This is the single dominant signal in
  the data and almost every downstream step operates on observed/expected, not observed.
- **A/B compartments** (`cooltools.eigs_cis`) — eigendecomposition of the
  distance-normalized correlation matrix, with the PC1 sign phased by GC content so that
  the A compartment is consistently the gene-rich one. Unphased eigenvectors flip sign
  arbitrarily between chromosomes.
- **Insulation profile and boundaries** (`cooltools.insulation`) at multiple window sizes.
  Reported at more than one window deliberately: TAD boundary calls are notoriously
  method- and scale-dependent, and presenting a single window as ground truth overstates
  what the data supports.

### Stage 3 — Integration (`05_integrate.py`)

The point of the whole exercise. For each ATAC peak within 1 Mb of an expressed gene TSS:

1. Look up the observed/expected contact frequency between the peak's bin and the TSS bin
   at 10 kb.
2. Build a **distance-matched null** by sampling non-peak bins at the same genomic
   separation. Without this the result is a restatement of the distance-decay curve — the
   most common way to get a spuriously "significant" 3D result.
3. Test whether peak→TSS pairs are enriched over the null, and whether enrichment scales
   with peak accessibility and with target-gene expression.
4. Separately, split peak→TSS pairs by whether they cross an insulation boundary, and
   compare contact enrichment within vs. across boundaries.

Expected result (and a check that the pipeline is behaving): accessible peaks contact
expressed promoters more than distance-matched background, and the effect is stronger for
pairs that do not cross a boundary.

---

## Running it

```bash
mamba env create -f environment.yml
conda activate chromatin-3d

bash scripts/00_download.sh
bash scripts/01_atac_prepare.sh ../data/atac/ENCFF962FMH.bam
python scripts/02_atac_qc.py
bash scripts/03_hic_pairs_demo.sh
python scripts/04_hic_features.py
python scripts/05_integrate.py
```

Compute: the whole thing runs on a laptop. Peak disk is ~40 GB (mostly the Hi-C matrix and
the genome FASTA), nearly all of which is deletable afterward; peak memory is a few GB,
since cooler reads ranges out of HDF5 rather than loading the matrix.

## Data

| Purpose | Source | Accession |
|---|---|---|
| Bulk ATAC-seq alignments, GM12878 | ENCODE | `ENCSR637XSC` / `ENCFF962FMH` |
| ATAC-seq peaks (reference) | ENCODE | `ENCFF748UZH` |

File accessions are pinned rather than resolved by output type: ENCODE stores several
files per type (one per replicate pair, and one per pipeline version), so type-based
selection is not reproducible across runs. The ATAC experiment has three isogenic
replicates and this analysis uses one — sufficient for QC and for the integration
analysis, which does not compare conditions.
| in situ Hi-C, GM12878 (`.mcool`) | 4DN | `4DNFIXP4QG5B` (set `4DNES3JX38V5`) |
| RNA-seq, GM12878 | ENCODE | set in `config.yaml` |

The 4DN portal requires login for direct download, but processed files are mirrored on the
public S3 bucket `4dn-open-data-public` — `00_download.sh` uses that path. Confirm the
accessions on the portal before a long run; 4DN and ENCODE both re-release files.

## Honest scope

- Public data only. I have not run these assays and have not designed a chromatin
  experiment at the bench.
- ATAC starts from ENCODE's filtered alignments, not from FASTQ — read alignment is
  deliberately out of scope (see Stage 1).
- The Hi-C stage starts from a published contact matrix, not from raw FASTQ. The
  `pairtools` stage demonstrates the read-level workflow on subsampled data.
- Single cell type, no biological replicates, so there is no differential analysis here.
  Comparative Hi-C has no settled framework equivalent to DESeq2, and I did not want to
  present feature-set overlap as if it were a statistical test.

## References

Buenrostro et al. 2013 (ATAC-seq) · Corces et al. 2017 (Omni-ATAC) · Rao et al. 2014
(in situ Hi-C) · Imakaev et al. 2012 (ICE) · Open2C `cooler`, `cooltools`, `pairtools`,
`bioframe` · ENCODE ATAC-seq pipeline specification · 4DN Hi-C processing pipeline
