#!/usr/bin/env python3
"""ATAC-seq QC: TSS enrichment, FRiP, fragment-size periodicity.

These three metrics replace the RNA-seq QC panel entirely. There is no
mitochondrial-percentage / genes-detected analogue that transfers -- ATAC
quality is judged on whether the insertion pattern still carries nucleosome
structure and still concentrates at promoters.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import pysam
import yaml


def load_tss(gtf: Path, chroms: list[str], protein_coding_only: bool = True) -> pd.DataFrame:
    """One TSS per gene, strand-aware, from a GENCODE GTF."""
    rows = []
    opener = __import__("gzip").open if str(gtf).endswith(".gz") else open
    with opener(gtf, "rt") as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if f[2] != "gene" or f[0] not in chroms:
                continue
            if protein_coding_only and 'gene_type "protein_coding"' not in f[8]:
                continue
            tss = int(f[3]) - 1 if f[6] == "+" else int(f[4]) - 1
            gid = f[8].split('gene_id "')[1].split('"')[0]
            rows.append((f[0], tss, f[6], gid))
    return pd.DataFrame(rows, columns=["chrom", "tss", "strand", "gene_id"])


def tss_enrichment(cutsites_tbx: str, tss: pd.DataFrame, flank: int,
                   core: int, flank_norm: int) -> tuple[float, np.ndarray]:
    """Mean insertion density in the TSS core over the density in the far flanks.

    ENCODE's hg38 threshold is >=5 (ideal >=7). The profile is strand-flipped so
    that minus-strand genes contribute in the same orientation.
    """
    tbx = pysam.TabixFile(cutsites_tbx)
    width = 2 * flank
    profile = np.zeros(width, dtype=np.float64)
    contigs = set(tbx.contigs)

    for chrom, pos, strand, _ in tss.itertuples(index=False):
        start, end = pos - flank, pos + flank
        if start < 0 or chrom not in contigs:
            continue
        try:
            hits = tbx.fetch(chrom, start, end)
        except ValueError:
            continue
        idx = [int(h.split("\t")[1]) - start for h in hits]
        idx = [i for i in idx if 0 <= i < width]
        if not idx:
            continue
        v = np.bincount(idx, minlength=width).astype(np.float64)
        profile += v[::-1] if strand == "-" else v

    # Background from the outer `flank_norm` bp at each end, where signal should
    # have decayed to genomic baseline.
    bg = np.concatenate([profile[:flank_norm], profile[-flank_norm:]]).mean()
    if bg <= 0:
        raise ValueError("Zero background in TSS flanks -- check chromosome naming.")
    norm = profile / bg
    score = float(norm[flank - core: flank + core].mean())
    return score, norm


def frip(bam: str, peaks: str) -> float:
    """Fraction of filtered reads in peaks. >0.2 is a reasonable bulk target."""
    import subprocess
    total = int(subprocess.run(["samtools", "view", "-c", bam],
                               capture_output=True, text=True, check=True).stdout)
    inpeak = subprocess.run(
        f"samtools view -c -L {peaks} {bam}", shell=True,
        capture_output=True, text=True, check=True)
    return int(inpeak.stdout) / total if total else float("nan")


def fragment_sizes(bam: str, max_size: int = 1000, max_reads: int = 5_000_000) -> np.ndarray:
    """Insert-size histogram from properly paired plus-strand reads."""
    counts = np.zeros(max_size + 1, dtype=np.int64)
    n = 0
    with pysam.AlignmentFile(bam) as af:
        for r in af:
            if not r.is_proper_pair or r.is_reverse or r.template_length <= 0:
                continue
            if r.template_length <= max_size:
                counts[r.template_length] += 1
            n += 1
            if n >= max_reads:
                break
    return counts


def nucleosome_signal(counts: np.ndarray) -> float:
    """Mono-nucleosomal (180-247 bp) over nucleosome-free (<100 bp) mass.

    A featureless smear -- no sub-100 bp peak, no ~200 bp shoulder -- means the
    library is over-transposed or the nuclei ruptured. No downstream analysis
    rescues that, which is why this is a pre-sequencing bench gate as well.
    """
    nfr = counts[1:100].sum()
    mono = counts[180:248].sum()
    return float(mono / nfr) if nfr else float("nan")


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--config", default="config/config.yaml")
    ap.add_argument("--sample", default="GM12878")
    args = ap.parse_args()

    cfg = yaml.safe_load(open(args.config))
    a, g = cfg["atac"], cfg["genome"]
    out = Path("results/atac")
    figs = Path("results/figures")
    figs.mkdir(parents=True, exist_ok=True)

    bam = str(out / f"{args.sample}.dedup.bam")
    peaks = str(out / f"{args.sample}_peaks.narrowPeak")
    cuts = str(out / f"{args.sample}.cutsites.bed.gz")

    tss = load_tss(Path(g["gtf"]), g["main_chroms"])
    print(f"{len(tss):,} protein-coding TSS")

    score, profile = tss_enrichment(cuts, tss, a["tss_flank"], a["tss_core"],
                                    a["tss_flank_norm"])
    f = frip(bam, peaks)
    sizes = fragment_sizes(bam)
    nuc = nucleosome_signal(sizes)
    n_peaks = sum(1 for _ in open(peaks))

    metrics = {
        "sample": args.sample,
        "tss_enrichment": round(score, 3),
        "tss_enrichment_pass": score >= 5,
        "frip": round(f, 4),
        "frip_pass": f >= 0.2,
        "nucleosome_signal": round(nuc, 3),
        "n_peaks": n_peaks,
    }
    Path("results/atac_qc.json").write_text(json.dumps(metrics, indent=2))
    print(json.dumps(metrics, indent=2))

    fig, ax = plt.subplots(1, 2, figsize=(11, 4))
    x = np.arange(-a["tss_flank"], a["tss_flank"])
    ax[0].plot(x, profile, lw=1)
    ax[0].axhline(1, ls=":", c="grey", lw=0.8)
    ax[0].set(xlabel="distance from TSS (bp)",
              ylabel="insertions / background",
              title=f"TSS enrichment = {score:.2f}")
    ax[1].plot(np.arange(len(sizes)), sizes, lw=1)
    for p in (200, 400, 600):
        ax[1].axvline(p, ls=":", c="grey", lw=0.8)
    ax[1].set(xlabel="fragment length (bp)", ylabel="count", yscale="log",
              title=f"nucleosome signal = {nuc:.2f}")
    fig.tight_layout()
    fig.savefig(figs / f"{args.sample}_atac_qc.png", dpi=150)


if __name__ == "__main__":
    main()
