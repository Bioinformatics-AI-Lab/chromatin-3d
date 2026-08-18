#!/usr/bin/env python3
"""Hi-C feature calling: distance decay, A/B compartments, insulation.

The organising fact of Hi-C analysis is that contact probability decays as a
power law over several orders of magnitude in genomic separation. That decay is
the dominant source of variance in the matrix, so nearly every feature caller
works on observed/expected rather than observed. There is no comparable
nuisance covariate in RNA-seq -- getting this step wrong produces results that
are really just a restatement of polymer physics.
"""
from __future__ import annotations

import argparse
import json
import warnings
from pathlib import Path

import bioframe
import cooler
import cooltools
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import yaml

# These come from pandas groupby calls inside cooltools/bioframe, not from this
# code, and there is nothing to fix on our side until those packages update.
warnings.filterwarnings("ignore", category=FutureWarning,
                        message=".*observed=False.*")


def jsonable(o):
    """Coerce numpy scalars to Python types.

    cooler and pandas return numpy int64/float64, which json.dumps rejects.
    Wrapping each value individually is easy to get wrong -- one missed field
    fails the whole write after the expensive work is already done.
    """
    if isinstance(o, (np.integer,)):
        return int(o)
    if isinstance(o, (np.floating,)):
        return float(o)
    if isinstance(o, (np.bool_,)):
        return bool(o)
    if isinstance(o, dict):
        return {str(k): jsonable(v) for k, v in o.items()}
    if isinstance(o, (list, tuple)):
        return [jsonable(v) for v in o]
    return o


def open_clr(mcool: str, resolution: int) -> cooler.Cooler:
    clr = cooler.Cooler(f"{mcool}::/resolutions/{resolution}")
    if "weight" not in clr.bins().columns:
        raise RuntimeError(
            "No 'weight' column -- the matrix is unbalanced. Run "
            "`cooler balance` first. ICE rescales rows and columns iteratively "
            "until every bin has equal marginal coverage; without it, "
            "mappability and restriction-site density dominate."
        )
    return clr


def make_view(clr: cooler.Cooler, chroms: list[str]) -> pd.DataFrame:
    sizes = clr.chromsizes
    keep = [c for c in chroms if c in sizes.index]
    return bioframe.make_viewframe({c: sizes[c] for c in keep})


def distance_decay(clr, view, nproc, figs: Path) -> pd.DataFrame:
    """P(s) curve. The slope near ~-1 over 10^5-10^7 bp is the polymer signature."""
    exp = cooltools.expected_cis(clr, view_df=view, smooth=True,
                                 aggregate_smoothed=True, nproc=nproc)
    agg = (exp.groupby("dist")["balanced.avg.smoothed.agg"].mean().dropna())
    s = agg.index.values * clr.binsize
    m = s > 0
    fig, ax = plt.subplots(figsize=(5, 4.5))
    ax.loglog(s[m], agg.values[m], lw=1.5)
    ax.set(xlabel="genomic separation (bp)", ylabel="P(contact)",
           title="distance decay")
    fig.tight_layout()
    fig.savefig(figs / "hic_distance_decay.png", dpi=150)
    return exp


def compartments(clr_lo, view, fasta: str, figs: Path) -> pd.DataFrame:
    """A/B compartments via eigendecomposition, phased by GC content.

    The sign of an eigenvector is arbitrary, so PC1 flips between chromosomes
    unless it is phased against a track that is known to track the A
    compartment. GC content is the standard choice; gene density also works.
    Reporting unphased eigenvectors is a common and quietly serious error.
    """
    bins = clr_lo.bins()[:]
    gc = bioframe.frac_gc(bins[["chrom", "start", "end"]],
                          bioframe.load_fasta(fasta))
    _, eigvecs = cooltools.eigs_cis(clr_lo, phasing_track=gc, view_df=view,
                                    n_eigs=3)
    e1 = eigvecs[["chrom", "start", "end", "E1"]].dropna()
    frac_a = float((e1["E1"] > 0).mean())

    fig, ax = plt.subplots(figsize=(11, 2.6))
    sub = e1[e1.chrom == view.chrom.iloc[0]]
    ax.fill_between(sub.start, sub.E1, 0, where=sub.E1 > 0, color="#B23A48", lw=0)
    ax.fill_between(sub.start, sub.E1, 0, where=sub.E1 <= 0, color="#1F3864", lw=0)
    ax.set(xlabel=f"{view.chrom.iloc[0]} position (bp)", ylabel="E1",
           title="A/B compartments (GC-phased)")
    fig.tight_layout()
    fig.savefig(figs / "hic_compartments.png", dpi=150)
    return e1.assign(_frac_a=frac_a)


def insulation(clr, view, windows: list[int], nproc, figs: Path) -> pd.DataFrame:
    """Insulation score at several window sizes.

    Deliberately multi-scale. TAD boundary calls are strongly method- and
    scale-dependent and concordance between callers is poor, so presenting a
    single window as ground truth claims more than the data supports.
    """
    ins = cooltools.insulation(clr, windows, view_df=view, nproc=nproc,
                               verbose=False)
    counts = {}
    for w in windows:
        col = f"is_boundary_{w}"
        counts[w] = int(ins[col].sum()) if col in ins else 0

    fig, ax = plt.subplots(figsize=(11, 3))
    sub = ins[ins.chrom == view.chrom.iloc[0]]
    for w in windows:
        c = f"log2_insulation_score_{w}"
        if c in sub:
            ax.plot(sub.start, sub[c], lw=0.8, label=f"{w // 1000} kb")
    ax.legend(frameon=False, fontsize=8)
    ax.set(xlabel=f"{view.chrom.iloc[0]} position (bp)",
           ylabel="log2 insulation", title="insulation profile")
    fig.tight_layout()
    fig.savefig(figs / "hic_insulation.png", dpi=150)
    return ins.assign(**{f"_n_boundaries_{w}": v for w, v in counts.items()})


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--config", default="config/config.yaml")
    ap.add_argument("--mcool", default="data/hic/GM12878.mcool")
    args = ap.parse_args()

    cfg = yaml.safe_load(open(args.config))
    h, g = cfg["hic"], cfg["genome"]
    nproc = cfg.get("threads", 4)
    out = Path("results/hic"); out.mkdir(parents=True, exist_ok=True)
    figs = Path("results/figures"); figs.mkdir(parents=True, exist_ok=True)

    clr = open_clr(args.mcool, h["resolution"])
    clr_lo = open_clr(args.mcool, h["compartment_resolution"])
    view = make_view(clr, g["main_chroms"])
    view_lo = make_view(clr_lo, g["main_chroms"])
    print(f"{clr.binsize:,} bp bins, {clr.info['nbins']:,} bins, "
          f"{clr.info['nnz']:,} non-zero pixels")

    exp = distance_decay(clr, view, nproc, figs)
    exp.to_parquet(out / "expected_cis.parquet")

    e1 = compartments(clr_lo, view_lo, g["fasta"], figs)
    e1.drop(columns=["_frac_a"]).to_parquet(out / "compartments.parquet")

    ins = insulation(clr, view, h["insulation_windows"], nproc, figs)
    ins[[c for c in ins.columns if not c.startswith("_")]].to_parquet(
        out / "insulation.parquet")

    summary = {
        "resolution": clr.binsize,
        "n_bins": int(clr.info["nbins"]),
        "nnz_pixels": int(clr.info["nnz"]),
        "frac_A_compartment": round(float(e1["_frac_a"].iloc[0]), 3),
        "n_boundaries": {str(w): int(ins[f"_n_boundaries_{w}"].iloc[0])
                         for w in h["insulation_windows"]},
    }
    summary = jsonable(summary)
    Path("results/hic_qc.json").write_text(json.dumps(summary, indent=2))
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
