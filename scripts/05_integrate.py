#!/usr/bin/env python3
"""Do accessible peaks contact expressed promoters more than chance?

This is the part that makes the project an analysis rather than a pipeline run.

The trap: contact frequency falls off steeply with genomic separation, and ATAC
peaks cluster near promoters, so peak-TSS pairs are systematically *closer*
than randomly chosen pairs. Comparing raw contact frequency between the two
therefore recovers the distance-decay curve and calls it biology. Every claim
below is made against a **distance-matched null** built at the same separation
in the same chromosome, which is the only version of the comparison that means
anything.

Secondary question: does an insulation boundary between a peak and a promoter
reduce their contact? If the boundaries called in 04_hic_features.py are real,
it should.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import cooler
import cooltools
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import yaml
from scipy import stats


def load_peaks(path: str, chroms: list[str]) -> pd.DataFrame:
    cols = ["chrom", "start", "end", "name", "score", "strand",
            "signal", "pval", "qval", "summit"]
    df = pd.read_csv(path, sep="\t", header=None, names=cols)
    df = df[df.chrom.isin(chroms)].copy()
    df["point"] = df.start + df.summit          # summit, not midpoint
    return df[["chrom", "start", "end", "point", "signal"]]


def load_expressed_tss(gtf: str, rna_tsv: str | None, chroms: list[str],
                       tpm_min: float) -> pd.DataFrame:
    """Protein-coding TSS, restricted to expressed genes when RNA is available."""
    import gzip
    rows = []
    opener = gzip.open if gtf.endswith(".gz") else open
    with opener(gtf, "rt") as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if f[2] != "gene" or f[0] not in chroms:
                continue
            if 'gene_type "protein_coding"' not in f[8]:
                continue
            tss = int(f[3]) - 1 if f[6] == "+" else int(f[4]) - 1
            gid = f[8].split('gene_id "')[1].split('"')[0]
            rows.append((f[0], tss, gid.split(".")[0]))
    tss = pd.DataFrame(rows, columns=["chrom", "tss", "gene_id"])

    if rna_tsv:
        rna = pd.read_csv(rna_tsv, sep="\t")
        col = next(c for c in rna.columns if c.upper() == "TPM")
        idc = next(c for c in rna.columns if "gene_id" in c.lower())
        rna[idc] = rna[idc].astype(str).str.split(".").str[0]
        rna = rna.groupby(idc)[col].max()
        tss["tpm"] = tss.gene_id.map(rna)
        tss = tss[tss.tpm >= tpm_min]
    else:
        tss["tpm"] = np.nan
        print("No RNA quantification configured -- using all protein-coding TSS. "
              "Set rna.encode_quantification in config.yaml to filter.")
    return tss


def oe_lookup(clr, exp_by_chrom: dict, chrom: str, b1: int, b2: int,
              mat_cache: dict) -> float:
    """Observed/expected contact between two bins of the same chromosome."""
    if chrom not in mat_cache:
        mat_cache[chrom] = clr.matrix(balance=True, sparse=True).fetch(chrom).tocsr()
    m = mat_cache[chrom]
    if b1 >= m.shape[0] or b2 >= m.shape[0]:
        return np.nan
    obs = m[b1, b2]
    e = exp_by_chrom[chrom].get(abs(b2 - b1), np.nan)
    if not np.isfinite(e) or e <= 0:
        return np.nan
    return float(obs) / e


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--config", default="config/config.yaml")
    ap.add_argument("--mcool", default="data/hic/GM12878.mcool")
    ap.add_argument("--peaks", default="results/atac/GM12878_peaks.narrowPeak")
    args = ap.parse_args()

    cfg = yaml.safe_load(open(args.config))
    g, h, r, ic = cfg["genome"], cfg["hic"], cfg["rna"], cfg["integrate"]
    rng = np.random.default_rng(ic["random_seed"])
    res = h["resolution"]
    figs = Path("results/figures"); figs.mkdir(parents=True, exist_ok=True)

    clr = cooler.Cooler(f"{args.mcool}::/resolutions/{res}")
    peaks = load_peaks(args.peaks, g["main_chroms"])
    tss = load_expressed_tss(g["gtf"], r.get("encode_quantification"),
                             g["main_chroms"], r["expressed_tpm_min"])
    print(f"{len(peaks):,} peaks, {len(tss):,} TSS")

    exp = pd.read_parquet("results/hic/expected_cis.parquet")
    exp_by_chrom = {
        c: dict(zip(d["dist"], d["balanced.avg.smoothed.agg"]))
        for c, d in exp.assign(chrom=exp.region1.str.split(":").str[0]).groupby("chrom")
    }
    ins = pd.read_parquet("results/hic/insulation.parquet")
    bcol = f"is_boundary_{h['insulation_windows'][1]}"
    boundaries = {c: np.sort(d.loc[d[bcol].fillna(False), "start"].values)
                  for c, d in ins.groupby("chrom", observed=True)} if bcol in ins else {}

    offsets = (clr.bins()[:].reset_index().groupby("chrom", observed=True)["index"]
               .min().to_dict())
    mat_cache: dict = {}
    records = []

    for chrom, ptab in peaks.groupby("chrom"):
        ttab = tss[tss.chrom == chrom]
        if ttab.empty or chrom not in offsets:
            continue
        off = offsets[chrom]
        tpos = ttab.tss.values
        for point, signal in zip(ptab.point.values, ptab.signal.values):
            d = np.abs(tpos - point)
            sel = (d >= ic["min_peak_tss_distance"]) & (d <= ic["max_peak_tss_distance"])
            if not sel.any():
                continue
            for t, tpm in zip(tpos[sel], ttab.tpm.values[sel]):
                b1, b2 = point // res, t // res
                if b1 == b2:
                    continue
                oe = oe_lookup(clr, exp_by_chrom, chrom, b1, b2, mat_cache)
                if not np.isfinite(oe):
                    continue
                sep = abs(b2 - b1)

                # Distance-matched null: same chromosome, same bin separation,
                # anchored at a random position rather than at a peak.
                null = []
                nb = mat_cache[chrom].shape[0]
                for _ in range(ic["n_null_per_pair"]):
                    a = int(rng.integers(0, max(nb - sep - 1, 1)))
                    v = oe_lookup(clr, exp_by_chrom, chrom, a, a + sep, mat_cache)
                    if np.isfinite(v):
                        null.append(v)
                if not null:
                    continue

                lo, hi = sorted((min(point, t), max(point, t)))
                bnd = boundaries.get(chrom, np.array([]))
                crosses = bool(((bnd > lo) & (bnd < hi)).any())

                records.append({
                    "chrom": chrom, "peak": point, "tss": t, "tpm": tpm,
                    "peak_signal": signal, "separation_bp": sep * res,
                    "oe": oe, "oe_null": float(np.mean(null)),
                    "log2_vs_null": float(np.log2(oe / np.mean(null)))
                    if np.mean(null) > 0 else np.nan,
                    "crosses_boundary": crosses,
                })

    df = pd.DataFrame(records)

    # A pair with zero observed contact gives log2(0) = -inf. Dropping those
    # silently would bias the enrichment estimate upward, since they are exactly
    # the pairs with the least contact. Count them and report the count next to
    # the statistic computed on the finite remainder.
    df["log2_vs_null"] = df["log2_vs_null"].replace([np.inf, -np.inf], np.nan)
    n_zero = int(df["log2_vs_null"].isna().sum())
    df.to_parquet("results/peak_tss_contacts.parquet")
    full = df
    df = df.dropna(subset=["log2_vs_null"])
    print(f"{len(df):,} peak-TSS pairs with finite enrichment; "
          f"{n_zero:,} zero-contact pairs excluded from log-ratio statistics")

    # Wilcoxon runs on raw O/E, so it is unaffected by the log and uses every
    # pair, including the zero-contact ones.
    w = stats.wilcoxon(full.oe, full.oe_null)
    within = df.loc[~df.crosses_boundary, "log2_vs_null"]
    across = df.loc[df.crosses_boundary, "log2_vs_null"]
    mw = stats.mannwhitneyu(within, across, alternative="greater")
    rho = stats.spearmanr(df.peak_signal, df.log2_vs_null)

    summary = {
        "n_pairs_finite": int(len(df)),
        "n_zero_contact_excluded": n_zero,
        "median_log2_enrichment_vs_distance_matched_null":
            round(float(df.log2_vs_null.median()), 4),
        "wilcoxon_p": float(w.pvalue),
        "median_within_boundary": round(float(within.median()), 4),
        "median_across_boundary": round(float(across.median()), 4),
        "boundary_mannwhitney_p": float(mw.pvalue),
        "spearman_peak_signal_vs_enrichment": {
            "rho": round(float(rho.statistic), 4), "p": float(rho.pvalue),
            "note": "n is large enough that a negligible rho is highly "
                    "significant; the effect size is what matters"},
    }
    Path("results/integration.json").write_text(json.dumps(summary, indent=2))
    print(json.dumps(summary, indent=2))

    fig, ax = plt.subplots(1, 2, figsize=(10, 4))
    ax[0].hist(df.log2_vs_null, bins=60, color="#1F3864")
    ax[0].axvline(0, c="grey", ls=":")
    ax[0].set(xlabel="log2(O/E vs distance-matched null)", ylabel="peak-TSS pairs",
              title="contact enrichment")
    ax[1].boxplot([within.values, across.values],
                  tick_labels=[f"within\n(n={len(within):,})",
                               f"across\n(n={len(across):,})"],
                  showfliers=False)
    ax[1].axhline(0, c="grey", ls=":")
    ax[1].set(ylabel="log2 enrichment", title="insulation boundary between peak and TSS")
    fig.tight_layout()
    fig.savefig(figs / "peak_tss_enrichment.png", dpi=150)


if __name__ == "__main__":
    main()
