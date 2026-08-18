#!/usr/bin/env python3
"""Recompute integration statistics and figures from peak_tss_contacts.parquet.

The original run wrote the parquet and the JSON, then failed while plotting:
pairs with zero observed contact give log2(0) = -inf, and matplotlib cannot
compute bin edges over an infinite range.

That -inf is not only a plotting problem. Dropping those pairs silently biases
the enrichment estimate upward, because zero-contact pairs are exactly the ones
with the least contact. This script counts them, reports the count, and
computes the log-ratio statistics on the finite remainder.
"""
from __future__ import annotations

import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats


def main() -> None:
    figs = Path("results/figures")
    figs.mkdir(parents=True, exist_ok=True)

    df = pd.read_parquet("results/peak_tss_contacts.parquet")
    n_total = len(df)

    df["log2_vs_null"] = df["log2_vs_null"].replace([np.inf, -np.inf], np.nan)
    n_zero = int(df["log2_vs_null"].isna().sum())
    df = df.dropna(subset=["log2_vs_null"])
    print(f"{n_total:,} pairs total; {n_zero:,} with zero observed contact "
          f"excluded from log-ratio statistics; {len(df):,} finite")

    # Wilcoxon runs on the raw O/E values, so it is unaffected by the log and
    # uses every pair including the zero-contact ones.
    full = pd.read_parquet("results/peak_tss_contacts.parquet")
    w = stats.wilcoxon(full.oe, full.oe_null)

    within = df.loc[~df.crosses_boundary, "log2_vs_null"]
    across = df.loc[df.crosses_boundary, "log2_vs_null"]
    mw = stats.mannwhitneyu(within, across, alternative="greater")
    rho = stats.spearmanr(df.peak_signal, df.log2_vs_null)

    summary = {
        "n_pairs_total": int(n_total),
        "n_zero_contact_excluded": n_zero,
        "n_pairs_finite": int(len(df)),
        "median_log2_enrichment_vs_distance_matched_null":
            round(float(df.log2_vs_null.median()), 4),
        "wilcoxon_p_on_raw_oe": float(w.pvalue),
        "median_within_boundary": round(float(within.median()), 4),
        "median_across_boundary": round(float(across.median()), 4),
        "n_within_boundary": int(len(within)),
        "n_across_boundary": int(len(across)),
        "boundary_mannwhitney_p": float(mw.pvalue),
        "spearman_peak_signal_vs_enrichment": {
            "rho": round(float(rho.statistic), 4),
            "p": float(rho.pvalue),
            "note": ("n is large enough that a negligible rho is highly "
                     "significant; the effect size is what matters"),
        },
    }
    Path("results/integration.json").write_text(json.dumps(summary, indent=2))
    print(json.dumps(summary, indent=2))

    fig, ax = plt.subplots(1, 2, figsize=(10, 4))
    ax[0].hist(df.log2_vs_null, bins=60, color="#1F3864")
    ax[0].axvline(0, c="grey", ls=":")
    ax[0].axvline(df.log2_vs_null.median(), c="#B23A48", ls="-", lw=1.2,
                  label=f"median {df.log2_vs_null.median():.2f}")
    ax[0].legend(frameon=False, fontsize=8)
    ax[0].set(xlabel="log2(O/E vs distance-matched null)",
              ylabel="peak-TSS pairs", title="contact enrichment")

    ax[1].boxplot([within.values, across.values],
                  tick_labels=[f"within\n(n={len(within):,})",
                               f"across\n(n={len(across):,})"],
                  showfliers=False)
    ax[1].axhline(0, c="grey", ls=":")
    ax[1].set(ylabel="log2 enrichment",
              title="insulation boundary between peak and TSS")
    fig.tight_layout()
    fig.savefig(figs / "peak_tss_enrichment.png", dpi=150)
    print("wrote results/figures/peak_tss_enrichment.png")


if __name__ == "__main__":
    main()
