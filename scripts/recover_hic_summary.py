#!/usr/bin/env python3
"""Rebuild results/hic_qc.json from the parquet outputs of 04_hic_features.py.

Exists because the original run wrote every parquet and figure and then failed
on the JSON summary (numpy int64 is not JSON-serializable). The expensive work
was already on disk, so recomputing the summary is cheap. Safe to delete once
04_hic_features.py has been re-run end to end.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import cooler
import numpy as np
import pandas as pd
import yaml


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--config", default="config/config.yaml")
    ap.add_argument("--mcool", default="data/hic/GM12878.mcool")
    args = ap.parse_args()

    cfg = yaml.safe_load(open(args.config))
    h = cfg["hic"]
    out = Path("results/hic")

    clr = cooler.Cooler(f"{args.mcool}::/resolutions/{h['resolution']}")

    e1 = pd.read_parquet(out / "compartments.parquet")
    frac_a = float((e1["E1"] > 0).mean())

    ins = pd.read_parquet(out / "insulation.parquet")
    n_boundaries = {}
    for w in h["insulation_windows"]:
        col = f"is_boundary_{w}"
        n_boundaries[str(w)] = int(ins[col].fillna(False).sum()) if col in ins else 0

    summary = {
        "resolution": int(clr.binsize),
        "n_bins": int(clr.info["nbins"]),
        "nnz_pixels": int(clr.info["nnz"]),
        "frac_A_compartment": round(frac_a, 3),
        "n_boundaries": n_boundaries,
    }
    Path("results/hic_qc.json").write_text(json.dumps(summary, indent=2))
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
