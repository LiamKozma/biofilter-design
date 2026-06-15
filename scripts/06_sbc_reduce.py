#!/usr/bin/env python3
"""Aggregate SBC array outputs and test rank uniformity.

Gathers the per-task rank files written by 05_sbc.py, then for each parameter
runs a chi-square goodness-of-fit test against the uniform distribution expected
under correct inference. A small p-value flags a parameter whose posterior is
biased or mis-calibrated.

    python scripts/06_sbc_reduce.py configs/sapelo2_full.yaml
"""
from __future__ import annotations

import glob
import json
import sys

import numpy as np
from scipy import stats

from _common import RESULTS, load_config


def main(cfg_path):
    cfg = load_config(cfg_path)
    name = cfg["run_name"]

    files = sorted(glob.glob(str(RESULTS / f"sbc_{name}_task*.json")))
    if not files:
        files = sorted(glob.glob(str(RESULTS / f"sbc_{name}.json")))
    if not files:
        print("no SBC outputs found")
        return

    ranks, names = [], None
    for f in files:
        d = json.loads(open(f).read())
        names = d["param_names"]
        ranks.extend(d["ranks"])
    ranks = np.array(ranks)            # (n_datasets, n_params), values in [0,1]
    n = ranks.shape[0]
    print(f"aggregated {n} SBC datasets, {ranks.shape[1]} parameters\n")

    n_bins = min(10, max(3, n // 5))
    report = {}
    print(f"{'parameter':18s} {'chi2':>8s} {'p_unif':>8s}  verdict")
    for j, nm in enumerate(names):
        counts, _ = np.histogram(ranks[:, j], bins=n_bins, range=(0, 1))
        expected = n / n_bins
        chi2 = float(np.sum((counts - expected) ** 2 / expected))
        p = float(stats.chi2.sf(chi2, n_bins - 1))
        verdict = "ok" if p > 0.05 else "CHECK"
        print(f"{nm:18s} {chi2:8.2f} {p:8.3f}  {verdict}")
        report[nm] = {"chi2": chi2, "p_uniform": p, "counts": counts.tolist()}

    out = {"n_datasets": n, "n_bins": n_bins, "per_param": report}
    (RESULTS / f"sbc_reduce_{name}.json").write_text(json.dumps(out, indent=2))
    print(f"\n-> results/sbc_reduce_{name}.json")


if __name__ == "__main__":
    main(sys.argv[1] if len(sys.argv) > 1 else "configs/local_smoke.yaml")
