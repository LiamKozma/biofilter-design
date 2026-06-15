#!/usr/bin/env python3
"""Stage E: posterior-predictive design under uncertainty.

Propagates the calibrated parameter posterior (Stage C) through the industrial
scale-up: for each posterior draw, find the bed length (and hence diameter and
volume, under the 3:1 rule) that achieves 95% hexanal removal. The output is a
credible interval on the required volume -- the honest version of the report's
single 3500 m^3 figure.

    python scripts/04_design_uq.py configs/local_smoke.yaml
"""
from __future__ import annotations

import json
import sys

import numpy as np

from _common import RESULTS, load_config
from biofilter.design import DesignDrawer
from biofilter.pool import get_pool


def main(cfg_path):
    cfg = load_config(cfg_path)
    name = cfg["run_name"]

    chain = np.load(RESULTS / f"calib_chain_{name}.npy")
    # Thin to a manageable, decorrelated set of draws for the design sweep.
    n_draw = min(2000, chain.shape[0])
    idx = np.linspace(0, chain.shape[0] - 1, n_draw).astype(int)
    draws = chain[idx]

    drawer = DesignDrawer(cfg["compounds"], target_compound="Hexanal")
    with get_pool(cfg["calibration"]["backend"]) as pool:
        if not pool.is_master():
            return
        res = np.array(pool.map(drawer, list(draws)))

    L, D, V = res[:, 0], res[:, 1], res[:, 2]
    ok = np.isfinite(V)
    V = V[ok]
    summary = {
        "run_name": name,
        "n_draws": int(ok.sum()),
        "frac_feasible": float(ok.mean()),
        "volume_m3": {
            "median": float(np.median(V)),
            "ci95": [float(np.percentile(V, 2.5)), float(np.percentile(V, 97.5))],
            "ci50": [float(np.percentile(V, 25)), float(np.percentile(V, 75))],
        },
        "length_m_median": float(np.median(L[ok])),
        "diameter_m_median": float(np.median(D[ok])),
        "report_point_estimate_m3": 3500.0,
    }
    (RESULTS / f"design_uq_{name}.json").write_text(json.dumps(summary, indent=2))
    print(json.dumps(summary["volume_m3"], indent=2))
    print(f"feasible draws: {ok.mean():.1%}  -> results/design_uq_{name}.json")


if __name__ == "__main__":
    main(sys.argv[1] if len(sys.argv) > 1 else "configs/local_smoke.yaml")
