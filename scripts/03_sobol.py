#!/usr/bin/env python3
"""Stage D: variance-based global sensitivity analysis (Sobol indices).

Saltelli-samples the mechanistic parameter space and evaluates the industrial-
scale hexanal removal efficiency for every sample, then computes first-order and
total Sobol indices. This answers the design question the original report could
not: of all the transport and kinetic parameters, which ones actually control
removal -- i.e. where pilot measurements would most reduce design risk.

Total forward-model evaluations = n_base * (2D + 2); with the production config
this is ~4e5 coupled BVP solves, an embarrassingly parallel sweep ideal for HPC.

    python scripts/03_sobol.py configs/local_smoke.yaml
"""
from __future__ import annotations

import json
import sys

import numpy as np
from SALib.analyze import sobol as sobol_analyze
from SALib.sample import sobol as sobol_sample

from _common import RESULTS, load_config
from biofilter.design import RemovalAtLength
from biofilter.pool import get_pool

# Parameter ranges in log10 space (bounds bracket the priors / literature).
PROBLEM = {
    "num_vars": 7,
    "names": ["log10_Rmax", "log10_Ks", "log10_Df", "log10_Lf",
              "log10_a", "log10_Dax", "log10_H"],
    "bounds": [
        [-0.5, 1.5],    # Rmax
        [-0.5, 1.5],    # Ks
        [-10.0, -8.0],  # Df
        [-5.0, -3.0],   # Lf
        [2.5, 4.0],     # a
        [-5.0, -3.0],   # Dax
        [0.5, 2.0],     # H
    ],
}

# Evaluate removal at a fixed reference industrial geometry so the index reflects
# parameter influence, not geometry feedback.
REF_LENGTH = 34.23   # m, the report's design length


def main(cfg_path):
    cfg = load_config(cfg_path)
    scfg = cfg["sobol"]
    X = sobol_sample.sample(PROBLEM, scfg["n_base"], calc_second_order=False)
    print(f"Sobol: {X.shape[0]} forward solves ({PROBLEM['num_vars']} params)")

    evaluator = RemovalAtLength(REF_LENGTH)
    with get_pool(scfg["backend"]) as pool:
        if not pool.is_master():
            return
        Y = np.array(pool.map(evaluator, list(X)))

    Si = sobol_analyze.analyze(PROBLEM, Y, calc_second_order=False, print_to_console=False)
    out = {
        "names": PROBLEM["names"],
        "S1": Si["S1"].tolist(),
        "S1_conf": Si["S1_conf"].tolist(),
        "ST": Si["ST"].tolist(),
        "ST_conf": Si["ST_conf"].tolist(),
        "n_eval": int(X.shape[0]),
        "ref_length_m": REF_LENGTH,
    }
    name = cfg["run_name"]
    (RESULTS / f"sobol_{name}.json").write_text(json.dumps(out, indent=2))
    print("\nTotal-order Sobol indices (ST):")
    for n, st, s1 in sorted(zip(PROBLEM["names"], Si["ST"], Si["S1"]),
                            key=lambda t: -t[1]):
        print(f"  {n:12s} ST={st:.3f}  S1={s1:.3f}")
    print(f"-> results/sobol_{name}.json")


if __name__ == "__main__":
    main(sys.argv[1] if len(sys.argv) > 1 else "configs/local_smoke.yaml")
