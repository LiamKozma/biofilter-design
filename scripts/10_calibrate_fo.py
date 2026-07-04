#!/usr/bin/env python3
"""v2 Stage C: Bayesian calibration of the identifiable first-order model.

Infers one pseudo-first-order rate constant ``k1`` (1/s) and one observation
noise ``sigma`` (ppmv) per compound.  The forward model is analytic, so this runs
in seconds on a single core -- no MPI, no HDF5 checkpoint.  Writes the flattened
posterior chain and a convergence summary in the same shape as Stage C v1.

    python scripts/10_calibrate_fo.py configs/first_order.yaml
"""
from __future__ import annotations

import json
import sys
import time

import emcee
import numpy as np

from _common import RESULTS, load_config
from biofilter import data, firstorder as FO


def main(cfg_path):
    cfg = load_config(cfg_path)
    ccfg = cfg["calibration"]
    rng = np.random.default_rng(cfg["seed"])
    name = cfg["run_name"]

    df = data.load_profiles()
    fo = FO.build_fo_data(df, cfg["compounds"])
    target = FO.LogPosterior(fo)
    ndim = fo.n_params
    nwalk = ccfg["n_walkers"]

    p0 = np.array([FO.prior_sample(rng, fo) for _ in range(nwalk)])

    t0 = time.time()
    sampler = emcee.EnsembleSampler(nwalk, ndim, target)
    sampler.run_mcmc(p0, ccfg["n_steps"], progress=True)
    elapsed = time.time() - t0

    chain = sampler.get_chain(discard=ccfg["n_burn"], flat=True)
    try:
        tau = sampler.get_autocorr_time(quiet=True)
        ess = float((ccfg["n_steps"] - ccfg["n_burn"]) * nwalk / np.nanmax(tau))
    except Exception:
        ess = float("nan")
    acc = float(np.mean(sampler.acceptance_fraction))

    np.save(RESULTS / f"calib_chain_{name}.npy", chain)

    # Also report k1 in linear units and, for the design compound, the implied
    # required industrial volume so the summary is directly interpretable.
    lin = {}
    for i, nm in enumerate(fo.param_names):
        col = np.exp(chain[:, i])
        lin[nm.replace("log_", "")] = {
            "median": float(np.median(col)),
            "ci95": [float(np.percentile(col, 2.5)), float(np.percentile(col, 97.5))],
        }

    summary = {
        "run_name": name,
        "model": "pseudo_first_order",
        "compounds": fo.compounds,
        "param_names": fo.param_names,
        "n_walkers": nwalk,
        "n_steps": ccfg["n_steps"],
        "n_burn": ccfg["n_burn"],
        "elapsed_s": elapsed,
        "acceptance_fraction": acc,
        "ess": ess,
        "posterior_median": {n: float(np.median(chain[:, i]))
                             for i, n in enumerate(fo.param_names)},
        "posterior_ci95": {n: [float(np.percentile(chain[:, i], 2.5)),
                               float(np.percentile(chain[:, i], 97.5))]
                           for i, n in enumerate(fo.param_names)},
        "linear_units": lin,
    }
    (RESULTS / f"calib_summary_{name}.json").write_text(json.dumps(summary, indent=2))
    print(f"\nacc={acc:.2f} ESS~{ess:.0f} elapsed={elapsed:.1f}s"
          f" -> results/calib_chain_{name}.npy")
    for c in fo.compounds:
        k = lin[f"k1[{c}]"]
        print(f"  k1[{c:8s}] = {k['median']:.4f} 1/s  "
              f"[{k['ci95'][0]:.4f}, {k['ci95'][1]:.4f}]")


if __name__ == "__main__":
    main(sys.argv[1] if len(sys.argv) > 1 else "configs/first_order.yaml")
