#!/usr/bin/env python3
"""Stage C: PDE-in-the-loop Bayesian calibration of the two-scale model.

Runs an affine-invariant ensemble sampler (emcee) whose every walker move
triggers a coupled biofilm+column BVP solve per compound and run. Parallelised
over walkers with multiprocessing locally or MPI on Sapelo2. Writes the flattened
posterior chain and a convergence summary.

    python scripts/02_calibrate.py configs/local_smoke.yaml
    srun python scripts/02_calibrate.py configs/sapelo2_full.yaml   # on Sapelo2
"""
from __future__ import annotations

import json
import sys
import time

import emcee
import numpy as np

from _common import RESULTS, load_config, operating_uL
from biofilter import data, likelihood as L
from biofilter.pool import get_pool


def main(cfg_path):
    cfg = load_config(cfg_path)
    ccfg = cfg["calibration"]
    rng = np.random.default_rng(cfg["seed"])

    u, length = operating_uL()
    df = data.load_profiles()
    cd = L.build_calibration_data(df, cfg["compounds"], length, u)
    target = L.LogPosterior(cd)          # picklable across mp/MPI workers
    ndim = cd.n_params
    nwalk = ccfg["n_walkers"]

    # Initialise walkers from the prior, lightly jittered.
    p0 = np.array([L.prior_sample(rng, cd) for _ in range(nwalk)])

    with get_pool(ccfg["backend"]) as pool:
        if not pool.is_master():
            return
        sampler = emcee.EnsembleSampler(nwalk, ndim, target, pool=pool)
        t0 = time.time()
        sampler.run_mcmc(p0, ccfg["n_steps"], progress=True)
        elapsed = time.time() - t0

    chain = sampler.get_chain(discard=ccfg["n_burn"], flat=True)
    try:
        tau = sampler.get_autocorr_time(quiet=True)
        ess = float((ccfg["n_steps"] - ccfg["n_burn"]) * nwalk / np.nanmax(tau))
    except Exception:
        tau, ess = np.full(ndim, np.nan), float("nan")
    acc = float(np.mean(sampler.acceptance_fraction))

    name = cfg["run_name"]
    np.save(RESULTS / f"calib_chain_{name}.npy", chain)
    summary = {
        "run_name": name,
        "compounds": cd.compounds,
        "param_names": cd.param_names,
        "n_walkers": nwalk,
        "n_steps": ccfg["n_steps"],
        "n_burn": ccfg["n_burn"],
        "elapsed_s": elapsed,
        "acceptance_fraction": acc,
        "ess": ess,
        "posterior_median": {n: float(np.median(chain[:, i]))
                             for i, n in enumerate(cd.param_names)},
        "posterior_ci95": {n: [float(np.percentile(chain[:, i], 2.5)),
                               float(np.percentile(chain[:, i], 97.5))]
                           for i, n in enumerate(cd.param_names)},
    }
    (RESULTS / f"calib_summary_{name}.json").write_text(json.dumps(summary, indent=2))
    print(f"\nacc={acc:.2f} ESS~{ess:.0f} elapsed={elapsed:.1f}s "
          f"-> results/calib_chain_{name}.npy")


if __name__ == "__main__":
    main(sys.argv[1] if len(sys.argv) > 1 else "configs/local_smoke.yaml")
