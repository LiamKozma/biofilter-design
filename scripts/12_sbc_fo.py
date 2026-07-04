#!/usr/bin/env python3
"""v2 Stage F+G: Simulation-Based Calibration for the first-order model.

Proves the v2 inference is statistically calibrated.  For each of ``n_datasets``:
draw ``theta_true`` from the prior, simulate a bench dataset at it, refit with
emcee, and record the rank of each true parameter within its posterior.  Under
correct inference the ranks are uniform; a per-parameter chi-square tests that.

Because the forward model is analytic this runs locally in minutes across the
available cores -- no SLURM array.  Combines the v1 Stage F (array) and Stage G
(reduce) into one script and writes ``sbc_reduce_<name>.json`` directly.

    python scripts/12_sbc_fo.py configs/first_order.yaml
"""
from __future__ import annotations

import json
import sys
import time
from multiprocessing import Pool

import emcee
import numpy as np
from scipy import stats

from _common import RESULTS, load_config
from biofilter import data, firstorder as FO


CFG = None          # populated per-worker via initializer
FO_TEMPLATE = None  # design matrix (Cin, t_s) reused for every synthetic dataset


def _init(cfg, fo_template):
    global CFG, FO_TEMPLATE
    CFG, FO_TEMPLATE = cfg, fo_template


def one_dataset(task: int):
    """Run one SBC replicate; return the per-parameter ranks in [0,1]."""
    scfg = CFG["sbc"]
    rng = np.random.default_rng(CFG["seed"] + 1000 + task)
    theta_true = FO.prior_sample(rng, FO_TEMPLATE)
    sim = FO.simulate(theta_true, FO_TEMPLATE, rng)

    target = FO.LogPosterior(sim)
    ndim = sim.n_params
    nwalk = scfg["inner_walkers"]
    p0 = np.array([FO.prior_sample(rng, sim) for _ in range(nwalk)])
    sampler = emcee.EnsembleSampler(nwalk, ndim, target)
    sampler.run_mcmc(p0, scfg["inner_steps"], progress=False)
    chain = sampler.get_chain(discard=scfg["inner_burn"],
                              thin=scfg.get("inner_thin", 1), flat=True)

    ranks = [float(np.mean(chain[:, j] < theta_true[j])) for j in range(ndim)]
    return ranks


def main(cfg_path):
    cfg = load_config(cfg_path)
    name = cfg["run_name"]
    n = cfg["sbc"]["n_datasets"]

    df = data.load_profiles()
    fo_template = FO.build_fo_data(df, cfg["compounds"])
    names = fo_template.param_names

    t0 = time.time()
    with Pool(initializer=_init, initargs=(cfg, fo_template)) as pool:
        all_ranks = pool.map(one_dataset, range(n))
    ranks = np.array(all_ranks)
    elapsed = time.time() - t0
    print(f"aggregated {ranks.shape[0]} SBC datasets, {ranks.shape[1]} parameters "
          f"in {elapsed:.0f}s\n")

    n_bins = min(10, max(3, n // 20))
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

    out = {"run_name": name, "model": "pseudo_first_order",
           "n_datasets": int(n), "n_bins": int(n_bins), "per_param": report}
    (RESULTS / f"sbc_reduce_{name}.json").write_text(json.dumps(out, indent=2))
    print(f"\n-> results/sbc_reduce_{name}.json")


if __name__ == "__main__":
    main(sys.argv[1] if len(sys.argv) > 1 else "configs/first_order.yaml")
