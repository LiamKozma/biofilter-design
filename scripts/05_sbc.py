#!/usr/bin/env python3
"""Stage F: Simulation-Based Calibration (SBC) -- validating the inference.

A Bayesian result is only trustworthy if the inference machinery is itself
calibrated. SBC tests exactly that: draw parameters from the prior, simulate a
dataset, run the *entire* calibration, and record the rank of each true value
within its posterior. If the sampler and model are correct, those ranks are
uniformly distributed. Systematic deviations (∪ or ∩ shaped rank histograms)
expose bias or mis-estimated uncertainty.

This is the workload that genuinely requires HPC: each of the N datasets is a
full PDE-in-the-loop MCMC. Run as a SLURM job array -- one task per dataset --
then aggregate the per-task rank files.

    # one dataset (array task):
    python scripts/05_sbc.py configs/sapelo2_full.yaml --task $SLURM_ARRAY_TASK_ID
    # all datasets locally (smoke):
    python scripts/05_sbc.py configs/local_smoke.yaml
"""
from __future__ import annotations

import argparse
import json

import emcee
import numpy as np

from _common import RESULTS, load_config, operating_uL
from biofilter import data, likelihood as L
from biofilter.simulator import forward_profile
from biofilter.design import params_from_theta
from biofilter.pool import get_pool


def simulate_dataset(theta_true, template: L.CalibrationData, rng) -> L.CalibrationData:
    """Generate synthetic observations from ``theta_true`` on the real design."""
    sigma = np.exp(theta_true[5])
    per = {}
    for comp in template.compounds:
        p = params_from_theta(theta_true, template.compounds, target_compound=comp)
        runs = []
        for Cin, pos, _y in template.per_compound[comp]:
            pred, ok = forward_profile(p, Cin, pos, template.length, template.u)
            if not ok:
                pred = np.full_like(np.asarray(pos, float), Cin)
            noisy = np.clip(pred + rng.normal(0, sigma, size=len(pred)), 0.0, None)
            runs.append((Cin, np.asarray(pos), noisy))
        per[comp] = runs
    return L.CalibrationData(template.compounds, template.length, template.u, per)


def run_one(idx, cfg, template, pool):
    scfg = cfg["sbc"]
    rng = np.random.default_rng(cfg["seed"] + 1000 + idx)
    theta_true = L.prior_sample(rng, template)
    cd = simulate_dataset(theta_true, template, rng)

    ndim = cd.n_params
    nwalk = scfg["inner_walkers"]
    p0 = np.array([L.prior_sample(rng, cd) for _ in range(nwalk)])
    target = L.LogPosterior(cd)                  # picklable, uses the task's cores
    sampler = emcee.EnsembleSampler(nwalk, ndim, target, pool=pool)
    sampler.run_mcmc(p0, scfg["inner_steps"], progress=False)
    chain = sampler.get_chain(discard=scfg["inner_burn"], flat=True)

    # SBC rank: fraction of posterior draws below the true value, per parameter.
    ranks = np.mean(chain < theta_true[None, :], axis=0)
    return ranks.tolist(), theta_true.tolist()


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("config")
    ap.add_argument("--task", type=int, default=None,
                    help="SLURM array index: run only dataset #task")
    args = ap.parse_args()
    cfg = load_config(args.config)
    name = cfg["run_name"]

    u, length = operating_uL()
    df = data.load_profiles()
    template = L.build_calibration_data(df, cfg["compounds"], length, u)

    n = cfg["sbc"]["n_datasets"]
    indices = [args.task] if args.task is not None else list(range(n))

    out_ranks, out_true = [], []
    with get_pool(cfg["sbc"]["backend"]) as pool:
        if not pool.is_master():
            return
        for i in indices:
            ranks, true = run_one(i, cfg, template, pool)
            out_ranks.append(ranks)
            out_true.append(true)
            print(f"dataset {i}: ranks={np.round(ranks, 2)}")

    suffix = f"_task{args.task}" if args.task is not None else ""
    payload = {"param_names": template.param_names, "ranks": out_ranks, "theta_true": out_true}
    (RESULTS / f"sbc_{name}{suffix}.json").write_text(json.dumps(payload, indent=2))
    print(f"-> results/sbc_{name}{suffix}.json")


if __name__ == "__main__":
    main()
