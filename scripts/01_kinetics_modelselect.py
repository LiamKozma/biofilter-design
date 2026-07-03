#!/usr/bin/env python3
"""Stage A: Bayesian reaction-order selection and rate-limiting determination.

For each aldehyde, fit the zeroth/first/second-order integrated rate laws with a
hierarchical model (shared k, run-specific C0) via emcee, score them with WAIC
and PSIS-LOO, and report the posterior rate constant. The rate-limiting compound
is the one with the slowest posterior rate at its inlet concentration, reported
with a credible interval.

Laptop-scale by design -- this stage frames the problem; the HPC load lives in
the mechanistic calibration (Stage C onward).
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import arviz as az
import emcee
import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
from biofilter import data, kinetics  # noqa: E402

ROOT = Path(__file__).resolve().parents[1]
RESULTS = ROOT / "results"
RESULTS.mkdir(exist_ok=True)

ALDEHYDES = ["3-MB", "2-MB", "Hexanal"]
N_WALKERS = 32
N_STEPS = 4000
N_BURN = 1000
RNG = np.random.default_rng(20220509)


def fit_order(order, t, y, ridx, n_runs, C0_obs):
    ndim = 2 + n_runs
    log_prob = kinetics.make_logprob(order, t, y, ridx, n_runs, C0_obs)
    # Initialise walkers near a sensible point.
    p0 = np.column_stack(
        [
            RNG.normal(np.log(1.0), 0.3, N_WALKERS),                 # log k
            RNG.normal(np.log(1.0), 0.3, N_WALKERS),                 # log sigma
            *[RNG.normal(C0_obs[i], 0.1 * C0_obs[i] + 1e-3, N_WALKERS)
              for i in range(n_runs)],                                # C0_r
        ]
    )
    sampler = emcee.EnsembleSampler(N_WALKERS, ndim, log_prob)
    sampler.run_mcmc(p0, N_STEPS, progress=False)
    chain = sampler.get_chain(discard=N_BURN, flat=True)
    tau = sampler.get_autocorr_time(quiet=True)
    ess = (N_STEPS - N_BURN) * N_WALKERS / np.nanmax(tau)

    loglik = kinetics.pointwise_loglik(order, t, y, ridx, chain)
    idata = az.from_dict(
        posterior={"k": np.exp(chain[None, :, 0])},
        log_likelihood={"y": loglik[None, :, :]},
    )
    loo = az.loo(idata)                       # PSIS-LOO; elpd, higher is better
    k_samp = np.exp(chain[:, 0])
    return {
        "order": order,
        "k_median": float(np.median(k_samp)),
        "k_ci95": [float(np.percentile(k_samp, 2.5)), float(np.percentile(k_samp, 97.5))],
        "elpd_loo": float(loo.elpd_loo),
        "p_loo": float(loo.p_loo),
        "ess": float(ess),
        "_k_samples": k_samp,
    }


def main():
    df = data.load_profiles()
    summary = {}
    rate_at_inlet = {}  # for rate-limiting comparison, with uncertainty

    for comp in ALDEHYDES:
        t, y, ridx, n_runs, C0_obs, runs = data.compound_arrays(df, comp)
        print(f"\n=== {comp} | {len(y)} obs across {n_runs} runs ===")
        fits = {}
        for order in kinetics.ORDERS:
            try:
                fits[order] = fit_order(order, t, y, ridx, n_runs, C0_obs)
            except Exception as exc:  # noqa: BLE001
                print(f"  {order}: FAILED ({exc})")
                continue
            f = fits[order]
            print(f"  {order:7s}: k={f['k_median']:.4g} "
                  f"[{f['k_ci95'][0]:.3g},{f['k_ci95'][1]:.3g}]  "
                  f"elpd_LOO={f['elpd_loo']:.2f}  p_loo={f['p_loo']:.1f}  ESS~{f['ess']:.0f}")

        # Higher elpd_loo = better expected predictive accuracy.
        best = max(fits.values(), key=lambda f: f["elpd_loo"])
        print(f"  -> best order by LOO: {best['order']}")

        # Rate at the mean inlet concentration, propagating k uncertainty.
        C0_mean = float(np.mean(C0_obs))
        ks = best["_k_samples"]
        if best["order"] == "zeroth":
            rate = ks
        elif best["order"] == "first":
            rate = ks * C0_mean
        else:
            rate = ks * C0_mean**2
        rate_at_inlet[comp] = rate

        summary[comp] = {
            "best_order": best["order"],
            "C0_mean_ppmv": C0_mean,
            "orders": {o: {k: v for k, v in f.items() if not k.startswith("_")}
                       for o, f in fits.items()},
            "rate_at_inlet_median_ppmv_per_min": float(np.median(rate)),
            "rate_at_inlet_ci95": [float(np.percentile(rate, 2.5)),
                                   float(np.percentile(rate, 97.5))],
        }

    # Rate-limiting = slowest posterior rate. Report P(compound is slowest).
    comps = list(rate_at_inlet)
    stacked = np.vstack([rate_at_inlet[c] for c in comps])  # (n_comp, n_draws) ragged
    n = min(arr.shape[0] for arr in rate_at_inlet.values())
    stacked = np.vstack([rate_at_inlet[c][:n] for c in comps])
    slowest = np.argmin(stacked, axis=0)
    p_slowest = {comps[i]: float(np.mean(slowest == i)) for i in range(len(comps))}
    print("\n=== Rate-limiting probability P(slowest at inlet) ===")
    for c, p in sorted(p_slowest.items(), key=lambda kv: -kv[1]):
        print(f"  {c:8s}: {p:.3f}")

    summary["_rate_limiting"] = {
        "P_slowest": p_slowest,
        "verdict": max(p_slowest, key=p_slowest.get),
    }
    out = RESULTS / "stageA_kinetics.json"
    out.write_text(json.dumps(summary, indent=2))
    print(f"\nwrote {out.relative_to(ROOT)}")


if __name__ == "__main__":
    main()
