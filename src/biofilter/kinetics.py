"""Stage A -- Bayesian integral-method kinetics and model selection.

The original project assigned each aldehyde a reaction order by linearising the
integrated rate law, fitting OLS, and picking the order with the highest R^2.
That procedure is fragile here: only 3-5 points per run, the linearising
transforms (``ln(C0/C)``, ``1/C``) are undefined once a concentration hits zero,
and R^2 is not a valid criterion for comparing models with different response
variables.

This module replaces it with a hierarchical Bayesian treatment. For a compound
measured over ``R`` runs, the rate constant ``k`` is shared (one chemistry) while
each run gets its own inlet concentration ``C0_r``. The three candidate integrated
rate laws are

    zeroth : C(t) = max(C0 - k t, 0)
    first  : C(t) = C0 exp(-k t)
    second : C(t) = C0 / (1 + k C0 t)

with ``t`` the empty-bed residence time to each port. Models are compared by
WAIC / PSIS-LOO (via arviz) -- predictive criteria that are valid across models
because the response variable (the observed concentration) is the same for all
three. The posterior on ``k`` then carries the uncertainty the R^2 method threw
away, which is what lets us state *with a credible interval* which aldehyde is
rate-limiting.
"""
from __future__ import annotations

import numpy as np

ORDERS = ("zeroth", "first", "second")


def predict(order: str, t, C0, k):
    """Integrated rate law: concentration at residence time(s) ``t``."""
    t = np.asarray(t, dtype=float)
    if order == "zeroth":
        return np.maximum(C0 - k * t, 0.0)
    if order == "first":
        return C0 * np.exp(-k * t)
    if order == "second":
        return C0 / (1.0 + k * C0 * t)
    raise ValueError(f"unknown order {order!r}")


def _unpack(theta, n_runs):
    """theta = [log_k, log_sigma, C0_1, ..., C0_R]."""
    log_k = theta[0]
    log_sigma = theta[1]
    C0 = theta[2:]
    return np.exp(log_k), np.exp(log_sigma), C0


def make_logprob(order, t, y, run_idx, n_runs, C0_obs):
    """Build a log-posterior closure for one compound under one reaction order.

    ``t``, ``y`` are 1-D arrays of residence time and concentration; ``run_idx``
    maps each observation to its run (0..n_runs-1). ``C0_obs`` are the measured
    inlet concentrations, used to centre the ``C0`` priors.
    """
    t = np.asarray(t, float)
    y = np.asarray(y, float)
    run_idx = np.asarray(run_idx, int)
    C0_obs = np.asarray(C0_obs, float)
    log_k_mu, log_k_sd = np.log(1.0), 3.0          # weakly informative, k ~ O(1) /min

    def log_prob(theta):
        k, sigma, C0 = _unpack(theta, n_runs)
        if not np.isfinite(k) or k <= 0 or sigma <= 0:
            return -np.inf
        if np.any(C0 <= 0):
            return -np.inf
        mu = predict(order, t, C0[run_idx], k)
        # Gaussian likelihood with a small floor; concentrations in ppmv.
        resid = (y - mu) / sigma
        ll = -0.5 * np.sum(resid**2) - len(y) * np.log(sigma)
        # Priors.
        lp = -0.5 * ((np.log(k) - log_k_mu) / log_k_sd) ** 2
        lp += -0.5 * ((np.log(sigma) - np.log(1.0)) / 1.5) ** 2
        # C0 priors centred on the observed inlet, generous SD.
        lp += np.sum(-0.5 * ((C0 - C0_obs) / (0.5 * C0_obs + 1e-9)) ** 2)
        return ll + lp

    return log_prob


def pointwise_loglik(order, t, y, run_idx, samples):
    """Per-observation log-likelihood for every posterior draw.

    Returns an array of shape ``(n_draws, n_obs)`` suitable for
    ``arviz.from_dict(log_likelihood=...)`` -> WAIC / LOO.
    """
    t = np.asarray(t, float)
    y = np.asarray(y, float)
    run_idx = np.asarray(run_idx, int)
    out = np.empty((samples.shape[0], len(y)))
    for i, theta in enumerate(samples):
        k, sigma, C0 = _unpack(theta, C0_count := samples.shape[1] - 2)
        mu = predict(order, t, C0[run_idx], k)
        out[i] = -0.5 * ((y - mu) / sigma) ** 2 - 0.5 * np.log(2 * np.pi) - np.log(sigma)
    return out
