"""Stage C -- Bayesian calibration of the two-scale mechanistic model.

The likelihood embeds the full coupled PDE forward model: every evaluation of the
log-posterior solves the biofilm BVP family and the column ADR BVP for each
compound and run, then compares predicted port concentrations to the bench-scale
measurements. This PDE-in-the-loop structure is what makes the inference
genuinely expensive -- and what a single lumped rate constant can never deliver:
posterior distributions over *physical* transport and kinetic parameters, with
their correlations and non-identifiabilities laid bare.

Model structure (hierarchical):
    shared transport     : Df, Lf, a, Dax, H     (one biofilm/bed geometry)
    per-compound kinetics : Rmax_c, Ks_c          (distinct chemistry per VOC)
    observation noise     : sigma                 (ppmv, Gaussian)

All strictly-positive parameters are sampled in log space. Priors are weakly
informative and centred on literature-scale values for gas-phase biofilters.
Inlet concentrations are taken per run from the measured inlet port.
"""
from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np

from .column import ColumnParams
from .simulator import forward_profile


# ---- parameter layout -------------------------------------------------------
# theta = [logDf, logLf, log_a, logDax, logH, logSigma,
#          logRmax_1, logKs_1, logRmax_2, logKs_2, ...]
SHARED_NAMES = ["log_Df", "log_Lf", "log_a", "log_Dax", "log_H", "log_sigma"]


@dataclass
class CalibrationData:
    """Bundled observations for one or more compounds across runs."""
    compounds: list[str]
    length: float
    u: float
    # per compound: list of (Cin, positions, y) over runs
    per_compound: dict = field(default_factory=dict)

    @property
    def n_params(self) -> int:
        return len(SHARED_NAMES) + 2 * len(self.compounds)

    @property
    def param_names(self) -> list[str]:
        names = list(SHARED_NAMES)
        for c in self.compounds:
            names += [f"log_Rmax[{c}]", f"log_Ks[{c}]"]
        return names


# Weakly-informative priors as (mu, sd) on the log parameter. Centres reflect
# typical gas-phase biofilter scales (SI units).
PRIOR = {
    "log_Df": (np.log(1e-9), 1.0),     # effective diffusivity, m^2/s
    "log_Lf": (np.log(1e-4), 1.0),     # film thickness, m (~100 um)
    "log_a": (np.log(2e3), 1.5),       # specific area, 1/m
    "log_Dax": (np.log(1e-4), 1.5),    # axial dispersion, m^2/s
    "log_H": (np.log(20.0), 1.0),      # dimensionless Henry (gas/liquid)
    "log_sigma": (np.log(1.0), 1.0),   # obs noise, ppmv
    "log_Rmax": (np.log(2.0), 2.0),    # per-compound max rate
    "log_Ks": (np.log(5.0), 2.0),      # per-compound half-saturation
}


def _split(theta, n_comp):
    shared = theta[: len(SHARED_NAMES)]
    kin = theta[len(SHARED_NAMES):].reshape(n_comp, 2)
    return shared, kin


def log_prior(theta, n_comp):
    shared, kin = _split(theta, n_comp)
    lp = 0.0
    for val, name in zip(shared, SHARED_NAMES):
        mu, sd = PRIOR[name]
        lp += -0.5 * ((val - mu) / sd) ** 2
    for rmax, ks in kin:
        mu, sd = PRIOR["log_Rmax"]
        lp += -0.5 * ((rmax - mu) / sd) ** 2
        mu, sd = PRIOR["log_Ks"]
        lp += -0.5 * ((ks - mu) / sd) ** 2
    return lp


def _params_for(shared, rmax_log, ks_log) -> ColumnParams:
    log_Df, log_Lf, log_a, log_Dax, log_H, _ = shared
    return ColumnParams(
        Rmax=np.exp(rmax_log),
        Ks=np.exp(ks_log),
        Df=np.exp(log_Df),
        Lf=np.exp(log_Lf),
        a=np.exp(log_a),
        Dax=np.exp(log_Dax),
        H=np.exp(log_H),
    )


def log_likelihood(theta, data: CalibrationData):
    shared, kin = _split(theta, len(data.compounds))
    sigma = np.exp(shared[5])
    if not np.isfinite(sigma) or sigma <= 0:
        return -np.inf
    ll = 0.0
    n = 0
    for ci, comp in enumerate(data.compounds):
        params = _params_for(shared, kin[ci, 0], kin[ci, 1])
        for Cin, positions, y in data.per_compound[comp]:
            pred, ok = forward_profile(params, Cin, positions, data.length, data.u)
            if not ok or np.any(~np.isfinite(pred)):
                return -np.inf
            resid = (np.asarray(y) - pred) / sigma
            ll += -0.5 * np.sum(resid**2) - len(y) * np.log(sigma)
            n += len(y)
    return ll


def log_posterior(theta, data: CalibrationData):
    lp = log_prior(theta, len(data.compounds))
    if not np.isfinite(lp):
        return -np.inf
    ll = log_likelihood(theta, data)
    if not np.isfinite(ll):
        return -np.inf
    return lp + ll


class LogPosterior:
    """Picklable log-posterior target.

    Defined in the library (not a driver ``__main__``) so it survives the pickle
    round-trip that ``multiprocessing`` (spawn/forkserver) and ``schwimmbad``'s
    ``MPIPool`` perform when shipping the target to worker ranks. Carries its own
    data, so no module-global state is relied upon.
    """

    def __init__(self, data: CalibrationData):
        self.data = data

    def __call__(self, theta):
        return log_posterior(theta, self.data)


def prior_sample(rng, data: CalibrationData):
    """Draw one parameter vector from the prior (for init / SBC)."""
    theta = np.empty(data.n_params)
    for i, name in enumerate(SHARED_NAMES):
        mu, sd = PRIOR[name]
        theta[i] = rng.normal(mu, sd)
    j = len(SHARED_NAMES)
    for _ in data.compounds:
        mu, sd = PRIOR["log_Rmax"]
        theta[j] = rng.normal(mu, sd); j += 1
        mu, sd = PRIOR["log_Ks"]
        theta[j] = rng.normal(mu, sd); j += 1
    return theta


def build_calibration_data(df, compounds, length, u):
    """Assemble a CalibrationData from a tidy profiles DataFrame."""
    from .data import compound_arrays  # local import to avoid cycle

    per = {}
    for comp in compounds:
        t, y, ridx, n_runs, C0_obs, runs = compound_arrays(df, comp)
        runs_list = []
        for r in range(n_runs):
            mask = ridx == r
            pos = df[df["compound"] == comp].sort_values(["run_date", "position_m"])
            # rebuild positions for this run from the same ordering
            sub = df[(df["compound"] == comp) & (df["run_date"] == runs[r])]
            sub = sub.sort_values("position_m")
            positions = sub["position_m"].to_numpy()
            yy = sub["conc_ppmv"].to_numpy()
            Cin = float(sub[sub["position_m"] == 0.0]["conc_ppmv"].iloc[0])
            # predict/compare at non-inlet ports (inlet is the boundary condition)
            keep = positions > 0.0
            runs_list.append((Cin, positions[keep], yy[keep]))
        per[comp] = runs_list
    return CalibrationData(compounds=list(compounds), length=length, u=u, per_compound=per)
