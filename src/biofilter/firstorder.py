"""v2 reanalysis -- identifiable pseudo-first-order model.

The full two-scale Monod model (``likelihood.py``) is *statistically calibrated*
(it passes SBC) yet returns design credible intervals spanning three orders of
magnitude.  The Deep-Research diagnosis (see ``deepresearchresults.txt``) traces
this to **structural non-identifiability**: the bench inlet concentrations sit far
below the half-saturation constant, so the biofilter operates entirely in the
pseudo-first-order regime where only the lumped ratio ``Rmax/Ks`` -- an effective
first-order rate -- is identifiable, never ``Rmax`` and ``Ks`` separately.  The
free Henry constant and the ``Lf``/``a`` product add further confounding, and the
"over-performing bed" (near-complete removal in the first third of the column)
leaves the slow-kinetics tail of the posterior unbounded.

This module infers **what the data actually constrains**: a single macroscopic
first-order rate constant ``k1`` (1/s) per compound, plus one observation-noise
scale.  With axial dispersion dropped (bed Peclet number >> 100, so plug flow)
the column forward model is analytic,

    Cg(t) = Cin * exp(-k1 * t),

so calibration is milliseconds/eval instead of a nested BVP solve -- and the
posterior on ``k1`` is tight and physically meaningful.

The mechanistic transport/thermodynamic constants are *fixed* to literature
values rather than inferred; they are retained here only to interpret an
identified ``k1`` (e.g. sanity-check the implied Thiele regime), never as free
knobs.  Separating ``Rmax`` from ``Ks`` would require new experiments that push
the column out of the starvation regime (higher inlet load / shorter EBCT /
transient step tests) -- a Discussion-section recommendation, not more compute.
"""
from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np

# --- Fixed literature constants (SI; dimensionless Henry = gas/liquid) --------
# Henry's law constants are tabulated physical chemistry, not tunable: fixing
# them removes a whole dimension of phase-partitioning confounding.  Values are
# dimensionless gas-over-liquid ratios at ~25 C from the Deep-Research synthesis.
HENRY = {
    "Hexanal": 0.0087,
    "3-MB": 0.0005,     # 3-methyl-1-butanol (isoamyl alcohol)
    "2-MB": 0.0005,     # 2-methyl-1-butanol
}
DF_FILM = 6.0e-10       # effective diffusivity in wet biofilm, m^2/s
# (Lf ~ 5e-5 m, a ~ 800 1/m are the tightened-prior central values; used only
#  for post-hoc interpretation of k1, not inferred here.)


# ---- parameter layout -------------------------------------------------------
# theta = [log_k1_c for each compound, log_sigma_c for each compound]
# Per-compound sigma: the alcohols (16-40 ppmv) and hexanal (2-7 ppmv) live on
# very different concentration scales, so one shared absolute noise would let the
# alcohols swamp hexanal's k1 -- the parameter the design actually depends on.
@dataclass
class FirstOrderData:
    """Pooled bench observations for the first-order reanalysis.

    ``per_compound[c]`` is a list of ``(Cin, t_s, y)`` tuples, one per run, with
    residence time in *seconds* and the inlet port excluded (it is the boundary
    condition).  Runs are pooled: a single ``k1`` per compound, with run-to-run
    scatter absorbed by the per-compound observation noise ``sigma``.
    """
    compounds: list[str]
    per_compound: dict = field(default_factory=dict)

    @property
    def n_params(self) -> int:
        return 2 * len(self.compounds)

    @property
    def param_names(self) -> list[str]:
        return ([f"log_k1[{c}]" for c in self.compounds]
                + [f"log_sigma[{c}]" for c in self.compounds])


# Weakly-informative priors as (mu, sd) on the log parameter.  log_k1 is centred
# on 0.05 1/s and spans ~0.007-0.37 1/s at +-2sd, covering the literature range
# for C4-C6 alcohols/aldehydes (~0.01-0.15 1/s).
PRIOR = {
    "log_k1": (np.log(0.05), 1.0),     # per-compound first-order rate, 1/s
    "log_sigma": (np.log(1.0), 1.0),   # per-compound obs noise, ppmv
}


def forward(k1: float, Cin: float, t_s: np.ndarray) -> np.ndarray:
    """Plug-flow first-order column: Cg(t) = Cin * exp(-k1 * t)."""
    return Cin * np.exp(-k1 * np.asarray(t_s))


def log_prior(theta, n_comp):
    mu_k, sd_k = PRIOR["log_k1"]
    mu_s, sd_s = PRIOR["log_sigma"]
    lp = 0.0
    for v in theta[:n_comp]:
        lp += -0.5 * ((v - mu_k) / sd_k) ** 2
    for v in theta[n_comp:2 * n_comp]:
        lp += -0.5 * ((v - mu_s) / sd_s) ** 2
    return lp


def log_likelihood(theta, data: FirstOrderData):
    n_comp = len(data.compounds)
    ll = 0.0
    for ci, comp in enumerate(data.compounds):
        k1 = np.exp(theta[ci])
        sigma = np.exp(theta[n_comp + ci])
        if not np.isfinite(sigma) or sigma <= 0:
            return -np.inf
        for Cin, t_s, y in data.per_compound[comp]:
            pred = forward(k1, Cin, t_s)
            resid = (np.asarray(y) - pred) / sigma
            ll += -0.5 * np.sum(resid**2) - len(y) * np.log(sigma)
    return ll


def log_posterior(theta, data: FirstOrderData):
    lp = log_prior(theta, len(data.compounds))
    if not np.isfinite(lp):
        return -np.inf
    ll = log_likelihood(theta, data)
    if not np.isfinite(ll):
        return -np.inf
    return lp + ll


class LogPosterior:
    """Picklable log-posterior target (carries its own data)."""

    def __init__(self, data: FirstOrderData):
        self.data = data

    def __call__(self, theta):
        return log_posterior(theta, self.data)


def prior_sample(rng, data: FirstOrderData):
    """Draw one parameter vector from the prior (for init / SBC)."""
    n_comp = len(data.compounds)
    theta = np.empty(data.n_params)
    mu_k, sd_k = PRIOR["log_k1"]
    mu_s, sd_s = PRIOR["log_sigma"]
    theta[:n_comp] = rng.normal(mu_k, sd_k, size=n_comp)
    theta[n_comp:2 * n_comp] = rng.normal(mu_s, sd_s, size=n_comp)
    return theta


def build_fo_data(df, compounds) -> FirstOrderData:
    """Assemble FirstOrderData from a tidy profiles DataFrame (runs pooled)."""
    from .data import compound_arrays  # local import to avoid cycle

    per = {}
    for comp in compounds:
        _, _, _, n_runs, _, runs = compound_arrays(df, comp)
        runs_list = []
        for r in range(n_runs):
            sub = df[(df["compound"] == comp) & (df["run_date"] == runs[r])]
            sub = sub.sort_values("position_m")
            positions = sub["position_m"].to_numpy()
            t_min = sub["residence_time_min"].to_numpy()
            yy = sub["conc_ppmv"].to_numpy()
            Cin = float(sub[sub["position_m"] == 0.0]["conc_ppmv"].iloc[0])
            keep = positions > 0.0     # inlet is the boundary condition
            runs_list.append((Cin, t_min[keep] * 60.0, yy[keep]))
        per[comp] = runs_list
    return FirstOrderData(compounds=list(compounds), per_compound=per)


def simulate(theta, data: FirstOrderData, rng) -> FirstOrderData:
    """Generate a synthetic dataset at ``theta`` (for SBC), reusing the design
    matrix (Cin, t_s) of ``data`` and adding Gaussian noise.

    The noise is *not* clipped at zero: SBC is only valid when the data-generating
    process is exactly the model the likelihood assumes (plain Gaussian, see
    ``log_likelihood``).  Clipping would fold negative draws upward and bias the
    inferred noise scale -- most where predictions sit near zero (hexanal, and the
    per-compound sigmas), which is precisely what a clipped simulate() flagged as
    spurious SBC failures.  (Detection-limit censoring of the *real* zeros is a
    separate real-data modelling refinement; it does not affect the k1 posteriors
    that drive the design.)
    """
    n_comp = len(data.compounds)
    per = {}
    for ci, comp in enumerate(data.compounds):
        k1 = np.exp(theta[ci])
        sigma = np.exp(theta[n_comp + ci])
        runs_list = []
        for Cin, t_s, y in data.per_compound[comp]:
            pred = forward(k1, Cin, t_s)
            noisy = pred + rng.normal(0.0, sigma, size=pred.shape)
            runs_list.append((Cin, t_s, noisy))
        per[comp] = runs_list
    return FirstOrderData(compounds=list(data.compounds), per_compound=per)


# --- Design under uncertainty (analytic quantile-based sizing) ---------------
Q_FULL_M3_MIN = 2831.68
Q_FULL_M3_S = Q_FULL_M3_MIN / 60.0
HEXANAL_CIN_PPMV = 2.0          # not used in sizing (first-order is linear) but
TARGET_REMOVAL = 0.95           # kept for provenance with the original duty spec
LENGTH_TO_DIAM = 3.0


def required_ebrt(k1: float, target: float = TARGET_REMOVAL) -> float:
    """Empty-bed residence time (s) for ``target`` fractional removal.

    First order: 1 - exp(-k1 * EBRT) = target  =>  EBRT = -ln(1-target)/k1.
    """
    return -np.log(1.0 - target) / k1


def volume_from_ebrt(ebrt_s: float) -> float:
    """Industrial bed volume (m^3) = Q * EBRT."""
    return Q_FULL_M3_S * ebrt_s


def geometry_from_volume(V: float):
    """Under L = 3D: V = pi*L^3/36  ->  (L, D)."""
    L = (36.0 * V / np.pi) ** (1.0 / 3.0)
    return L, L / LENGTH_TO_DIAM


def required_volume(k1: float, target: float = TARGET_REMOVAL) -> float:
    return volume_from_ebrt(required_ebrt(k1, target))
