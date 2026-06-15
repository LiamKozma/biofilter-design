"""Stage E -- scale-up and design under uncertainty.

The original report scaled the bench column to an industrial duty
(Q = 100,000 ft^3/min = 2831.68 m^3/min) by integrating a single first-order
rate constant in a plug-flow model, obtaining one number: V approx 3500 m^3,
giving D = 11.41 m, L = 34.23 m under BIOREM's 3:1 length-to-diameter rule.

Here we replace that point estimate with the *posterior predictive design*. For
each draw of the calibrated mechanistic parameters we solve the full two-scale
model at industrial scale and find the bed length that achieves the target
removal, honouring the same geometric constraint. The result is a credible
interval on the required volume -- an honest statement of how much the sparse
bench data actually constrains the design, which a deterministic calculation
cannot provide.
"""
from __future__ import annotations

import numpy as np
from scipy.optimize import brentq

from .column import ColumnParams
from .simulator import forward_profile

# Industrial duty (report Section 3).
Q_FULL_M3_MIN = 2831.68
Q_FULL_M3_S = Q_FULL_M3_MIN / 60.0
HEXANAL_CIN_PPMV = 2.0          # influent hexanal (Table 2 sample mean)
TARGET_REMOVAL = 0.95
LENGTH_TO_DIAM = 3.0            # l <= 3 D structural rule (BIOREM)


def geometry_from_length(length: float):
    """Under l = 3 D: diameter, cross-sectional area, superficial velocity."""
    D = length / LENGTH_TO_DIAM
    A = np.pi * D**2 / 4.0
    u = Q_FULL_M3_S / A
    V = A * length
    return D, A, u, V


def removal_at_length(params: ColumnParams, length: float, Cin: float = HEXANAL_CIN_PPMV) -> float:
    D, A, u, V = geometry_from_length(length)
    pred, ok = forward_profile(params, Cin, [length], length, u)
    if not ok:
        return np.nan
    return 1.0 - pred[0] / Cin


def required_length(
    params: ColumnParams,
    target: float = TARGET_REMOVAL,
    Cin: float = HEXANAL_CIN_PPMV,
    bracket=(1.0, 200.0),
):
    """Bed length achieving ``target`` removal; returns (L, D, V) or NaNs.

    Note: increasing L also widens D (3:1 rule) and so *raises* superficial
    velocity, shortening residence time -- the trade-off the deterministic
    plug-flow sizing ignores. We root-find on the net removal.
    """
    def f(L):
        e = removal_at_length(params, L, Cin) - target
        # brentq cannot tolerate NaN; map solver failures to a large negative
        # residual (treated as "removal far below target" so the root search
        # walks toward longer beds rather than aborting).
        return e if np.isfinite(e) else -target

    lo, hi = bracket
    flo, fhi = f(lo), f(hi)
    if flo * fhi > 0:
        return np.nan, np.nan, np.nan
    try:
        L = brentq(f, lo, hi, xtol=1e-2, rtol=1e-4, maxiter=100)
    except (ValueError, RuntimeError):
        return np.nan, np.nan, np.nan
    D, A, u, V = geometry_from_length(L)
    return L, D, V


class RemovalAtLength:
    """Picklable Sobol target: log10 parameter vector -> removal at a fixed bed.

    Input row order: [log10 Rmax, Ks, Df, Lf, a, Dax, H]. Returns the industrial-
    scale hexanal removal at ``ref_length`` (geometry held fixed so the index
    isolates parameter influence). Lives in the library so it pickles across MPI
    workers.
    """

    def __init__(self, ref_length: float, Cin: float = HEXANAL_CIN_PPMV):
        self.ref_length = ref_length
        self.Cin = Cin

    def __call__(self, log10_row):
        rmax, ks, df, lf, a, dax, h = (10.0 ** np.asarray(log10_row))
        p = ColumnParams(Rmax=rmax, Ks=ks, Df=df, Lf=lf, a=a, Dax=dax, H=h)
        e = removal_at_length(p, self.ref_length, self.Cin)
        return e if np.isfinite(e) else 0.0


class DesignDrawer:
    """Picklable design-UQ target: posterior theta -> (L, D, V) at 95% removal."""

    def __init__(self, compounds, target_compound="Hexanal", target=TARGET_REMOVAL):
        self.compounds = list(compounds)
        self.target_compound = target_compound
        self.target = target

    def __call__(self, theta):
        p = params_from_theta(theta, self.compounds, self.target_compound)
        return required_length(p, target=self.target)


def params_from_theta(theta, compounds, target_compound="Hexanal"):
    """Build ColumnParams for one compound from a calibration theta vector."""
    from .likelihood import SHARED_NAMES

    shared = theta[: len(SHARED_NAMES)]
    ci = compounds.index(target_compound)
    rmax_log = theta[len(SHARED_NAMES) + 2 * ci]
    ks_log = theta[len(SHARED_NAMES) + 2 * ci + 1]
    log_Df, log_Lf, log_a, log_Dax, log_H, _ = shared
    return ColumnParams(
        Rmax=np.exp(rmax_log), Ks=np.exp(ks_log), Df=np.exp(log_Df),
        Lf=np.exp(log_Lf), a=np.exp(log_a), Dax=np.exp(log_Dax), H=np.exp(log_H),
    )
