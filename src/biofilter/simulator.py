"""Coupled two-scale forward model: parameters -> predicted port concentrations.

One forward evaluation builds the effectiveness-factor interpolant for the given
biofilm parameters (one nonlinear BVP family), then solves the column ADR BVP and
samples the predicted gas concentration at the measurement ports. This is the
unit of work that Bayesian calibration (Stage C), global sensitivity (Stage D),
and design-under-uncertainty (Stage E) call hundreds of thousands of times.
"""
from __future__ import annotations

import numpy as np

from .biofilm import effectiveness_table
from .column import ColumnParams, solve_column


def forward_profile(
    params: ColumnParams,
    Cin: float,
    positions,
    length: float,
    u: float,
    n_mesh: int = 41,
    eta_points: int = 16,
):
    """Predicted gas concentration (ppmv) at each position in ``positions``.

    Returns ``(pred, ok)`` where ``ok`` is the solver success flag. The eta table
    is built once and reused for the single column solve.
    """
    positions = np.asarray(positions, float)
    sb_max = max(5.0, (Cin / params.H) / params.Ks * 1.5)
    eta = effectiveness_table(params.phi, params.Ks, sb_max=sb_max, n_points=eta_points)
    x, Cg, sol = solve_column(params, Cin, length, u, n_mesh=n_mesh, eta=eta)
    if not sol.success:
        return np.full_like(positions, np.nan), False
    pred = np.interp(positions, x, Cg)
    return pred, True
