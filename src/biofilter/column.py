"""Column scale: advection-dispersion-reaction BVP for the packed bed.

Zooming out from the biofilm, the gas-phase concentration ``Cg(x)`` along the
bed obeys advection at superficial velocity ``u``, axial dispersion ``Dax``, and
a volumetric sink supplied by the biofilm model:

    Dax Cg'' - u Cg' - R(Cg) = 0,   0 < x < L

    R(Cg) = a * Lf * eta(Ss) * Rmax * Ss / (Ks + Ss),   Ss = Cg / H

``a`` is the specific biofilm surface area (m^2 per m^3 reactor), ``Lf`` the film
thickness, so ``a*Lf`` is the biofilm volume fraction. ``H`` is a dimensionless
Henry coefficient partitioning the gas concentration into the liquid film the
microbes actually see (``Ss``). The effectiveness factor ``eta(Ss)`` makes the
local rate honest about diffusion limitation -- without it the column model
assumes every cell in every film is fully fed and over-predicts removal.

Boundary conditions are the Danckwerts pair for a closed-closed vessel:

    inlet  (x=0):  u Cg - Dax Cg' = u Cin
    outlet (x=L):  Cg' = 0

The coupling ``eta(Ss)`` is the whole point of the two-scale treatment: a single
lumped first-order rate constant cannot tell you *why* a column under-performs or
which knob (surface area vs. residence time) recovers it.
"""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from scipy.integrate import solve_bvp

from .biofilm import EffectivenessInterpolant, effectiveness_table, thiele_modulus


@dataclass
class ColumnParams:
    """Physical parameters of one biofilter configuration.

    All in SI. ``mu_max``, ``Xf``, ``Y`` enter only through ``Rmax``; the
    calibration parameterises ``Rmax`` directly to avoid an unidentifiable
    product.
    """
    Rmax: float      # max volumetric biofilm reaction rate, ppmv-equiv / s
    Ks: float        # Monod half-saturation (liquid-film units), ppmv-equiv
    Df: float        # effective substrate diffusivity in film, m^2/s
    Lf: float        # biofilm thickness, m
    a: float         # specific biofilm surface area, 1/m  (m^2 / m^3)
    Dax: float       # axial dispersion coefficient, m^2/s
    H: float         # dimensionless Henry partition (gas/liquid)

    @property
    def phi(self) -> float:
        return thiele_modulus(self.Lf, self.Rmax, self.Df, self.Ks)


def _build_eta(params: ColumnParams, sb_max: float, n_points: int) -> EffectivenessInterpolant:
    return effectiveness_table(params.phi, params.Ks, sb_max=sb_max, n_points=n_points)


def solve_column(
    params: ColumnParams,
    Cin: float,
    length: float,
    u: float,
    n_mesh: int = 41,
    eta_points: int = 18,
    eta: EffectivenessInterpolant | None = None,
):
    """Solve the steady ADR BVP. Returns ``(x, Cg, sol)``.

    ``u`` is the superficial gas velocity (m/s); ``Cin`` the inlet gas
    concentration (ppmv); ``length`` the packed-bed depth (m). ``eta`` may be a
    pre-built interpolant (reused across calls with identical ``phi, Ks``).
    """
    if eta is None:
        sb_max = max(5.0, (Cin / params.H) / params.Ks * 1.5)
        eta = _build_eta(params, sb_max=sb_max, n_points=eta_points)

    inv_H = 1.0 / params.H
    aLf_Rmax = params.a * params.Lf * params.Rmax

    def sink(Cg):
        Cg = np.clip(Cg, 0.0, None)
        Ss = Cg * inv_H
        return aLf_Rmax * eta(Ss) * Ss / (params.Ks + Ss)

    def rhs(x, y):
        # y[0] = Cg, y[1] = Cg'
        Cg, Cgp = y
        Cgpp = (u * Cgp + sink(Cg)) / params.Dax
        return np.vstack((Cgp, Cgpp))

    def bc(ya, yb):
        # Danckwerts inlet, zero-gradient outlet.
        return np.array([u * ya[0] - params.Dax * ya[1] - u * Cin, yb[1]])

    x = np.linspace(0.0, length, n_mesh)
    # Initial guess: gentle exponential decay.
    Cg0 = Cin * np.exp(-2.0 * x / length)
    Cgp0 = np.gradient(Cg0, x)
    y0 = np.vstack((Cg0, Cgp0))

    sol = solve_bvp(rhs, bc, x, y0, max_nodes=20000, tol=1e-5)
    return sol.x, sol.y[0], sol


def removal_efficiency(params: ColumnParams, Cin: float, length: float, u: float, **kw) -> float:
    """Inlet-to-outlet fractional removal E = 1 - Cg(L)/Cg(0)."""
    x, Cg, sol = solve_column(params, Cin, length, u, **kw)
    if not sol.success:
        return np.nan
    return float(1.0 - Cg[-1] / Cg[0])
