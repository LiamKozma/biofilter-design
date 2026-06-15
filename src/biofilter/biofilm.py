"""Biofilm scale: reaction-diffusion BVP and the effectiveness factor.

Inside a flat biofilm of thickness ``Lf`` the dissolved-substrate concentration
``S(z)`` obeys a balance between molecular diffusion and Monod uptake:

    Df S''(z) = Rmax * S / (Ks + S),     Rmax = mu_max * Xf / Y

with a no-flux condition at the support (``S'(0) = 0``) and the bulk interface
concentration imposed at the outer face (``S(Lf) = Sb``).

Non-dimensionalising with ``s = S/Ks`` and ``zeta = z/Lf`` collapses every
parameter into a single group -- the Thiele modulus ``phi`` -- so the entire
family of profiles is a one-parameter problem in ``(phi, sb)``:

    s''(zeta) = phi^2 * s / (1 + s),     s'(0) = 0,  s(1) = sb = Sb/Ks
    phi = Lf * sqrt(Rmax / (Df * Ks))

The effectiveness factor is the fraction of the film's intrinsic capacity that
is actually realised once diffusion limitation is accounted for:

    eta = (diffusive flux into film) / (Lf * surface reaction rate)
        = s'(1) * (1 + sb) / (phi^2 * sb)     in (0, 1].

``phi << 1`` -> reaction-limited, full penetration, ``eta -> 1``.
``phi >> 1`` -> diffusion-limited, only an outer shell works, ``eta -> 1/phi``
(first-order limit) and the interior is dead weight.

The expensive operation is solving the nonlinear BVP. Calibration and
sensitivity analysis evaluate ``eta`` thousands of times per forward model, so
:func:`effectiveness_table` solves the BVP on a grid of bulk concentrations once
and returns a fast monotone interpolant.
"""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from scipy.integrate import solve_bvp
from scipy.interpolate import PchipInterpolator


def thiele_modulus(Lf: float, Rmax: float, Df: float, Ks: float) -> float:
    """Thiele modulus phi = Lf * sqrt(Rmax / (Df * Ks))."""
    return Lf * np.sqrt(Rmax / (Df * Ks))


def _bvp_rhs(zeta, y, phi2):
    # y[0] = s, y[1] = s'
    s = y[0]
    return np.vstack((y[1], phi2 * s / (1.0 + s)))


def _bvp_bc(ya, yb, sb):
    # s'(0) = 0  and  s(1) = sb
    return np.array([ya[1], yb[0] - sb])


def effectiveness_factor(phi: float, sb: float, n_mesh: int = 41) -> float:
    """Effectiveness factor for a flat Monod biofilm at Thiele modulus ``phi``.

    ``sb`` is the dimensionless bulk concentration ``Sb/Ks``. Returns a value in
    (0, 1]. Falls back to the closed-form first-order result when ``sb`` is tiny
    (the Monod term linearises and the BVP becomes stiff/degenerate).
    """
    if sb <= 0.0:
        return 1.0
    phi2 = phi * phi
    # First-order limit (sb << 1): eta = tanh(phi)/phi, analytic and robust.
    if sb < 1e-3 or phi < 1e-6:
        return float(np.tanh(phi) / phi) if phi > 1e-6 else 1.0

    zeta = np.linspace(0.0, 1.0, n_mesh)
    # Initial guess: linear from a reduced core value up to sb.
    s0 = np.maximum(sb * (0.3 + 0.7 * zeta), 1e-6)
    sp0 = np.gradient(s0, zeta)
    y0 = np.vstack((s0, sp0))

    sol = solve_bvp(
        lambda z, y: _bvp_rhs(z, y, phi2),
        lambda ya, yb: _bvp_bc(ya, yb, sb),
        zeta,
        y0,
        max_nodes=4000,
        tol=1e-6,
    )
    if not sol.success:
        # Conservative fallback: first-order estimate.
        return float(np.tanh(phi) / phi)

    s_prime_1 = sol.y[1, -1]
    eta = s_prime_1 * (1.0 + sb) / (phi2 * sb)
    return float(np.clip(eta, 1e-6, 1.0))


@dataclass
class EffectivenessInterpolant:
    """Fast eta(Sb_gas) callable for a fixed parameter set.

    Built once per forward model. ``Sb`` here is the *interface liquid*
    concentration the biofilm sees; the column module supplies it after Henry
    partitioning of the gas-phase concentration.
    """
    phi: float
    Ks: float
    _interp: PchipInterpolator
    _sb_max: float

    def __call__(self, Sb):
        Sb = np.asarray(Sb, dtype=float)
        sb = np.clip(Sb, 0.0, self._sb_max * self.Ks) / self.Ks
        return self._interp(sb)


def effectiveness_table(
    phi: float,
    Ks: float,
    sb_max: float = 50.0,
    n_points: int = 24,
) -> EffectivenessInterpolant:
    """Tabulate eta over dimensionless bulk concentration ``sb in [0, sb_max]``.

    Uses log-spaced nodes (resolution is needed at low ``sb`` where eta varies
    fastest) plus the endpoints, then a shape-preserving PCHIP interpolant so the
    column solver gets a smooth, monotone, cheap eta(Sb).
    """
    sb_nodes = np.concatenate(([0.0], np.logspace(-3, np.log10(sb_max), n_points)))
    eta_nodes = np.array([effectiveness_factor(phi, sb) for sb in sb_nodes])
    interp = PchipInterpolator(sb_nodes, eta_nodes, extrapolate=True)
    return EffectivenessInterpolant(phi=phi, Ks=Ks, _interp=interp, _sb_max=sb_max)
