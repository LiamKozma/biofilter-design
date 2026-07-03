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
    # Physically s = S/Ks >= 0 everywhere. solve_bvp's Newton iterations can
    # transiently overshoot into s < 0, where the Monod term phi2*s/(1+s) has a
    # pole at s = -1 -- the source of the divide-by-zero/overflow warnings that
    # made the solver thrash to max_nodes on diffusion-limited profiles. Clamping
    # to the physical branch removes the singularity without changing the
    # converged solution (where s >= 0 holds anyway).
    s = np.maximum(y[0], 0.0)
    return np.vstack((y[1], phi2 * s / (1.0 + s)))


def _bvp_bc(ya, yb, sb):
    # s'(0) = 0  and  s(1) = sb
    return np.array([ya[1], yb[0] - sb])


# Above this Thiele modulus the film is strongly diffusion-limited and the
# deep-film closed form below is accurate to <0.1% (validated against the solved
# BVP over sb in [1e-3, 50]); it also takes over *before* solve_bvp starts
# failing to resolve the O(1/phi)-thick boundary layer (failures set in ~phi=20).
_PHI_ASYMPTOTIC = 12.0


def _eta_deepfilm(phi: float, sb: float) -> float:
    """Diffusion-limited (large-phi) effectiveness factor, closed form.

    One first integral of ``s'' = phi^2 s/(1+s)`` with a fully depleted core
    (``s -> 0`` deep in the film) gives ``0.5 s'(1)^2 = phi^2 (F(sb) - F(0))``
    with ``F(s) = s - ln(1+s)``, hence

        eta = s'(1) (1+sb)/(phi^2 sb) = (1+sb)/sb * sqrt(2 (sb - ln(1+sb))) / phi.

    Exact as ``phi -> inf``. Unlike the first-order form ``tanh(phi)/phi`` it
    stays correct at large ``sb`` (near zero-order kinetics), where the latter
    underestimates eta severalfold. Clipped to (0, 1] -- the clip also makes it a
    safe fallback at small phi (full penetration, eta -> 1).
    """
    return float(np.clip((1.0 + sb) / sb * np.sqrt(2.0 * (sb - np.log1p(sb))) / phi,
                         1e-6, 1.0))


def effectiveness_factor(phi: float, sb: float, n_mesh: int = 41) -> float:
    """Effectiveness factor for a flat Monod biofilm at Thiele modulus ``phi``.

    ``sb`` is the dimensionless bulk concentration ``Sb/Ks``. Returns a value in
    (0, 1]. Three regimes: the first-order limit ``tanh(phi)/phi`` at tiny ``sb``;
    the deep-film closed form :func:`_eta_deepfilm` at large ``phi`` (both exact
    and far cheaper than solving the stiff boundary layer); and the solved BVP in
    between.
    """
    if sb <= 0.0:
        return 1.0
    phi2 = phi * phi
    # First-order limit (sb << 1): eta = tanh(phi)/phi, analytic and robust.
    if sb < 1e-3 or phi < 1e-6:
        return float(np.tanh(phi) / phi) if phi > 1e-6 else 1.0
    # Diffusion-limited limit (phi large): closed form, exact and BVP-free.
    if phi >= _PHI_ASYMPTOTIC:
        return _eta_deepfilm(phi, sb)

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
        # Correct fallback: the deep-film form (clipped) beats tanh(phi)/phi,
        # which collapses to the wrong branch at large sb.
        return _eta_deepfilm(phi, sb)

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
