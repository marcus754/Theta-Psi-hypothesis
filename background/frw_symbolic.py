# -*- coding: utf-8 -*-
"""
Symbolic helper utilities for the minisuperspace reduction with fields
θ_eff(t) and R(t), where R is the reduced amplitude representative of the
complex energy sector Ψ = R e^{iφ}.

The module keeps the effective homogeneous Lagrangian in reduced amplitude form

    L_eff = 3/2 (θ̇/θ)^2 + 1/2 Ṙ^2 - 1/2 V0 R^2 + 1/2 γ^2 θ̇^2

and exposes convenience functions for ρ, p, уравнения движения и преобразование
в каноническую переменную χ = √3 ln θ_eff, чтобы легче видеть физику без
искусственных полюсов.

Notes
-----
- `theta_t` in this file is the FRW scalar reduction `θ_eff(t)` of the
  fundamental tensor clock sector `Θ_{μν}`.
- `R_t` is the reduced homogeneous amplitude representative `R` of the
  complex energy sector `Ψ = R e^{iφ}`; the optional phase is tracked
  separately by the background `Q`-charge.
- This file is a minisuperspace helper, not the full tensor ontology.
"""
from __future__ import annotations

import math
from dataclasses import dataclass
from typing import Dict, Tuple

import sympy as sp

SQRT3 = math.sqrt(3.0)
sqrt3_sym = sp.sqrt(3)

__all__ = [
    "lagrangian_expr",
    "energy_density_expr",
    "pressure_expr",
    "eom_theta_expr",
    "eom_R_expr",
    "lagrangian",
    "lagrangian_chi",
    "energy_density",
    "pressure",
    "canonical_momenta",
    "chi_from_theta",
    "chi_dot_from_theta",
    "theta_from_chi",
    "theta_dot_from_chi",
    "chi_from_theta_eff",
    "chi_dot_from_theta_eff",
    "theta_eff_from_chi",
    "theta_eff_dot_from_chi",
    "theta_eff_from_x",
    "theta_eff_dot_from_x",
]

# ---------------------------------------------------------------------------
# Symbolic definitions

t = sp.symbols("t", real=True)
theta_t = sp.Function("theta")(t)
R_t = sp.Function("R")(t)
theta_dot_t = sp.diff(theta_t, t)
R_dot_t = sp.diff(R_t, t)
gamma_sym, V0_sym = sp.symbols("gamma V0", real=True)

lagrangian_expr = (
    sp.Rational(3, 2) * (theta_dot_t / theta_t) ** 2
    + sp.Rational(1, 2) * R_dot_t**2
    - sp.Rational(1, 2) * V0_sym * R_t**2
    + sp.Rational(1, 2) * gamma_sym**2 * theta_dot_t**2
)

# Canonical variable definitions
chi_t = sp.Function("chi")(t)
chi_dot_t = sp.diff(chi_t, t)
theta_from_chi_expr = sp.exp(chi_t / sqrt3_sym)
theta_dot_from_chi_expr = sp.simplify(sp.diff(theta_from_chi_expr, t))

lagrangian_chi_expr = sp.simplify(
    sp.Rational(1, 2) * chi_dot_t**2
    + sp.Rational(1, 2) * R_dot_t**2
    - sp.Rational(1, 2) * V0_sym * R_t**2
    + sp.Rational(1, 2) * gamma_sym**2 * theta_dot_from_chi_expr**2
)

# Minisuperspace energy density and pressure (up to factor a^3)
energy_density_expr = sp.simplify(
    (sp.diff(lagrangian_expr, theta_dot_t) * theta_dot_t)
    + (sp.diff(lagrangian_expr, R_dot_t) * R_dot_t)
    - lagrangian_expr
)
pressure_expr = sp.simplify(lagrangian_expr)


def _euler_lagrange(L, q):
    """Convenience wrapper for Euler–Lagrange equation."""
    return sp.simplify(sp.diff(sp.diff(L, sp.diff(q, t)), t) - sp.diff(L, q))


eom_theta_expr = _euler_lagrange(lagrangian_expr, theta_t)
eom_R_expr = _euler_lagrange(lagrangian_expr, R_t)

# ---------------------------------------------------------------------------
# Numerical helpers (lambdified functions)

theta_sym, thetadot_sym, R_sym, Rdot_sym = sp.symbols(
    "theta thetadot R Rdot", real=True
)

lagrangian_simple = lagrangian_expr.subs(
    {theta_t: theta_sym, theta_dot_t: thetadot_sym, R_t: R_sym, R_dot_t: Rdot_sym}
)
energy_density_simple = energy_density_expr.subs(
    {theta_t: theta_sym, theta_dot_t: thetadot_sym, R_t: R_sym, R_dot_t: Rdot_sym}
)
pressure_simple = pressure_expr.subs(
    {theta_t: theta_sym, theta_dot_t: thetadot_sym, R_t: R_sym, R_dot_t: Rdot_sym}
)

chi_sym, chidot_sym = sp.symbols("chi chidot", real=True)
theta_from_chi_simple = sp.exp(chi_sym / sqrt3_sym)
theta_dot_from_chi_simple = (theta_from_chi_simple / sqrt3_sym) * chidot_sym
lagrangian_chi_simple = (
    sp.Rational(1, 2) * chidot_sym**2
    + sp.Rational(1, 2) * Rdot_sym**2
    - sp.Rational(1, 2) * V0_sym * R_sym**2
    + sp.Rational(1, 2) * gamma_sym**2 * theta_dot_from_chi_simple**2
)

_lagrangian_func = sp.lambdify(
    (theta_sym, thetadot_sym, R_sym, Rdot_sym, gamma_sym, V0_sym),
    lagrangian_simple,
    modules="math",
)
_rho_func = sp.lambdify(
    (theta_sym, thetadot_sym, R_sym, Rdot_sym, gamma_sym, V0_sym),
    energy_density_simple,
    modules="math",
)
_p_func = sp.lambdify(
    (theta_sym, thetadot_sym, R_sym, Rdot_sym, gamma_sym, V0_sym),
    pressure_simple,
    modules="math",
)

_lagrangian_chi_func = sp.lambdify(
    (chi_sym, chidot_sym, R_sym, Rdot_sym, gamma_sym, V0_sym),
    lagrangian_chi_simple,
    modules="math",
)

_pi_theta_func = sp.lambdify(
    (theta_sym, thetadot_sym, gamma_sym),
    sp.diff(lagrangian_simple, thetadot_sym),
    modules="math",
)
_pi_R_func = sp.lambdify(
    (R_sym, Rdot_sym),
    sp.diff(lagrangian_simple, Rdot_sym),
    modules="math",
)


# ---------------------------------------------------------------------------
# Public numerical API

def lagrangian(
    theta: float,
    thetadot: float,
    R: float,
    Rdot: float,
    gamma: float,
    V0: float,
) -> float:
    """Evaluate the minisuperspace Lagrangian L_eff."""
    return float(_lagrangian_func(theta, thetadot, R, Rdot, gamma, V0))


def lagrangian_chi(
    chi: float,
    chidot: float,
    R: float,
    Rdot: float,
    gamma: float,
    V0: float,
) -> float:
    """Lagrangian expressed in the canonical variable χ = √3 ln θ."""
    return float(
        _lagrangian_chi_func(chi, chidot, R, Rdot, gamma, V0)
    )


def energy_density(
    theta: float,
    thetadot: float,
    R: float,
    Rdot: float,
    gamma: float,
    V0: float,
) -> float:
    """Energy density ρ (without the overall a^3 factor)."""
    return float(_rho_func(theta, thetadot, R, Rdot, gamma, V0))


def pressure(
    theta: float,
    thetadot: float,
    R: float,
    Rdot: float,
    gamma: float,
    V0: float,
) -> float:
    """Pressure p = L_eff (again up to a^3)."""
    return float(_p_func(theta, thetadot, R, Rdot, gamma, V0))


@dataclass
class CanonicalMomenta:
    pi_theta: float
    pi_R: float


def canonical_momenta(
    theta: float,
    thetadot: float,
    R: float,
    Rdot: float,
    gamma: float,
) -> CanonicalMomenta:
    """
    Return canonical momenta (π_θ, π_Ψ) corresponding to L_eff.

    π_θ = ∂L/∂θ̇,  π_R = ∂L/∂Ṙ.
    """
    pi_theta = float(_pi_theta_func(theta, thetadot, gamma))
    pi_R = float(_pi_R_func(R, Rdot))
    return CanonicalMomenta(pi_theta=pi_theta, pi_R=pi_R)


def chi_from_theta(theta: float) -> float:
    """Canonical logarithmic variable; θ must stay positive."""
    if theta <= 0.0:
        raise ValueError("theta must be positive to define χ = √3 ln θ")
    return SQRT3 * math.log(theta)


def chi_dot_from_theta(theta: float, thetadot: float) -> float:
    if theta <= 0.0:
        raise ValueError("theta must be positive to define χ̇ = √3 θ̇/θ")
    return SQRT3 * (thetadot / theta)


def theta_from_chi(chi: float) -> float:
    """Inverse of χ; returns the positive θ that generated χ."""
    return math.exp(chi / SQRT3)


def theta_dot_from_chi(chi: float, chidot: float) -> float:
    """Derivative θ̇ expressed through χ and its derivative."""
    return theta_from_chi(chi) * chidot / SQRT3


def chi_from_theta_eff(theta_eff: float) -> float:
    """Alias for χ(θ_eff)."""
    return chi_from_theta(theta_eff)


def chi_dot_from_theta_eff(theta_eff: float, theta_eff_dot: float) -> float:
    """Alias for χ̇(θ_eff, θ̇_eff)."""
    return chi_dot_from_theta(theta_eff, theta_eff_dot)


def theta_eff_from_chi(chi: float) -> float:
    """Alias for θ_eff(χ)."""
    return theta_from_chi(chi)


def theta_eff_dot_from_chi(chi: float, chidot: float) -> float:
    """Alias for θ̇_eff(χ, χ̇)."""
    return theta_dot_from_chi(chi, chidot)


# ---------------------------------------------------------------------------
# Utility: equation of motion in the x-variable


def theta_from_x(x: float, gamma: float) -> float:
    return (SQRT3 / gamma) * math.tanh(x)


def theta_dot_from_x(x: float, xdot: float, gamma: float) -> float:
    return (SQRT3 / gamma) * (1.0 / math.cosh(x) ** 2) * xdot


def theta_eff_from_x(x: float, gamma: float) -> float:
    """Alias for θ_eff(x, γ)."""
    return theta_from_x(x, gamma)


def theta_eff_dot_from_x(x: float, xdot: float, gamma: float) -> float:
    """Alias for θ̇_eff(x, ẋ, γ)."""
    return theta_dot_from_x(x, xdot, gamma)


def x_eom_rhs(x: float, xdot: float, H: float, gamma: float) -> float:
    """
    Right-hand side for ẍ in the regularized variable x = artanh(γ θ / √3).

    Derived from θ̈ + 3H θ̇ = θ̇² / [θ (1 - γ² θ² / 3)].
    """
    if abs(x) < 1e-9:
        # series expansion: coth x ≈ 1/x + x/3
        coth_x = (1.0 / max(x, 1e-12)) + x / 3.0
    else:
        coth_x = 1.0 / math.tanh(x)
    return -3.0 * H * xdot + (xdot * xdot) * (2.0 * math.tanh(x) + coth_x)


# ---------------------------------------------------------------------------
# CLI helper

def _pretty(expr) -> str:
    return sp.pretty(expr, use_unicode=False)


if __name__ == "__main__":
    print("L_eff(theta) =")
    print(_pretty(lagrangian_expr))
    print("\nEnergy density rho =")
    print(_pretty(energy_density_expr))
    print("\nPressure p =")
    print(_pretty(pressure_expr))
    print("\nEuler-Lagrange for theta:")
    print(_pretty(eom_theta_expr))
    print("\nEuler-Lagrange for R:")
    print(_pretty(eom_R_expr))
