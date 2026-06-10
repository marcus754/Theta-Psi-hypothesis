# -*- coding: utf-8 -*-
"""
Canonical symbolic ledger for the covariant Θ–Ψ action.

This module formalizes the canonical reduced action in one place and exposes
the exact same blocks that appear in the docs and tests.

Canonical reduced action:

    S_red = ∫ √(-g) [ EH + L_Θ + L_Ψ - U_red + M_*^4 F(J_refr) ] d^4x

with
    theta_eff^2 = Θ_{μν} u^μ u^ν,
    Φ_eff = ln(theta_eff / theta_0),
    J_refr = Φ_eff * I_grad,
    I_grad = |∇_⊥ θ|^2 + w_Ψ |∇_⊥ Ψ|^2.

In the full interpretation Ψ is a complex order parameter `Ψ = R e^{iφ}`.
This ledger keeps the reduced amplitude representative and the FRW phase charge
handled at the background level.

The module keeps this in symbolic form so the downstream docs, tests and
diagnostics stay aligned.
"""
from __future__ import annotations

from typing import Dict

import sympy as sp


t = sp.symbols("t", real=True)
gamma, V0 = sp.symbols("gamma V0", nonnegative=True, real=True)
m_theta2, lam_theta, lam_theta_psi = sp.symbols("m_theta2 lam_theta lam_theta_psi", real=True)
Mstar4 = sp.symbols("Mstar4", positive=True, real=True)
theta0 = sp.symbols("theta_0", positive=True, real=True)


def theta_eff_projection() -> Dict[str, sp.Expr]:
    """Return symbolic theta_eff projection from Θ_{μν} onto u^μ."""
    Theta = sp.MatrixSymbol("Theta", 4, 4)
    u = sp.MatrixSymbol("u", 4, 1)
    theta_sq = sum(Theta[i, j] * u[i, 0] * u[j, 0] for i in range(4) for j in range(4))
    theta_eff = sp.sqrt(theta_sq)
    return {"Theta": Theta, "u": u, "theta_eff_sq": theta_sq, "theta_eff": theta_eff}


def phi_eff_canonical(theta_eff: sp.Expr, theta_0: sp.Expr | None = None) -> sp.Expr:
    """Canonical bridge scalar Φ_eff = ln(theta_eff/theta_0)."""
    th0 = theta0 if theta_0 is None else theta_0
    return sp.log(theta_eff / th0)


def reduced_potential(theta_eff: sp.Expr, R: sp.Expr) -> sp.Expr:
    """Return the canonical reduced potential."""
    return (
        sp.Rational(1, 2) * m_theta2 * theta_eff**2
        + sp.Rational(1, 4) * lam_theta * theta_eff**4
        + sp.Rational(1, 2) * V0 * R**2
        + sp.Rational(1, 2) * lam_theta_psi * theta_eff**2 * R**2
    )


def canonical_action_ledger() -> Dict[str, sp.Expr]:
    """Return the current canonical action blocks with explicit labels."""
    proj = theta_eff_projection()
    theta_eff = proj["theta_eff"]
    R = sp.Symbol("R", nonnegative=True, real=True)
    phi = phi_eff_canonical(theta_eff)
    i_grad = sp.Symbol("I_grad", nonnegative=True, real=True)
    j_refr = sp.Symbol("J_refr", nonnegative=True, real=True)
    l_theta = sp.Rational(3, 2) * sp.Symbol("grad_ln_theta_sq", nonnegative=True, real=True) + sp.Rational(1, 2) * gamma**2 * sp.Symbol("grad_theta_sq", nonnegative=True, real=True)
    l_R = sp.Rational(1, 2) * sp.Symbol("grad_R_sq", nonnegative=True, real=True)
    u_red = reduced_potential(theta_eff, R)
    l_refr = Mstar4 * sp.Function("F")(j_refr)
    s_red = sp.Symbol("sqrt(-g)") * (sp.Symbol("EH") + sp.Symbol("L_Theta") + sp.Symbol("L_R") - u_red + l_refr)
    return {
        **proj,
        "R": R,
        "Phi_eff": phi,
        "I_grad": i_grad,
        "J_refr": j_refr,
        "L_Theta": l_theta,
        "L_R": l_R,
        "U_red": u_red,
        "L_refr": l_refr,
        "S_red": s_red,
    }


def full_covariant_action_ledger() -> Dict[str, sp.Expr]:
    """Return the full covariant completion built from the reduced canon."""
    red = canonical_action_ledger()
    g_hat = sp.MatrixSymbol("g_hat", 4, 4)
    matter = sp.Symbol("matter")
    s_matter = sp.Symbol("S_m[g_hat,matter]")
    s_cov = red["S_red"] + s_matter
    return {
        **red,
        "g_hat": g_hat,
        "matter": matter,
        "S_matter": s_matter,
        "S_cov": s_cov,
    }


def frw_reduced_lagrangian() -> Dict[str, sp.Expr]:
    """Return the current homogeneous FRW reduced Lagrangian and pieces."""
    theta = sp.Function("theta")(t)
    R = sp.Function("R")(t)
    theta_dot = sp.diff(theta, t)
    R_dot = sp.diff(R, t)

    k_theta = sp.Rational(3, 2) * (theta_dot / theta) ** 2 + sp.Rational(1, 2) * gamma**2 * theta_dot**2
    k_R = sp.Rational(1, 2) * R_dot**2
    u_red = sp.Rational(1, 2) * V0 * R**2
    l_eff = sp.simplify(k_theta + k_R - u_red)
    return {"theta": theta, "R": R, "K_theta": k_theta, "K_R": k_R, "U_frw": u_red, "L_eff": l_eff}


def stationary_refractive_sector() -> Dict[str, sp.Expr]:
    """Return the symbolic stationary J_refr sector."""
    r = sp.symbols("r", positive=True, real=True)
    theta = sp.Function("theta")(r)
    R = sp.Function("R")(r)
    theta_p = sp.diff(theta, r)
    R_p = sp.diff(R, r)
    wR = sp.symbols("w_R", positive=True, real=True)
    i_grad = theta_p**2 + wR * R_p**2
    phi = phi_eff_canonical(theta)
    j = sp.simplify(phi * i_grad)
    F = sp.Function("F")
    l_refr = Mstar4 * F(j)
    n = sp.Function("n")(j)
    return {
        "r": r,
        "theta": theta,
        "R": R,
        "I_grad": i_grad,
        "Phi_eff": phi,
        "J_refr": j,
        "n_of_J": n,
        "L_refr": l_refr,
    }


def reduced_stress_energy_blocks() -> Dict[str, sp.Expr]:
    """Return formal reduced stress-energy building blocks."""
    # Formal covariant placeholders; these are not expanded in a chosen metric.
    grad_ln_theta_sq = sp.Symbol("grad_ln_theta_sq", real=True)
    grad_theta_sq = sp.Symbol("grad_theta_sq", real=True)
    grad_psi_sq = sp.Symbol("grad_psi_sq", real=True)
    i_grad = sp.Symbol("I_grad", nonnegative=True, real=True)
    j_refr = sp.Symbol("J_refr", nonnegative=True, real=True)
    T_theta = 3 * grad_ln_theta_sq + gamma**2 * grad_theta_sq
    T_R = grad_psi_sq
    T_refr = Mstar4 * sp.Function("F")(j_refr)
    return {
        "I_grad": i_grad,
        "J_refr": j_refr,
        "T_theta_scalar_block": T_theta,
        "T_R_scalar_block": T_R,
        "T_refr_scalar_block": T_refr,
    }


def minisuperspace_eom_ledger() -> Dict[str, sp.Expr]:
    """Return the FRW reduced equations of motion and source blocks."""
    theta = sp.Function("theta")(t)
    R = sp.Function("R")(t)
    theta_dot = sp.diff(theta, t)
    theta_ddot = sp.diff(theta, t, 2)
    R_dot = sp.diff(R, t)
    R_ddot = sp.diff(R, t, 2)

    e_theta = sp.simplify((sp.Integer(3) / theta**2 + gamma**2) * theta_ddot - sp.Integer(3) * theta_dot**2 / theta**3)
    e_R = sp.simplify(R_ddot + V0 * R)
    rho_red = sp.simplify(
        sp.Rational(3, 2) * (theta_dot / theta) ** 2
        + sp.Rational(1, 2) * R_dot**2
        + sp.Rational(1, 2) * V0 * R**2
        + sp.Rational(1, 2) * gamma**2 * theta_dot**2
    )
    p_red = sp.simplify(
        sp.Rational(3, 2) * (theta_dot / theta) ** 2
        + sp.Rational(1, 2) * R_dot**2
        - sp.Rational(1, 2) * V0 * R**2
        + sp.Rational(1, 2) * gamma**2 * theta_dot**2
    )
    t00_red = rho_red
    return {
        "theta": theta,
        "R": R,
        "E_theta": e_theta,
        "E_R": e_R,
        "rho_red": rho_red,
        "p_red": p_red,
        "T00_red": t00_red,
    }


def covariant_action_summary() -> Dict[str, sp.Expr]:
    """Return the canonical symbolic summary of the reduced covariant action."""
    return {
        **full_covariant_action_ledger(),
        **reduced_stress_energy_blocks(),
        **minisuperspace_eom_ledger(),
        **frw_reduced_lagrangian(),
        **stationary_refractive_sector(),
        "status": sp.Symbol("canonical_covariant"),
    }
