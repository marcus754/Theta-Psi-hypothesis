# -*- coding: utf-8 -*-
"""Derivation ledgers from the current Θ–Ψ action.

The purpose of this module is deliberately narrow: separate quantities that
follow from the written action from quantities that still require an extra
operator, boundary condition, or sourced solution. It is an anti-fit layer, not
another fitting interface.
"""
from __future__ import annotations

from typing import Dict

import sympy as sp


def adm_minimal_quadratic_derivation() -> Dict[str, object]:
    """Return the ADM quadratic blocks implied by the current minimal action."""
    theta, R, gamma, V0, Q, a = sp.symbols(
        "theta R gamma V0 Q a", positive=True, real=True
    )

    Ktheta = sp.Integer(3) / theta**2 + gamma**2
    K = sp.diag(Ktheta, sp.Integer(1), R**2)
    G = sp.diag(Ktheta, sp.Integer(1), R**2)
    M = sp.diag(
        sp.Integer(0),
        V0 + sp.Integer(3) * Q**2 / (a**6 * R**4),
        sp.Integer(0),
    )
    zero_mixings = {
        "epsK_theta_R": sp.Integer(0),
        "epsK_theta_phi": sp.Integer(0),
        "epsK_R_phi": sp.Integer(0),
        "epsG_theta_R": sp.Integer(0),
        "epsG_theta_phi": sp.Integer(0),
        "epsG_R_phi": sp.Integer(0),
        "alpha_k2_theta": sp.Integer(0),
        "alpha_k2_R": sp.Integer(0),
        "alpha_k2_phi": sp.Integer(0),
    }
    return {
        "fields": ("delta_theta", "delta_R", "delta_phi"),
        "K_from_action": K,
        "G_from_action": G,
        "M_from_action": M,
        "derived_zero_mixings": zero_mixings,
        "cs2_from_action": (sp.Integer(1), sp.Integer(1), sp.Integer(1)),
        "requires_extra_operators_for_nonzero_mixings": sp.Integer(1),
        "status": sp.Symbol("adm_minimal_quadratic_derived_from_action"),
    }


def refractive_function_derivation_status() -> Dict[str, object]:
    """Return what the current action derives about F(J_refr)."""
    J = sp.Symbol("J_refr", nonnegative=True, real=True)
    F = sp.Function("F")
    return {
        "action_term": sp.Symbol("Mstar4") * F(J),
        "derived_argument": J,
        "derived_constraints": (
            "F(0)=0",
            "F'(J)>=0",
            "finite response required for no-horizon branch",
        ),
        "not_derived_function_choice": "F_v1(J)=1-exp(-J) after normalization; J_star is not physical",
        "reason": "the current action names F(J_refr) but does not specify its functional form; an explicit scale can be absorbed into J_refr",
        "status": sp.Symbol("refractive_function_not_yet_derived"),
    }


def refractive_semigroup_principle() -> Dict[str, object]:
    """Return the closure principle that makes the exponential form natural."""
    S = sp.Function("S")
    J1 = sp.Symbol("J1", nonnegative=True, real=True)
    J2 = sp.Symbol("J2", nonnegative=True, real=True)
    J = sp.Symbol("J", nonnegative=True, real=True)
    J_star = sp.Symbol("J_star", positive=True, real=True)
    semigroup_law = sp.Eq(S(J1 + J2), S(J1) * S(J2))
    differential_law = sp.Eq(sp.diff(S(J), J), -S(J) / J_star)
    normalized_solution = sp.Eq(S(J), sp.exp(-J / J_star))
    F_solution = sp.Eq(1 - S(J), 1 - sp.exp(-J / J_star))
    return {
        "residual": "S(J) = 1 - F(J)",
        "semigroup_law": semigroup_law,
        "differential_law": differential_law,
        "normalized_solution": normalized_solution,
        "F_solution": F_solution,
        "interpretation": (
            "independent increments act multiplicatively on the remaining response"
        ),
        "status": sp.Symbol("refractive_semigroup_principle"),
    }


def strongfield_bridge_derivation_status() -> Dict[str, object]:
    """Return the derivation status of compact-object bridge scales."""
    return {
        "derived_requirement": "solve sourced static Euler-Lagrange equations",
        "not_derived_scales": (
            "theta_c_scale",
            "m2_scale",
            "lam_scale",
            "powerlaw2d_coefficients",
        ),
        "allowed_before_derivation": "diagnostic scans only",
        "forbidden_as_canon": sp.Integer(1),
        "status": sp.Symbol("strongfield_bridge_not_yet_derived"),
    }


def derivation_status_ledger() -> Dict[str, object]:
    """Return the current action-derivation status for risky knobs."""
    return {
        "adm": adm_minimal_quadratic_derivation(),
        "refractive": refractive_function_derivation_status(),
        "refractive_semigroup": refractive_semigroup_principle(),
        "strongfield_bridge": strongfield_bridge_derivation_status(),
        "status": sp.Symbol("action_derivation_status"),
    }
