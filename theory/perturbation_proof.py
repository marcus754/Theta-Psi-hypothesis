# -*- coding: utf-8 -*-
"""Canonical perturbation ledger for the Θ–Ψ completion.

This module collects the full covariant action ledger and the ADM perturbation
ledger in one place so the repo can inspect a single proof object:

1. full covariant completion `S_cov = S_red + S_m[g_hat, matter]`;
2. quadratic perturbation sector `(K, G, M)`;
3. strict direct-vs-Schur comparison for the reduced ADM completion;
4. FRW/stationary health summary from `theory.twofield_stability`.
"""
from __future__ import annotations

from typing import Dict, Any

import sympy as sp

from theory.adm_symbolic import effective_matrices_sigma_strict_compare
from theory.covariant_action import full_covariant_action_ledger
from theory.twofield_stability import stability_ledger


def _matrix_close(A: sp.Matrix, B: sp.Matrix, tol: float = 1e-9) -> bool:
    if A.shape != B.shape:
        return False
    for i in range(A.shape[0]):
        for j in range(A.shape[1]):
            diff = sp.N(A[i, j] - B[i, j])
            if abs(complex(diff)) > tol:
                return False
    return True


def perturbation_proof_ledger(
    theta: float,
    R: float,
    a: float,
    k: float,
    *,
    gamma: float,
    Q: float,
    V0: float,
    H: float,
    theta_dot: float,
    R_dot: float,
    Hdot: float,
) -> Dict[str, Any]:
    """Return the full perturbation proof ledger for the current canon."""
    action = full_covariant_action_ledger()
    direct, schur = effective_matrices_sigma_strict_compare(
        theta,
        R,
        a,
        k,
        gamma=gamma,
        Q=Q,
        V0=V0,
        H=H,
        theta_dot=theta_dot,
        R_dot=R_dot,
        Hdot=Hdot,
    )
    (K_direct, G_direct, M_direct) = direct
    (K_schur, G_schur, M_schur) = schur
    health = stability_ledger(
        theta=theta,
        gamma=gamma,
        theta_prime=theta_dot,
        psi_prime=R_dot,
        mtheta2=0.0,
        V0=V0,
    )
    routes_match = bool(
        _matrix_close(K_direct, K_schur)
        and _matrix_close(G_direct, G_schur)
        and _matrix_close(M_direct, M_schur)
    )
    proof_closed = bool(routes_match and health["healthy"])
    return {
        **action,
        "K_direct": K_direct,
        "G_direct": G_direct,
        "M_direct": M_direct,
        "K_schur": K_schur,
        "G_schur": G_schur,
        "M_schur": M_schur,
        "routes_match": routes_match,
        "health": health,
        "proof_closed": proof_closed,
        "status": sp.Symbol("closed_ledger"),
    }
