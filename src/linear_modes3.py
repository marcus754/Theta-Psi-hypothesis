# -*- coding: utf-8 -*-
"""
3×3 linear modes scaffold for (δθ, δR, δφ) with complex Ψ = R e^{iφ}.

Kinetic matrix K and gradient matrix G are modeled as diagonal by default with
physically motivated weights:
  Kθθ = 3/θ^2 + γ^2,  KRR = 1,  Kφφ = R^2
  Gθθ = cθ^2 * Kθθ,   GRR = cR^2, Gφφ = cφ^2 * R^2

Optional small mixings epsK_ij, epsG_ij are supported for future studies.

This module is a stepping stone before full ADM-based derivation of K,G,M.
"""
from __future__ import annotations

from typing import Dict, Tuple
import math
import sympy as sp


def kinetic_matrix_3(
    theta: float,
    gamma: float,
    R: float,
    *,
    epsK_theta_R: float = 0.0,
    epsK_theta_phi: float = 0.0,
    epsK_R_phi: float = 0.0,
):
    if theta == 0.0:
        raise ValueError("theta=0 is singular; work in x-variable to avoid exact zero.")
    Kth = 3.0 / (theta * theta) + gamma * gamma
    KRR = 1.0
    Kpp = R * R
    return (
        (Kth, epsK_theta_R, epsK_theta_phi),
        (epsK_theta_R, KRR, epsK_R_phi),
        (epsK_theta_phi, epsK_R_phi, Kpp),
    )


def gradient_matrix_3(
    theta: float,
    gamma: float,
    R: float,
    *,
    ctheta2: float = 1.0,
    cR2: float = 1.0,
    cphi2: float = 1.0,
    epsG_theta_R: float = 0.0,
    epsG_theta_phi: float = 0.0,
    epsG_R_phi: float = 0.0,
):
    Kth = 3.0 / (theta * theta) + gamma * gamma
    Gth = ctheta2 * Kth
    GRR = cR2
    Gpp = cphi2 * (R * R)
    return (
        (Gth, epsG_theta_R, epsG_theta_phi),
        (epsG_theta_R, GRR, epsG_R_phi),
        (epsG_theta_phi, epsG_R_phi, Gpp),
    )


def _to_sympy(M):
    return sp.Matrix([[M[i][j] for j in range(3)] for i in range(3)])


def cs2_eigs_3(
    theta: float,
    gamma: float,
    R: float,
    *,
    ctheta2: float = 1.0,
    cR2: float = 1.0,
    cphi2: float = 1.0,
    epsK_theta_R: float = 0.0,
    epsK_theta_phi: float = 0.0,
    epsK_R_phi: float = 0.0,
    epsG_theta_R: float = 0.0,
    epsG_theta_phi: float = 0.0,
    epsG_R_phi: float = 0.0,
):
    K = kinetic_matrix_3(theta, gamma, R, epsK_theta_R=epsK_theta_R, epsK_theta_phi=epsK_theta_phi, epsK_R_phi=epsK_R_phi)
    G = gradient_matrix_3(theta, gamma, R, ctheta2=ctheta2, cR2=cR2, cphi2=cphi2,
                          epsG_theta_R=epsG_theta_R, epsG_theta_phi=epsG_theta_phi, epsG_R_phi=epsG_R_phi)
    Ksym = _to_sympy(K)
    Gsym = _to_sympy(G)
    Asym = Ksym.LUsolve(Gsym)  # K^{-1} G (solve is better conditioned than explicit inverse)
    eigs = [complex(ev) for ev in Asym.eigenvals().keys()]
    # Return real parts (imaginary parts should be ~0 for symmetric K,G)
    return tuple(float(ev.real) for ev in eigs)


def _is_pd_sym(M) -> bool:
    # Sylvester via principal minors with sympy
    Ms = _to_sympy(M)
    for i in range(1, 4):
        if Ms[:i, :i].det() <= 0:
            return False
    return True


def stability_3(
    theta: float,
    gamma: float,
    R: float,
    *,
    ctheta2: float = 1.0,
    cR2: float = 1.0,
    cphi2: float = 1.0,
    epsK_theta_R: float = 0.0,
    epsK_theta_phi: float = 0.0,
    epsK_R_phi: float = 0.0,
    epsG_theta_R: float = 0.0,
    epsG_theta_phi: float = 0.0,
    epsG_R_phi: float = 0.0,
) -> Dict[str, object]:
    K = kinetic_matrix_3(theta, gamma, R, epsK_theta_R=epsK_theta_R, epsK_theta_phi=epsK_theta_phi, epsK_R_phi=epsK_R_phi)
    G = gradient_matrix_3(theta, gamma, R, ctheta2=ctheta2, cR2=cR2, cphi2=cphi2,
                          epsG_theta_R=epsG_theta_R, epsG_theta_phi=epsG_theta_phi, epsG_R_phi=epsG_R_phi)
    no_ghost = _is_pd_sym(K)
    eigs = cs2_eigs_3(theta, gamma, R, ctheta2=ctheta2, cR2=cR2, cphi2=cphi2,
                      epsK_theta_R=epsK_theta_R, epsK_theta_phi=epsK_theta_phi, epsK_R_phi=epsK_R_phi,
                      epsG_theta_R=epsG_theta_R, epsG_theta_phi=epsG_theta_phi, epsG_R_phi=epsG_R_phi)
    no_grad = all(e > 0.0 for e in eigs)
    return {"no_ghost": no_ghost, "no_gradient_instability": no_grad, "cs2_eigs": eigs, "K": K, "G": G}

