# -*- coding: utf-8 -*-
"""Superluminality veto using K and G from the linear-mode scaffold.

We flag a configuration as superluminal when any high-k sound speed squared
from the eigenvalues of K^{-1} G exceeds 1 by more than a small tolerance.

Parameters `H` and `V0` are kept for forward compatibility with a future
ADM-based implementation where gradient coefficients may depend on background
curvature and potential; they are not used directly in this placeholder.
Optionally support small kinetic/gradient mixing via `epsK`, `epsG`.
"""
from __future__ import annotations

from src.linear_modes import cs2_eigs, kinetic_matrix, gradient_matrix
from src.linear_modes_adm import cs2_eigs_adm, stability_finite_k_adm


def has_superluminal(
    theta: float,
    H: float,
    gamma: float,
    V0: float,
    ctheta2: float = 1.0,
    cpsi2: float = 1.0,
    *,
    epsK: float = 0.0,
    epsG: float = 0.0,
    tol: float = 1e-9,
) -> bool:
    """Return True if any scalar mode has c_s^2 > 1 + tol.

    Uses eigenvalues of K^{-1} G computed by the current linear-mode scaffold.
    Includes optional small mixing (epsK, epsG). A tiny tolerance `tol` avoids
    false positives due to round-off near exactly luminal propagation.
    """
    # Construct matrices (kept here for potential sanity checks in the future)
    _ = kinetic_matrix(theta, gamma, epsK)
    _ = gradient_matrix(theta, gamma, ctheta2, cpsi2, epsG)
    c1, c2 = cs2_eigs(theta, gamma, ctheta2, cpsi2, epsK, epsG)
    return (c1 > 1.0 + tol) or (c2 > 1.0 + tol)


def superluminal_details(
    theta: float,
    H: float,
    gamma: float,
    V0: float,
    ctheta2: float = 1.0,
    cpsi2: float = 1.0,
    *,
    epsK: float = 0.0,
    epsG: float = 0.0,
    tol: float = 1e-9,
):
    """Return a dict with eigenvalues and flag for superluminality.

    Useful for debugging/reporting in scans.
    """
    c1, c2 = cs2_eigs(theta, gamma, ctheta2, cpsi2, epsK, epsG)
    return {
        "cs2_eigs": (c1, c2),
        "superluminal": (c1 > 1.0 + tol) or (c2 > 1.0 + tol),
        "tol": tol,
        "epsK": epsK,
        "epsG": epsG,
    }


def has_superluminal_adm(
    theta: float,
    R: float,
    a: float,
    k: float,
    gamma: float,
    Q: float,
    V0: float,
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
    tol: float = 1e-9,
    alpha_k2_theta: float = 0.0,
    alpha_k2_R: float = 0.0,
    alpha_k2_phi: float = 0.0,
    # gravitational mass corrections
    gravM: bool = False,
    H_list = None,
    theta_dot_list = None,
    R_dot_list = None,
    Hdot_list = None,
    adm_strict: bool = False,
) -> bool:
    """ADM-flavored superluminality check with 3×3 modes.

    Returns True if any eigenmode has c_s^2 > 1 + tol.
    In the current approximation the eigenvalues are high-k limits and do not
    depend explicitly on k (we accept k for forward compatibility).
    """
    # If strict mode and background derivatives are provided via lists, use the first element
    H = None
    thd = None
    Rd = None
    Hdot = None
    if adm_strict and H_list is not None and len(H_list):
        H = H_list[0]
    if adm_strict and theta_dot_list is not None and len(theta_dot_list):
        thd = theta_dot_list[0]
    if adm_strict and R_dot_list is not None and len(R_dot_list):
        Rd = R_dot_list[0]
    if adm_strict and Hdot_list is not None and len(Hdot_list):
        Hdot = Hdot_list[0]
    c = cs2_eigs_adm(
        theta, R, a, k, gamma, Q, V0,
        ctheta2=ctheta2, cR2=cR2, cphi2=cphi2,
        epsK_theta_R=epsK_theta_R, epsK_theta_phi=epsK_theta_phi, epsK_R_phi=epsK_R_phi,
        epsG_theta_R=epsG_theta_R, epsG_theta_phi=epsG_theta_phi, epsG_R_phi=epsG_R_phi,
        alpha_k2_theta=alpha_k2_theta, alpha_k2_R=alpha_k2_R, alpha_k2_phi=alpha_k2_phi,
        elim=('sigma_strict' if adm_strict else None),
        H=H, theta_dot=thd, R_dot=Rd, Hdot=Hdot,
    )
    return any(ci > 1.0 + tol for ci in c)


def finite_k_stable_along_traj(
    theta_list,
    R_list,
    a_list,
    gamma: float,
    Q: float,
    V0: float,
    *,
    ctheta2: float = 1.0,
    cR2: float = 1.0,
    cphi2: float = 1.0,
    kmin: float = 1e-2,
    kmax: float = 1.0,
    nk: int = 5,
    epsK_theta_R: float = 0.0,
    epsK_theta_phi: float = 0.0,
    epsK_R_phi: float = 0.0,
    epsG_theta_R: float = 0.0,
    epsG_theta_phi: float = 0.0,
    epsG_R_phi: float = 0.0,
    alpha_k2_theta: float = 0.0,
    alpha_k2_R: float = 0.0,
    alpha_k2_phi: float = 0.0,
    # gravitational mass corrections
    gravM: bool = False,
    H_list = None,
    theta_dot_list = None,
    R_dot_list = None,
    Hdot_list = None,
    adm_strict: bool = False,
) -> bool:
    """Check ω^2(k) > 0 for all modes on a small k-grid along a trajectory."""
    import numpy as np
    ks = np.geomspace(kmin, kmax, max(nk, 2))
    for i, (th, R, a) in enumerate(zip(theta_list, R_list, a_list)):
        for k in ks:
            Hi = None if H_list is None or i >= len(H_list) else H_list[i]
            thd = None if theta_dot_list is None or i >= len(theta_dot_list) else theta_dot_list[i]
            Rd = None if R_dot_list is None or i >= len(R_dot_list) else R_dot_list[i]
            Hdi = None if Hdot_list is None or i >= len(Hdot_list) else Hdot_list[i]
            st = stability_finite_k_adm(
                th,
                R,
                a,
                float(k),
                gamma,
                Q,
                V0,
                ctheta2=ctheta2,
                cR2=cR2,
                cphi2=cphi2,
                epsK_theta_R=epsK_theta_R,
                epsK_theta_phi=epsK_theta_phi,
                epsK_R_phi=epsK_R_phi,
                epsG_theta_R=epsG_theta_R,
                epsG_theta_phi=epsG_theta_phi,
                epsG_R_phi=epsG_R_phi,
                alpha_k2_theta=alpha_k2_theta,
                alpha_k2_R=alpha_k2_R,
                alpha_k2_phi=alpha_k2_phi,
                elim=('sigma_strict' if adm_strict else None),
                gravM=gravM,
                H=Hi,
                theta_dot=thd,
                R_dot=Rd,
                Hdot=Hdi,
            )
            if not (st["no_ghost"] and st["no_tachyon"]):
                return False
    return True
