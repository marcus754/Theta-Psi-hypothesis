# -*- coding: utf-8 -*-
"""
ADM-based linear modes (3×3) for (δθ, δR, δφ) — numerically robust version.

Builds K(a,k), G(a,k), M(a,k) via theory.adm_symbolic and computes
generalized eigenvalues using a Cholesky-based reduction:
  G v = λ K v  ⇒  (L^{-1} G L^{-T}) w = λ w,  with K = L L^T.

Provides high-k sound speeds (eigs of K^{-1}G) and finite-k dispersions
ω^2 (eigs of K^{-1}[G k^2/a^2 + M]). This replaces prior SymPy eigenvalue
calls with fast, deterministic NumPy routines.
"""
from __future__ import annotations

from typing import Dict, Tuple
import numpy as np

from theory.adm_symbolic import (
    kinetic_matrix,
    gradient_matrix,
    mass_matrix,
    effective_matrices_sym,
    mass_matrix_sigma_total,
    effective_matrices_sigma_strict,
    effective_matrices_sigma_strict_schur,
)


def _to_ndarray(M) -> np.ndarray:
    return np.array(M, dtype=float)


def _chol_pd(M: np.ndarray) -> Tuple[np.ndarray, bool]:
    """Try Cholesky; return (L, ok) with K = L L^T."""
    try:
        L = np.linalg.cholesky(M)
        return L, True
    except np.linalg.LinAlgError:
        return np.empty((0, 0)), False


def _gen_eigs(K: np.ndarray, B: np.ndarray) -> Tuple[float, ...]:
    """Generalized eigenvalues of (B, K) using Cholesky reduction."""
    L, ok = _chol_pd(K)
    if not ok:
        # Fallback: try symmetric solve via eig(K^{-1}B)
        A = np.linalg.solve(K, B)
        vals = np.linalg.eigvals(A)
        vals = np.real(vals)
        return tuple(sorted(float(v) for v in vals))
    # Symmetrize numerically to reduce round-off
    Linv = np.linalg.inv(L)
    S = Linv @ B @ Linv.T
    S = 0.5 * (S + S.T)
    vals = np.linalg.eigvalsh(S)
    return tuple(float(v) for v in vals)


def _augment_with_sigma(
    Kd: np.ndarray,
    Gd: np.ndarray,
    Md: np.ndarray | None = None,
    *,
    include_sigma: bool = False,
    csigma2: float = 1.0,
    msigma2: float = 0.0,
):
    if not include_sigma:
        return Kd, Gd, Md
    n = int(Kd.shape[0])
    Kx = np.zeros((n + 1, n + 1), dtype=float)
    Gx = np.zeros((n + 1, n + 1), dtype=float)
    Kx[:n, :n] = Kd
    Gx[:n, :n] = Gd
    Kx[n, n] = 1.0
    Gx[n, n] = float(csigma2)
    if Md is None:
        return Kx, Gx, None
    Mx = np.zeros((n + 1, n + 1), dtype=float)
    Mx[:n, :n] = Md
    Mx[n, n] = float(msigma2)
    return Kx, Gx, Mx


def cs2_eigs_adm(
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
    alpha_k2_theta: float = 0.0,
    alpha_k2_R: float = 0.0,
    alpha_k2_phi: float = 0.0,
    elim: dict | str | bool | None = None,
    gravM: bool = False,
    H: float | None = None,
    theta_dot: float | None = None,
    R_dot: float | None = None,
    Hdot: float | None = None,
    include_sigma: bool = False,
    csigma2: float = 1.0,
    msigma2: float = 0.0,
):
    elim_kwargs = None
    # Interpret elim flag: None → no elimination; dict → explicit couplings; 'sigma_strict' → strict σ-model
    if isinstance(elim, str) and elim == 'sigma_strict':
        Ke, Ge, _ = effective_matrices_sigma_strict_schur(
            theta, R, a, k,
            gamma=gamma, Q=Q, V0=V0,
            H=H,
            theta_dot=theta_dot,
            R_dot=R_dot,
            Hdot=Hdot,
        )
        Kd = _to_ndarray(Ke)
        Gd = _to_ndarray(Ge)
    elif isinstance(elim, dict):
        elim_kwargs = elim
    else:
        elim_kwargs = None

    if elim_kwargs is None:
        K = kinetic_matrix(
            theta,
            R,
            a,
            k,
            gamma,
            Q,
            V0,
            epsK_theta_R=epsK_theta_R,
            epsK_theta_phi=epsK_theta_phi,
            epsK_R_phi=epsK_R_phi,
        )
        G = gradient_matrix(
            theta,
            R,
            a,
            k,
            gamma,
            Q,
            V0,
            ctheta2=ctheta2,
            cR2=cR2,
            cphi2=cphi2,
            epsG_theta_R=epsG_theta_R,
            epsG_theta_phi=epsG_theta_phi,
            epsG_R_phi=epsG_R_phi,
            alpha_k2_theta=alpha_k2_theta,
            alpha_k2_R=alpha_k2_R,
            alpha_k2_phi=alpha_k2_phi
        )
        Kd = _to_ndarray(K)
        Gd = _to_ndarray(G)
    else:
        Ke, Ge, _ = effective_matrices_sym(
            theta,
            R,
            a,
            k,
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
            **elim_kwargs,
        )
        Kd = _to_ndarray(Ke)
        Gd = _to_ndarray(Ge)
    Kd, Gd, _ = _augment_with_sigma(
        Kd,
        Gd,
        None,
        include_sigma=include_sigma,
        csigma2=csigma2,
        msigma2=msigma2,
    )
    eigs = _gen_eigs(Kd, Gd)
    return eigs


def stability_adm(
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
    alpha_k2_theta: float = 0.0,
    alpha_k2_R: float = 0.0,
    alpha_k2_phi: float = 0.0,
    elim: dict | str | bool | None = None,
    H: float | None = None,
    theta_dot: float | None = None,
    R_dot: float | None = None,
    Hdot: float | None = None,
    include_sigma: bool = False,
    csigma2: float = 1.0,
    msigma2: float = 0.0,
) -> Dict[str, object]:
    # Interpret elim flag: None → no elimination; dict → explicit couplings; 'sigma_strict' → strict σ-model
    if isinstance(elim, str) and elim == 'sigma_strict':
        Ke, Ge, _ = effective_matrices_sigma_strict_schur(
            theta, R, a, k,
            gamma=gamma, Q=Q, V0=V0,
            H=(0.0 if H is None else H),
            theta_dot=(0.0 if theta_dot is None else theta_dot),
            R_dot=(0.0 if R_dot is None else R_dot),
            Hdot=(0.0 if Hdot is None else Hdot),
        )
        Kd = _to_ndarray(Ke)
        G = Ge
    elif isinstance(elim, dict):
        elim_kwargs = elim
    else:
        elim_kwargs = None

    if elim_kwargs is None:
        K = kinetic_matrix(
            theta,
            R,
            a,
            k,
            gamma,
            Q,
            V0,
            epsK_theta_R=epsK_theta_R,
            epsK_theta_phi=epsK_theta_phi,
            epsK_R_phi=epsK_R_phi,
        )
        G = gradient_matrix(
            theta,
            R,
            a,
            k,
            gamma,
            Q,
            V0,
            ctheta2=ctheta2,
            cR2=cR2,
            cphi2=cphi2,
            epsG_theta_R=epsG_theta_R,
            epsG_theta_phi=epsG_theta_phi,
            epsG_R_phi=epsG_R_phi,
            alpha_k2_theta=alpha_k2_theta,
            alpha_k2_R=alpha_k2_R,
            alpha_k2_phi=alpha_k2_phi
        )
        Kd = _to_ndarray(K)
    else:
        Ke, G, _ = effective_matrices_sym(
            theta,
            R,
            a,
            k,
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
            **elim_kwargs,
        )
        Kd = _to_ndarray(Ke)
    Gd = _to_ndarray(G)
    Kd, Gd, _ = _augment_with_sigma(
        Kd,
        Gd,
        None,
        include_sigma=include_sigma,
        csigma2=csigma2,
        msigma2=msigma2,
    )
    # Ghost-free if K is SPD (Cholesky succeeds)
    _, no_ghost = _chol_pd(Kd)
    eigs = cs2_eigs_adm(
        theta,
        R,
        a,
        k,
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
        elim=elim_kwargs,
        H=H,
        theta_dot=theta_dot,
        R_dot=R_dot,
        Hdot=Hdot,
        include_sigma=include_sigma,
        csigma2=csigma2,
        msigma2=msigma2,
    )
    # Allow tiny negatives from numerical round-off
    no_grad = all(e > -1e-12 for e in eigs)
    if elim_kwargs is None:
        return {"no_ghost": no_ghost, "no_gradient_instability": no_grad, "cs2_eigs": eigs, "K": Kd, "G": Gd}
    else:
        Ke, Ge, _ = effective_matrices_sym(
            theta,
            R,
            a,
            k,
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
            **elim_kwargs,
        )
        return {"no_ghost": no_ghost, "no_gradient_instability": no_grad, "cs2_eigs": eigs, "K": Kd, "G": Gd}


def dispersion_eigs_adm(
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
    alpha_k2_theta: float = 0.0,
    alpha_k2_R: float = 0.0,
    alpha_k2_phi: float = 0.0,
    elim: dict | str | bool | None = None,
    gravM: bool = False,
    H: float | None = None,
    theta_dot: float | None = None,
    R_dot: float | None = None,
    Hdot: float | None = None,
    include_sigma: bool = False,
    csigma2: float = 1.0,
    msigma2: float = 0.0,
):
    # Interpret elim flag: None → no elimination; dict → explicit couplings; 'sigma_strict' → strict σ-model
    prepared = False
    if isinstance(elim, str) and elim == 'sigma_strict':
        Ke, Ge, Me = effective_matrices_sigma_strict_schur(
            theta, R, a, k,
            gamma=gamma, Q=Q, V0=V0,
            H=(0.0 if H is None else H),
            theta_dot=(0.0 if theta_dot is None else theta_dot),
            R_dot=(0.0 if R_dot is None else R_dot),
            Hdot=(0.0 if Hdot is None else Hdot),
        )
        Kd = _to_ndarray(Ke)
        block = _to_ndarray(Ge) * ((k * k) / (a * a) if a != 0.0 else 0.0) + _to_ndarray(Me)
        elim_kwargs = None
        prepared = True
    elif isinstance(elim, dict):
        elim_kwargs = elim
    else:
        elim_kwargs = None

    """Eigenvalues of K^{-1} [ G (k^2/a^2) + M ] ≡ ω^2 (finite-k dispersion)."""
    if elim_kwargs is None and not prepared:
        K = kinetic_matrix(
            theta,
            R,
            a,
            k,
            gamma,
            Q,
            V0,
            epsK_theta_R=epsK_theta_R,
            epsK_theta_phi=epsK_theta_phi,
            epsK_R_phi=epsK_R_phi,
        )
        G = gradient_matrix(
            theta,
            R,
            a,
            k,
            gamma,
            Q,
            V0,
            ctheta2=ctheta2,
            cR2=cR2,
            cphi2=cphi2,
            epsG_theta_R=epsG_theta_R,
            epsG_theta_phi=epsG_theta_phi,
            epsG_R_phi=epsG_R_phi,
            alpha_k2_theta=alpha_k2_theta,
            alpha_k2_R=alpha_k2_R,
            alpha_k2_phi=alpha_k2_phi
        )
        if gravM:
            if any(v is None for v in (H, theta_dot, R_dot, Hdot)):
                raise ValueError("gravM=True requires H, theta_dot, R_dot, Hdot")
            M = mass_matrix_sigma_total(theta, R, a, gamma=gamma, Q=Q, V0=V0, H=H, theta_dot=theta_dot, R_dot=R_dot, Hdot=Hdot)
        else:
            M = mass_matrix(theta, R, a, k, gamma, Q, V0)
        Kd = _to_ndarray(K)
        block = _to_ndarray(G) * ((k * k) / (a * a) if a != 0.0 else 0.0) + _to_ndarray(M)
    elif elim_kwargs is not None:
        Ke, Ge, Me = effective_matrices_sym(
            theta,
            R,
            a,
            k,
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
            **elim_kwargs,
        )
        Kd = _to_ndarray(Ke)
        block = _to_ndarray(Ge) * ((k * k) / (a * a) if a != 0.0 else 0.0) + _to_ndarray(Me)
    Kd, _, block_aug = _augment_with_sigma(
        Kd,
        np.zeros_like(Kd),
        block,
        include_sigma=include_sigma,
        csigma2=csigma2,
        msigma2=msigma2,
    )
    block = block_aug if block_aug is not None else block
    eigs = _gen_eigs(Kd, block)
    return eigs


def stability_finite_k_adm(
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
    alpha_k2_theta: float = 0.0,
    alpha_k2_R: float = 0.0,
    alpha_k2_phi: float = 0.0,
    elim: dict | str | bool | None = None,
    gravM: bool = False,
    H: float | None = None,
    theta_dot: float | None = None,
    R_dot: float | None = None,
    Hdot: float | None = None,
    include_sigma: bool = False,
    csigma2: float = 1.0,
    msigma2: float = 0.0,
):
    # Interpret elim flag: None → no elimination; True/'default' → default couplings
    if elim is True or (isinstance(elim, str) and str(elim).lower() == 'default'):
        elim_kwargs = default_elimination_couplings(gamma)
    elif isinstance(elim, dict):
        elim_kwargs = elim
    else:
        elim_kwargs = None

    """Finite-k stability: K ≻ 0 and ω^2 > 0 for all modes at given k."""
    if elim_kwargs is None:
        K = kinetic_matrix(
            theta,
            R,
            a,
            k,
            gamma,
            Q,
            V0,
            epsK_theta_R=epsK_theta_R,
            epsK_theta_phi=epsK_theta_phi,
            epsK_R_phi=epsK_R_phi,
        )
        Kd = _to_ndarray(K)
        w2 = dispersion_eigs_adm(
            theta,
            R,
            a,
            k,
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
            gravM=gravM,
            H=H,
            theta_dot=theta_dot,
            R_dot=R_dot,
            Hdot=Hdot,
            include_sigma=include_sigma,
            csigma2=csigma2,
            msigma2=msigma2,
        )
    else:
        Ke, _, _ = effective_matrices_sym(
            theta,
            R,
            a,
            k,
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
            **elim_kwargs,
        )
        Kd = _to_ndarray(Ke)
        w2 = dispersion_eigs_adm(
            theta,
            R,
            a,
            k,
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
            elim=elim_kwargs,
            include_sigma=include_sigma,
            csigma2=csigma2,
            msigma2=msigma2,
        )
    if include_sigma:
        Kd, _, _ = _augment_with_sigma(
            Kd,
            np.zeros_like(Kd),
            None,
            include_sigma=True,
            csigma2=csigma2,
            msigma2=msigma2,
        )
    _, no_ghost = _chol_pd(Kd)
    # Allow tiny negatives from numerical round-off
    no_tachyon = all(val > -1e-12 for val in w2)
    return {"no_ghost": no_ghost, "no_tachyon": no_tachyon, "w2_eigs": w2}
