# -*- coding: utf-8 -*-
"""
Unified health evaluation along a background trajectory.

Aggregates the existing vetoes and stability checks into a single function
suited for scans and reports. Designed to work with outputs produced by
frw_background.integrate_background[_rk45] and integrate_background_lnN.
"""
from __future__ import annotations

from typing import Dict, Any, Iterable, Callable, Optional

from .cutoff import cutoff_ok as _cutoff_ok
from .ppn_local import ppn_proxy_ok as _ppn_proxy_ok
from .early_universe import bbn_safe as _bbn_safe
from .superluminal import (
    has_superluminal as _has_superluminal_2x2,
    has_superluminal_adm as _has_superluminal_adm,
    finite_k_stable_along_traj as _finite_k_stable_along_traj,
)
from src.linear_modes import stability as _stability_2x2
from .no_horizon import no_horizon_from_theta as _no_horizon_from_theta


def _interp_linear(x: Iterable[float], y: Iterable[float]) -> Callable[[float], float]:
    xs = list(x); ys = list(y)
    assert len(xs) == len(ys) and len(xs) >= 2
    # Ensure ascending by x
    pairs = sorted(zip(xs, ys), key=lambda p: p[0])
    xs = [p[0] for p in pairs]; ys = [p[1] for p in pairs]

    def f(xq: float) -> float:
        if xq <= xs[0]:
            return ys[0]
        if xq >= xs[-1]:
            return ys[-1]
        # binary search
        lo, hi = 0, len(xs) - 1
        while hi - lo > 1:
            mid = (lo + hi) // 2
            if xs[mid] <= xq:
                lo = mid
            else:
                hi = mid
        t = (xq - xs[lo]) / (xs[hi] - xs[lo])
        return (1.0 - t) * ys[lo] + t * ys[hi]

    return f


def evaluate_health(
    traj: Dict[str, Iterable[float]],
    *,
    gamma: float,
    V0: float,
    Q: float = 0.0,
    # Linear modes baseline
    ctheta2: float = 1.0,
    cpsi2: float = 1.0,
    # ADM options
    use_adm: bool = False,
    cR2: float = 1.0,
    cphi2: float = 1.0,
    epsKtr: float = 0.0,
    epsKtphi: float = 0.0,
    epsKRphi: float = 0.0,
    epsGtr: float = 0.0,
    epsGtphi: float = 0.0,
    epsGRphi: float = 0.0,
    alpha_k2_theta: float = 0.0,
    alpha_k2_R: float = 0.0,
    alpha_k2_phi: float = 0.0,
    finite_k: bool = False,
    gravM: bool = False,
    kmin: float = 1e-2,
    kmax: float = 1.0,
    nk: int = 5,
    # Cutoff and local constraints
    cutoff: float = 1e-1,
    ppn_eps_max: float = 1e-2,
    # First-principles mode: ignore heuristic vetoes (cutoff, local PPN veto)
    strict_first_principles: bool = True,
    # Early universe safety
    check_bbn: bool = False,
    # Optical-index finite-branch check (optional)
    check_n_bounds: bool = False,
    theta_scale: float = 1.0,
    check_no_horizon: bool = False,
) -> Dict[str, Any]:
    """
    Return dict of flags and summary diagnostics.

    Expected keys in traj: 'a','H','theta','theta_dot','psi'.
    """
    theta = list(traj.get("theta", []))
    theta_dot = list(traj.get("theta_dot", []))
    a = list(traj.get("a", []))
    H = list(traj.get("H", []))
    R = list(traj.get("psi", []))

    # Linear stability (2×2 scaffold) along trajectory
    ok_lin = True
    for th in theta:
        st = _stability_2x2(th, gamma, ctheta2, cpsi2)
        if not (st["no_ghost"] and st["no_gradient_instability"]):
            ok_lin = False
            break

    # Superluminality
    ok_sub = True
    if use_adm:
        for idx, (th, r_i, a_i) in enumerate(zip(theta, R, a)):
            if _has_superluminal_adm(
                th, r_i, a_i, 1.0, gamma, Q, V0,
                ctheta2=ctheta2, cR2=cR2, cphi2=cphi2,
                epsK_theta_R=epsKtr, epsK_theta_phi=epsKtphi, epsK_R_phi=epsKRphi,
                epsG_theta_R=epsGtr, epsG_theta_phi=epsGtphi, epsG_R_phi=epsGRphi,
                alpha_k2_theta=alpha_k2_theta, alpha_k2_R=alpha_k2_R, alpha_k2_phi=alpha_k2_phi,
                gravM=gravM,
                H_list=[list(traj.get("H", []))[idx]] if len(list(traj.get("H", [])))>idx else None,
                theta_dot_list=[list(traj.get("theta_dot", []))[idx]] if len(list(traj.get("theta_dot", [])))>idx else None,
                R_dot_list=[list(traj.get("psi_dot", []))[idx]] if len(list(traj.get("psi_dot", [])))>idx else None,
                Hdot_list=[list(traj.get("Hdot", []))[idx]] if len(list(traj.get("Hdot", [])))>idx else None,
                adm_strict=True,
            ):
                ok_sub = False
                break
        # In ADM mode, enforce finite-k stability by default unless explicitly disabled
        if use_adm and not finite_k:
            finite_k = True
        if ok_sub and finite_k:
            ok_sub = _finite_k_stable_along_traj(
                theta, R, a, gamma, Q, V0,
                ctheta2=ctheta2, cR2=cR2, cphi2=cphi2,
                kmin=kmin, kmax=kmax, nk=nk,
                epsK_theta_R=epsKtr, epsK_theta_phi=epsKtphi, epsK_R_phi=epsKRphi,
                epsG_theta_R=epsGtr, epsG_theta_phi=epsGtphi, epsG_R_phi=epsGRphi,
                alpha_k2_theta=alpha_k2_theta, alpha_k2_R=alpha_k2_R, alpha_k2_phi=alpha_k2_phi,
                gravM=gravM,
                H_list=list(traj.get("H", [])),
                theta_dot_list=list(traj.get("theta_dot", [])),
                R_dot_list=list(traj.get("psi_dot", [])),
                Hdot_list=list(traj.get("Hdot", [])),
                adm_strict=True,
            )
    else:
        for th, Hi in zip(theta, H):
            if _has_superluminal_2x2(th, Hi, gamma, V0, ctheta2, cpsi2):
                ok_sub = False
                break

    # EFT cutoff and PPN: optionally bypass in first-principles mode
    if strict_first_principles:
        ok_cutoff = True
        ok_ppn = True
    else:
        ok_cutoff = _cutoff_ok(theta, theta_dot, gamma, cutoff)
        ok_ppn = all(_ppn_proxy_ok(th, gamma, eps_max=ppn_eps_max) for th in theta)

    # BBN monotonicity (optional; interpolate H(a) from trajectory)
    ok_bbn = True
    if check_bbn and len(a) >= 2:
        H_of_a = _interp_linear(a, H)
        ok_bbn = _bbn_safe(H_of_a)

    ok_n_bounds = True
    if check_n_bounds:
        ok_n_bounds = _no_horizon_from_theta(theta, theta_scale=theta_scale)["ok_no_horizon"]

    ok_no_horizon = True
    if check_no_horizon:
        ok_no_horizon = bool(_no_horizon_from_theta(theta, theta_scale=theta_scale)["ok_no_horizon"])

    ok = ok_lin and ok_sub and ok_cutoff and ok_ppn and ok_bbn and ok_n_bounds and ok_no_horizon
    return {
        "ok": ok,
        "ok_linear": ok_lin,
        "ok_sub": ok_sub,
        "ok_cutoff": ok_cutoff,
        "ok_ppn": ok_ppn,
        "ok_bbn": ok_bbn,
        "ok_n_bounds": ok_n_bounds,
        "ok_no_horizon": ok_no_horizon,
    }
