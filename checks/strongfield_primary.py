# -*- coding: utf-8 -*-
"""
Primary strong-field summary for the Θ–Ψ compact branch.

The point of this module is to keep the theory-facing strong-field claim
minimal: finite surface redshift, finite delay, and no true horizon on the
physical strong-field branch.
"""
from __future__ import annotations

from dataclasses import dataclass

from checks.compact_star import CompactStarParams, solve_compact_star_profile


@dataclass(frozen=True)
class StrongFieldPrimarySummary:
    theta_c: float
    m2: float
    lam: float
    n_peak: float
    z_surface: float
    delay_proxy: float
    ok_no_horizon: bool
    classification: str


def summarize_strongfield_primary(
    *,
    theta_c: float,
    m2: float = 1.0,
    lam: float = 0.1,
    r_max: float = 5.0,
    dr: float = 0.02,
    theta_scale: float = 1.0,
    psi_grad_weight: float = 1.0,
    theta_psi_coupling: float = 1.0,
    refractive_strength: float = 1.0,
    refractive_j_star: float = 1.0,
    kappa_refr: float = 1.0,
) -> StrongFieldPrimarySummary:
    """
    Return the compact-branch strong-field summary used as the primary check.
    """
    params = CompactStarParams(
        theta_c=float(theta_c),
        m2=float(m2),
        lam=float(lam),
        r_max=float(r_max),
        dr=float(dr),
        theta_scale=float(theta_scale),
        psi_grad_weight=float(psi_grad_weight),
        theta_psi_coupling=float(theta_psi_coupling),
        refractive_strength=float(refractive_strength),
        refractive_j_star=float(refractive_j_star),
        kappa_refr=float(kappa_refr),
    )
    prof = solve_compact_star_profile(params)
    ok_no_horizon = bool(prof["ok_no_horizon"])
    return StrongFieldPrimarySummary(
        theta_c=float(theta_c),
        m2=float(m2),
        lam=float(lam),
        n_peak=float(prof["n_peak"]),
        z_surface=float(prof["z_surface"]),
        delay_proxy=float(prof["delay_proxy"]),
        ok_no_horizon=ok_no_horizon,
        classification="compact_no_horizon" if ok_no_horizon else "pathological_or_horizon",
    )
