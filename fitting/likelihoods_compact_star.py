# -*- coding: utf-8 -*-
"""
Likelihood for stationary compact-star strong-field observables.
"""
from __future__ import annotations

from dataclasses import dataclass
from typing import Dict

from fitting.core_api import CompactStarParams, solve_compact_star_profile


@dataclass
class CompactStarTargets:
    z_surface: float
    sigma_z_surface: float
    delay_proxy: float
    sigma_delay_proxy: float


def compact_star_metrics_model(
    *,
    theta_c: float,
    m2: float,
    lam: float,
    r_max: float = 20.0,
    dr: float = 0.02,
    theta_scale: float = 1.0,
    psi_grad_weight: float = 0.5,
    theta_psi_coupling: float = 0.25,
    refractive_strength: float = 0.5,
    refractive_j_star: float = 1.0,
    kappa_refr: float = 1.0,
) -> Dict[str, float]:
    p = CompactStarParams(
        theta_c=theta_c,
        m2=m2,
        lam=lam,
        r_max=r_max,
        dr=dr,
        theta_scale=theta_scale,
        psi_grad_weight=psi_grad_weight,
        theta_psi_coupling=theta_psi_coupling,
        refractive_strength=refractive_strength,
        refractive_j_star=refractive_j_star,
        kappa_refr=kappa_refr,
    )
    out = solve_compact_star_profile(p)
    return {
        "z_surface": float(out["z_surface"]),
        "delay_proxy": float(out["delay_proxy"]),
        "n_peak": float(out["n_peak"]),
    }


def loglike_compact_star(metrics: Dict[str, float], targets: CompactStarTargets) -> float:
    chi2 = (
        ((metrics["z_surface"] - targets.z_surface) / targets.sigma_z_surface) ** 2
        + ((metrics["delay_proxy"] - targets.delay_proxy) / targets.sigma_delay_proxy) ** 2
    )
    return -0.5 * chi2


def default_compact_star_targets(
    *,
    theta_c: float = 1.5,
    m2: float = 1.0,
    lam: float = 0.1,
    r_max: float = 20.0,
    dr: float = 0.02,
    theta_scale: float = 1.0,
    psi_grad_weight: float = 0.5,
    theta_psi_coupling: float = 0.25,
    refractive_strength: float = 0.5,
    refractive_j_star: float = 1.0,
    kappa_refr: float = 1.0,
) -> CompactStarTargets:
    m = compact_star_metrics_model(
        theta_c=theta_c,
        m2=m2,
        lam=lam,
        r_max=r_max,
        dr=dr,
        theta_scale=theta_scale,
        psi_grad_weight=psi_grad_weight,
        theta_psi_coupling=theta_psi_coupling,
        refractive_strength=refractive_strength,
        refractive_j_star=refractive_j_star,
        kappa_refr=kappa_refr,
    )
    return CompactStarTargets(
        z_surface=m["z_surface"],
        sigma_z_surface=max(1.0, 0.5 * abs(m["z_surface"])),
        delay_proxy=m["delay_proxy"],
        sigma_delay_proxy=max(1.0, 0.5 * abs(m["delay_proxy"])),
    )
