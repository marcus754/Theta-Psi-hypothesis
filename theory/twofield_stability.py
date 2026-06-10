# -*- coding: utf-8 -*-
"""
Minimal stability helpers for the reduced two-field Θ–Ψ sector.

These helpers provide one compact place to check whether the current reduced
canon passes the health conditions needed before fit:

1. FRW homogeneous sector:
   K = diag(3/theta^2 + gamma^2, 1)
   G = diag(c_theta^2 K_theta, c_psi^2)

2. Stationary radial strong-field sector:
   local quasilinear 2x2 operator derived from
       L_stat = 1/2 I_grad - U(theta, psi) + alpha F(J_refr)
"""
from __future__ import annotations

import math
from typing import Dict

from theory.optical_metric import phi_eff_from_theta_eff
from theory.refractive_sector import f_v1_prime_with_scale, f_v1_double_prime_with_scale


def frw_stability_conditions(
    theta: float,
    gamma: float,
    *,
    ctheta2: float = 1.0,
    cpsi2: float = 1.0,
    mtheta2: float = 0.0,
    V0: float = 0.0,
) -> Dict[str, float | bool]:
    """Return minimal FRW no-ghost / no-gradient / no-tachyon diagnostics."""
    th = float(theta)
    if abs(th) < 1e-12:
        raise ValueError("theta too close to zero for reduced FRW stability check")
    k_theta = 3.0 / (th * th) + float(gamma) * float(gamma)
    k_psi = 1.0
    g_theta = float(ctheta2) * k_theta
    g_psi = float(cpsi2)
    no_ghost = (k_theta > 0.0) and (k_psi > 0.0)
    no_gradient_instability = (g_theta > 0.0) and (g_psi > 0.0)
    no_tachyon = (float(mtheta2) >= 0.0) and (float(V0) >= 0.0)
    return {
        "k_theta": k_theta,
        "k_psi": k_psi,
        "g_theta": g_theta,
        "g_psi": g_psi,
        "no_ghost": no_ghost,
        "no_gradient_instability": no_gradient_instability,
        "no_tachyon": no_tachyon,
        "no_extra_modes": bool(no_ghost and no_gradient_instability),
    }


def stationary_operator_coefficients(
    *,
    theta: float,
    theta_prime: float,
    psi_prime: float,
    psi_grad_weight: float = 1.0,
    refractive_strength: float = 1.0,
    refractive_j_star: float = 1.0,
    theta_scale: float = 1.0,
) -> Dict[str, float]:
    """Return local coefficients of the quasilinear radial operator."""
    wpsi = max(float(psi_grad_weight), 1e-12)
    alpha = float(refractive_strength)
    phi_eff = phi_eff_from_theta_eff(theta, theta_scale=theta_scale)
    eps = 1e-12
    scale = max(float(theta_scale), eps)
    x = float(theta) / scale
    phi_prime = (x / max(math.sqrt(x * x + eps * eps), eps)) / scale
    i_grad = max(0.0, float(theta_prime) ** 2 + wpsi * float(psi_prime) ** 2)
    j_refr = phi_eff * i_grad
    fp = f_v1_prime_with_scale(j_refr, j_star=refractive_j_star)
    fpp = f_v1_double_prime_with_scale(j_refr, j_star=refractive_j_star)
    a = 1.0 + 2.0 * alpha * phi_eff * fp
    a_theta = 2.0 * alpha * phi_prime * (fp + phi_eff * i_grad * fpp)
    a_theta_prime = 4.0 * alpha * phi_eff * phi_eff * fpp * float(theta_prime)
    a_psi_prime = 4.0 * alpha * phi_eff * phi_eff * fpp * wpsi * float(psi_prime)
    m11 = a + a_theta_prime * float(theta_prime)
    m12 = a_psi_prime * float(theta_prime)
    m21 = a_theta_prime * wpsi * float(psi_prime)
    m22 = wpsi * (a + a_psi_prime * float(psi_prime))
    det = m11 * m22 - m12 * m21
    return {
        "phi_eff": phi_eff,
        "phi_prime": phi_prime,
        "i_grad": i_grad,
        "j_refr": j_refr,
        "fprime": fp,
        "fdouble": fpp,
        "a": a,
        "a_theta": a_theta,
        "a_theta_prime": a_theta_prime,
        "a_psi_prime": a_psi_prime,
        "m11": m11,
        "m12": m12,
        "m21": m21,
        "m22": m22,
        "det": det,
        "wpsi": wpsi,
    }


def stationary_stability_conditions(**kwargs) -> Dict[str, float | bool]:
    """Return minimal ellipticity / positivity diagnostics for the radial sector."""
    d = stationary_operator_coefficients(**kwargs)
    d["positive_refr_gain"] = bool(d["a"] > 0.0)
    d["elliptic"] = bool(d["m11"] > 0.0 and d["m22"] > 0.0 and d["det"] > 0.0)
    d["no_extra_modes"] = bool(d["positive_refr_gain"] and d["elliptic"])
    return d


def stability_ledger(
    theta: float,
    gamma: float,
    *,
    theta_prime: float = 0.0,
    psi_prime: float = 0.0,
    ctheta2: float = 1.0,
    cpsi2: float = 1.0,
    mtheta2: float = 0.0,
    V0: float = 0.0,
    psi_grad_weight: float = 1.0,
    refractive_strength: float = 1.0,
    refractive_j_star: float = 1.0,
    theta_scale: float = 1.0,
) -> Dict[str, float | bool | Dict[str, float | bool]]:
    """Return a combined FRW/stationary health ledger.

    The ledger is intentionally explicit so the repo can check one boolean:
    the reduced sector is acceptable only if both background and stationary
    diagnostics are healthy.
    """
    frw = frw_stability_conditions(
        theta=theta,
        gamma=gamma,
        ctheta2=ctheta2,
        cpsi2=cpsi2,
        mtheta2=mtheta2,
        V0=V0,
    )
    stationary = stationary_stability_conditions(
        theta=theta,
        theta_prime=theta_prime,
        psi_prime=psi_prime,
        psi_grad_weight=psi_grad_weight,
        refractive_strength=refractive_strength,
        refractive_j_star=refractive_j_star,
        theta_scale=theta_scale,
    )
    no_ghost = bool(frw["no_ghost"])
    no_gradient = bool(frw["no_gradient_instability"])
    no_tachyon = bool(frw["no_tachyon"])
    no_extra_modes = bool(frw["no_extra_modes"] and stationary["no_extra_modes"])
    healthy = bool(no_ghost and no_gradient and no_tachyon and no_extra_modes)
    return {
        "frw": frw,
        "stationary": stationary,
        "no_ghost": no_ghost,
        "no_gradient_instability": no_gradient,
        "no_tachyon": no_tachyon,
        "no_extra_modes": no_extra_modes,
        "healthy": healthy,
    }
