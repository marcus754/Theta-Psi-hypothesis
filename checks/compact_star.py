# -*- coding: utf-8 -*-
"""
Stationary spherically-symmetric compact-star closure for the Θ–Ψ branch.

We solve a minimal variational radial system for theta(r) and psi(r) built from

    L_stat = 1/2 * I_grad - U(theta, psi) + alpha * F(J_refr),
    I_grad = theta'^2 + w_psi * psi'^2,
    J_refr = Phi_eff(theta) * I_grad,

with the stationary potential

    U(theta, psi)
      = 1/2*m2*theta^2 + 1/4*lam*theta^4
      - 1/2*m2*psi^2 + 1/2*c_tp*lam*theta^2*psi^2.

The Euler-Lagrange equations are written in quasilinear 2x2 form and solved
locally for (theta'', psi''). This keeps the public API stable while moving
the strong-field closure away from the old one-field n(theta) ansatz.
"""
from __future__ import annotations

from dataclasses import dataclass
import math
from typing import Dict, List, Iterable

from theory.optical_metric import refractive_index_from_refractive_invariant, phi_eff_from_theta_eff
from theory.refractive_sector import (
    f_v1,
    f_v1_with_scale,
    f_v1_prime_with_scale,
    f_v1_double_prime_with_scale,
)


@dataclass
class CompactStarParams:
    theta_c: float = 2.0
    m2: float = 1.0
    lam: float = 0.1
    r_max: float = 20.0
    dr: float = 0.02
    theta_scale: float = 1.0
    psi_core_scale: float = 1.0
    psi_grad_weight: float = 1.0
    theta_psi_coupling: float = 1.0
    refractive_strength: float = 1.0
    refractive_j_star: float = 1.0
    kappa_refr: float = 1.0


def _phi_eff_prime(theta: float, theta_scale: float) -> float:
    eps = 1e-12
    scale = max(float(theta_scale), eps)
    x = float(theta) / scale
    denom = math.sqrt(x * x + eps * eps)
    return (x / max(denom, eps)) / scale


def _potential_partials(theta: float, psi: float, m2: float, lam: float, theta_psi_coupling: float) -> tuple[float, float]:
    """Return dU/dtheta and dU/dpsi for the stationary two-field potential."""
    u_theta = (
        m2 * theta
        + lam * theta * theta * theta
        + theta_psi_coupling * lam * theta * psi * psi
    )
    u_psi = -m2 * psi + theta_psi_coupling * lam * theta * theta * psi
    return u_theta, u_psi


def _rhs_coupled(
    r: float,
    theta: float,
    theta_prime: float,
    psi: float,
    psi_prime: float,
    m2: float,
    lam: float,
    theta_psi_coupling: float,
    psi_grad_weight: float,
    refractive_strength: float,
    refractive_j_star: float,
    theta_scale: float,
) -> tuple[float, float, float, float]:
    # theta' = theta_prime, psi' = psi_prime
    rr = r if r > 1e-8 else 1e-8
    dtheta = theta_prime
    dpsi = psi_prime
    wpsi = max(float(psi_grad_weight), 1e-12)
    alpha = float(refractive_strength)
    i_grad = max(0.0, theta_prime * theta_prime + wpsi * psi_prime * psi_prime)
    phi_eff = phi_eff_from_theta_eff(theta, theta_scale=theta_scale)
    j_refr = phi_eff * i_grad
    fprime = f_v1_prime_with_scale(j_refr, j_star=refractive_j_star)
    fdouble = f_v1_double_prime_with_scale(j_refr, j_star=refractive_j_star)
    phi_prime = _phi_eff_prime(theta, theta_scale)
    a = max(1.0 + 2.0 * alpha * phi_eff * fprime, 1e-8)
    a_theta = 2.0 * alpha * phi_prime * (fprime + phi_eff * i_grad * fdouble)
    a_theta_prime = 4.0 * alpha * phi_eff * phi_eff * fdouble * theta_prime
    a_psi_prime = 4.0 * alpha * phi_eff * phi_eff * fdouble * wpsi * psi_prime
    u_theta, u_psi = _potential_partials(theta, psi, m2, lam, theta_psi_coupling)

    m11 = a + a_theta_prime * theta_prime
    m12 = a_psi_prime * theta_prime
    m21 = a_theta_prime * wpsi * psi_prime
    m22 = wpsi * (a + a_psi_prime * psi_prime)

    rhs1 = -(a_theta * theta_prime * theta_prime + 2.0 * a * theta_prime / rr + u_theta - alpha * fprime * phi_prime * i_grad)
    rhs2 = -(a_theta * wpsi * theta_prime * psi_prime + 2.0 * a * wpsi * psi_prime / rr + u_psi)

    det = m11 * m22 - m12 * m21
    if abs(det) < 1e-10:
        # Fallback to the diagonal weak-coupling limit if the local system degenerates.
        dtheta_prime = rhs1 / max(m11, 1e-8)
        dpsi_prime = rhs2 / max(m22, 1e-8)
    else:
        dtheta_prime = (rhs1 * m22 - rhs2 * m12) / det
        dpsi_prime = (m11 * rhs2 - m21 * rhs1) / det
    return dtheta, dtheta_prime, dpsi, dpsi_prime

def solve_compact_star_profile(p: CompactStarParams) -> Dict[str, List[float] | float | bool]:
    nsteps = max(2, int(p.r_max / p.dr))
    r: List[float] = [0.0]
    theta: List[float] = [p.theta_c]
    theta_prime: List[float] = [0.0]
    psi: List[float] = [max(float(p.psi_core_scale) * max(abs(float(p.theta_c)), float(p.theta_scale)), 1e-9)]
    psi_prime: List[float] = [0.0]

    for i in range(nsteps):
        ri = r[-1]
        ti = theta[-1]
        thi = theta_prime[-1]
        pi = psi[-1]
        ppi = psi_prime[-1]
        h = p.dr

        k1_t, k1_thp, k1_p, k1_pp = _rhs_coupled(
            ri, ti, thi, pi, ppi, p.m2, p.lam, p.theta_psi_coupling,
            p.psi_grad_weight, p.refractive_strength, p.refractive_j_star, p.theta_scale,
        )
        k2_t, k2_thp, k2_p, k2_pp = _rhs_coupled(
            ri + 0.5 * h,
            ti + 0.5 * h * k1_t,
            thi + 0.5 * h * k1_thp,
            pi + 0.5 * h * k1_p,
            ppi + 0.5 * h * k1_pp,
            p.m2,
            p.lam,
            p.theta_psi_coupling,
            p.psi_grad_weight,
            p.refractive_strength,
            p.refractive_j_star,
            p.theta_scale,
        )
        k3_t, k3_thp, k3_p, k3_pp = _rhs_coupled(
            ri + 0.5 * h,
            ti + 0.5 * h * k2_t,
            thi + 0.5 * h * k2_thp,
            pi + 0.5 * h * k2_p,
            ppi + 0.5 * h * k2_pp,
            p.m2,
            p.lam,
            p.theta_psi_coupling,
            p.psi_grad_weight,
            p.refractive_strength,
            p.refractive_j_star,
            p.theta_scale,
        )
        k4_t, k4_thp, k4_p, k4_pp = _rhs_coupled(
            ri + h,
            ti + h * k3_t,
            thi + h * k3_thp,
            pi + h * k3_p,
            ppi + h * k3_pp,
            p.m2,
            p.lam,
            p.theta_psi_coupling,
            p.psi_grad_weight,
            p.refractive_strength,
            p.refractive_j_star,
            p.theta_scale,
        )

        t_next = ti + (h / 6.0) * (k1_t + 2.0 * k2_t + 2.0 * k3_t + k4_t)
        thp_next = thi + (h / 6.0) * (k1_thp + 2.0 * k2_thp + 2.0 * k3_thp + k4_thp)
        p_next = pi + (h / 6.0) * (k1_p + 2.0 * k2_p + 2.0 * k3_p + k4_p)
        pp_next = ppi + (h / 6.0) * (k1_pp + 2.0 * k2_pp + 2.0 * k3_pp + k4_pp)
        r_next = ri + h

        r.append(r_next)
        theta.append(t_next)
        theta_prime.append(thp_next)
        psi.append(max(p_next, 0.0))
        psi_prime.append(pp_next)

    phi_eff = [phi_eff_from_theta_eff(ti, theta_scale=p.theta_scale) for ti in theta]
    i_grad = [
        max(0.0, vi * vi + float(p.psi_grad_weight) * wi * wi)
        for vi, wi in zip(theta_prime, psi_prime)
    ]
    j_refr = [ph * ig for ph, ig in zip(phi_eff, i_grad)]
    f_refr = [f_v1_with_scale(jv, j_star=p.refractive_j_star) for jv in j_refr]
    nvals = [
        refractive_index_from_refractive_invariant(jv, kappa_refr=p.kappa_refr)
        for jv in j_refr
    ]
    n_peak = max(nvals) if nvals else 1.0
    z_surface = n_peak - 1.0
    # Trapezoid integral of excess optical path
    delay = 0.0
    for i in range(1, len(r)):
        delay += 0.5 * ((nvals[i - 1] - 1.0) + (nvals[i] - 1.0)) * (r[i] - r[i - 1])

    return {
        "r": r,
        "theta": theta,
        "theta_prime": theta_prime,
        "psi": psi,
        "psi_prime": psi_prime,
        "phi_eff": phi_eff,
        "i_grad": i_grad,
        "j_refr": j_refr,
        "f_refr": f_refr,
        "n": nvals,
        "z_surface": z_surface,
        "delay_proxy": delay,
        "n_peak": n_peak,
        "ok_no_horizon": bool(math.isfinite(n_peak) and n_peak >= 1.0),
    }


def scan_compact_star_grid(
    theta_c_values: Iterable[float],
    m2_values: Iterable[float],
    lam_values: Iterable[float],
    *,
    r_max: float = 20.0,
    dr: float = 0.02,
    theta_scale: float = 1.0,
) -> List[dict]:
    rows: List[dict] = []
    for tc in theta_c_values:
        for m2 in m2_values:
            for lam in lam_values:
                p = CompactStarParams(
                    theta_c=float(tc),
                    m2=float(m2),
                    lam=float(lam),
                    r_max=r_max,
                    dr=dr,
                    theta_scale=theta_scale,
                )
                out = solve_compact_star_profile(p)
                rows.append(
                    {
                        "theta_c": tc,
                        "m2": m2,
                        "lam": lam,
                        "z_surface": float(out["z_surface"]),
                        "delay_proxy": float(out["delay_proxy"]),
                        "n_peak": float(out["n_peak"]),
                        "ok_no_horizon": bool(out["ok_no_horizon"]),
                    }
                )
    return rows
