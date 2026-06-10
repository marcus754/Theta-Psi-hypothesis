# -*- coding: utf-8 -*-
"""
Minimal static vacuum profile for the reduced time field.

This module provides a first local-gravity layer that is closer to the reduced
Theta-Psi core than the compact-star closure.  The working variable is

    chi = sqrt(3) * ln(theta),

for which the homogeneous minisuperspace kinetic term becomes canonical.  In a
static, spherically-symmetric exterior region with no local source terms, the
minimal reduced equation is taken to be Laplace:

    chi'' + (2/r) chi' = 0.

Its regular exterior solution is

    chi(r) = chi_inf + source_strength / r,

which yields a monotone theta(r) and, via n(theta), finite redshift and delay
for the diagnostic optical-index profile.

This is deliberately only an exterior vacuum scaffold: it does not yet solve
the full sourced local field equations for matter.
"""
from __future__ import annotations

from dataclasses import dataclass
from typing import Dict, List
import math

from theory.optical_metric import diagnostic_refractive_index_from_theta


SQRT3 = math.sqrt(3.0)


@dataclass
class StaticVacuumParams:
    theta_inf: float = 1.0
    source_strength: float = 1.0
    r_min: float = 1.0
    r_max: float = 100.0
    dr: float = 0.1
    theta_scale: float = 1.0
    kappa_wf: float = 2.0


def chi_of_r(r: float, *, theta_inf: float, source_strength: float) -> float:
    if r <= 0.0:
        raise ValueError("r must be > 0")
    if theta_inf <= 0.0:
        raise ValueError("theta_inf must be > 0")
    chi_inf = SQRT3 * math.log(theta_inf)
    return chi_inf + float(source_strength) / float(r)


def theta_of_r(r: float, *, theta_inf: float, source_strength: float) -> float:
    return math.exp(chi_of_r(r, theta_inf=theta_inf, source_strength=source_strength) / SQRT3)


def build_static_vacuum_profile(p: StaticVacuumParams) -> Dict[str, List[float] | float | bool]:
    if p.r_min <= 0.0:
        raise ValueError("r_min must be > 0")
    if p.r_max <= p.r_min:
        raise ValueError("r_max must be > r_min")
    if p.dr <= 0.0:
        raise ValueError("dr must be > 0")

    nsteps = max(2, int(math.ceil((p.r_max - p.r_min) / p.dr)))
    r: List[float] = [p.r_min + i * p.dr for i in range(nsteps + 1)]
    if r[-1] < p.r_max:
        r.append(p.r_max)

    theta = [
        theta_of_r(ri, theta_inf=p.theta_inf, source_strength=p.source_strength)
        for ri in r
    ]
    chi = [SQRT3 * math.log(max(th, 1e-300)) for th in theta]
    nvals = [
        diagnostic_refractive_index_from_theta(
            th,
            theta_scale=p.theta_scale,
            kappa_wf=p.kappa_wf,
            profile="asinh",
        )
        for th in theta
    ]

    # Exterior optical excess delay relative to n=1 vacuum.
    delay = 0.0
    for i in range(1, len(r)):
        delay += 0.5 * ((nvals[i - 1] - 1.0) + (nvals[i] - 1.0)) * (r[i] - r[i - 1])

    z_edge = (nvals[0] / max(nvals[-1], 1e-30)) - 1.0
    monotone = all(theta[i] >= theta[i + 1] for i in range(len(theta) - 1))

    return {
        "r": r,
        "chi": chi,
        "theta": theta,
        "n": nvals,
        "z_edge": z_edge,
        "delay_proxy": delay,
        "monotone_theta": monotone,
        "ok_no_horizon": True,
    }
