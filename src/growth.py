# -*- coding: utf-8 -*-
"""
Linear growth factor D(a) in the quasi-static sub-horizon approximation.

Equation in N = ln a:
    D'' + (2 + d ln H / d ln a) D' - (3/2) Ω_m(a) μ(a) D = 0

with μ(a)=1 in GR. This module uses the background Omegas to compute H(a),
E2(a) and Ω_m(a), and integrates the ODE over N.
"""
from __future__ import annotations

from typing import List, Callable, Tuple
import math

from background.components import Omegas, E2


def _dlnH_dln_a(a: float, om: Omegas) -> float:
    e2 = E2(a, om)
    # d ln H / d ln a = (1/2) d ln E2 / d ln a
    # d ln E2 / d ln a = (a / E2) dE2/da
    dE2_da = -4.0 * om.Omega_r * (a ** -5) - 3.0 * om.Omega_m * (a ** -4)
    return 0.5 * (a / e2) * dE2_da


def _Omega_m_of_a(a: float, om: Omegas) -> float:
    e2 = E2(a, om)
    return (om.Omega_m * (a ** -3)) / e2


def growth_factor(
    a_min: float,
    a_max: float,
    *,
    omegas: Omegas,
    nsteps: int = 200,
    mu_of_a: Callable[[float], float] | None = None,
) -> Tuple[List[float], List[float]]:
    """
    Return (a_grid, D_norm) with D(a=1)=1 normalization.

    Integrates the growth ODE in N=ln a with simple RK4.
    """
    assert 0.0 < a_min < a_max
    mu = mu_of_a or (lambda a: 1.0)
    N0 = math.log(a_min)
    N1 = math.log(a_max)
    h = (N1 - N0) / max(nsteps, 1)

    # State: (D, D') as functions of N
    # Start deep in matter era: use growing-mode-like initial condition
    D = a_min  # proportional to a at early times (rough proxy)
    Dp = D     # dD/dN ≈ D
    a = a_min

    a_list: List[float] = []
    D_list: List[float] = []

    def rhs(a: float, D: float, Dp: float) -> Tuple[float, float]:
        dlnH = _dlnH_dln_a(a, omegas)
        Om = _Omega_m_of_a(a, omegas)
        Dpp = - (2.0 + dlnH) * Dp + 1.5 * Om * mu(a) * D
        return Dp, Dpp

    def rk4_step(a: float, D: float, Dp: float, h: float):
        # Evolve in N; update a multiplicatively
        k1_D, k1_Dp = rhs(a, D, Dp)
        a2 = math.exp(math.log(a) + 0.5 * h)
        k2_D, k2_Dp = rhs(a2, D + 0.5 * h * k1_D, Dp + 0.5 * h * k1_Dp)
        k3_D, k3_Dp = rhs(a2, D + 0.5 * h * k2_D, Dp + 0.5 * h * k2_Dp)
        a4 = math.exp(math.log(a) + h)
        k4_D, k4_Dp = rhs(a4, D + h * k3_D, Dp + h * k3_Dp)
        Dn = D + (h / 6.0) * (k1_D + 2*k2_D + 2*k3_D + k4_D)
        Dpn = Dp + (h / 6.0) * (k1_Dp + 2*k2_Dp + 2*k3_Dp + k4_Dp)
        an = a4
        return an, Dn, Dpn

    for _ in range(nsteps + 1):
        a_list.append(a)
        D_list.append(D)
        a, D, Dp = rk4_step(a, D, Dp, h)

    # Normalize so that D(a=1) = 1 if a=1 is in the grid; otherwise use last
    if a_list[-1] != 1.0 and (a_list[0] - 1.0) * (a_list[-1] - 1.0) < 0:
        # Interpolate to a=1
        # Find bracketing indices
        for i in range(len(a_list) - 1):
            if a_list[i] <= 1.0 <= a_list[i+1]:
                t = (1.0 - a_list[i]) / (a_list[i+1] - a_list[i])
                D1 = (1.0 - t) * D_list[i] + t * D_list[i+1]
                scale = D1 if D1 != 0.0 else 1.0
                break
    else:
        scale = D_list[-1] if D_list[-1] != 0.0 else 1.0

    D_norm = [Di/scale for Di in D_list]
    return a_list, D_norm

