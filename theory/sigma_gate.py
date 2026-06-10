# -*- coding: utf-8 -*-
"""
Exploratory legacy gate for the old sigma-layer.

Kept only as a compatibility helper for historical scans and tests.
"""
from __future__ import annotations

from typing import Iterable
import math


def gate_sigma(sigma: float) -> float:
    """G(sigma) = sigma^2 / (1 + sigma^2) in [0, 1)."""
    s2 = float(sigma) * float(sigma)
    return s2 / (1.0 + s2)


def i_grad_1d(
    theta_profile: Iterable[float],
    psi_profile: Iterable[float],
    *,
    dr: float,
    l_theta: float = 1.0,
    l_psi: float = 1.0,
) -> float:
    """
    Simple 1D proxy for I_grad = l_theta^2 |grad theta|^2 + l_psi^2 |grad psi|^2.

    Uses central differences on a radial/1D grid with spacing dr.
    """
    th = list(theta_profile)
    ps = list(psi_profile)
    n = min(len(th), len(ps))
    if n < 3:
        return 0.0
    if dr <= 0.0:
        raise ValueError("dr must be > 0")
    lt2 = float(l_theta) * float(l_theta)
    lp2 = float(l_psi) * float(l_psi)
    s = 0.0
    for i in range(1, n - 1):
        dth = (th[i + 1] - th[i - 1]) / (2.0 * dr)
        dps = (ps[i + 1] - ps[i - 1]) / (2.0 * dr)
        s += lt2 * dth * dth + lp2 * dps * dps
    return s / float(n - 2)


def sigma_eom_source(i_grad: float, lam: float) -> float:
    """Source term in sigma equation: + lam * I_grad."""
    return float(lam) * float(i_grad)
