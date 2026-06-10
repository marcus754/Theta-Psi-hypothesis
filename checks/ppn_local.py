# -*- coding: utf-8 -*-
"""
Simple local PPN-style veto.

The effective coupling ε ≃ 1 / (γ² θ² + 3) should stay below a small number
so that local tests of gravity remain satisfied.  This is *not* a rigorous PPN
calculation, only a quick veto to avoid obviously large couplings during
parameter scans.
"""
from __future__ import annotations


def ppn_proxy_ok(theta: float, gamma: float, eps_max: float = 1e-2) -> bool:
    try:
        eps = 1.0 / (gamma * gamma * theta * theta + 3.0)
    except ZeroDivisionError:
        return False
    return abs(eps) < eps_max
