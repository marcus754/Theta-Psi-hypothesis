# -*- coding: utf-8 -*-
"""Numerical continuity check for FRW background trajectories.

Checks residual of:
    d rho_tot / dt + 3 H (rho_tot + p_tot) = 0,
with rho_tot = 3 H^2 and p_tot reconstructed from Hdot.
"""
from __future__ import annotations

from dataclasses import dataclass
from typing import Dict, List
import math


@dataclass
class ContinuityReport:
    max_abs_residual: float
    max_rel_residual: float
    mean_rel_residual: float
    n_samples: int


def continuity_report_from_background(track: Dict[str, List[float]]) -> ContinuityReport:
    t = list(track.get("t", []))
    H = list(track.get("H", []))
    Hdot = list(track.get("Hdot", []))
    n = min(len(t), len(H), len(Hdot))
    if n < 3:
        return ContinuityReport(max_abs_residual=float("nan"), max_rel_residual=float("nan"), mean_rel_residual=float("nan"), n_samples=n)

    max_abs = 0.0
    max_rel = 0.0
    s_rel = 0.0
    m = 0
    for i in range(1, n - 1):
        dt = t[i + 1] - t[i - 1]
        if abs(dt) < 1e-18:
            continue
        rho_prev = 3.0 * H[i - 1] * H[i - 1]
        rho_next = 3.0 * H[i + 1] * H[i + 1]
        drho_dt = (rho_next - rho_prev) / dt

        rho_i = 3.0 * H[i] * H[i]
        # Raychaudhuri: Hdot = -1/2 (rho+p) => p = -2 Hdot - rho
        p_i = -2.0 * Hdot[i] - rho_i
        cont = drho_dt + 3.0 * H[i] * (rho_i + p_i)

        scale = abs(drho_dt) + abs(3.0 * H[i] * (rho_i + p_i)) + 1e-30
        rel = abs(cont) / scale
        if not (math.isfinite(cont) and math.isfinite(rel)):
            continue

        max_abs = max(max_abs, abs(cont))
        max_rel = max(max_rel, rel)
        s_rel += rel
        m += 1

    mean_rel = s_rel / max(m, 1)
    return ContinuityReport(max_abs_residual=max_abs, max_rel_residual=max_rel, mean_rel_residual=mean_rel, n_samples=m)
