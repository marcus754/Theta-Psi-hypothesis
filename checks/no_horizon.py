# -*- coding: utf-8 -*-
"""No-horizon criteria for the strong-field branch of the Θ–Ψ model."""
from __future__ import annotations

import math
from typing import Iterable, Dict

from theory.optical_metric import canonical_refractive_index_from_theta


def no_horizon_from_n_values(n_values: Iterable[float], *, eps: float = 1e-12) -> bool:
    """
    Return True if all sampled n are finite, positive and bounded away from zero.
    """
    vals = list(n_values)
    if not vals:
        return False
    for nv in vals:
        if not math.isfinite(nv):
            return False
        if nv <= eps:
            return False
    return True


def no_horizon_from_theta(
    theta_values: Iterable[float],
    *,
    theta_scale: float = 1.0,
) -> Dict[str, float | bool]:
    """
    Evaluate no-horizon condition via the optical index n(theta).

    The physical statement is that the branch stays finite and does not reach
    a true horizon.
    """
    n_vals = [
        canonical_refractive_index_from_theta(th, theta_scale=theta_scale)
        for th in theta_values
    ]
    ok = no_horizon_from_n_values(n_vals)
    return {
        "ok_no_horizon": ok,
        "n_min": min(n_vals) if n_vals else float("nan"),
        "n_peak": max(n_vals) if n_vals else float("nan"),
        "n_peak_sample": max(n_vals) if n_vals else float("nan"),
    }
