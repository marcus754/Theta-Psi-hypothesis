# -*- coding: utf-8 -*-
"""
Event-level forward likelihood from explicit transfer.

Empirical mode is restricted to ring observables. Legacy echo-like proxies are
not accepted here as event-level data.
"""
from __future__ import annotations

from typing import Dict, List, Tuple

from fitting.forward_transfer import (
    SourceSpec,
    ring_diameter_uas_from_profile,
)


def evaluate_forward_events(
    profile: Dict[str, List[float] | float | bool],
    events: List[dict],
    *,
    ring_barrier_r_code: float = 3.0,
    ring_band_index: int = 2,
    ring_n_vac_tol: float = 0.01,
    ring_r_min_code: float = 2.0,
    echo_barrier_r_code: float = 3.0,
    echo_log_enhancement_coeff: float = 0.6366197723675814,
) -> List[Tuple[str, float, float, float]]:
    """
    Returns list of tuples: (name, predicted, target, sigma)

    Legacy `echo_*` kwargs are accepted for compatibility with older callers
    but ignored in empirical ring-only mode.
    """
    out = []
    for e in events:
        obs = str(e["observable"])
        if obs == "ring_uas":
            pred = ring_diameter_uas_from_profile(
                profile,
                source=SourceSpec(mass_msun=float(e["mass_msun"]), distance_m=float(e["distance_m"])),
                ring_barrier_r_code=ring_barrier_r_code,
                ring_band_index=ring_band_index,
                ring_n_vac_tol=ring_n_vac_tol,
                ring_r_min_code=ring_r_min_code,
            )
        else:
            raise ValueError(
                f"Unsupported event observable in empirical mode: {obs}. "
                "Allowed observable: ring_uas"
            )
        out.append((str(e["name"]), float(pred), float(e["target"]), float(e["sigma"])))
    return out


def loglike_forward_events(
    profile: Dict[str, List[float] | float | bool],
    events: List[dict],
    **kwargs,
) -> float:
    chi2 = 0.0
    for _, pred, target, sigma in evaluate_forward_events(profile, events, **kwargs):
        chi2 += ((pred - target) / sigma) ** 2
    return -0.5 * chi2
