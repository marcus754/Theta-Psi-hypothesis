# -*- coding: utf-8 -*-
"""
Forward observables layer for the compact-star strong-field sector.

This module maps model-internal compact metrics (z_surface, delay, n_peak)
to observational observables (EHT ring size, GW echo delay scale).
"""
from __future__ import annotations

from typing import Dict


def calibrate_v1_coefficients(
    *,
    reference_metrics: Dict[str, float],
    target_ring_m87_uas: float,
    target_ring_sgra_uas: float,
    target_echo_delay_ms: float | None = None,
) -> Dict[str, float]:
    """
    Stage-1 calibration: infer k_* from a reference compact-star point.
    """
    z = max(0.0, float(reference_metrics["z_surface"]))
    d = max(1e-18, float(reference_metrics["delay_proxy"]))
    root = (1.0 + z) ** 0.5
    return {
        "k_ring_m87": float(target_ring_m87_uas) / max(root, 1e-18),
        "k_ring_sgra": float(target_ring_sgra_uas) / max(root, 1e-18),
        "k_echo_delay": (float(target_echo_delay_ms) / d) if target_echo_delay_ms is not None else 1.0,
    }


def forward_observables_from_compact_metrics(
    metrics: Dict[str, float],
    *,
    k_ring_m87: float = 3.0,
    k_ring_sgra: float = 3.6,
    k_echo_delay: float = 1.0,
) -> Dict[str, float]:
    """
    Map compact metrics -> observational observables.

    Mapping (explicit, version-v1):
      ring_m87_uas  = k_ring_m87 * sqrt(1 + z_surface)
      ring_sgra_uas = k_ring_sgra * sqrt(1 + z_surface)
      echo_delay_ms = k_echo_delay * delay
    """
    z = max(0.0, float(metrics["z_surface"]))
    d = max(0.0, float(metrics["delay_proxy"]))
    root = (1.0 + z) ** 0.5
    return {
        "ring_m87_uas": float(k_ring_m87) * root,
        "ring_sgra_uas": float(k_ring_sgra) * root,
        "echo_delay_ms": float(k_echo_delay) * d,
    }
