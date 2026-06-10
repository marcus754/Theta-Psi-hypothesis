# -*- coding: utf-8 -*-
"""
Likelihood for forward strong-field observables (EHT/GW proxies).
"""
from __future__ import annotations

from dataclasses import dataclass
from typing import Dict


@dataclass
class ForwardTargets:
    ring_m87_uas: float
    sigma_ring_m87_uas: float
    ring_sgra_uas: float
    sigma_ring_sgra_uas: float
    echo_delay_ms: float | None = None
    sigma_echo_delay_ms: float | None = None


def loglike_forward(obs: Dict[str, float], t: ForwardTargets) -> float:
    chi2 = (
        ((obs["ring_m87_uas"] - t.ring_m87_uas) / t.sigma_ring_m87_uas) ** 2
        + ((obs["ring_sgra_uas"] - t.ring_sgra_uas) / t.sigma_ring_sgra_uas) ** 2
    )
    if t.echo_delay_ms is not None and t.sigma_echo_delay_ms is not None and "echo_delay_ms" in obs:
        chi2 += ((obs["echo_delay_ms"] - t.echo_delay_ms) / t.sigma_echo_delay_ms) ** 2
    return -0.5 * chi2
