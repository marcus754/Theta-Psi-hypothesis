# -*- coding: utf-8 -*-
"""
Atomic-scale weak-field consistency checks for the Θ–Ψ setup.

This module does not try to derive atomic spectra. It only checks that the
local, weak-field limit does not inject any spurious refractive or PPN shift
at atomic scales. That is the correct consistency requirement for the
relativistic chemistry behind the color of gold and similar effects.
"""
from __future__ import annotations

from dataclasses import dataclass

from checks.ppn_full import compute_ppn_params
from theory.optical_metric import canonical_refractive_index_from_theta


@dataclass(frozen=True)
class AtomicRelativisticConsistency:
    theta_eff: float
    theta_scale: float
    n_value: float
    delta_n: float
    gamma_ppn: float
    beta_ppn: float
    ok_n: bool
    ok_ppn: bool
    ok: bool


def atomic_relativistic_consistency(
    theta_eff: float = 1e-30,
    *,
    theta_scale: float = 1.0,
    n_tol: float = 1e-12,
    ppn_tol: float = 1e-15,
) -> AtomicRelativisticConsistency:
    """
    Check that the atomic/chemical regime stays in the standard local limit.

    The Θ–Ψ sector is allowed to influence strong-field propagation, but at
    atomic scales it must reduce to:
    - unit optical response: n -> 1,
    - local PPN baseline: gamma = beta = 1.
    """
    n_value = canonical_refractive_index_from_theta(theta_eff, theta_scale=theta_scale)
    ppn = compute_ppn_params(gamma_model=1.0, theta0=max(theta_scale, 1e-30))
    delta_n = abs(n_value - 1.0)
    ok_n = delta_n <= n_tol
    ok_ppn = abs(ppn.gamma - 1.0) <= ppn_tol and abs(ppn.beta - 1.0) <= ppn_tol
    return AtomicRelativisticConsistency(
        theta_eff=float(theta_eff),
        theta_scale=float(theta_scale),
        n_value=float(n_value),
        delta_n=float(delta_n),
        gamma_ppn=float(ppn.gamma),
        beta_ppn=float(ppn.beta),
        ok_n=ok_n,
        ok_ppn=ok_ppn,
        ok=ok_n and ok_ppn,
    )
