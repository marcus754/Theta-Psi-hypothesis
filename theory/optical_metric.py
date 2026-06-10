# -*- coding: utf-8 -*-
"""
Optical metric utilities for the Θ–Ψ hypothesis.

Baseline: a single universal index n(θ) induces an effective (optical) metric
    g̃_{μν} = η_{μν} + (n(θ)^2 − 1) u_μ u_ν,
where u_μ is a unit timelike 4-vector constructed from θ (local clock) and
η_{μν} = diag(1, -1, -1, -1) is the Minkowski metric in local inertial frame.

This module exposes:
- optical_metric(n, u): construct g̃_{μν} from a scalar n and timelike u.
- normalize_timelike(u): return u normalized w.r.t. η so that η(u,u)=+1.
- tensor_speed_optical_metric(): return c_T^2=1 as the low-energy EH response.

See docs/14_optical_metric_emergent.md and docs/12_tensor_speed.md.
"""
from __future__ import annotations

import math
import warnings
from typing import Iterable, Sequence, Tuple


ETA = ( (1.0, 0.0, 0.0, 0.0),
        (0.0,-1.0, 0.0, 0.0),
        (0.0, 0.0,-1.0, 0.0),
        (0.0, 0.0, 0.0,-1.0) )


def _eta_contract(u: Sequence[float], v: Sequence[float]) -> float:
    return (
        ETA[0][0]*u[0]*v[0] + ETA[1][1]*u[1]*v[1]
        + ETA[2][2]*u[2]*v[2] + ETA[3][3]*u[3]*v[3]
    )


def normalize_timelike(u: Sequence[float]) -> Tuple[float, float, float, float]:
    """Normalize u so that η(u,u)=+1 (timelike unit) if possible."""
    s = _eta_contract(u, u)
    if s <= 0.0:
        raise ValueError("Vector is not timelike w.r.t. η")
    scale = (s) ** 0.5
    return (u[0]/scale, u[1]/scale, u[2]/scale, u[3]/scale)


def optical_metric(n: float, u: Sequence[float]) -> Tuple[Tuple[float, ...], ...]:
    """Return g̃_{μν} = η_{μν} + (n^2-1) u_μ u_ν.

    u is treated as covariant components in the local inertial frame.
    """
    a = float(n)*float(n) - 1.0
    uu = (
        (u[0]*u[0], u[0]*u[1], u[0]*u[2], u[0]*u[3]),
        (u[1]*u[0], u[1]*u[1], u[1]*u[2], u[1]*u[3]),
        (u[2]*u[0], u[2]*u[1], u[2]*u[2], u[2]*u[3]),
        (u[3]*u[0], u[3]*u[1], u[3]*u[2], u[3]*u[3]),
    )
    return tuple(
        tuple(ETA[i][j] + a*uu[i][j] for j in range(4)) for i in range(4)
    )


def _phi_eff_from_theta(theta: float, theta_scale: float) -> float:
    """
    Smooth nonnegative weak-field proxy built from theta.

    Canonical definition: Phi_eff = ln(theta / theta_scale).
    Numerical regularization ensures it is smooth and vanishes at theta = theta_scale.
    Interpreted as |Phi_eff| for matching n >= 1.
    """
    eps = 1e-12
    # Normalizing theta by theta_scale so that x=1 is the background (vacuum).
    x = float(theta) / max(float(theta_scale), eps)
    # The potential Phi_eff is the log-deviation from the background.
    val = math.log(max(x, eps))
    # We return the absolute value to ensure n >= 1 even for negative fluctuations.
    return math.sqrt(val * val + eps * eps) - eps


CANONICAL_WEAK_FIELD_KAPPA = 2.0
CANONICAL_REFRACTIVE_PROFILE = "asinh"
DIAGNOSTIC_REFRACTIVE_PROFILES = ("linear", "exp", "tanh2")


def _refractive_index_from_theta_profile(
    theta: float,
    *,
    theta_scale: float = 1.0,
    kappa_wf: float = 2.0,
    profile: str = CANONICAL_REFRACTIVE_PROFILE,
) -> float:
    """Internal profile family. Noncanonical profiles are diagnostics only."""
    if theta_scale <= 0.0:
        raise ValueError("theta_scale must be > 0")
    phi_eff = _phi_eff_from_theta(theta, theta_scale)
    prof = str(profile).strip().lower()
    if prof == "linear":
        return 1.0 + float(kappa_wf) * phi_eff
    if prof == "asinh":
        return 1.0 + float(kappa_wf) * math.asinh(phi_eff)
    if prof == "exp":
        return 1.0 + float(kappa_wf) * (1.0 - math.exp(-phi_eff))
    if prof == "tanh2":
        s = math.tanh(float(theta) / float(theta_scale))
        return 1.0 + float(kappa_wf) * (s * s)
    raise ValueError(f"Unknown profile for refractive_index_from_theta: {profile}")


def canonical_refractive_index_from_theta(
    theta: float,
    *,
    theta_scale: float = 1.0,
) -> float:
    """Canonical weak-field index fixed before data-facing fits."""
    return _refractive_index_from_theta_profile(
        theta,
        theta_scale=theta_scale,
        kappa_wf=CANONICAL_WEAK_FIELD_KAPPA,
        profile=CANONICAL_REFRACTIVE_PROFILE,
    )


def diagnostic_refractive_index_from_theta(
    theta: float,
    *,
    theta_scale: float = 1.0,
    kappa_wf: float = 2.0,
    profile: str = CANONICAL_REFRACTIVE_PROFILE,
) -> float:
    """Diagnostic index family for scans; not a canonical prediction."""
    return _refractive_index_from_theta_profile(
        theta,
        theta_scale=theta_scale,
        kappa_wf=kappa_wf,
        profile=profile,
    )


def refractive_index_from_theta(
    theta: float,
    *,
    theta_scale: float = 1.0,
    kappa_wf: float = CANONICAL_WEAK_FIELD_KAPPA,
    profile: str | None = None,
) -> float:
    """Backward-compatible wrapper around the canonical weak-field index.

    Passing a profile or noncanonical kappa is diagnostic-only and emits a
    warning so data-facing code cannot silently treat it as canonical.
    """
    if profile is None and float(kappa_wf) == CANONICAL_WEAK_FIELD_KAPPA:
        return canonical_refractive_index_from_theta(theta, theta_scale=theta_scale)
    prof = CANONICAL_REFRACTIVE_PROFILE if profile is None else str(profile)
    warnings.warn(
        "profile/kappa_wf overrides are diagnostic-only; use "
        "canonical_refractive_index_from_theta for canonical predictions",
        DeprecationWarning,
        stacklevel=2,
    )
    return diagnostic_refractive_index_from_theta(
        theta,
        theta_scale=theta_scale,
        kappa_wf=kappa_wf,
        profile=prof,
    )


def refractive_index_from_refractive_invariant(
    j_refr: float,
    *,
    kappa_refr: float = 1.0,
) -> float:
    """Refractive index from J_refr without a physical ceiling parameter."""
    j = max(0.0, float(j_refr))
    return 1.0 + float(kappa_refr) * math.asinh(j)


def phi_eff_from_theta_eff(theta_eff: float, theta_scale: float = 1.0) -> float:
    """Public helper for the reduced weak-field proxy Φ_eff(θ_eff)."""
    if theta_scale <= 0.0:
        raise ValueError("theta_scale must be > 0")
    return _phi_eff_from_theta(theta_eff, theta_scale)


def redshift_factor(
    theta_emit: float,
    theta_obs: float,
    *,
    theta_scale: float = 1.0,
) -> float:
    """Return 1+z = n_emit / n_obs for the finite n(theta) profile."""
    n_emit = refractive_index_from_theta(theta_emit, theta_scale=theta_scale)
    n_obs = refractive_index_from_theta(theta_obs, theta_scale=theta_scale)
    return n_emit / n_obs


def redshift_factor_theta_eff(
    theta_eff_emit: float,
    theta_eff_obs: float,
    *,
    theta_scale: float = 1.0,
) -> float:
    """Alias for redshift_factor with explicit θ_eff naming."""
    return redshift_factor(
        theta_eff_emit,
        theta_eff_obs,
        theta_scale=theta_scale,
    )


def optical_time_of_flight(
    dl: Iterable[float],
    theta_path: Iterable[float],
    *,
    theta_scale: float = 1.0,
) -> float:
    """Compute t = Σ n(theta_i) * dl_i for a discretized path."""
    dls = list(dl)
    thetas = list(theta_path)
    if len(dls) != len(thetas):
        raise ValueError("dl and theta_path must have equal length")
    return sum(
        refractive_index_from_theta(th, theta_scale=theta_scale) * dli
        for dli, th in zip(dls, thetas)
    )


def optical_time_of_flight_theta_eff(
    dl: Iterable[float],
    theta_eff_path: Iterable[float],
    *,
    theta_scale: float = 1.0,
) -> float:
    """Alias for optical_time_of_flight with explicit θ_eff naming."""
    return optical_time_of_flight(
        dl,
        theta_eff_path,
        theta_scale=theta_scale,
    )


def n_finite_ok(
    theta_values: Iterable[float],
    *,
    theta_scale: float = 1.0,
) -> bool:
    """Validate that n(theta) stays finite and above unity on a set of values."""
    for th in theta_values:
        nval = refractive_index_from_theta(th, theta_scale=theta_scale)
        if not math.isfinite(nval) or nval < 1.0:
            return False
    return True


def n_finite_ok_theta_eff(
    theta_eff_values: Iterable[float],
    *,
    theta_scale: float = 1.0,
) -> bool:
    """Alias for n_finite_ok with explicit θ_eff naming."""
    return n_finite_ok(theta_eff_values, theta_scale=theta_scale)


def tensor_speed_optical_metric() -> float:
    """Low-energy EH response on g̃ yields luminal tensor waves: c_T^2 = 1."""
    return 1.0
