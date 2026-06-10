# -*- coding: utf-8 -*-
"""
Compact-object likelihood helpers for finite no-horizon observables.
"""
from __future__ import annotations

from dataclasses import dataclass
import math
from typing import Dict

from fitting.core_api import BackgroundParams, Omegas, integrate_background_lnN
from checks.compact_report import compute_compact_metrics


@dataclass
class CompactTargets:
    # mode="gaussian": symmetric point-target fit for each metric.
    # mode="constraints": one-sided soft barriers (min / max thresholds).
    mode: str = "gaussian"
    z_surface: float | None = None
    sigma_z_surface: float | None = None
    delay_proxy: float | None = None
    sigma_delay_proxy: float | None = None
    n_peak: float | None = None
    sigma_n_peak: float | None = None
    z_surface_min: float | None = None
    sigma_z_surface_min: float | None = None
    delay_proxy_min: float | None = None
    sigma_delay_proxy_min: float | None = None
    n_peak_max: float | None = None
    sigma_n_peak_max: float | None = None


def compact_metrics_model(
    *,
    gamma: float,
    V0: float,
    Q: float = 0.0,
    omegas: Omegas | None = None,
    theta_scale: float = 1.0,
    a0: float = 1e-3,
    x0: float = 0.1,
    xdot0: float = 0.0,
    psi0: float = 0.05,
    psidot0: float = 0.0,
    dN: float = 0.02,
) -> Dict[str, float]:
    if omegas is None:
        omegas = Omegas()
    p = BackgroundParams(gamma=gamma, V0=V0, Q=Q, omegas=omegas)
    out = integrate_background_lnN(
        a0, x0, xdot0, psi0, psidot0, N1=0.0, p=p, hN=dN
    )
    return compute_compact_metrics(out["theta"], theta_scale=theta_scale, n_obs=1.0)


def loglike_compact(metrics: Dict[str, float], targets: CompactTargets) -> float:
    if targets.mode == "gaussian":
        if (
            targets.z_surface is None
            or targets.sigma_z_surface is None
            or targets.delay_proxy is None
            or targets.sigma_delay_proxy is None
            or targets.n_peak is None
            or targets.sigma_n_peak is None
        ):
            raise ValueError("gaussian compact targets require z/delay/n_peak with sigmas")
        terms = (
            ((metrics["z_surface"] - targets.z_surface) / targets.sigma_z_surface) ** 2,
            ((metrics["delay_proxy"] - targets.delay_proxy) / targets.sigma_delay_proxy) ** 2,
            ((metrics["n_peak"] - targets.n_peak) / targets.sigma_n_peak) ** 2,
        )
        return -0.5 * sum(terms)

    if targets.mode == "constraints":
        chi2 = 0.0
        if (
            targets.z_surface_min is not None
            and targets.sigma_z_surface_min is not None
            and metrics["z_surface"] < targets.z_surface_min
        ):
            chi2 += ((targets.z_surface_min - metrics["z_surface"]) / targets.sigma_z_surface_min) ** 2
        if (
            targets.delay_proxy_min is not None
            and targets.sigma_delay_proxy_min is not None
            and metrics["delay_proxy"] < targets.delay_proxy_min
        ):
            chi2 += ((targets.delay_proxy_min - metrics["delay_proxy"]) / targets.sigma_delay_proxy_min) ** 2
        if (
            targets.n_peak_max is not None
            and targets.sigma_n_peak_max is not None
            and metrics["n_peak"] > targets.n_peak_max
        ):
            chi2 += ((metrics["n_peak"] - targets.n_peak_max) / targets.sigma_n_peak_max) ** 2
        return -0.5 * chi2

    raise ValueError(f"Unknown compact target mode: {targets.mode}")


def default_compact_targets(
    *,
    omegas: Omegas | None = None,
    theta_scale: float = 1.0,
) -> CompactTargets:
    """
    Build a broad demo target around a fiducial point to enable combined runs
    without external compact-object datasets.
    """
    m = compact_metrics_model(
        gamma=1.0,
        V0=0.1,
        Q=0.0,
        omegas=omegas,
        theta_scale=theta_scale,
    )
    # Weak one-sided constraints: compact object should look very dark
    # (large redshift and delay), while avoiding hard point targets.
    return CompactTargets(
        mode="constraints",
        z_surface_min=max(1.0, 0.5 * abs(m["z_surface"])),
        sigma_z_surface_min=max(1.0, 0.5 * abs(m["z_surface"])),
        delay_proxy_min=max(1.0, 0.5 * abs(m["delay_proxy"])),
        sigma_delay_proxy_min=max(1.0, 0.5 * abs(m["delay_proxy"])),
        n_peak_max=None,
        sigma_n_peak_max=None,
    )
