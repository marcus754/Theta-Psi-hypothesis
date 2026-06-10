# -*- coding: utf-8 -*-
"""
Stage-2 explicit transfer: compact-star profile -> strong-field observables.

Empirical mode uses ring observables. Echo-like delay helpers are retained only
as compatibility diagnostics and are not part of the falsifiability contour.
"""
from __future__ import annotations

from dataclasses import dataclass
from typing import Dict, List
import math


UAS_PER_RAD = 206264806247.09634
G = 6.67430e-11
C = 299792458.0
M_SUN = 1.98847e30


@dataclass
class SourceSpec:
    mass_msun: float
    distance_m: float


M87_SPEC = SourceSpec(mass_msun=6.5e9, distance_m=16.8e6 * 3.085677581491367e16)
SGRA_SPEC = SourceSpec(mass_msun=4.297e6, distance_m=8.277e3 * 3.085677581491367e16)


def _trapz(y: List[float], x: List[float]) -> float:
    s = 0.0
    for i in range(1, len(x)):
        s += 0.5 * (y[i - 1] + y[i]) * (x[i] - x[i - 1])
    return s


def _interp_profile_n(profile: Dict[str, List[float] | float | bool], rq: float) -> float:
    r = [float(v) for v in profile["r"]]  # type: ignore[index]
    n = [float(v) for v in profile["n"]]  # type: ignore[index]
    if not r or not n or len(r) != len(n):
        return 1.0
    if rq <= r[0]:
        return n[0]
    for i in range(1, len(r)):
        if rq <= r[i]:
            x0, x1 = r[i - 1], r[i]
            y0, y1 = n[i - 1], n[i]
            t = (rq - x0) / max(x1 - x0, 1e-18)
            return (1.0 - t) * y0 + t * y1
    return n[-1]


def _local_minima_indices(y: List[float]) -> List[int]:
    idx: List[int] = []
    for i in range(1, len(y) - 1):
        if y[i] <= y[i - 1] and y[i] < y[i + 1]:
            idx.append(i)
    return idx


def critical_impact_parameter_code(
    profile: Dict[str, List[float] | float | bool],
    *,
    ring_barrier_r_code: float = 3.0,
    ring_band_index: int = 2,
    ring_n_vac_tol: float = 0.01,
    ring_r_min_code: float = 2.0,
) -> float:
    r = [float(v) for v in profile["r"]]  # type: ignore[index]
    n = [float(v) for v in profile["n"]]  # type: ignore[index]
    if not r or not n or len(r) != len(n):
        nbar = max(1.0, _interp_profile_n(profile, ring_barrier_r_code))
        return float(ring_barrier_r_code) * nbar

    b = [ri * ni for ri, ni in zip(r, n)]
    mins = _local_minima_indices(b)
    cand = [
        i
        for i in mins
        if r[i] >= float(ring_r_min_code) and n[i] <= 1.0 + float(ring_n_vac_tol)
    ]
    if cand:
        k = max(1, int(ring_band_index))
        j = cand[min(k - 1, len(cand) - 1)]
        return b[j]

    nbar = max(1.0, _interp_profile_n(profile, ring_barrier_r_code))
    return float(ring_barrier_r_code) * nbar


def ring_diameter_uas_from_profile(
    profile: Dict[str, List[float] | float | bool],
    *,
    source: SourceSpec,
    ring_barrier_r_code: float = 3.0,
    ring_band_index: int = 2,
    ring_n_vac_tol: float = 0.01,
    ring_r_min_code: float = 2.0,
) -> float:
    """
    Angular ring diameter from critical impact parameter.

    Assumption: code r is in units r_g = GM/c^2.
    """
    b_code = max(
        0.0,
        critical_impact_parameter_code(
            profile,
            ring_barrier_r_code=ring_barrier_r_code,
            ring_band_index=ring_band_index,
            ring_n_vac_tol=ring_n_vac_tol,
            ring_r_min_code=ring_r_min_code,
        ),
    )
    rg = G * (source.mass_msun * M_SUN) / (C * C)
    b_m = b_code * rg
    theta_rad = 2.0 * b_m / max(source.distance_m, 1e-30)
    return theta_rad * UAS_PER_RAD


def echo_delay_ms_from_profile(
    profile: Dict[str, List[float] | float | bool],
    *,
    mass_msun: float = 60.0,
    barrier_r_code: float = 3.0,
    include_log_enhancement: bool = True,
    log_enhancement_coeff: float = 2.0 / math.pi,
) -> float:
    """
    Echo delay scale from explicit optical travel time.

    tau_code ~ 2 * ∫_0^{r_barrier} n(r) dr  (in units r_g/c),
    tau_ms = tau_code * (GM/c^3) * 1e3.
    """
    r = [float(v) for v in profile["r"]]  # type: ignore[index]
    n = [float(v) for v in profile["n"]]  # type: ignore[index]
    if not r or not n or len(r) != len(n):
        return 0.0

    rr: List[float] = []
    nn: List[float] = []
    for ri, ni in zip(r, n):
        if ri <= barrier_r_code:
            rr.append(ri)
            nn.append(ni)
        else:
            break
    if len(rr) < 2:
        rr = r[:2]
        nn = n[:2]
    tau_code = 2.0 * _trapz(nn, rr)
    if include_log_enhancement:
        z = max(0.0, float(profile.get("z_surface", 0.0)))
        tau_code *= (1.0 + float(log_enhancement_coeff) * math.log1p(z))
    t_g = G * (mass_msun * M_SUN) / (C ** 3)
    return tau_code * t_g * 1e3


def explicit_forward_observables(
    profile: Dict[str, List[float] | float | bool],
    *,
    m87: SourceSpec = M87_SPEC,
    sgra: SourceSpec = SGRA_SPEC,
    ring_barrier_r_code: float = 3.0,
    ring_band_index: int = 2,
    ring_n_vac_tol: float = 0.01,
    ring_r_min_code: float = 2.0,
    echo_mass_msun: float = 60.0,
    echo_barrier_r_code: float = 3.0,
    echo_log_enhancement_coeff: float = 2.0 / math.pi,
) -> Dict[str, float]:
    return {
        "ring_m87_uas": ring_diameter_uas_from_profile(
            profile,
            source=m87,
            ring_barrier_r_code=ring_barrier_r_code,
            ring_band_index=ring_band_index,
            ring_n_vac_tol=ring_n_vac_tol,
            ring_r_min_code=ring_r_min_code,
        ),
        "ring_sgra_uas": ring_diameter_uas_from_profile(
            profile,
            source=sgra,
            ring_barrier_r_code=ring_barrier_r_code,
            ring_band_index=ring_band_index,
            ring_n_vac_tol=ring_n_vac_tol,
            ring_r_min_code=ring_r_min_code,
        ),
        "echo_delay_ms": echo_delay_ms_from_profile(
            profile,
            mass_msun=echo_mass_msun,
            barrier_r_code=echo_barrier_r_code,
            include_log_enhancement=True,
            log_enhancement_coeff=echo_log_enhancement_coeff,
        ),
    }
