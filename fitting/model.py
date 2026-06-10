# -*- coding: utf-8 -*-
"""
Model predictions for background observables.

Computes H(z) in dimensionless units by integrating the FRW background with
fields using the ln a integrator, then sampling the supplied z grid.
"""
from __future__ import annotations

from bisect import bisect_left
from typing import List, Tuple
import math

from fitting.core_api import BackgroundParams as BGParams, Omegas, integrate_background_lnN


def Hz_dimless(
    z_points: List[float],
    *,
    gamma: float,
    V0: float,
    Q: float = 0.0,
    omegas: Omegas | None = None,
    a_min: float = 1e-3,
    x0: float = 0.1,
    xdot0: float = 0.0,
    R0: float = 0.05,
    Rdot0: float = 0.0,
    dN: float = 0.01,
    **_legacy_kwargs,
) -> List[float]:
    """
    Return H(z) at requested redshifts in dimensionless units.

    Integrates background in segments from a_min to each target scale factor.
    Assumes monotonically increasing a; z_points may be unsorted (will be
    evaluated in ascending a internally and re-ordered back).
    """
    if omegas is None:
        omegas = Omegas()
    if not z_points:
        return []
    # Integrate once on a fixed ln(a) grid, then interpolate H(a) for all z.
    p = BGParams(gamma=gamma, V0=V0, Q=Q, omegas=omegas)
    amin_eff = max(1e-8, min(a_min, min(1.0 / (1.0 + z) for z in z_points)))
    track = integrate_background_lnN(
        amin_eff,
        x0,
        xdot0,
        R0,
        Rdot0,
        N1=0.0,
        p=p,
        hN=dN,
    )
    a_grid = track["a"]
    h_grid = track["H"]

    def h_at_a(a_t: float) -> float:
        if a_t <= a_grid[0]:
            return float(h_grid[0])
        if a_t >= a_grid[-1]:
            return float(h_grid[-1])
        j = bisect_left(a_grid, a_t)
        if j <= 0:
            return float(h_grid[0])
        if j >= len(a_grid):
            return float(h_grid[-1])
        a0 = a_grid[j - 1]
        a1 = a_grid[j]
        h0 = h_grid[j - 1]
        h1 = h_grid[j]
        t = (a_t - a0) / max(a1 - a0, 1e-18)
        return float(h0 + t * (h1 - h0))

    return [h_at_a(1.0 / (1.0 + z)) for z in z_points]


def _cumtrapz(y: List[float], x: List[float]) -> List[float]:
    out = [0.0] * len(x)
    s = 0.0
    for i in range(1, len(x)):
        dx = x[i] - x[i - 1]
        s += 0.5 * (y[i] + y[i - 1]) * dx
        out[i] = s
    return out


def lumdist_dimless(
    z_points: List[float],
    *,
    gamma: float,
    V0: float,
    Q: float = 0.0,
    omegas: Omegas | None = None,
    a_min: float = 1e-3,
    dN: float = 0.01,
    nint: int = 200,
    **_legacy_kwargs,
) -> List[float]:
    """
    Luminosity distance in dimensionless units for flat FRW: D_L = (1+z) ∫_0^z dz'/H(z').
    Uses trapezoidal rule over `nint` points per target z.
    """
    if omegas is None:
        omegas = Omegas()
    if not z_points:
        return []
    zmax = max(z_points)
    n = max(nint, 1)
    zs = [zmax * i / n for i in range(n + 1)]
    Hs = Hz_dimless(
        zs,
        gamma=gamma,
        V0=V0,
        Q=Q,
        omegas=omegas,
        a_min=a_min,
        dN=dN,
    )
    invH = [1.0 / max(h, 1e-18) for h in Hs]
    dc_grid = _cumtrapz(invH, zs)

    def dc_at_z(z_t: float) -> float:
        if z_t <= zs[0]:
            return dc_grid[0]
        if z_t >= zs[-1]:
            return dc_grid[-1]
        j = bisect_left(zs, z_t)
        if j <= 0:
            return dc_grid[0]
        z0 = zs[j - 1]
        z1 = zs[j]
        d0 = dc_grid[j - 1]
        d1 = dc_grid[j]
        t = (z_t - z0) / max(z1 - z0, 1e-18)
        return d0 + t * (d1 - d0)

    return [(1.0 + zt) * dc_at_z(zt) for zt in z_points]


def mu_dimless(
    z_points: List[float],
    *,
    gamma: float,
    V0: float,
    Q: float = 0.0,
    omegas: Omegas | None = None,
    a_min: float = 1e-3,
    dN: float = 0.01,
    nint: int = 200,
    mu0: float = 0.0,
    **_legacy_kwargs,
) -> List[float]:
    """
    Dimensionless distance-modulus helper: μ = μ0 + 5 log10 D_L.
    μ0 absorbs unit conversion; here default μ0=0 for demo datasets.
    """
    Dl = lumdist_dimless(
        z_points,
        gamma=gamma,
        V0=V0,
        Q=Q,
        omegas=omegas,
        a_min=a_min,
        dN=dN,
        nint=nint,
    )
    return [mu0 + 5.0 * math.log10(max(d, 1e-18)) for d in Dl]


def DV_dimless(
    z_points: List[float],
    *,
    gamma: float,
    V0: float,
    Q: float = 0.0,
    omegas: Omegas | None = None,
    a_min: float = 1e-3,
    dN: float = 0.01,
    nint: int = 200,
    **_legacy_kwargs,
) -> List[float]:
    """
    BAO volume distance: D_V(z) = [ (1+z)^2 D_A(z)^2 z/H(z) ]^{1/3}, with c=1.
    D_A = D_C/(1+z), D_C = ∫ dz/H.
    """
    if omegas is None:
        omegas = Omegas()
    Dl = lumdist_dimless(
        z_points,
        gamma=gamma,
        V0=V0,
        Q=Q,
        omegas=omegas,
        a_min=a_min,
        dN=dN,
        nint=nint,
    )
    # D_A = D_L/(1+z)^2; D_C = (1+z)D_A = D_L/(1+z)
    Dc = [dl/max(1.0+z,1e-18) for dl,z in zip(Dl, z_points)]
    Hz = Hz_dimless(
        z_points,
        gamma=gamma,
        V0=V0,
        Q=Q,
        omegas=omegas,
        a_min=a_min,
        dN=dN,
    )
    out: List[float] = []
    for z, dc, h in zip(z_points, Dc, Hz):
        val = (( (1.0+z)*(1.0+z) * (dc/(1.0+z))**2 * z/max(h,1e-18) )) ** (1.0/3.0)
        out.append(val)
    return out
