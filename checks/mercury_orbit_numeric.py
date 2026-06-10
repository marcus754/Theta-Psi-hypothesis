# -*- coding: utf-8 -*-
"""
Numerical Mercury perihelion precession from the weak-field orbit equation.

We integrate the 1PN-like Binet equation written for the Θ–Ψ weak-field
effective metric (through PPN gamma/beta):

    u'' + u = 1/p + kappa * u^2,

where
    u = 1/r,
    p = a (1 - e^2),
    kappa = (2 + 2*gamma - beta) * GM/c^2.

For gamma=beta=1 this reduces to GR:
    u'' + u = 1/p + 3(GM/c^2) u^2.
"""
from __future__ import annotations

from dataclasses import dataclass
import math
from typing import List

from .mercury_precession import ARCSEC_PER_RAD, MercuryOrbit


@dataclass
class NumericPrecessionResult:
    delta_per_orbit_rad: float
    arcsec_per_century: float
    perihelion_count: int


def _rk4_step(phi: float, y1: float, y2: float, h: float, *, inv_p: float, kappa: float) -> tuple[float, float]:
    # y1=u, y2=du/dphi
    def f1(_phi: float, _y1: float, _y2: float) -> float:
        return _y2

    def f2(_phi: float, _y1: float, _y2: float) -> float:
        return inv_p + kappa * (_y1 * _y1) - _y1

    k1y1 = h * f1(phi, y1, y2)
    k1y2 = h * f2(phi, y1, y2)
    k2y1 = h * f1(phi + 0.5 * h, y1 + 0.5 * k1y1, y2 + 0.5 * k1y2)
    k2y2 = h * f2(phi + 0.5 * h, y1 + 0.5 * k1y1, y2 + 0.5 * k1y2)
    k3y1 = h * f1(phi + 0.5 * h, y1 + 0.5 * k2y1, y2 + 0.5 * k2y2)
    k3y2 = h * f2(phi + 0.5 * h, y1 + 0.5 * k2y1, y2 + 0.5 * k2y2)
    k4y1 = h * f1(phi + h, y1 + k3y1, y2 + k3y2)
    k4y2 = h * f2(phi + h, y1 + k3y1, y2 + k3y2)

    y1n = y1 + (k1y1 + 2.0 * k2y1 + 2.0 * k3y1 + k4y1) / 6.0
    y2n = y2 + (k1y2 + 2.0 * k2y2 + 2.0 * k3y2 + k4y2) / 6.0
    return y1n, y2n


def _find_local_max_positions(phi: List[float], u: List[float]) -> List[float]:
    out: List[float] = []
    for i in range(1, len(u) - 1):
        if u[i] >= u[i - 1] and u[i] > u[i + 1]:
            # Quadratic interpolation around (i-1, i, i+1)
            x0, x1, x2 = phi[i - 1], phi[i], phi[i + 1]
            y0, y1, y2 = u[i - 1], u[i], u[i + 1]
            # vertex shift from x1 for equally spaced points:
            denom = (y0 - 2.0 * y1 + y2)
            if abs(denom) > 1e-20:
                dx = 0.5 * (y0 - y2) / denom * (x2 - x1)
                out.append(x1 + dx)
            else:
                out.append(x1)
    return out


def mercury_precession_numeric(
    *,
    gamma_ppn: float,
    beta_ppn: float,
    orbit: MercuryOrbit | None = None,
    gm_sun_m3_s2: float = 1.32712440018e20,
    c_m_s: float = 299792458.0,
    n_orbits: int = 12,
    steps_per_orbit: int = 24000,
) -> NumericPrecessionResult:
    """Integrate orbit equation and return anomalous precession."""
    if orbit is None:
        orbit = MercuryOrbit()
    a = float(orbit.semimajor_axis_m)
    e = float(orbit.eccentricity)
    p = a * (1.0 - e * e)
    inv_p = 1.0 / p
    kappa = (2.0 + 2.0 * float(gamma_ppn) - float(beta_ppn)) * (float(gm_sun_m3_s2) / (c_m_s * c_m_s))

    # Start at perihelion: r_p=a(1-e), so u0=1/r_p and du/dphi=0
    u0 = 1.0 / (a * (1.0 - e))
    up0 = 0.0

    hphi = (2.0 * math.pi) / float(steps_per_orbit)
    n_steps = int(max(4, n_orbits * steps_per_orbit))

    phi_vals = [0.0]
    u_vals = [u0]
    phi = 0.0
    u = u0
    up = up0
    for _ in range(n_steps):
        u, up = _rk4_step(phi, u, up, hphi, inv_p=inv_p, kappa=kappa)
        phi += hphi
        phi_vals.append(phi)
        u_vals.append(u)

    maxima = _find_local_max_positions(phi_vals, u_vals)
    if len(maxima) < 3:
        raise RuntimeError("Not enough detected perihelia; increase integration resolution")

    # Use mean shift between successive perihelia
    shifts = []
    for i in range(1, len(maxima)):
        dphi = maxima[i] - maxima[i - 1]
        shifts.append(dphi - 2.0 * math.pi)
    delta_per_orbit = sum(shifts) / len(shifts)

    orbits_per_century = (100.0 * 365.25) / float(orbit.period_days)
    arcsec_per_century = delta_per_orbit * ARCSEC_PER_RAD * orbits_per_century
    return NumericPrecessionResult(
        delta_per_orbit_rad=delta_per_orbit,
        arcsec_per_century=arcsec_per_century,
        perihelion_count=len(maxima),
    )
