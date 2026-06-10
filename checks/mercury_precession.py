# -*- coding: utf-8 -*-
"""
Mercury perihelion-precession check for the Θ–Ψ model.

Uses the standard weak-field PPN expression for the relativistic anomalous
advance per orbit:

    Δϖ = ((2 + 2γ - β) / 3) * (6π GM / (a (1 - e^2) c^2)).

For γ=β=1 this reduces to the GR value.
"""
from __future__ import annotations

from dataclasses import dataclass
import math


ARCSEC_PER_RAD = 206264.80624709636


@dataclass
class MercuryOrbit:
    # SI units
    semimajor_axis_m: float = 5.7909175e10
    eccentricity: float = 0.205630
    period_days: float = 87.9691


@dataclass
class MercuryObservation:
    # Observed anomalous (relativistic) advance, arcsec/century
    arcsec_per_century: float = 42.98
    sigma_arcsec_per_century: float = 0.04


@dataclass
class MercuryPrecessionResult:
    delta_per_orbit_rad: float
    arcsec_per_century: float
    pull_sigma: float | None


def anomalous_precession_ppn(
    *,
    gamma_ppn: float,
    beta_ppn: float,
    orbit: MercuryOrbit | None = None,
    gm_sun_m3_s2: float = 1.32712440018e20,
    c_m_s: float = 299792458.0,
) -> float:
    """Return anomalous perihelion advance in arcsec/century."""
    if orbit is None:
        orbit = MercuryOrbit()
    a = float(orbit.semimajor_axis_m)
    e = float(orbit.eccentricity)
    p_days = float(orbit.period_days)
    if not (a > 0.0 and 0.0 <= e < 1.0 and p_days > 0.0):
        raise ValueError("Invalid orbital parameters")

    pref = (2.0 + 2.0 * float(gamma_ppn) - float(beta_ppn)) / 3.0
    delta_rad = pref * (6.0 * math.pi * float(gm_sun_m3_s2)) / (a * (1.0 - e * e) * c_m_s * c_m_s)

    orbits_per_century = (100.0 * 365.25) / p_days
    return delta_rad * ARCSEC_PER_RAD * orbits_per_century


def mercury_precession_result(
    *,
    gamma_ppn: float,
    beta_ppn: float,
    orbit: MercuryOrbit | None = None,
    obs: MercuryObservation | None = None,
) -> MercuryPrecessionResult:
    """Return predicted arcsec/century and optional pull vs observation."""
    if orbit is None:
        orbit = MercuryOrbit()
    val = anomalous_precession_ppn(gamma_ppn=gamma_ppn, beta_ppn=beta_ppn, orbit=orbit)
    # Recover per-orbit rad from arcsec/century conversion.
    orbits_per_century = (100.0 * 365.25) / float(orbit.period_days)
    delta_per_orbit_rad = val / (ARCSEC_PER_RAD * orbits_per_century)

    pull = None
    if obs is not None and obs.sigma_arcsec_per_century > 0.0:
        pull = (val - obs.arcsec_per_century) / obs.sigma_arcsec_per_century

    return MercuryPrecessionResult(
        delta_per_orbit_rad=delta_per_orbit_rad,
        arcsec_per_century=val,
        pull_sigma=pull,
    )
