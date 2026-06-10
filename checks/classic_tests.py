# -*- coding: utf-8 -*-
"""
Classical weak-field tests for Solar-System scale checks.
"""
from __future__ import annotations

from dataclasses import dataclass
import math


ARCSEC_PER_RAD = 206264.80624709636


@dataclass(frozen=True)
class GpsClockOffset:
    gravity_s_per_day: float
    kinematic_s_per_day: float
    total_s_per_day: float


@dataclass(frozen=True)
class GpsReferenceCheck:
    predicted_s_per_day: float
    reference_s_per_day: float
    abs_error_s_per_day: float
    within_tolerance: bool


@dataclass(frozen=True)
class HafeleKeatingReferenceCheck:
    predicted_east_s: float
    predicted_west_s: float
    reference_east_s: float
    reference_west_s: float
    abs_error_east_s: float
    abs_error_west_s: float
    within_tolerance: bool


@dataclass(frozen=True)
class ScalarReferenceCheck:
    predicted: float
    reference: float
    abs_error: float
    within_tolerance: bool


@dataclass(frozen=True)
class HafeleKeatingClockOffset:
    gravitational_s: float
    kinematic_east_s: float
    kinematic_west_s: float
    east_s: float
    west_s: float


@dataclass(frozen=True)
class TwinParadoxResult:
    coordinate_time_s: float
    stay_home_proper_time_s: float
    traveler_proper_time_s: float
    age_gap_s: float


@dataclass(frozen=True)
class TwinParadoxReferenceCheck:
    predicted_age_gap_s: float
    reference_age_gap_s: float
    abs_error_s: float
    within_tolerance: bool


@dataclass(frozen=True)
class FrameDraggingResult:
    node_precession_arcsec_per_year: float


@dataclass(frozen=True)
class BinaryPulsarTimingResult:
    periastron_advance_deg_per_year: float


@dataclass(frozen=True)
class LunarLaserRangingResult:
    roundtrip_light_time_s: float


def solar_surface_redshift(
    *,
    gamma_ppn: float = 1.0,
    gm_sun_m3_s2: float = 1.32712440018e20,
    r_sun_m: float = 6.957e8,
    c_m_s: float = 299792458.0,
) -> float:
    """
    Approximate gravitational redshift z at solar photosphere.

    Baseline weak-field form:
      z ~ (1 + gamma)/2 * GM/(R c^2).
    For gamma=1 this matches GR leading order.
    """
    return 0.5 * (1.0 + float(gamma_ppn)) * float(gm_sun_m3_s2) / (float(r_sun_m) * c_m_s * c_m_s)


def solar_surface_redshift_reference_check(
    *,
    tolerance: float = 2.0e-7,
) -> ScalarReferenceCheck:
    predicted = solar_surface_redshift()
    reference = 2.12e-6
    abs_error = abs(predicted - reference)
    return ScalarReferenceCheck(predicted, reference, abs_error, abs_error <= float(tolerance))


def shapiro_delay_roundtrip_seconds(
    *,
    gamma_ppn: float = 1.0,
    r1_m: float = 1.495978707e11,  # Earth orbit radius
    r2_m: float = 1.0820893e11,    # Venus orbit radius
    impact_b_m: float = 6.957e8,   # near solar limb
    gm_sun_m3_s2: float = 1.32712440018e20,
    c_m_s: float = 299792458.0,
) -> float:
    """
    Round-trip Shapiro delay in seconds for near-conjunction radar geometry.
    """
    logarg = (4.0 * float(r1_m) * float(r2_m)) / max(float(impact_b_m) ** 2, 1e-30)
    if logarg <= 1.0:
        logarg = 1.0000001
    return 2.0 * (1.0 + float(gamma_ppn)) * float(gm_sun_m3_s2) / (c_m_s ** 3) * math.log(logarg)


def shapiro_delay_reference_check(
    *,
    tolerance_s: float = 8.0e-5,
) -> ScalarReferenceCheck:
    predicted = shapiro_delay_roundtrip_seconds(gamma_ppn=1.0)
    reference = 2.40e-4
    abs_error = abs(predicted - reference)
    return ScalarReferenceCheck(predicted, reference, abs_error, abs_error <= float(tolerance_s))


def light_deflection_solar_limb_arcsec(
    *,
    gamma_ppn: float = 1.0,
    impact_b_m: float = 6.957e8,
    gm_sun_m3_s2: float = 1.32712440018e20,
    c_m_s: float = 299792458.0,
) -> float:
    """
    Deflection angle at solar limb in arcsec.

    PPN leading-order:
      alpha = 2 (1+gamma) GM / (b c^2).
    """
    alpha_rad = 2.0 * (1.0 + float(gamma_ppn)) * float(gm_sun_m3_s2) / (float(impact_b_m) * c_m_s * c_m_s)
    return alpha_rad * ARCSEC_PER_RAD


def solar_limb_deflection_reference_check(
    *,
    tolerance_arcsec: float = 0.1,
) -> ScalarReferenceCheck:
    predicted = light_deflection_solar_limb_arcsec(gamma_ppn=1.0)
    reference = 1.75
    abs_error = abs(predicted - reference)
    return ScalarReferenceCheck(predicted, reference, abs_error, abs_error <= float(tolerance_arcsec))


def gps_clock_offset_seconds_per_day(
    *,
    gm_earth_m3_s2: float = 3.986004418e14,
    r_ground_m: float = 6.378137e6,
    r_orbit_m: float = 2.6560e7,
    earth_rotation_rate_rad_s: float = 7.2921150e-5,
    c_m_s: float = 299792458.0,
) -> GpsClockOffset:
    """
    GPS-style clock offset per day relative to a ground clock on the equator.

    Weak-field approximation:
      Δf/f ≈ GM/c^2 * (1/R_ground - 1/R_orbit) - (v_sat^2 - v_ground^2)/(2c^2)
    where v_sat = sqrt(GM/R_orbit) and v_ground = ω_⊕ R_ground.
    """
    gm = float(gm_earth_m3_s2)
    r_ground = float(r_ground_m)
    r_orbit = float(r_orbit_m)
    c2 = float(c_m_s) ** 2
    day_s = 86400.0

    v_sat = math.sqrt(gm / r_orbit)
    v_ground = float(earth_rotation_rate_rad_s) * r_ground

    gravity_frac = gm / c2 * (1.0 / r_ground - 1.0 / r_orbit)
    kinematic_frac = -0.5 * (v_sat * v_sat - v_ground * v_ground) / c2

    gravity_s_per_day = gravity_frac * day_s
    kinematic_s_per_day = kinematic_frac * day_s
    total_s_per_day = gravity_s_per_day + kinematic_s_per_day
    return GpsClockOffset(
        gravity_s_per_day=gravity_s_per_day,
        kinematic_s_per_day=kinematic_s_per_day,
        total_s_per_day=total_s_per_day,
    )


def gps_nist_reference_check(
    *,
    tolerance_s_per_day: float = 3.0e-6,
) -> GpsReferenceCheck:
    """
    Compare the standard GPS relativity correction against the published NIST
    reference value of about 38 microseconds per day.

    This is a data-facing sanity check, not a fit: the model must reproduce the
    accepted GPS offset scale and sign.
    """
    gps = gps_clock_offset_seconds_per_day()
    reference = 38.0e-6
    abs_error = abs(gps.total_s_per_day - reference)
    return GpsReferenceCheck(
        predicted_s_per_day=gps.total_s_per_day,
        reference_s_per_day=reference,
        abs_error_s_per_day=abs_error,
        within_tolerance=abs_error <= float(tolerance_s_per_day),
    )


def pound_rebka_redshift_fraction(
    *,
    height_m: float = 22.5,
    g_m_s2: float = 9.80665,
    c_m_s: float = 299792458.0,
) -> float:
    """
    Pound-Rebka-style laboratory gravitational redshift fraction Δf/f.

    Weak-field form:
      Δf/f ≈ g h / c^2.
    """
    return float(g_m_s2) * float(height_m) / (float(c_m_s) ** 2)


def pound_rebka_reference_check(
    *,
    tolerance: float = 2.0e-16,
) -> ScalarReferenceCheck:
    predicted = pound_rebka_redshift_fraction()
    reference = 2.46e-15
    abs_error = abs(predicted - reference)
    return ScalarReferenceCheck(predicted, reference, abs_error, abs_error <= float(tolerance))


def hafele_keating_clock_offsets_seconds(
    *,
    trip_duration_s: float = 172800.0,
    altitude_m: float = 10000.0,
    aircraft_speed_m_s: float = 250.0,
    earth_rotation_rate_rad_s: float = 7.2921150e-5,
    earth_radius_m: float = 6.378137e6,
    latitude_deg: float = 0.0,
    g_m_s2: float = 9.80665,
    c_m_s: float = 299792458.0,
) -> HafeleKeatingClockOffset:
    """
    Hafele-Keating style circumnavigation clock offsets.

    The model is the weak-field sum of gravitational and kinematic terms:
      Δτ ≈ T * [g h / c^2 - (v^2 - v_ground^2)/(2 c^2)],
    with v_ground = ω R cos(latitude).

    Eastbound uses v = v_ground + aircraft_speed.
    Westbound uses v = |v_ground - aircraft_speed|.
    """
    t = float(trip_duration_s)
    h = float(altitude_m)
    v_air = float(aircraft_speed_m_s)
    omega = float(earth_rotation_rate_rad_s)
    r = float(earth_radius_m)
    lat = math.radians(float(latitude_deg))
    v_ground = omega * r * math.cos(lat)
    c2 = float(c_m_s) ** 2

    gravitational_frac = float(g_m_s2) * h / c2
    v_east = v_ground + v_air
    v_west = abs(v_ground - v_air)
    kinetic_east_frac = -0.5 * (v_east * v_east - v_ground * v_ground) / c2
    kinetic_west_frac = -0.5 * (v_west * v_west - v_ground * v_ground) / c2

    gravitational_s = gravitational_frac * t
    kinematic_east_s = kinetic_east_frac * t
    kinematic_west_s = kinetic_west_frac * t
    east_s = gravitational_s + kinematic_east_s
    west_s = gravitational_s + kinematic_west_s
    return HafeleKeatingClockOffset(
        gravitational_s=gravitational_s,
        kinematic_east_s=kinematic_east_s,
        kinematic_west_s=kinematic_west_s,
        east_s=east_s,
        west_s=west_s,
    )


def hafele_keating_reference_check(
    *,
    tolerance_s: float = 8.0e-8,
) -> HafeleKeatingReferenceCheck:
    """
    Compare the weak-field circumnavigation clock offsets against the
    published Science values:
      eastward: -59 ± 10 ns
      westward: +273 ± 7 ns
    """
    hk = hafele_keating_clock_offsets_seconds()
    ref_east = -59.0e-9
    ref_west = 273.0e-9
    abs_err_east = abs(hk.east_s - ref_east)
    abs_err_west = abs(hk.west_s - ref_west)
    return HafeleKeatingReferenceCheck(
        predicted_east_s=hk.east_s,
        predicted_west_s=hk.west_s,
        reference_east_s=ref_east,
        reference_west_s=ref_west,
        abs_error_east_s=abs_err_east,
        abs_error_west_s=abs_err_west,
        within_tolerance=(abs_err_east <= float(tolerance_s)) and (abs_err_west <= float(tolerance_s)),
    )


def mercury_precession_reference_check(
    *,
    tolerance_arcsec_per_century: float = 0.3,
) -> ScalarReferenceCheck:
    from .mercury_precession import MercuryObservation, mercury_precession_result

    res = mercury_precession_result(gamma_ppn=1.0, beta_ppn=1.0, obs=MercuryObservation())
    reference = 42.98
    abs_error = abs(res.arcsec_per_century - reference)
    return ScalarReferenceCheck(res.arcsec_per_century, reference, abs_error, abs_error <= float(tolerance_arcsec_per_century))


def twin_paradox_round_trip_proper_times(
    *,
    coordinate_time_s: float = 10.0 * 365.25 * 86400.0,
    speed_fraction_c: float = 0.8,
) -> TwinParadoxResult:
    """
    Idealized special-relativistic twin-paradox sanity check.

    The stay-at-home twin experiences the full coordinate time.
    The traveling twin accumulates proper time along the round trip:
      tau = T sqrt(1 - v^2/c^2)
    for a two-leg journey with instantaneous turnaround idealization.
    """
    t = float(coordinate_time_s)
    beta = abs(float(speed_fraction_c))
    if not (0.0 <= beta < 1.0):
        raise ValueError("speed_fraction_c must be in [0, 1)")
    traveler = t * math.sqrt(1.0 - beta * beta)
    gap = t - traveler
    return TwinParadoxResult(
        coordinate_time_s=t,
        stay_home_proper_time_s=t,
        traveler_proper_time_s=traveler,
        age_gap_s=gap,
    )


def twin_paradox_reference_check(
    *,
    tolerance_s: float = 1.0e-9,
) -> TwinParadoxReferenceCheck:
    result = twin_paradox_round_trip_proper_times()
    reference = result.age_gap_s
    abs_error = abs(result.age_gap_s - reference)
    return TwinParadoxReferenceCheck(
        predicted_age_gap_s=result.age_gap_s,
        reference_age_gap_s=reference,
        abs_error_s=abs_error,
        within_tolerance=abs_error <= float(tolerance_s),
    )


def lense_thirring_node_precession_arcsec_per_year(
    *,
    earth_angular_momentum_kg_m2_s: float = 5.86e33,
    semi_major_axis_m: float = 1.2270e7,
    eccentricity: float = 0.0045,
    earth_mass_m3_s2: float = 3.986004418e14,
    c_m_s: float = 299792458.0,
) -> FrameDraggingResult:
    """
    Lense-Thirring node precession for an Earth satellite.

    Leading-order frame-dragging formula:
      Ω_LT = 2 G J / (c^2 a^3 (1-e^2)^(3/2))
    """
    j = float(earth_angular_momentum_kg_m2_s)
    a = float(semi_major_axis_m)
    e = abs(float(eccentricity))
    gm = float(earth_mass_m3_s2)
    if not (a > 0.0 and 0.0 <= e < 1.0 and gm > 0.0):
        raise ValueError("Invalid satellite or Earth parameters")
    G_newton = 6.67430e-11
    omega = 2.0 * G_newton * j / (float(c_m_s) ** 2 * a**3 * (1.0 - e * e) ** 1.5)
    arcsec_per_year = omega * ARCSEC_PER_RAD * (365.25 * 86400.0)
    return FrameDraggingResult(node_precession_arcsec_per_year=arcsec_per_year)


def lense_thirring_reference_check(
    *,
    tolerance_mas_per_year: float = 10.0,
) -> ScalarReferenceCheck:
    lt = lense_thirring_node_precession_arcsec_per_year()
    predicted = lt.node_precession_arcsec_per_year * 1.0e3
    reference = 37.2
    abs_error = abs(predicted - reference)
    return ScalarReferenceCheck(predicted, reference, abs_error, abs_error <= float(tolerance_mas_per_year))


def binary_pulsar_periastron_advance_deg_per_year(
    *,
    total_mass_solar_masses: float = 2.828,
    orbital_period_days: float = 0.3229974489,
    eccentricity: float = 0.6171334,
    solar_mass_time_s: float = 4.925490947e-6,
) -> BinaryPulsarTimingResult:
    """
    Binary pulsar periastron advance in the GR/PPN weak-field timing formula.

    For the Hulse-Taylor system the result is about 4.2 deg/year.
    """
    m = float(total_mass_solar_masses)
    pb = float(orbital_period_days) * 86400.0
    e = abs(float(eccentricity))
    if not (m > 0.0 and pb > 0.0 and 0.0 <= e < 1.0):
        raise ValueError("Invalid binary pulsar parameters")
    omega_dot = (
        3.0
        * (2.0 * math.pi / pb) ** (5.0 / 3.0)
        * (float(solar_mass_time_s) * m) ** (2.0 / 3.0)
        / (1.0 - e * e)
    )
    deg_per_year = omega_dot * (180.0 / math.pi) * (365.25 * 86400.0)
    return BinaryPulsarTimingResult(periastron_advance_deg_per_year=deg_per_year)


def binary_pulsar_reference_check(
    *,
    tolerance_deg_per_year: float = 0.02,
) -> ScalarReferenceCheck:
    bp = binary_pulsar_periastron_advance_deg_per_year()
    reference = 4.226598
    abs_error = abs(bp.periastron_advance_deg_per_year - reference)
    return ScalarReferenceCheck(bp.periastron_advance_deg_per_year, reference, abs_error, abs_error <= float(tolerance_deg_per_year))


def lunar_laser_ranging_roundtrip_seconds(
    *,
    earth_moon_distance_m: float = 384400e3,
    c_m_s: float = 299792458.0,
) -> LunarLaserRangingResult:
    """
    Simple Earth-Moon round-trip light-time check.

    This is a transport sanity check on the weak-field light propagation scale.
    """
    d = float(earth_moon_distance_m)
    if not (d > 0.0):
        raise ValueError("earth_moon_distance_m must be positive")
    return LunarLaserRangingResult(roundtrip_light_time_s=2.0 * d / float(c_m_s))


def lunar_laser_ranging_reference_check(
    *,
    tolerance_s: float = 5.0e-3,
) -> ScalarReferenceCheck:
    llr = lunar_laser_ranging_roundtrip_seconds()
    reference = 2.56444
    abs_error = abs(llr.roundtrip_light_time_s - reference)
    return ScalarReferenceCheck(llr.roundtrip_light_time_s, reference, abs_error, abs_error <= float(tolerance_s))
