# -*- coding: utf-8 -*-
"""
Unified prediction suite for Θ–Ψ weak-field tests.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Sequence

from checks.ppn_full import compute_ppn_params
from checks.mercury_precession import MercuryObservation, mercury_precession_result
from checks.mercury_orbit_numeric import mercury_precession_numeric
from checks.classic_tests import (
    solar_surface_redshift,
    solar_surface_redshift_reference_check,
    shapiro_delay_roundtrip_seconds,
    shapiro_delay_reference_check,
    light_deflection_solar_limb_arcsec,
    solar_limb_deflection_reference_check,
    gps_clock_offset_seconds_per_day,
    gps_nist_reference_check,
    hafele_keating_reference_check,
    pound_rebka_reference_check,
    pound_rebka_redshift_fraction,
    hafele_keating_clock_offsets_seconds,
    twin_paradox_round_trip_proper_times,
    twin_paradox_reference_check,
    lense_thirring_node_precession_arcsec_per_year,
    lense_thirring_reference_check,
    binary_pulsar_periastron_advance_deg_per_year,
    binary_pulsar_reference_check,
    lunar_laser_ranging_roundtrip_seconds,
    lunar_laser_ranging_reference_check,
    mercury_precession_reference_check,
)


def parse_args(argv: Sequence[str] | None = None) -> argparse.Namespace:
    p = argparse.ArgumentParser(description="Run unified Θ–Ψ prediction suite")
    p.add_argument("--gamma_model", type=float, default=1.0)
    p.add_argument("--theta0", type=float, default=0.0)
    p.add_argument("--mercury_obs", type=float, default=42.98)
    p.add_argument("--mercury_sigma", type=float, default=0.04)
    p.add_argument("--n_orbits", type=int, default=12)
    p.add_argument("--steps_per_orbit", type=int, default=18000)
    p.add_argument("--output_json", type=str, default="results/prediction_suite.json")
    p.add_argument("--output_md", type=str, default="results/prediction_suite.md")
    return p.parse_args(argv)


def main(argv: Sequence[str] | None = None) -> int:
    args = parse_args(argv)
    ppn = compute_ppn_params(gamma_model=args.gamma_model, theta0=args.theta0)

    obs = MercuryObservation(args.mercury_obs, args.mercury_sigma)
    mer_ana = mercury_precession_result(gamma_ppn=ppn.gamma, beta_ppn=ppn.beta, obs=obs)
    mer_num = mercury_precession_numeric(
        gamma_ppn=ppn.gamma,
        beta_ppn=ppn.beta,
        n_orbits=args.n_orbits,
        steps_per_orbit=args.steps_per_orbit,
    )
    mer_num_pull = (mer_num.arcsec_per_century - obs.arcsec_per_century) / obs.sigma_arcsec_per_century

    z_sun = solar_surface_redshift(gamma_ppn=ppn.gamma)
    z_sun_ref = solar_surface_redshift_reference_check()
    shapiro = shapiro_delay_roundtrip_seconds(gamma_ppn=ppn.gamma)
    shapiro_ref = shapiro_delay_reference_check()
    alpha = light_deflection_solar_limb_arcsec(gamma_ppn=ppn.gamma)
    alpha_ref = solar_limb_deflection_reference_check()
    gps = gps_clock_offset_seconds_per_day()
    gps_ref = gps_nist_reference_check()
    mer_ref = mercury_precession_reference_check()
    pound_rebka = pound_rebka_redshift_fraction()
    pr_ref = pound_rebka_reference_check()
    hk = hafele_keating_clock_offsets_seconds()
    hk_ref = hafele_keating_reference_check()
    twin = twin_paradox_round_trip_proper_times()
    twin_ref = twin_paradox_reference_check()
    lt = lense_thirring_node_precession_arcsec_per_year()
    lt_ref = lense_thirring_reference_check()
    bp = binary_pulsar_periastron_advance_deg_per_year()
    bp_ref = binary_pulsar_reference_check()
    llr = lunar_laser_ranging_roundtrip_seconds()
    llr_ref = lunar_laser_ranging_reference_check()

    out = {
        "ppn": {"gamma": ppn.gamma, "beta": ppn.beta},
        "mercury_analytic_arcsec_per_century": mer_ana.arcsec_per_century,
        "mercury_analytic_pull_sigma": mer_ana.pull_sigma,
        "mercury_numeric_arcsec_per_century": mer_num.arcsec_per_century,
        "mercury_numeric_pull_sigma": mer_num_pull,
        "solar_redshift_z": z_sun,
        "solar_redshift_reference_z": z_sun_ref.reference,
        "solar_redshift_abs_error": z_sun_ref.abs_error,
        "solar_redshift_within_tolerance": z_sun_ref.within_tolerance,
        "shapiro_roundtrip_s": shapiro,
        "shapiro_reference_s": shapiro_ref.reference,
        "shapiro_abs_error_s": shapiro_ref.abs_error,
        "shapiro_within_tolerance": shapiro_ref.within_tolerance,
        "solar_limb_deflection_arcsec": alpha,
        "solar_limb_reference_arcsec": alpha_ref.reference,
        "solar_limb_abs_error_arcsec": alpha_ref.abs_error,
        "solar_limb_within_tolerance": alpha_ref.within_tolerance,
        "gps_clock_offset_s_per_day": gps.total_s_per_day,
        "gps_nist_reference_s_per_day": gps_ref.reference_s_per_day,
        "gps_nist_abs_error_s_per_day": gps_ref.abs_error_s_per_day,
        "gps_nist_within_tolerance": gps_ref.within_tolerance,
        "mercury_reference_arcsec_per_century": mer_ref.reference,
        "mercury_abs_error_arcsec_per_century": mer_ref.abs_error,
        "mercury_within_tolerance": mer_ref.within_tolerance,
        "pound_rebka_reference": pr_ref.reference,
        "pound_rebka_abs_error": pr_ref.abs_error,
        "pound_rebka_within_tolerance": pr_ref.within_tolerance,
        "hafele_keating_reference_east_s": hk_ref.reference_east_s,
        "hafele_keating_reference_west_s": hk_ref.reference_west_s,
        "hafele_keating_abs_error_east_s": hk_ref.abs_error_east_s,
        "hafele_keating_abs_error_west_s": hk_ref.abs_error_west_s,
        "hafele_keating_within_tolerance": hk_ref.within_tolerance,
        "pound_rebka_redshift": pound_rebka,
        "hafele_keating_east_s": hk.east_s,
        "hafele_keating_west_s": hk.west_s,
        "twin_paradox_age_gap_s": twin.age_gap_s,
        "twin_paradox_reference_age_gap_s": twin_ref.reference_age_gap_s,
        "twin_paradox_abs_error_s": twin_ref.abs_error_s,
        "twin_paradox_within_tolerance": twin_ref.within_tolerance,
        "lense_thirring_node_precession_arcsec_per_year": lt.node_precession_arcsec_per_year,
        "lense_thirring_reference_mas_per_year": lt_ref.reference,
        "lense_thirring_abs_error_mas_per_year": lt_ref.abs_error,
        "lense_thirring_within_tolerance": lt_ref.within_tolerance,
        "binary_pulsar_periastron_advance_deg_per_year": bp.periastron_advance_deg_per_year,
        "binary_pulsar_reference_deg_per_year": bp_ref.reference,
        "binary_pulsar_abs_error_deg_per_year": bp_ref.abs_error,
        "binary_pulsar_within_tolerance": bp_ref.within_tolerance,
        "lunar_laser_ranging_roundtrip_s": llr.roundtrip_light_time_s,
        "lunar_laser_ranging_reference_s": llr_ref.reference,
        "lunar_laser_ranging_abs_error_s": llr_ref.abs_error,
        "lunar_laser_ranging_within_tolerance": llr_ref.within_tolerance,
    }

    out_json = Path(args.output_json)
    out_json.parent.mkdir(parents=True, exist_ok=True)
    out_json.write_text(json.dumps(out, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")

    lines = []
    lines.append("# Prediction Suite\n\n")
    lines.append(f"- PPN: gamma={ppn.gamma:.6g}, beta={ppn.beta:.6g}\n")
    lines.append("\n## Mercury\n")
    lines.append(f"- Analytic: {mer_ana.arcsec_per_century:.6f} arcsec/century (pull={mer_ana.pull_sigma:.3f} sigma)\n")
    lines.append(f"- Numeric: {mer_num.arcsec_per_century:.6f} arcsec/century (pull={mer_num_pull:.3f} sigma)\n")
    lines.append("\n## Other Classical Tests\n")
    lines.append(f"- Solar redshift: {z_sun:.9g} (ref={z_sun_ref.reference:.9g}, err={z_sun_ref.abs_error:.3g}, ok={z_sun_ref.within_tolerance})\n")
    lines.append(f"- Shapiro roundtrip delay: {shapiro:.9g} s (ref={shapiro_ref.reference:.9g}, err={shapiro_ref.abs_error:.3g}, ok={shapiro_ref.within_tolerance})\n")
    lines.append(f"- Solar-limb deflection: {alpha:.9g} arcsec (ref={alpha_ref.reference:.9g}, err={alpha_ref.abs_error:.3g}, ok={alpha_ref.within_tolerance})\n")
    lines.append(f"- GPS clock offset: {gps.total_s_per_day:.9g} s/day (ref={gps_ref.reference_s_per_day:.9g}, err={gps_ref.abs_error_s_per_day:.3g}, ok={gps_ref.within_tolerance})\n")
    lines.append(f"- Pound-Rebka redshift: {pound_rebka:.9g} (ref={pr_ref.reference:.9g}, err={pr_ref.abs_error:.3g}, ok={pr_ref.within_tolerance})\n")
    lines.append(f"- Hafele-Keating east/west: {hk.east_s:.9g} s / {hk.west_s:.9g} s (ref={hk_ref.reference_east_s:.9g} s / {hk_ref.reference_west_s:.9g} s, ok={hk_ref.within_tolerance})\n")
    lines.append(f"- Twin-paradox age gap (10 years at 0.8c): {twin.age_gap_s:.9g} s (ref={twin_ref.reference_age_gap_s:.9g} s, err={twin_ref.abs_error_s:.3g}, ok={twin_ref.within_tolerance})\n")
    lines.append(f"- Lense-Thirring node precession: {lt.node_precession_arcsec_per_year:.9g} arcsec/year (ref={lt_ref.reference:.9g} mas/year, err={lt_ref.abs_error:.3g}, ok={lt_ref.within_tolerance})\n")
    lines.append(f"- Binary-pulsar periastron advance: {bp.periastron_advance_deg_per_year:.9g} deg/year (ref={bp_ref.reference:.9g}, err={bp_ref.abs_error:.3g}, ok={bp_ref.within_tolerance})\n")
    lines.append(f"- Lunar laser ranging roundtrip: {llr.roundtrip_light_time_s:.9g} s (ref={llr_ref.reference:.9g}, err={llr_ref.abs_error:.3g}, ok={llr_ref.within_tolerance})\n")

    out_md = Path(args.output_md)
    out_md.parent.mkdir(parents=True, exist_ok=True)
    out_md.write_text("".join(lines), encoding="utf-8")
    print(f"Wrote prediction suite JSON to {args.output_json}")
    print(f"Wrote prediction suite report to {args.output_md}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
