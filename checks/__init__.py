"""Physical checks and derived observables for Θ–Ψ.

This package belongs to the kernel-side validation layer:
- weak-field tests;
- strong-field profile diagnostics;
- health / causality / stability vetoes.

It is intentionally separate from the observational policy layer in
`observational/`, where pass/fail protocols and model comparisons live.
"""

from .cutoff import cutoff_ok
from .superluminal import has_superluminal, has_superluminal_adm, finite_k_stable_along_traj
from .early_universe import bbn_safe
from .ppn_local import ppn_proxy_ok
from .ppn_full import compute_ppn_params
from .atomic_relativistic import AtomicRelativisticConsistency, atomic_relativistic_consistency
from .no_horizon import no_horizon_from_n_values, no_horizon_from_theta
from .compact_star import CompactStarParams, solve_compact_star_profile, scan_compact_star_grid
from .strongfield_primary import StrongFieldPrimarySummary, summarize_strongfield_primary
from .strongfield_reference import StrongFieldEventReference, StrongFieldReferenceSummary, strongfield_public_reference_check
from .static_vacuum import StaticVacuumParams, build_static_vacuum_profile, theta_of_r
from .mercury_precession import (
    MercuryOrbit,
    MercuryObservation,
    MercuryPrecessionResult,
    anomalous_precession_ppn,
    mercury_precession_result,
)
from .mercury_orbit_numeric import NumericPrecessionResult, mercury_precession_numeric
from .classic_tests import (
    GpsClockOffset,
    GpsReferenceCheck,
    HafeleKeatingReferenceCheck,
    ScalarReferenceCheck,
    TwinParadoxReferenceCheck,
    HafeleKeatingClockOffset,
    TwinParadoxResult,
    FrameDraggingResult,
    BinaryPulsarTimingResult,
    LunarLaserRangingResult,
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

__all__ = [
    "cutoff_ok",
    "has_superluminal",
    "has_superluminal_adm",
    "bbn_safe",
    "ppn_proxy_ok",
    "compute_ppn_params",
    "AtomicRelativisticConsistency",
    "atomic_relativistic_consistency",
    "finite_k_stable_along_traj",
    "no_horizon_from_n_values",
    "no_horizon_from_theta",
    "CompactStarParams",
    "solve_compact_star_profile",
    "scan_compact_star_grid",
    "StrongFieldPrimarySummary",
    "summarize_strongfield_primary",
    "StrongFieldEventReference",
    "StrongFieldReferenceSummary",
    "strongfield_public_reference_check",
    "StaticVacuumParams",
    "build_static_vacuum_profile",
    "theta_of_r",
    "MercuryOrbit",
    "MercuryObservation",
    "MercuryPrecessionResult",
    "anomalous_precession_ppn",
    "mercury_precession_result",
    "NumericPrecessionResult",
    "mercury_precession_numeric",
    "solar_surface_redshift",
    "shapiro_delay_roundtrip_seconds",
    "light_deflection_solar_limb_arcsec",
    "GpsClockOffset",
    "GpsReferenceCheck",
    "HafeleKeatingReferenceCheck",
    "ScalarReferenceCheck",
    "TwinParadoxReferenceCheck",
    "HafeleKeatingClockOffset",
    "TwinParadoxResult",
    "FrameDraggingResult",
    "BinaryPulsarTimingResult",
    "LunarLaserRangingResult",
    "gps_clock_offset_seconds_per_day",
    "gps_nist_reference_check",
    "hafele_keating_reference_check",
    "pound_rebka_reference_check",
    "pound_rebka_redshift_fraction",
    "hafele_keating_clock_offsets_seconds",
    "twin_paradox_round_trip_proper_times",
    "twin_paradox_reference_check",
    "lense_thirring_node_precession_arcsec_per_year",
    "lense_thirring_reference_check",
    "binary_pulsar_periastron_advance_deg_per_year",
    "binary_pulsar_reference_check",
    "lunar_laser_ranging_roundtrip_seconds",
    "lunar_laser_ranging_reference_check",
    "mercury_precession_reference_check",
    "solar_surface_redshift_reference_check",
    "shapiro_delay_reference_check",
    "solar_limb_deflection_reference_check",
]
