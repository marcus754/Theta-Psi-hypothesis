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


def test_solar_redshift_order_of_magnitude():
    z = solar_surface_redshift(gamma_ppn=1.0)
    assert 1e-6 < z < 5e-6


def test_solar_redshift_matches_reference_scale():
    ref = solar_surface_redshift_reference_check()
    assert ref.within_tolerance


def test_shapiro_delay_order_of_magnitude():
    dt = shapiro_delay_roundtrip_seconds(gamma_ppn=1.0)
    # Near-conjunction round-trip is of order ~100 microseconds.
    assert 1e-5 < dt < 5e-4


def test_shapiro_delay_matches_reference_scale():
    ref = shapiro_delay_reference_check()
    assert ref.within_tolerance


def test_solar_limb_deflection_near_gr_value():
    a = light_deflection_solar_limb_arcsec(gamma_ppn=1.0)
    assert 1.0 < a < 2.5


def test_solar_limb_deflection_matches_reference_scale():
    ref = solar_limb_deflection_reference_check()
    assert ref.within_tolerance


def test_gps_clock_offset_is_positive_and_microseconds_per_day():
    gps = gps_clock_offset_seconds_per_day()
    assert gps.gravity_s_per_day > 0.0
    assert gps.kinematic_s_per_day < 0.0
    assert 3e-5 < gps.total_s_per_day < 5e-5


def test_gps_matches_nist_reference_scale():
    ref = gps_nist_reference_check()
    assert ref.within_tolerance
    assert ref.abs_error_s_per_day < 3e-6


def test_pound_rebka_redshift_is_tiny_and_positive():
    z = pound_rebka_redshift_fraction()
    assert 2e-15 < z < 3e-15


def test_pound_rebka_matches_reference_scale():
    ref = pound_rebka_reference_check()
    assert ref.within_tolerance


def test_hafele_keating_has_opposite_signs_for_east_and_west_routes():
    hk = hafele_keating_clock_offsets_seconds()
    assert hk.east_s < 0.0
    assert hk.west_s > 0.0
    assert abs(hk.east_s) < 1e-6
    assert abs(hk.west_s) < 1e-6


def test_hafele_keating_matches_reference_scale():
    ref = hafele_keating_reference_check()
    assert ref.within_tolerance
    assert ref.abs_error_east_s < 5e-8
    assert ref.abs_error_west_s < 8e-8


def test_twin_paradox_traveler_ages_less():
    twin = twin_paradox_round_trip_proper_times()
    assert twin.stay_home_proper_time_s > twin.traveler_proper_time_s
    assert twin.age_gap_s > 0.0


def test_twin_paradox_matches_reference_scale():
    ref = twin_paradox_reference_check()
    assert ref.within_tolerance


def test_lense_thirring_is_small_and_positive():
    lt = lense_thirring_node_precession_arcsec_per_year()
    assert 1e-2 < lt.node_precession_arcsec_per_year < 1e-1


def test_lense_thirring_matches_reference_scale():
    ref = lense_thirring_reference_check()
    assert ref.within_tolerance


def test_binary_pulsar_periastron_advance_is_order_degrees_per_year():
    bp = binary_pulsar_periastron_advance_deg_per_year()
    assert 3.0 < bp.periastron_advance_deg_per_year < 5.5


def test_binary_pulsar_matches_reference_scale():
    ref = binary_pulsar_reference_check()
    assert ref.within_tolerance


def test_lunar_laser_ranging_roundtrip_is_seconds_scale():
    llr = lunar_laser_ranging_roundtrip_seconds()
    assert 2.4 < llr.roundtrip_light_time_s < 2.7


def test_lunar_laser_ranging_matches_reference_scale():
    ref = lunar_laser_ranging_reference_check()
    assert ref.within_tolerance


def test_mercury_precession_matches_reference_scale():
    ref = mercury_precession_reference_check()
    assert ref.within_tolerance
