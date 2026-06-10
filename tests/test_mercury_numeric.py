from checks.mercury_orbit_numeric import mercury_precession_numeric
from checks.mercury_precession import anomalous_precession_ppn


def test_numeric_mercury_precession_near_gr():
    num = mercury_precession_numeric(gamma_ppn=1.0, beta_ppn=1.0, n_orbits=10, steps_per_orbit=12000)
    assert 42.0 < num.arcsec_per_century < 44.0


def test_numeric_matches_ppn_formula_reasonably():
    num = mercury_precession_numeric(gamma_ppn=1.0, beta_ppn=1.0, n_orbits=12, steps_per_orbit=16000)
    ana = anomalous_precession_ppn(gamma_ppn=1.0, beta_ppn=1.0)
    # Numerical integration should agree at sub-arcsec level with this resolution.
    assert abs(num.arcsec_per_century - ana) < 0.5
