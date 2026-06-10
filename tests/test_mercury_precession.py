from checks.mercury_precession import anomalous_precession_ppn, mercury_precession_result, MercuryObservation


def test_mercury_precession_gr_like_value():
    val = anomalous_precession_ppn(gamma_ppn=1.0, beta_ppn=1.0)
    # Canonical GR anomalous advance is about 43 arcsec/century.
    assert 42.0 < val < 44.0


def test_mercury_precession_pull_is_finite():
    obs = MercuryObservation(arcsec_per_century=42.98, sigma_arcsec_per_century=0.04)
    res = mercury_precession_result(gamma_ppn=1.0, beta_ppn=1.0, obs=obs)
    assert res.pull_sigma is not None
    assert abs(res.pull_sigma) < 5.0
