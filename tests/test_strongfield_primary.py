from checks.strongfield_primary import summarize_strongfield_primary


def test_strongfield_primary_returns_finite_no_horizon_summary():
    res = summarize_strongfield_primary(theta_c=1.5)
    assert res.ok_no_horizon is True
    assert res.classification == "compact_no_horizon"
    assert res.n_peak > 1.0
    assert res.z_surface >= 0.0
    assert res.delay_proxy >= 0.0


def test_strongfield_primary_grows_with_compactness():
    lo = summarize_strongfield_primary(theta_c=1.0)
    hi = summarize_strongfield_primary(theta_c=3.0)
    assert hi.z_surface > lo.z_surface
    assert hi.delay_proxy > lo.delay_proxy
    assert hi.n_peak > lo.n_peak
