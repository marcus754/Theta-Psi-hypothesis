from checks.atomic_relativistic import atomic_relativistic_consistency


def test_atomic_relativistic_consistency_stays_in_local_limit():
    res = atomic_relativistic_consistency(theta_eff=1e-30, theta_scale=1.0)
    assert res.ok
    assert res.ok_n
    assert res.ok_ppn
    assert abs(res.n_value - 1.0) < 1e-12
    assert abs(res.gamma_ppn - 1.0) < 1e-15
    assert abs(res.beta_ppn - 1.0) < 1e-15
