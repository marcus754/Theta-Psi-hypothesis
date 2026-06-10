from checks.no_horizon import no_horizon_from_n_values, no_horizon_from_theta


def test_no_horizon_from_n_values_basic():
    assert no_horizon_from_n_values([1.0, 2.0, 10.0])
    assert not no_horizon_from_n_values([])
    assert not no_horizon_from_n_values([1.0, float("inf")])


def test_no_horizon_from_theta_finite_branch():
    theta = [0.0, 0.1, 1.0, 2.0]
    res = no_horizon_from_theta(theta, theta_scale=1.0)
    assert res["ok_no_horizon"]
    assert res["n_min"] >= 1.0
    assert res["n_peak"] >= res["n_min"]
