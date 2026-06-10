from checks.ppn_full import compute_ppn_params


def test_ppn_full_defaults_to_gr_values():
    # In the baseline (EH gravity), expect GR values regardless of model gamma.
    res = compute_ppn_params(gamma_model=1.0, theta0=0.1)
    assert abs(res.gamma - 1.0) < 1e-15
    assert abs(res.beta - 1.0) < 1e-15

