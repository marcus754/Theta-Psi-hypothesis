from theory.optical_metric import (
    CANONICAL_REFRACTIVE_PROFILE,
    CANONICAL_WEAK_FIELD_KAPPA,
    DIAGNOSTIC_REFRACTIVE_PROFILES,
    canonical_refractive_index_from_theta,
    diagnostic_refractive_index_from_theta,
    optical_metric,
    normalize_timelike,
    tensor_speed_optical_metric,
    refractive_index_from_theta,
    redshift_factor,
    optical_time_of_flight,
    n_finite_ok,
)


def test_optical_metric_local_frame():
    # Local inertial frame, u=(1,0,0,0)
    u = (1.0, 0.0, 0.0, 0.0)
    n = 1.3
    g = optical_metric(n, u)
    # g00 = n^2, gij = -delta_ij
    assert abs(g[0][0] - n*n) < 1e-15
    for i in range(1, 4):
        for j in range(1, 4):
            if i == j:
                assert abs(g[i][j] + 1.0) < 1e-15
            else:
                assert abs(g[i][j]) < 1e-15


def test_normalize_timelike_and_cT():
    # Non-unit timelike vector should be normalized
    u = (2.0, 0.0, 0.0, 0.0)
    un = normalize_timelike(u)
    assert abs(un[0] - 1.0) < 1e-15
    # Tensor speed from optical metric baseline is luminal
    assert tensor_speed_optical_metric() == 1.0


def test_n_and_redshift_are_finite():
    # Vacuum background: theta = theta_scale should yield n = 1.0
    n_bg = refractive_index_from_theta(1.0, theta_scale=1.0)
    # Large theta should yield larger n
    n_big = refractive_index_from_theta(1e6, theta_scale=1.0)
    
    assert abs(n_bg - 1.0) < 1e-12
    assert n_big > 1.0
    
    # Redshift between background and compact region (large theta)
    zfac = redshift_factor(10.0, 1.0, theta_scale=1.0)
    assert zfac > 1.0


def test_canonical_refractive_index_has_no_profile_knob():
    # On background (theta = theta_scale), n should be 1.0
    nbg = canonical_refractive_index_from_theta(1.0, theta_scale=1.0)
    # Above background, n should be > 1.0
    n_higher = canonical_refractive_index_from_theta(2.0, theta_scale=1.0)
    
    assert CANONICAL_WEAK_FIELD_KAPPA == 2.0
    assert CANONICAL_REFRACTIVE_PROFILE == "asinh"
    assert DIAGNOSTIC_REFRACTIVE_PROFILES == ("linear", "exp", "tanh2")
    assert abs(nbg - 1.0) < 1e-12
    assert n_higher > nbg


def test_profile_overrides_are_diagnostic_only():
    import pytest

    with pytest.warns(DeprecationWarning):
        legacy = refractive_index_from_theta(1.0, theta_scale=1.0, profile="linear")
    diagnostic = diagnostic_refractive_index_from_theta(1.0, theta_scale=1.0, profile="linear")
    assert legacy == diagnostic


def test_optical_time_and_bounds_check():
    dl = [1.0, 2.0, 3.0]
    theta = [0.0, 1.0, 2.0]
    t = optical_time_of_flight(dl, theta, theta_scale=1.0)
    assert t > 0.0
    assert n_finite_ok(theta, theta_scale=1.0)
