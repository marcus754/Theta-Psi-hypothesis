from checks.static_vacuum import StaticVacuumParams, build_static_vacuum_profile, theta_of_r


def test_theta_profile_decreases_outward_for_positive_source_strength():
    th_near = theta_of_r(2.0, theta_inf=1.0, source_strength=3.0)
    th_far = theta_of_r(20.0, theta_inf=1.0, source_strength=3.0)
    assert th_near > th_far


def test_static_vacuum_profile_is_monotone_and_finite():
    prof = build_static_vacuum_profile(
        StaticVacuumParams(
            theta_inf=1.0,
            source_strength=2.0,
            r_min=2.0,
            r_max=40.0,
            dr=0.2,
            theta_scale=1.0,
        )
    )
    assert prof["monotone_theta"]
    assert prof["ok_no_horizon"]
    assert prof["z_edge"] >= 0.0
    assert prof["delay_proxy"] >= 0.0
