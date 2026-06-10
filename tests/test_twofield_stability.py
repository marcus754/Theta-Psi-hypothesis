from theory.twofield_stability import frw_stability_conditions, stationary_stability_conditions, stability_ledger


def test_frw_twofield_stability_basic():
    d = frw_stability_conditions(theta=1.2, gamma=0.8, ctheta2=0.9, cpsi2=1.0, mtheta2=0.1, V0=0.2)
    assert d["no_ghost"]
    assert d["no_gradient_instability"]
    assert d["no_tachyon"]


def test_stationary_twofield_operator_is_elliptic_in_nominal_regime():
    d = stationary_stability_conditions(
        theta=1.5,
        theta_prime=0.08,
        psi_prime=0.05,
        psi_grad_weight=1.0,
        refractive_strength=1.0,
        refractive_j_star=1.0,
        theta_scale=1.0,
    )
    assert d["positive_refr_gain"]
    assert d["elliptic"]
    assert d["det"] > 0.0


def test_stability_ledger_is_healthy_in_nominal_regime():
    d = stability_ledger(
        theta=1.2,
        gamma=0.8,
        theta_prime=0.08,
        psi_prime=0.05,
        ctheta2=0.9,
        cpsi2=1.0,
        mtheta2=0.1,
        V0=0.2,
        psi_grad_weight=1.0,
        refractive_strength=1.0,
        refractive_j_star=1.0,
        theta_scale=1.0,
    )
    assert d["healthy"]
    assert d["no_ghost"]
    assert d["no_gradient_instability"]
    assert d["no_tachyon"]
    assert d["no_extra_modes"]
