from checks.compact_star import CompactStarParams, solve_compact_star_profile
from fitting.forward_transfer import explicit_forward_observables


def test_explicit_forward_transfer_returns_positive_observables():
    p = CompactStarParams(theta_c=0.5, m2=1.0, lam=0.001)
    prof = solve_compact_star_profile(p)
    obs = explicit_forward_observables(prof)
    assert obs["ring_m87_uas"] > 0.0
    assert obs["ring_sgra_uas"] > 0.0
    assert obs["echo_delay_ms"] > 0.0
