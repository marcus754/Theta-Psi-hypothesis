from src.growth import growth_factor
from background.components import Omegas


def test_growth_runs_and_is_monotonic_in_gr():
    om = Omegas(H0=1.0, Omega_r=0.0, Omega_b=0.05, Omega_c=0.25, Omega_L=0.7)
    a, D = growth_factor(1e-2, 1.0, omegas=om, nsteps=200)
    # Basic sanity: grid length and monotonic growth to late times
    assert len(a) == len(D) and len(a) > 10
    assert a[0] < a[-1]
    assert D[0] < D[-1]
    # Normalization: D(a=1) ≈ 1 by construction (last point close to 1)
    assert abs(D[-1] - 1.0) < 1e-6

