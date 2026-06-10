import math

from background.components import Omegas, rho_r, rho_m, rho_L
import background.frw_background as bg


def test_radiation_dominance_small_a():
    om = Omegas(H0=1e-4, Omega_r=1e-2, Omega_b=0.15, Omega_c=0.15, Omega_L=0.69)
    a = 1e-3
    assert rho_r(a, om) > rho_m(a, om)
    assert rho_r(a, om) > rho_L(a, om)


def test_background_demo_increases_a():
    demo = bg.demo()
    assert demo["a_end"] > 1e-3
    assert demo["H_end"] > 0.0
