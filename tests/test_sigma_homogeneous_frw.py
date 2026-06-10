import math

import background.frw_background as bg
from background.components import Omegas


def test_j_refr_stays_zero_on_homogeneous_frw():
    p = bg.Params(gamma=1.0, V0=0.1, Q=0.0, omegas=Omegas(H0=1e-4))
    out = bg.integrate_background_lnN(
        1e-3,
        0.1,
        0.0,
        0.05,
        0.0,
        N1=math.log(1e-3) + 0.05,
        p=p,
        hN=0.01,
    )
    assert all(abs(j) < 1e-15 for j in out["J_refr"])
    assert all(abs(g) < 1e-15 for g in out["I_grad"])

def test_background_output_exposes_two_field_refractive_track():
    p = bg.Params(gamma=0.5, V0=0.2, Q=0.0, omegas=Omegas(H0=1e-4))
    out = bg.integrate_background_lnN(1e-3, 0.1, 0.0, 0.05, 0.0, N1=math.log(1e-3) + 0.02, p=p, hN=0.01)
    assert "J_refr" in out
    assert "I_grad" in out
