import math

import background.frw_background as bg
from background.components import Omegas


def test_integrate_background_lnN_increases_a():
    p = bg.Params(gamma=1.0, V0=0.1, Q=0.0, omegas=Omegas(H0=1e-4))
    a0 = 1e-3
    x0, xdot0 = 0.1, 0.0
    R0, Rdot0 = 0.05, 0.0
    N0 = math.log(a0)
    out = bg.integrate_background_lnN(a0, x0, xdot0, R0, Rdot0, N1=N0 + 0.1, p=p, hN=0.01)
    assert math.isclose(out["a"][0], a0, rel_tol=0, abs_tol=1e-18)
    assert out["a"][-1] > a0
    assert out["H"][-1] > 0.0
