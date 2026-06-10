import math

import background.frw_background as bg


def test_phase_energy_positive_and_q_behavior():
    # Small but nonzero Q; start with moderate R so no singularity
    params = bg.Params(gamma=1.0, V0=0.1, Q=1e-8, omegas=__import__('background.components', fromlist=['Omegas']).Omegas(H0=1e-4))
    # (a, x, xdot, psi(=R), psidot)
    y0 = (1e-3, 0.1, 0.0, 0.05, 0.0)
    # Short adaptive run to ensure stability
    out = bg.integrate_background_rk45(0.0, y0, t1=0.01, p=params, h0=1e-4)
    # H stays positive and psi does not cross zero
    assert min(out["H"]) > 0.0
    assert min(out["psi"]) > 0.0
