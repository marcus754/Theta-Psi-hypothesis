import math

import background.frw_background as bg
from background.components import Omegas
from checks.superluminal import finite_k_stable_along_traj


def test_finite_k_stability_with_gravM_uses_background_derivatives():
    # Short trajectory to populate H, Hdot, theta_dot, psi_dot
    p = bg.Params(gamma=1.0, V0=0.1, Q=0.0, omegas=Omegas(H0=1e-4))
    y0 = (1e-3, 0.2, 0.0, 0.1, 0.0)
    out = bg.integrate_background_rk45(0.0, y0, t1=0.05, p=p, h0=1e-3)

    ok = finite_k_stable_along_traj(
        out["theta"], out["psi"], out["a"], p.gamma, p.Q, p.V0,
        ctheta2=0.9, cR2=0.9, cphi2=0.95,
        kmin=1e-2, kmax=1e-1, nk=3,
        gravM=True,
        H_list=out["H"],
        theta_dot_list=out["theta_dot"],
        R_dot_list=out["psi_dot"],
        Hdot_list=out["Hdot"],
    )
    assert ok

