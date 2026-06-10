import background.frw_background as bg
from background.components import Omegas
from checks.health import evaluate_health


def test_unified_health_ok_on_short_traj():
    p = bg.Params(gamma=1.0, V0=0.1, Q=0.0, omegas=Omegas(H0=1e-4))
    y0 = (1e-3, 0.2, 0.0, 0.1, 0.0)
    out = bg.integrate_background_rk45(0.0, y0, t1=0.05, p=p, h0=1e-3)
    res = evaluate_health(out, gamma=p.gamma, V0=p.V0, Q=p.Q, ctheta2=0.9, cpsi2=0.9,
                          cutoff=1e-1, ppn_eps_max=0.5, check_bbn=False)
    assert res["ok"]
    assert res["ok_linear"] and res["ok_sub"] and res["ok_cutoff"] and res["ok_ppn"]
    assert res["ok_n_bounds"]


def test_unified_health_n_bounds_check_enabled():
    p = bg.Params(gamma=1.0, V0=0.1, Q=0.0, omegas=Omegas(H0=1e-4))
    y0 = (1e-3, 0.2, 0.0, 0.1, 0.0)
    out = bg.integrate_background_rk45(0.0, y0, t1=0.05, p=p, h0=1e-3)
    res = evaluate_health(
        out,
        gamma=p.gamma,
        V0=p.V0,
        Q=p.Q,
        check_n_bounds=True,
        theta_scale=1.0,
        check_bbn=False,
    )
    assert res["ok_n_bounds"]
    assert res["ok_no_horizon"]


def test_unified_health_no_horizon_check_enabled():
    p = bg.Params(gamma=1.0, V0=0.1, Q=0.0, omegas=Omegas(H0=1e-4))
    y0 = (1e-3, 0.2, 0.0, 0.1, 0.0)
    out = bg.integrate_background_rk45(0.0, y0, t1=0.05, p=p, h0=1e-3)
    res = evaluate_health(
        out,
        gamma=p.gamma,
        V0=p.V0,
        Q=p.Q,
        check_n_bounds=True,
        check_no_horizon=True,
        theta_scale=1.0,
        check_bbn=False,
    )
    assert res["ok_no_horizon"]
