import background.frw_background as bg
from background.components import Omegas
from checks.health import evaluate_health


def test_unified_health_ok_adm_strict_on_short_traj():
    p = bg.Params(gamma=1.0, V0=0.1, Q=0.0, omegas=Omegas(H0=1e-4))
    y0 = (1e-3, 0.2, 0.0, 0.1, 0.0)
    out = bg.integrate_background_rk45(0.0, y0, t1=0.05, p=p, h0=1e-3)
    res = evaluate_health(
        out,
        gamma=p.gamma,
        V0=p.V0,
        Q=p.Q,
        ctheta2=0.95,
        cpsi2=0.95,
        use_adm=True,
        cR2=0.95,
        cphi2=0.95,
        # finite_k left as default; ADM path forces it on
        check_bbn=False,
    )
    assert res["ok"]
    assert res["ok_linear"] and res["ok_sub"]

