from checks.compact_star import CompactStarParams, solve_compact_star_profile, scan_compact_star_grid


def test_solve_compact_star_profile_finite_outputs():
    p = CompactStarParams(theta_c=1.5, m2=1.0, lam=0.1, r_max=5.0, dr=0.02, theta_scale=1.0)
    out = solve_compact_star_profile(p)
    assert len(out["r"]) > 10
    assert len(out["psi"]) == len(out["r"])
    assert len(out["i_grad"]) == len(out["r"])
    assert len(out["j_refr"]) == len(out["r"])
    assert len(out["f_refr"]) == len(out["r"])
    assert out["z_surface"] >= 0.0
    assert out["delay_proxy"] >= 0.0
    assert out["n_peak"] >= 1.0
    assert max(out["j_refr"]) >= 0.0
    assert max(out["f_refr"]) >= 0.0
    assert out["ok_no_horizon"] is True


def test_scan_compact_star_grid_rows():
    rows = scan_compact_star_grid([1.0, 1.5], [1.0], [0.1, 0.2], r_max=3.0, dr=0.05)
    assert len(rows) == 4
    assert all(r["ok_no_horizon"] for r in rows)
