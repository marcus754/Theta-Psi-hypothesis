from scripts import run_grid_scan as rgs


def test_scan_point_runs():
    cfg = rgs.ScanConfig(nsteps=20, h=0.01)
    res = rgs.scan_point(0.8, 0.05, cfg)
    required = {
        "gamma",
        "V0",
        "ok",
        "ok_linear",
        "ok_sub",
        "ok_cutoff",
        "ok_n_bounds",
        "ok_no_horizon",
        "H_end",
        "theta_end",
        "theta_ratio_max",
        "n_surface",
        "n_peak",
        "z_surface_proxy",
        "delay_proxy",
    }
    assert required.issubset(set(res.keys()))
    assert res["n_peak"] >= res["n_surface"] >= 1.0


def test_run_scan_collects():
    cfg = rgs.ScanConfig(nsteps=10, h=0.01)
    rows = rgs.run_scan([0.8], [0.05, 0.1], cfg)
    assert len(rows) == 2
    assert all(row["gamma"] == 0.8 for row in rows)
