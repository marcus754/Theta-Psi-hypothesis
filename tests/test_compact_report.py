from pathlib import Path

from checks.compact_report import compute_compact_metrics, write_compact_markdown


def test_compact_metrics_finite_and_bounded():
    theta = [0.0, 0.2, 0.5, 1.0]
    m = compute_compact_metrics(theta, theta_scale=1.0, n_obs=1.0)
    assert m["n_surface"] >= 1.0
    assert m["n_peak"] >= m["n_surface"]
    assert m["z_surface"] >= 0.0
    assert m["delay"] >= 0.0
    assert m["has_horizon"] is False


def test_compact_metrics_grows_with_theta():
    theta_small = [0.0, 0.2, 0.5, 1.0]
    theta_big = [0.0, 0.5, 1.0, 1.5]
    m_small = compute_compact_metrics(theta_small, theta_scale=1.0, n_obs=1.0)
    m_big = compute_compact_metrics(theta_big, theta_scale=1.0, n_obs=1.0)
    assert m_big["z_surface"] > m_small["z_surface"]
    assert m_big["delay"] > m_small["delay"]


def test_write_compact_markdown(tmp_path: Path):
    rows = [
        {
            "gamma": 1.0,
            "V0": 0.1,
            "ok": True,
            "ok_n_bounds": True,
            "z_surface": 12.0,
            "n_surface": 13.0,
            "n_peak": 15.0,
            "delay": 30.0,
        }
    ]
    out = tmp_path / "compact_report.md"
    write_compact_markdown(rows, str(out))
    text = out.read_text(encoding="utf-8")
    assert "Compact Regime Report" in text
    assert "z_surface" in text
