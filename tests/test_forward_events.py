from pathlib import Path

from checks.compact_star import CompactStarParams, solve_compact_star_profile
from fitting.convert_public import convert_forward_events
from fitting.data_io import load_forward_event_catalog_json
from fitting.likelihoods_forward_events import evaluate_forward_events, loglike_forward_events


def test_load_forward_event_catalog_json():
    ev = load_forward_event_catalog_json(Path("fitting/data/forward_event_catalog.json"))
    assert len(ev) >= 1
    assert ev[0]["name"]


def test_forward_events_likelihood_finite():
    ev = load_forward_event_catalog_json(Path("fitting/data/forward_event_catalog.json"))
    prof = solve_compact_star_profile(CompactStarParams(theta_c=0.5, m2=1.0, lam=0.001))
    vals = evaluate_forward_events(prof, ev)
    assert len(vals) == len(ev)
    ll = loglike_forward_events(prof, ev)
    assert ll < 0.0


def test_convert_forward_events_roundtrip(tmp_path):
    src = Path("fitting/data/forward_events_public_proxy.csv")
    out = tmp_path / "forward_event_catalog.json"
    convert_forward_events(src, out)
    ev = load_forward_event_catalog_json(out)
    assert len(ev) == 2
    assert any(e["observable"] == "ring_uas" for e in ev)
    assert not any(e["observable"] == "echo_delay_ms" for e in ev)
