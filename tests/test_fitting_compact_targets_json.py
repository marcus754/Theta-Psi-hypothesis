from pathlib import Path

from fitting.data_io import load_compact_targets_json, load_compact_star_targets_json, load_strong_field_bounds_json
from fitting.run_mcmc import parse_args
from fitting.data_io import load_forward_targets_json, load_forward_event_catalog_json


def test_load_compact_targets_json_demo():
    p = Path("fitting/data/compact_demo_targets.json")
    d = load_compact_targets_json(p)
    assert d["mode"] == "gaussian"
    assert d["z_surface"] > 0.0
    assert d["sigma_z_surface"] > 0.0
    assert d["delay_proxy"] > 0.0
    assert d["sigma_delay_proxy"] > 0.0
    assert d["n_peak"] > 1.0
    assert d["sigma_n_peak"] > 0.0


def test_load_compact_targets_json_constraints_observational_proxy():
    p = Path("fitting/data/compact_constraints_observational_proxy.json")
    d = load_compact_targets_json(p)
    assert d["mode"] == "constraints"
    assert d["z_surface_min"] > 0.0
    assert d["sigma_z_surface_min"] > 0.0
    assert d["delay_proxy_min"] > 0.0
    assert d["sigma_delay_proxy_min"] > 0.0


def test_parse_args_accepts_compact_targets_json():
    ns = parse_args([
        "--use_compact",
        "--compact_targets_json",
        "fitting/data/compact_demo_targets.json",
        "--priors_json",
        "fitting/data/canonical_priors.json",
        "--steps",
        "5",
    ])
    assert ns.use_compact
    assert ns.compact_targets_json.endswith("compact_demo_targets.json")
    assert ns.priors_json.endswith("canonical_priors.json")


def test_parse_args_accepts_compact_star_link_background():
    ns = parse_args([
        "--use_compact_star",
        "--compact_star_link_background",
        "--compact_star_psi_grad_weight",
        "0.7",
        "--compact_star_theta_psi_coupling",
        "0.35",
        "--compact_star_refractive_strength",
        "0.8",
        "--compact_star_refractive_j_star",
        "1.4",
        "--compact_star_kappa_refr",
        "1.2",
        "--steps",
        "5",
    ])
    assert ns.use_compact_star
    assert ns.compact_star_link_background
    assert ns.compact_star_psi_grad_weight == 0.7
    assert ns.compact_star_theta_psi_coupling == 0.35
    assert ns.compact_star_refractive_strength == 0.8
    assert ns.compact_star_refractive_j_star == 1.4
    assert ns.compact_star_kappa_refr == 1.2


def test_parse_args_accepts_use_forward():
    ns = parse_args([
        "--use_forward",
        "--forward_targets_json",
        "fitting/data/forward_observational_proxy.json",
        "--steps",
        "5",
    ])
    assert ns.use_forward
    assert ns.forward_targets_json.endswith("forward_observational_proxy.json")


def test_load_strong_field_bounds_registry():
    p = Path("fitting/data/strong_field_bounds_registry.json")
    d = load_strong_field_bounds_json(p)
    assert "entries" in d
    assert len(d["entries"]) >= 1
    assert d["entries"][0]["observed_on"] != ""


def test_load_forward_targets_json_demo():
    p = Path("fitting/data/forward_observational_proxy.json")
    d = load_forward_targets_json(p)
    assert d["ring_m87_uas"] > 0.0
    assert d["sigma_ring_m87_uas"] > 0.0


def test_parse_args_accepts_use_forward_events():
    ns = parse_args([
        "--use_forward_events",
        "--forward_event_catalog_json",
        "fitting/data/forward_event_catalog.json",
        "--steps",
        "5",
    ])
    assert ns.use_forward_events
    assert ns.forward_event_catalog_json.endswith("forward_event_catalog.json")


def test_load_forward_event_catalog_json_demo():
    ev = load_forward_event_catalog_json(Path("fitting/data/forward_event_catalog.json"))
    assert len(ev) >= 1


def test_load_compact_star_targets_json_demo():
    p = Path("fitting/data/compact_star_demo_targets.json")
    d = load_compact_star_targets_json(p)
    assert d["z_surface"] > 0.0
    assert d["sigma_z_surface"] > 0.0
    assert d["delay_proxy"] > 0.0
    assert d["sigma_delay_proxy"] > 0.0
