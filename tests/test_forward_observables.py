from pathlib import Path

from fitting.forward_observables import forward_observables_from_compact_metrics
from fitting.likelihoods_forward import ForwardTargets, loglike_forward
from fitting.data_io import load_forward_targets_json


def test_forward_mapping_is_positive():
    obs = forward_observables_from_compact_metrics({"z_surface": 100.0, "delay_proxy": 700.0, "n_peak": 101.0})
    assert obs["ring_m87_uas"] > 0.0
    assert obs["ring_sgra_uas"] > 0.0
    assert obs["echo_delay_ms"] > 0.0


def test_loglike_forward_prefers_target_match():
    t = ForwardTargets(
        ring_m87_uas=42.0,
        sigma_ring_m87_uas=3.0,
        ring_sgra_uas=51.8,
        sigma_ring_sgra_uas=2.3,
    )
    good = {"ring_m87_uas": 42.0, "ring_sgra_uas": 51.8}
    bad = {"ring_m87_uas": 20.0, "ring_sgra_uas": 20.0}
    assert loglike_forward(good, t) > loglike_forward(bad, t)


def test_load_forward_targets_json():
    d = load_forward_targets_json(Path("fitting/data/forward_observational_proxy.json"))
    assert d["ring_m87_uas"] > 0.0
    assert d["sigma_ring_m87_uas"] > 0.0
