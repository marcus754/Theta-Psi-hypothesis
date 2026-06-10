import inspect
from math import sqrt

import pytest

from theory.epr_algebroid import (
    chsh_value,
    experimental_profile_from_angles,
    epr_algebroid_ledger,
    foliation_order_probabilities,
    ghz_mermin_ledger,
    marginal_a,
    marginal_b,
    optimal_chsh_directions,
    pr_box_ledger,
    singlet_correlation,
    singlet_joint_probability,
    tsirelson_gate,
)


def test_singlet_correlation_rule_is_fixed_dot_product():
    assert singlet_correlation((1, 0, 0), (1, 0, 0)) == pytest.approx(-1.0)
    assert singlet_correlation((1, 0, 0), (0, 1, 0)) == pytest.approx(0.0)
    assert singlet_correlation((0, 1, 0), (0, -2, 0)) == pytest.approx(1.0)


def test_singlet_chsh_violates_bell_and_hits_tsirelson():
    dirs = optimal_chsh_directions()
    s = chsh_value(dirs["a0"], dirs["a1"], dirs["b0"], dirs["b1"])
    assert abs(s) > 2.0
    assert abs(s) == pytest.approx(2.0 * sqrt(2.0))


def test_joint_probabilities_are_normalized_and_nonnegative():
    a = (1, 0, 0)
    b = (1, 1, 0)
    probs = [singlet_joint_probability(oa, ob, a, b) for oa in (-1, 1) for ob in (-1, 1)]
    assert sum(probs) == pytest.approx(1.0)
    assert all(p >= 0.0 for p in probs)


def test_no_signalling_marginals_do_not_depend_on_remote_setting():
    a = (1, 0, 0)
    b0 = (1, 1, 0)
    b1 = (-1, 2, 0)
    assert marginal_a(1, a, b0) == pytest.approx(0.5)
    assert marginal_a(1, a, b1) == pytest.approx(0.5)

    b = (0, 1, 0)
    a0 = (1, 1, 0)
    a1 = (2, -1, 0)
    assert marginal_b(-1, a0, b) == pytest.approx(0.5)
    assert marginal_b(-1, a1, b) == pytest.approx(0.5)


def test_anchor_maps_correlations_to_zero_signal():
    d = epr_algebroid_ledger()
    anchor = d["anchor"]
    assert anchor["rho_correlation"] == 0
    assert "v_signal_sq <= 1" in str(anchor["causal_signal_condition"])
    assert d["interpretation"].name == "Psi_nonseparable_Theta_causal_anchor"


def test_no_free_correlation_function_or_fit_parameters():
    d = epr_algebroid_ledger()
    assert d["free_correlation_function"] == 0
    assert d["free_fit_parameters"] == 0
    assert d["correlation_rule"].name == "E(a,b) = - dot(a,b)"

    signature = inspect.signature(singlet_correlation)
    assert tuple(signature.parameters) == ("a", "b")


def test_pr_box_is_rejected_by_tsirelson_gate():
    d = pr_box_ledger()
    assert d["no_signalling"] == 1
    assert d["S_CHSH"] == pytest.approx(4.0)
    assert d["gate"]["reject_super_tsirelson"] is True
    assert d["accepted_by_theta_psi"] == 0
    assert "super_tsirelson" in d["rejection_reason"].name


def test_no_signalling_alone_is_not_sufficient():
    d = pr_box_ledger()
    assert d["no_signalling"] == 1
    assert d["gate"]["within_quantum_bound"] is False
    assert d["accepted_by_theta_psi"] == 0


def test_tsirelson_gate_accepts_quantum_and_rejects_superquantum():
    quantum = tsirelson_gate(2.0 * sqrt(2.0))
    assert quantum["violates_bell"] is True
    assert quantum["within_quantum_bound"] is True
    assert quantum["reject_super_tsirelson"] is False

    superquantum = tsirelson_gate(3.0)
    assert superquantum["violates_bell"] is True
    assert superquantum["within_quantum_bound"] is False
    assert superquantum["reject_super_tsirelson"] is True


def test_foliation_order_independence_for_singlet():
    probs = foliation_order_probabilities((1, 0, 0), (1, 1, 0))
    assert probs["A_before_B"] == probs["B_before_A"]
    assert sum(probs["A_before_B"].values()) == pytest.approx(1.0)


def test_ghz_mermin_rejects_global_hidden_valuation():
    d = ghz_mermin_ledger()
    assert d["quantum_product"] == -1
    assert d["noncontextual_hidden_valuation_product"] == 1
    assert d["global_hidden_valuation_exists"] == 0
    assert d["contextuality_required"] == 1


def test_no_fit_angle_profile_matches_optimal_chsh():
    import math

    d = experimental_profile_from_angles(
        0.0,
        math.pi / 2.0,
        math.pi / 4.0,
        -math.pi / 4.0,
    )
    assert d["free_fit_parameters"] == 0
    assert abs(d["S_CHSH"]) == pytest.approx(2.0 * sqrt(2.0))
    assert d["gate"]["within_quantum_bound"] is True


def test_quantum_layer_bridge_uses_single_psi_state():
    d = epr_algebroid_ledger()
    bridge = d["quantum_layer_bridge"]
    assert bridge["shared_object"].name == "Psi_state"
    assert bridge["bridge_chain"] == (
        "Psi_state_layer",
        "two_subsystem_entanglement",
        "singlet_projection",
        "fixed_CHSH_correlations",
        "Theta_anchor_no_signalling",
    )
