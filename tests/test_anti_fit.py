import sympy as sp

from theory.anti_fit import anti_fit_ledger, canonical_promotion_check


def test_anti_fit_ledger_keeps_diagnostics_out_of_canon():
    d = anti_fit_ledger()
    assert d["status"] == sp.Symbol("anti_fit_guardrails")
    assert d["canonical_free_profile_count"] == 0
    assert d["canonical_bridge_fit_parameters"] == 0
    assert d["canonical_fixed"]["weak_field_n_profile"] == "asinh"
    assert d["diagnostic_only"]["weak_field_n_profile_alternatives"] == (
        "linear",
        "exp",
        "tanh2",
    )
    assert "asinh" not in d["diagnostic_only"]["weak_field_n_profile_alternatives"]
    assert "legacy_explicit_profile_selector" in d["diagnostic_only"]
    assert "strongfield_bridge_calibration" in d["diagnostic_only"]
    assert "inverse_bao_required_multipliers" in d["diagnostic_only"]


def test_anti_fit_ledger_keeps_fixed_physics_explicit():
    d = anti_fit_ledger()
    fixed = d["canonical_fixed"]
    assert fixed["primary_ontology"] == ("Theta_mu_nu", "Psi")
    assert fixed["weak_field_n_slope"] == 2
    assert fixed["weak_field_n_profile"] == "asinh"
    assert fixed["adm_minimal_mixings_from_action"] == {
        "epsK*": 0,
        "epsG*": 0,
        "alpha_k2*": 0,
    }
    assert fixed["matter_response_ratio"] == "n^2 = Z_s/Z_t"
    assert fixed["epr_correlation_rule"] == "E(a,b) = -a dot b"


def test_anti_fit_forbids_common_patching_moves():
    d = anti_fit_ledger()
    forbidden = d["forbidden_moves"]
    assert forbidden["choose_profile_after_data"] == 1
    assert forbidden["retune_bridge_on_test"] == 1
    assert forbidden["promote_inverse_diagnostic_to_correction"] == 1
    assert forbidden["use_adm_eps_without_derivation"] == 1
    assert forbidden["add_sigma_as_physical_field_without_new_canon"] == 1


def test_canonical_promotion_requires_all_closure_conditions():
    rejected = canonical_promotion_check(
        derived_from_action=True,
        fixed_before_fit=True,
        no_test_set_selection=False,
        no_extra_light_modes=True,
        independent_observable_prediction=True,
    )
    assert rejected["accepted_for_canon"] is False
    assert rejected["rejection_reason"] == sp.Symbol("missing_promotion_requirement")

    accepted = canonical_promotion_check(
        derived_from_action=True,
        fixed_before_fit=True,
        no_test_set_selection=True,
        no_extra_light_modes=True,
        independent_observable_prediction=True,
    )
    assert accepted["accepted_for_canon"] is True
    assert accepted["rejection_reason"] == sp.Symbol("none")
