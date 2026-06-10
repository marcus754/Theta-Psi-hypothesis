import sympy as sp

from theory.covariant_action import canonical_action_ledger, full_covariant_action_ledger
from theory.micro_action import (
    micro_action_ledger,
    micro_action_summary,
    micro_derivation_ledger,
    micro_measure_derivation_ledger,
    propagation_law_closure_ledger,
    micro_response_derivation_ledger,
    micro_symbol_audit,
    micro_to_reduced_bridge_ledger,
    micro_to_covariant_bridge_ledger,
    micro_to_strongfield_primary_bridge,
    vacuum_offset_source_ledger,
)


def test_micro_action_ledger_contains_local_blocks():
    d = micro_action_ledger()
    for key in (
        "theta_eff",
        "chi",
        "R_micro",
        "R_micro_definition",
        "A_Theta",
        "B_Theta",
        "A_R",
        "B_R",
        "Z_t",
        "Z_s",
        "grad_ln_theta_sq",
        "L_Theta_micro",
        "L_R_micro",
        "L_matter_micro",
        "L_int_micro",
        "L_micro",
        "S_micro",
        "n_sq",
        "n",
        "weak_field_matching",
        "bridge_to_covariant",
    ):
        assert key in d
    assert "R_micro" in str(d["A_Theta"])
    assert "R_micro" in str(d["n_sq"])
    assert "Function(" not in str(d["V_Theta"])
    assert "Function(" not in str(d["V_R"])
    assert "Function(" not in str(d["V_matter"])
    assert "g_hat" not in str(d["S_micro"])


def test_micro_response_derivation_ledger_is_minimal_and_positive():
    d = micro_response_derivation_ledger()
    assert "axioms" in d
    assert "terms" in d
    assert d["axioms"]["locality"] == 1
    assert d["axioms"]["weak_field_normalization"] == 1
    assert d["axioms"]["even_quadratic_order"] == 1
    assert d["axioms"]["positive_semidefinite"] == 1
    assert d["axioms"]["rotational_invariance"] == 1
    assert d["axioms"]["minimal_independent_terms"] == 1
    terms = d["terms"]
    assert "delta_theta_sq" in terms
    assert "psi_sq" in terms
    assert "chi_sq" in terms
    assert "grad_theta_sq" in terms
    assert "grad_R_sq" in terms
    assert "grad_chi_sq" in terms
    assert "Function(" not in str(d["R_micro_definition"])
    assert "theta_eff" in str(d["R_micro_definition"])
    assert "R" in str(d["R_micro_definition"])
    assert "chi" in str(d["R_micro_definition"])


def test_micro_measure_derivation_ledger_is_derived_normalization():
    d = micro_measure_derivation_ledger()
    assert "dt" in d
    assert "d3x" in d
    assert "measure_dt_d3x" in d
    assert "axioms" in d
    assert d["axioms"]["preferred_slicing"] == 1
    assert d["axioms"]["orientation_preserving"] == 1
    assert d["axioms"]["unit_coordinate_normalization"] == 1
    assert d["axioms"]["no_extra_field_weight"] == 1
    assert d["axioms"]["positive_volume_element"] == 1
    assert d["measure_dt_d3x"] == d["product"]
    assert "dt" in str(d["product"])
    assert "d3x" in str(d["product"])
    assert "derived_micro_measure" in str(d["status"])


def test_micro_action_weak_field_matching_is_unit():
    d = micro_action_ledger()
    matching = d["weak_field_matching"]
    assert matching["R_micro"] == sp.Integer(0)
    assert matching["Z_t"] == 1
    assert matching["Z_s"] == 1
    assert matching["A_Theta"] == 1
    assert "gamma" in str(matching["B_Theta"])
    assert "R_micro" not in str(matching["B_Theta"])


def test_micro_action_summary_is_provisional_and_bridges_to_covariant():
    d = micro_action_summary()
    assert "S_micro" in d
    assert "derived_microphysics" in str(d["status"])
    assert "S_micro" in str(d["bridge_to_covariant"])
    assert "S_cov" in str(d["bridge_to_covariant"])


def test_micro_to_reduced_bridge_ledger_has_explicit_reduction():
    d = micro_to_reduced_bridge_ledger()
    canonical = canonical_action_ledger()
    assert "S_micro" in d
    assert "S_red_from_micro" in d
    assert "heavy_modes_integrated_out" in d
    assert "bridge_conditions" in d
    assert "weak_field_reduction" in d
    assert "reduction_targets" in d
    assert "R_micro_to_0" in d["bridge_conditions"]
    assert "R_micro" in str(d["bridge_conditions"]["match_n_sq"])
    weak = d["weak_field_reduction"]
    assert "grad_ln_theta_sq" in str(weak["L_Theta_micro_weak"])
    assert "grad_R_sq" in str(weak["L_R_micro_weak"])
    assert "chi_dot_sq" in str(weak["L_matter_micro_weak"])
    assert weak["n_sq_weak"] == 1
    assert weak["n_weak"] == 1
    targets = d["reduction_targets"]
    assert d["S_red_from_micro"] == canonical["S_red"]
    assert targets["S_red"] == canonical["S_red"]
    assert targets["L_Theta_target"] == canonical["L_Theta"]
    assert targets["L_R_target"] == canonical["L_R"]
    assert targets["U_red_target"] == canonical["U_red"]
    assert targets["L_refr_target"] == canonical["L_refr"]
    assert targets["R_micro_target"] == 0


def test_micro_to_covariant_bridge_ledger_adds_completion():
    d = micro_to_covariant_bridge_ledger()
    canonical_cov = full_covariant_action_ledger()
    assert "S_micro" in d
    assert "S_red_from_micro" in d
    assert "S_cov" in d
    assert "matter_completion" in d
    assert "weak_field_reduction" in d
    assert d["S_cov"] == canonical_cov["S_cov"]
    assert d["matter_completion"] == canonical_cov["S_matter"]


def test_micro_derivation_ledger_orders_the_reduction_steps():
    d = micro_derivation_ledger()
    canonical = canonical_action_ledger()
    covariant = full_covariant_action_ledger()
    assert "R_micro" in d
    assert "R_micro_definition" in d
    assert d["S_red_derived"] == canonical["S_red"]
    assert d["S_cov_derived"] == covariant["S_cov"]
    assert d["S_matter_completion"] == covariant["S_matter"]
    assert "propagation_law_closure" in d
    assert d["propagation_law_closure"]["status"] == sp.Symbol("propagation_law_closed")
    assert d["derivation_steps"] == [
        "S_micro",
        "derive_R_micro",
        "weak_field_projection",
        "integrate_heavy_modes",
        "assemble_S_red",
        "add_S_matter_completion",
        "close_propagation_law",
        "S_cov",
    ]
    assert "derived_microphysics_derivation" in str(d["status"])


def test_propagation_law_closure_pins_points_1_2_3():
    d = propagation_law_closure_ledger()
    assert d["status"] == sp.Symbol("propagation_law_closed")
    assert d["closure_chain"] == (
        "S_micro",
        "R_micro",
        "universal_Z_t_Z_s",
        "n_sq_equals_Z_s_over_Z_t",
        "g_hat_bookkeeping_only",
        "weak_field_observables",
    )

    bookkeeping = d["bookkeeping_conditions"]
    assert bookkeeping["primary_ontology"] == (
        sp.Symbol("Theta_mu_nu"),
        sp.Symbol("Psi"),
    )
    assert bookkeeping["g_hat_independent_dof"] == 0
    assert bookkeeping["g_hat_role"] == sp.Symbol("technical_propagation_record")

    universal = d["universal_response"]
    assert universal["species_independent"] == 1
    assert universal["Z_t_chi_a"] == universal["Z_t_chi_b"]
    assert universal["Z_s_chi_a"] == universal["Z_s_chi_b"]
    assert universal["n_sq_chi_a"] == universal["n_sq_chi_b"]
    assert universal["n_sq_chi_a"] == d["n_sq"]

    weak = d["weak_field_closure"]
    assert weak["Z_t_weak"] == 1
    assert weak["Z_s_weak"] == 1
    assert weak["n_sq_weak"] == 1
    assert weak["n_weak"] == 1
    assert weak["normalization_fixed_before_fit"] == 1
    assert "Phi_eff" not in str(weak["n_linear_observable_normalization"])
    assert "log" in str(weak["n_linear_observable_normalization"])

    obs = d["observable_map"]
    assert "Delta_ln_n" in str(obs["redshift"])
    assert "integral_n_dl" in str(obs["delay"])
    assert "grad_perp_ln_n" in str(obs["lensing"])
    assert "grad_ln_n" in str(obs["probe_acceleration"])


def test_vacuum_offset_source_ledger_decouples_absolute_offset():
    d = vacuum_offset_source_ledger()
    assert d["status"] == sp.Symbol("vacuum_offset_source_map")
    assert d["additive_constant_decouples"] == 1
    assert d["source_rule"] == sp.Symbol("source_depends_on_response_not_absolute_offset")
    assert set(d["offset_sources"]) == {"theta_eff", "R", "chi", "R_micro"}
    assert all(value == 0 for value in d["offset_sources"].values())
    assert "Lambda_vac_offset" in str(d["absolute_offset_action"])


def test_vacuum_offset_source_ledger_keeps_response_and_boundary_sources():
    d = vacuum_offset_source_ledger()
    assert d["response_source_allowed"] == 1
    assert d["boundary_source_allowed"] == 1
    assert d["response_sources"]["R_micro"] != 0
    assert d["boundary_sources"]["R_micro"] != 0
    assert "R_micro" in str(d["n"])
    assert "Delta_response" in str(d["response_observable_action"])
    assert "Delta_Gamma_boundary" in str(d["boundary_observable_action"])


def test_micro_to_strongfield_primary_bridge_is_single_consequence():
    d = micro_to_strongfield_primary_bridge()
    assert "strongfield_primary_observables" in d
    assert "strongfield_primary_consequence" in d
    assert "strongfield_primary_keys" in d
    assert d["strongfield_primary_keys"] == ("z_surface", "delay_proxy", "ok_no_horizon")
    obs = d["strongfield_primary_observables"]
    assert "R_micro" in str(obs["z_surface"])
    assert "R_micro" in str(obs["delay_proxy"])
    assert "R_micro" in str(obs["n_peak"])
    assert obs["ok_no_horizon"] == 1
    assert "J_refr_micro" in d
    assert "n_strong" in d
    assert "compact_no_horizon" in str(d["classification"])
    consequence = d["strongfield_primary_consequence"]
    assert consequence["finite_redshift"] == 1
    assert consequence["finite_delay"] == 1
    assert consequence["no_horizon"] == 1
    assert "derived_microphysics_strongfield_bridge" in str(d["status"])


def test_micro_symbol_audit_classifies_symbols():
    audit = micro_symbol_audit()
    assert audit["theta_eff"]["status"] == sp.Symbol("derived")
    assert audit["R"]["status"] == sp.Symbol("assumed")
    assert audit["chi"]["status"] == sp.Symbol("assumed")
    assert audit["R_micro"]["status"] == sp.Symbol("derived")
    assert audit["n_sq"]["status"] == sp.Symbol("derived")
    assert audit["S_red"]["status"] == sp.Symbol("derived")
    assert audit["S_cov"]["status"] == sp.Symbol("derived")
    assert audit["finite_redshift"]["status"] == sp.Symbol("diagnostic")
    assert audit["finite_delay"]["status"] == sp.Symbol("diagnostic")
    assert audit["ok_no_horizon"]["status"] == sp.Symbol("diagnostic")
    assert "projection" in str(audit["theta_eff"]["reason"])
    assert "minimal local quadratic deviation" in str(audit["R_micro"]["reason"])
    assert audit["measure_dt_d3x"]["status"] == sp.Symbol("derived")
    assert "preferred-slicing volume form" in str(audit["measure_dt_d3x"]["reason"])
