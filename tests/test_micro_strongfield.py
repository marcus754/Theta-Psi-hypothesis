import sympy as sp

from theory.micro_strongfield import (
    dynamic_collapse_closure_ledger,
    micro_strongfield_primary_ledger,
    micro_to_strongfield_primary_bridge,
    strongfield_branch_closure_ledger,
)


def test_micro_strongfield_primary_ledger_is_theory_only():
    d = micro_strongfield_primary_ledger()
    assert "S_micro" in d
    assert "R_micro" in d
    assert "J_refr_micro" in d
    assert "n_strong" in d
    assert "z_surface" in d
    assert "delay_proxy" in d
    assert d["ok_no_horizon"] == 1
    assert "compact_no_horizon" in str(d["classification"])
    assert "finite_redshift_finite_delay_no_horizon" in str(d["single_consequence"])
    assert "derived_microphysics_strongfield_theory" in str(d["status"])
    assert "benchmark_point" not in d


def test_micro_strongfield_primary_bridge_returns_single_consequence():
    d = micro_to_strongfield_primary_bridge()
    assert "strongfield_primary_observables" in d
    assert "strongfield_primary_consequence" in d
    assert d["strongfield_primary_keys"] == ("z_surface", "delay_proxy", "ok_no_horizon")
    obs = d["strongfield_primary_observables"]
    assert obs["z_surface"] == d["z_surface"]
    assert obs["delay_proxy"] == d["delay_proxy"]
    assert obs["n_peak"] == d["n_strong"]
    assert obs["ok_no_horizon"] == 1
    assert "J_refr_micro" in d
    assert "n_strong" in d
    consequence = d["strongfield_primary_consequence"]
    assert consequence["finite_redshift"] == 1
    assert consequence["finite_delay"] == 1
    assert consequence["no_horizon"] == 1
    assert "derived_microphysics_strongfield_bridge" in str(d["status"])


def test_strongfield_branch_closure_pins_horizonless_criteria():
    d = strongfield_branch_closure_ledger()
    assert d["status"] == sp.Symbol("strongfield_branch_closed")
    assert d["closure_scope"] == sp.Symbol("stationary_strongfield_branch_only")
    assert d["dynamic_collapse_closed"] == 0
    assert d["closure_chain"] == (
        "local_Theta_Psi_basis",
        "stationary_radial_profile",
        "finite_positive_theta_eff",
        "finite_n_profile",
        "finite_redshift_delay",
        "outgoing_characteristic",
        "forward_observables",
        "falsification_criteria",
    )

    criteria = d["horizonless_criteria"]
    assert criteria["theta_positive"] == 1
    assert criteria["n_finite"] == 1
    assert criteria["finite_redshift"] == 1
    assert criteria["finite_delay"] == 1
    assert criteria["outgoing_characteristic_exists"] == 1
    assert criteria["finite_energy"] == 1
    assert criteria["radial_operator_elliptic"] == 1
    assert criteria["regular_center"] == 1

    bc = d["stationary_boundary_conditions"]
    assert bc["theta_eff_prime_center"] == 0
    assert bc["R_prime_center"] == 0
    assert "theta_0" in str(bc["theta_eff_infinity"])
    assert "R_0" in str(bc["R_infinity"])

    obs = d["forward_observable_map"]
    assert obs["n_profile"] == d["n_strong"]
    assert obs["z_surface"] == d["z_surface"]
    assert obs["delay"] == d["delay_proxy"]
    assert "ray_map_from_n_profile" in str(obs["ring_or_image"])
    assert "transfer_function_from_finite_delay" in str(obs["echo_transfer"])

    falsify = d["falsification_criteria"]
    for key in (
        "theta_eff_reaches_zero",
        "n_diverges",
        "energy_diverges",
        "radial_operator_not_elliptic",
        "no_outgoing_characteristic",
        "forward_observables_incompatible",
    ):
        assert falsify[key] == sp.Symbol("reject_branch")


def test_dynamic_collapse_closure_contract_pins_evolution_criteria():
    d = dynamic_collapse_closure_ledger()
    assert d["status"] == sp.Symbol("dynamic_collapse_closed")
    assert d["closure_scope"] == sp.Symbol("dynamic_collapse_ledger_contract")
    assert d["dynamic_collapse_closed"] == 1
    assert d["requires_numeric_pde_implementation"] == 1
    assert d["closure_chain"] == (
        "finite_energy_initial_data",
        "evolve_Theta_Psi",
        "track_theta_min_and_n_max",
        "preserve_outgoing_characteristic",
        "avoid_blow_up",
        "settle_or_disperse",
        "reject_invalid_outcomes",
    )

    criteria = d["dynamic_criteria"]
    assert criteria["initial_finite_energy"] == 1
    assert criteria["theta_eff_positive_all_t"] == 1
    assert criteria["n_finite_all_t"] == 1
    assert criteria["energy_finite_all_t"] == 1
    assert criteria["no_blow_up_before_branch_decision"] == 1
    assert criteria["outgoing_characteristic_exists_all_t"] == 1
    assert criteria["hyperbolic_evolution_operator"] == 1
    assert criteria["stationary_or_dispersive_outcome"] == 1

    outcomes = d["branch_outcomes"]
    assert outcomes["settle_to_stationary_branch"] == sp.Symbol("allowed")
    assert outcomes["disperse_to_weak_field"] == sp.Symbol("allowed")
    for key in (
        "theta_eff_zero",
        "n_infinite",
        "energy_blow_up",
        "loss_of_hyperbolicity",
        "no_outgoing_characteristic",
    ):
        assert outcomes[key] == sp.Symbol("reject")

    monitored = d["monitored_quantities"]
    assert "min_theta_eff" in str(monitored["min_theta_eff_t"])
    assert "max_n" in str(monitored["max_n_t"])
    assert "E_total" in str(monitored["E_total_t"])
    assert "outgoing_characteristic_exists_t" in str(monitored["outgoing_characteristic_t"])
    assert "hyperbolicity_margin_positive_t" in str(monitored["hyperbolicity_margin_t"])
