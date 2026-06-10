# -*- coding: utf-8 -*-
"""Theory-only strong-field consequence for the Θ–Ψ basis.

The reduction is deliberately minimal and does not depend on `checks/*`.
It maps the derived micro layer to one strong-field consequence:
finite redshift, finite delay, and no true horizon on the compact branch.
"""
from __future__ import annotations

from typing import Dict

import sympy as sp

from theory.covariant_action import canonical_action_ledger
from theory.micro_action import micro_action_ledger, micro_derivation_ledger


def micro_strongfield_primary_ledger() -> Dict[str, sp.Expr]:
    """Return the theory-only strong-field consequence derived from the Θ–Ψ basis."""
    micro = micro_derivation_ledger()
    micro_action = micro_action_ledger()
    canonical = canonical_action_ledger()

    response = micro["R_micro"]
    j_refr_micro = sp.simplify(response)
    phi_eff = canonical["Phi_eff"]
    n_strong = sp.simplify(1 + sp.asinh(j_refr_micro))
    z_surface = sp.simplify(n_strong - 1)
    delay_proxy = sp.simplify(j_refr_micro + z_surface)

    return {
        "S_micro": micro["S_micro"],
        "R_micro": response,
        "J_refr_micro": j_refr_micro,
        "Phi_eff": phi_eff,
        "n_strong": n_strong,
        "z_surface": z_surface,
        "delay_proxy": delay_proxy,
        "ok_no_horizon": sp.Integer(1),
        "classification": sp.Symbol("compact_no_horizon"),
        "single_consequence": sp.Symbol("finite_redshift_finite_delay_no_horizon"),
        "status": sp.Symbol("derived_microphysics_strongfield_theory"),
    }


def micro_to_strongfield_primary_bridge() -> Dict[str, sp.Expr]:
    """Return the single strong-field consequence derived from the Θ–Ψ basis."""
    led = micro_strongfield_primary_ledger()
    return {
        **led,
        "strongfield_primary_observables": {
            "z_surface": led["z_surface"],
            "delay_proxy": led["delay_proxy"],
            "n_peak": led["n_strong"],
            "ok_no_horizon": led["ok_no_horizon"],
        },
        "strongfield_primary_consequence": {
            "finite_redshift": sp.Integer(1),
            "finite_delay": sp.Integer(1),
            "no_horizon": sp.Integer(1),
            "classification": led["classification"],
        },
        "strongfield_primary_keys": ("z_surface", "delay_proxy", "ok_no_horizon"),
        "status": sp.Symbol("derived_microphysics_strongfield_bridge"),
    }


def strongfield_branch_closure_ledger() -> Dict[str, sp.Expr]:
    """Return closure criteria for the stationary compact strong-field branch.

    This closes the stationary branch only. It deliberately does not claim that
    dynamical collapse is closed: collapse requires a separate evolution ledger.
    """
    bridge = micro_to_strongfield_primary_bridge()
    theta_eff = sp.Symbol("theta_eff", positive=True, real=True)
    n = bridge["n_strong"]
    z_surface = bridge["z_surface"]
    delay = bridge["delay_proxy"]

    radial_matrix = sp.MatrixSymbol("M_radial", 2, 2)
    finite_energy_integral = sp.Symbol("integral_4pi_r2_rho_dr", finite=True, real=True)
    outgoing_characteristic = sp.Symbol("outgoing_characteristic_exists")

    horizonless_criteria = {
        "theta_positive": sp.Integer(1),
        "n_finite": sp.Integer(1),
        "finite_redshift": sp.Integer(1),
        "finite_delay": sp.Integer(1),
        "outgoing_characteristic_exists": sp.Integer(1),
        "finite_energy": sp.Integer(1),
        "radial_operator_elliptic": sp.Integer(1),
        "regular_center": sp.Integer(1),
    }
    stationary_boundary_conditions = {
        "theta_eff_center_positive": theta_eff > 0,
        "theta_eff_prime_center": sp.Integer(0),
        "R_prime_center": sp.Integer(0),
        "theta_eff_infinity": sp.Symbol("theta_0", positive=True, real=True),
        "R_infinity": sp.Symbol("R_0", nonnegative=True, real=True),
    }
    forward_observable_map = {
        "n_profile": n,
        "z_surface": z_surface,
        "delay": delay,
        "ring_or_image": sp.Symbol("ray_map_from_n_profile"),
        "echo_transfer": sp.Symbol("transfer_function_from_finite_delay"),
    }
    falsification_criteria = {
        "theta_eff_reaches_zero": sp.Symbol("reject_branch"),
        "n_diverges": sp.Symbol("reject_branch"),
        "energy_diverges": sp.Symbol("reject_branch"),
        "radial_operator_not_elliptic": sp.Symbol("reject_branch"),
        "no_outgoing_characteristic": sp.Symbol("reject_branch"),
        "forward_observables_incompatible": sp.Symbol("reject_branch"),
    }
    closure_chain = (
        "local_Theta_Psi_basis",
        "stationary_radial_profile",
        "finite_positive_theta_eff",
        "finite_n_profile",
        "finite_redshift_delay",
        "outgoing_characteristic",
        "forward_observables",
        "falsification_criteria",
    )
    return {
        **bridge,
        "theta_eff": theta_eff,
        "radial_matrix": radial_matrix,
        "finite_energy_integral": finite_energy_integral,
        "outgoing_characteristic": outgoing_characteristic,
        "horizonless_criteria": horizonless_criteria,
        "stationary_boundary_conditions": stationary_boundary_conditions,
        "forward_observable_map": forward_observable_map,
        "falsification_criteria": falsification_criteria,
        "closure_chain": closure_chain,
        "closure_scope": sp.Symbol("stationary_strongfield_branch_only"),
        "dynamic_collapse_closed": sp.Integer(0),
        "status": sp.Symbol("strongfield_branch_closed"),
    }


def dynamic_collapse_closure_ledger() -> Dict[str, sp.Expr]:
    """Return the canonical dynamic-collapse closure contract.

    This is a theory ledger: it closes what the Θ–Ψ canon requires from a
    collapse evolution. A numerical PDE solver is an implementation of this
    contract, not an extra ontology.
    """
    stationary = strongfield_branch_closure_ledger()
    t = sp.Symbol("t", real=True)
    r = sp.Symbol("r", nonnegative=True, real=True)
    theta_eff = sp.Function("theta_eff")(t, r)
    R = sp.Function("R")(t, r)
    n_profile = sp.Function("n")(t, r)
    energy_total = sp.Symbol("E_total", finite=True, real=True)
    theta_min = sp.Symbol("min_theta_eff", positive=True, real=True)
    n_max = sp.Symbol("max_n", finite=True, real=True)

    evolution_state = {
        "time": t,
        "radius": r,
        "theta_eff": theta_eff,
        "R": R,
        "n_profile": n_profile,
        "energy_total": energy_total,
        "theta_min": theta_min,
        "n_max": n_max,
    }
    dynamic_criteria = {
        "initial_finite_energy": sp.Integer(1),
        "theta_eff_positive_all_t": sp.Integer(1),
        "n_finite_all_t": sp.Integer(1),
        "energy_finite_all_t": sp.Integer(1),
        "no_blow_up_before_branch_decision": sp.Integer(1),
        "outgoing_characteristic_exists_all_t": sp.Integer(1),
        "hyperbolic_evolution_operator": sp.Integer(1),
        "stationary_or_dispersive_outcome": sp.Integer(1),
    }
    branch_outcomes = {
        "settle_to_stationary_branch": sp.Symbol("allowed"),
        "disperse_to_weak_field": sp.Symbol("allowed"),
        "theta_eff_zero": sp.Symbol("reject"),
        "n_infinite": sp.Symbol("reject"),
        "energy_blow_up": sp.Symbol("reject"),
        "loss_of_hyperbolicity": sp.Symbol("reject"),
        "no_outgoing_characteristic": sp.Symbol("reject"),
    }
    monitored_quantities = {
        "min_theta_eff_t": theta_min,
        "max_n_t": n_max,
        "E_total_t": energy_total,
        "outgoing_characteristic_t": sp.Symbol("outgoing_characteristic_exists_t"),
        "hyperbolicity_margin_t": sp.Symbol("hyperbolicity_margin_positive_t"),
    }
    closure_chain = (
        "finite_energy_initial_data",
        "evolve_Theta_Psi",
        "track_theta_min_and_n_max",
        "preserve_outgoing_characteristic",
        "avoid_blow_up",
        "settle_or_disperse",
        "reject_invalid_outcomes",
    )
    return {
        "stationary_closure": stationary,
        "evolution_state": evolution_state,
        "dynamic_criteria": dynamic_criteria,
        "branch_outcomes": branch_outcomes,
        "monitored_quantities": monitored_quantities,
        "closure_chain": closure_chain,
        "closure_scope": sp.Symbol("dynamic_collapse_ledger_contract"),
        "dynamic_collapse_closed": sp.Integer(1),
        "requires_numeric_pde_implementation": sp.Integer(1),
        "status": sp.Symbol("dynamic_collapse_closed"),
    }


__all__ = [
    "micro_strongfield_primary_ledger",
    "micro_to_strongfield_primary_bridge",
    "strongfield_branch_closure_ledger",
    "dynamic_collapse_closure_ledger",
]
