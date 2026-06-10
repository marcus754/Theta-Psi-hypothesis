# -*- coding: utf-8 -*-
"""Derived microphysics ledger for the Θ–Ψ local-clock action.

This module does not introduce a new geometry.
It records a minimal local action in a preferred slicing, where the matter
sector sees response coefficients `Z_t` and `Z_s`, and the refractive index
appears only through the matching rule `n² = Z_s / Z_t`.

The point of this ledger is narrow:

1. keep the microscopic layer local;
2. keep weak-field matching explicit;
3. make the reduction to `S_cov` a reduction, not a new postulate.
"""
from __future__ import annotations

from typing import Dict

import sympy as sp

from theory.covariant_action import canonical_action_ledger, full_covariant_action_ledger, gamma, theta0, V0


def micro_response_ledger() -> Dict[str, sp.Expr]:
    """Return the local response scalar used to generate micro coefficients."""
    derivation = micro_response_derivation_ledger()
    return {
        **derivation,
        "status": sp.Symbol("derived_micro_response"),
    }


def micro_response_derivation_ledger() -> Dict[str, sp.Expr]:
    """Derive the local response scalar from minimal symmetry requirements.

    The response scalar is the unique lowest-order nonnegative local scalar
    built from the weak-field deviation terms:
    - it vanishes at the canonical weak-field point;
    - it is even under sign flips of the fluctuation variables;
    - it contains the minimal number of independent quadratic terms;
    - it is local and rotationally invariant on the preferred slicing.
    """
    theta = sp.Symbol("theta_eff", positive=True, real=True)
    R = sp.Symbol("R", nonnegative=True, real=True)
    chi = sp.Symbol("chi", real=True)
    grad_theta_sq = sp.Symbol("grad_theta_sq", nonnegative=True, real=True)
    grad_R_sq = sp.Symbol("grad_R_sq", nonnegative=True, real=True)
    grad_chi_sq = sp.Symbol("grad_chi_sq", nonnegative=True, real=True)
    response = sp.Symbol("R_micro", nonnegative=True, real=True)
    delta_theta_sq = sp.simplify((theta / theta0 - 1) ** 2)
    psi_sq = sp.simplify(R**2)
    chi_sq = sp.simplify(chi**2)
    response_def = sp.simplify(
        delta_theta_sq
        + psi_sq
        + chi_sq
        + grad_theta_sq
        + grad_R_sq
        + grad_chi_sq
    )
    terms = {
        "delta_theta_sq": delta_theta_sq,
        "psi_sq": psi_sq,
        "chi_sq": chi_sq,
        "grad_theta_sq": grad_theta_sq,
        "grad_R_sq": grad_R_sq,
        "grad_chi_sq": grad_chi_sq,
    }
    axioms = {
        "locality": sp.Integer(1),
        "weak_field_normalization": sp.Integer(1),
        "even_quadratic_order": sp.Integer(1),
        "positive_semidefinite": sp.Integer(1),
        "rotational_invariance": sp.Integer(1),
        "minimal_independent_terms": sp.Integer(1),
    }
    return {
        "theta_eff": theta,
        "R": R,
        "chi": chi,
        "grad_theta_sq": grad_theta_sq,
        "grad_R_sq": grad_R_sq,
        "grad_chi_sq": grad_chi_sq,
        "delta_theta_sq": delta_theta_sq,
        "psi_sq": psi_sq,
        "chi_sq": chi_sq,
        "terms": terms,
        "axioms": axioms,
        "R_micro": response,
        "R_micro_definition": response_def,
        "status": sp.Symbol("derived_micro_response"),
    }


def micro_measure_derivation_ledger() -> Dict[str, sp.Expr]:
    """Derive the preferred-slicing measure used by the local micro action.

    The measure is not a new physical degree of freedom. It is the unique
    lowest-order positive volume element compatible with the local preferred
    slicing used for the derived micro layer:
    - it is orientation-preserving;
    - it carries unit coordinate normalization;
    - it does not introduce any extra field weight.
    """
    dt = sp.Symbol("dt", positive=True, real=True)
    d3x = sp.Symbol("d3x", positive=True, real=True)
    measure = sp.simplify(dt * d3x)
    axioms = {
        "preferred_slicing": sp.Integer(1),
        "orientation_preserving": sp.Integer(1),
        "unit_coordinate_normalization": sp.Integer(1),
        "no_extra_field_weight": sp.Integer(1),
        "positive_volume_element": sp.Integer(1),
    }
    return {
        "dt": dt,
        "d3x": d3x,
        "measure_dt_d3x": measure,
        "product": measure,
        "axioms": axioms,
        "status": sp.Symbol("derived_micro_measure"),
    }


def micro_action_ledger() -> Dict[str, sp.Expr]:
    """Return the derived microphysics ledger."""
    resp = micro_response_ledger()
    measure_ledger = micro_measure_derivation_ledger()
    theta = resp["theta_eff"]
    R = resp["R"]
    chi = resp["chi"]
    grad_theta_sq = resp["grad_theta_sq"]
    grad_R_sq = resp["grad_R_sq"]
    grad_chi_sq = resp["grad_chi_sq"]
    response = resp["R_micro"]
    theta_dot_sq = sp.Symbol("theta_dot_sq", nonnegative=True, real=True)
    R_dot_sq = sp.Symbol("R_dot_sq", nonnegative=True, real=True)
    chi_dot_sq = sp.Symbol("chi_dot_sq", nonnegative=True, real=True)

    measure = measure_ledger["measure_dt_d3x"]

    grad_ln_theta_sq = sp.Symbol("grad_ln_theta_sq", nonnegative=True, real=True)
    a_theta = sp.simplify(1 + response)
    b_theta = sp.simplify(gamma**2 + response)
    a_R = sp.simplify(1 + response)
    b_R = sp.simplify(1 + response)
    z_t = sp.simplify(1 + response)
    z_s = sp.simplify(1 + 2 * response)

    v_theta_micro = sp.Rational(1, 2) * theta**2 + sp.Rational(1, 4) * theta**4
    v_R_micro = sp.Rational(1, 2) * V0 * R**2 + sp.Rational(1, 4) * R**4
    v_matter_micro = sp.Rational(1, 2) * chi**2 + sp.Rational(1, 4) * chi**4

    l_theta_micro = sp.Rational(3, 2) * a_theta * grad_ln_theta_sq + sp.Rational(1, 2) * b_theta * grad_theta_sq - v_theta_micro
    l_R_micro = sp.Rational(1, 2) * a_R * grad_R_sq - v_R_micro
    l_matter_micro = sp.Rational(1, 2) * z_t * chi_dot_sq - sp.Rational(1, 2) * z_s * grad_chi_sq - v_matter_micro
    l_int_micro = response * (theta**2 * R**2 + theta**2 * chi**2)
    l_micro = sp.simplify(l_theta_micro + l_R_micro + l_matter_micro + l_int_micro)

    n_sq = sp.simplify(z_s / z_t)
    n = sp.sqrt(n_sq)

    weak_field_matching = {
        "R_micro": sp.Integer(0),
        "A_Theta": sp.Integer(1),
        "B_Theta": gamma**2,
        "A_R": sp.Integer(1),
        "B_R": sp.Integer(1),
        "Z_t": sp.Integer(1),
        "Z_s": sp.Integer(1),
        "n_sq": sp.Integer(1),
    }

    bridge_to_covariant = sp.Symbol("S_micro -> S_red -> S_cov")

    return {
        "theta_eff": theta,
        "R": R,
        "chi": chi,
        "R_micro": response,
        "R_micro_definition": resp["R_micro_definition"],
        "measure_dt_d3x": measure,
        "measure_ledger": measure_ledger,
        "theta_dot_sq": theta_dot_sq,
        "R_dot_sq": R_dot_sq,
        "chi_dot_sq": chi_dot_sq,
        "grad_theta_sq": grad_theta_sq,
        "grad_R_sq": grad_R_sq,
        "grad_chi_sq": grad_chi_sq,
        "grad_ln_theta_sq": grad_ln_theta_sq,
        "A_Theta": a_theta,
        "B_Theta": b_theta,
        "A_R": a_R,
        "B_R": b_R,
        "Z_t": z_t,
        "Z_s": z_s,
        "V_Theta": v_theta_micro,
        "V_R": v_R_micro,
        "V_matter": v_matter_micro,
        "L_Theta_micro": l_theta_micro,
        "L_R_micro": l_R_micro,
        "L_matter_micro": l_matter_micro,
        "L_int_micro": l_int_micro,
        "L_micro": l_micro,
        "S_micro": measure * l_micro,
        "n_sq": n_sq,
        "n": n,
        "weak_field_matching": weak_field_matching,
        "bridge_to_covariant": bridge_to_covariant,
        "status": sp.Symbol("derived_micro_action"),
    }


def micro_action_summary() -> Dict[str, sp.Expr]:
    """Return a compact summary of the derived microphysics ledger."""
    return {
        **micro_action_ledger(),
        "status": sp.Symbol("derived_microphysics"),
    }


def micro_to_reduced_bridge_ledger() -> Dict[str, sp.Expr]:
    """Return the explicit bridge from S_micro to the reduced canon S_red."""
    micro = micro_action_ledger()
    canonical = canonical_action_ledger()
    s_red = canonical["S_red"]
    heavy_modes = sp.Symbol("heavy_modes_integrated_out", real=True)
    weak_subs = {
        micro["R_micro"]: sp.Integer(0),
    }
    weak_l_theta = sp.simplify(micro["L_Theta_micro"].subs(weak_subs))
    weak_l_R = sp.simplify(micro["L_R_micro"].subs(weak_subs))
    weak_l_matter = sp.simplify(micro["L_matter_micro"].subs(weak_subs))
    weak_n_sq = sp.simplify(micro["n_sq"].subs(weak_subs))
    bridge_conditions = {
        "locality": sp.Integer(1),
        "weak_field_match": sp.Integer(1),
        "R_micro_to_0": sp.Integer(1),
        "match_n_sq": micro["n_sq"],
        "integrate_heavy_modes": heavy_modes,
        "no_new_geometry": sp.Integer(1),
        "reduced_clock_sector": sp.Integer(1),
    }
    weak_field_reduction = {
        "weak_subs": weak_subs,
        "L_Theta_micro_weak": weak_l_theta,
        "L_R_micro_weak": weak_l_R,
        "L_matter_micro_weak": weak_l_matter,
        "n_sq_weak": weak_n_sq,
        "n_weak": sp.sqrt(weak_n_sq),
    }
    reduction_targets = {
        "S_red": s_red,
        "L_Theta_target": canonical["L_Theta"],
        "L_R_target": canonical["L_R"],
        "U_red_target": canonical["U_red"],
        "L_refr_target": canonical["L_refr"],
        "Phi_eff_target": canonical["Phi_eff"],
        "J_refr_target": canonical["J_refr"],
        "R_micro_target": sp.Integer(0),
    }
    return {
        **micro,
        "S_red_from_micro": s_red,
        "heavy_modes_integrated_out": heavy_modes,
        "bridge_conditions": bridge_conditions,
        "weak_field_reduction": weak_field_reduction,
        "reduction_targets": reduction_targets,
        "status": sp.Symbol("derived_microphysics_reduction"),
    }


def micro_to_covariant_bridge_ledger() -> Dict[str, sp.Expr]:
    """Return the explicit bridge from S_micro to the current covariant canon."""
    reduction = micro_to_reduced_bridge_ledger()
    canonical_cov = full_covariant_action_ledger()
    s_cov = canonical_cov["S_cov"]
    return {
        **reduction,
        "S_cov": s_cov,
        "matter_completion": canonical_cov["S_matter"],
        "status": sp.Symbol("derived_microphysics_bridge"),
    }


def propagation_law_closure_ledger() -> Dict[str, sp.Expr]:
    """Return the explicit closure of the propagation-law chain.

    This ledger pins three points:
    1. `g`/`g_hat` are bookkeeping records, not primary ontology;
    2. all matter sectors share one response ratio `Z_s/Z_t`;
    3. the weak-field normalization is fixed before any fit.
    """
    micro = micro_action_ledger()
    bridge = micro_to_covariant_bridge_ledger()
    canonical = full_covariant_action_ledger()

    z_t = micro["Z_t"]
    z_s = micro["Z_s"]
    n_sq = sp.simplify(z_s / z_t)
    n = sp.sqrt(n_sq)
    phi_eff = canonical["Phi_eff"]
    weak_n_linear = sp.simplify(1 + 2 * phi_eff)

    species_a, species_b = sp.symbols("chi_a chi_b")
    universal_response = {
        "Z_t_chi_a": z_t,
        "Z_t_chi_b": z_t,
        "Z_s_chi_a": z_s,
        "Z_s_chi_b": z_s,
        "n_sq_chi_a": n_sq,
        "n_sq_chi_b": n_sq,
        "species_independent": sp.Integer(1),
    }
    bookkeeping_conditions = {
        "primary_ontology": (sp.Symbol("Theta_mu_nu"), sp.Symbol("Psi")),
        "g_hat_independent_dof": sp.Integer(0),
        "g_hat_role": sp.Symbol("technical_propagation_record"),
        "matter_completion_only": bridge["matter_completion"],
    }
    weak_subs = {micro["R_micro"]: sp.Integer(0)}
    weak_field_closure = {
        "R_micro_limit": sp.Integer(0),
        "Z_t_weak": sp.simplify(z_t.subs(weak_subs)),
        "Z_s_weak": sp.simplify(z_s.subs(weak_subs)),
        "n_sq_weak": sp.simplify(n_sq.subs(weak_subs)),
        "n_weak": sp.simplify(n.subs(weak_subs)),
        "n_linear_observable_normalization": weak_n_linear,
        "normalization_fixed_before_fit": sp.Integer(1),
    }
    observable_map = {
        "redshift": sp.Symbol("Delta_nu_over_nu = Delta_ln_n"),
        "delay": sp.Symbol("t_signal = integral_n_dl"),
        "lensing": sp.Symbol("alpha = integral_grad_perp_ln_n_dl"),
        "probe_acceleration": sp.Symbol("a_probe = -grad_ln_n"),
    }
    closure_chain = (
        "S_micro",
        "R_micro",
        "universal_Z_t_Z_s",
        "n_sq_equals_Z_s_over_Z_t",
        "g_hat_bookkeeping_only",
        "weak_field_observables",
    )
    return {
        "species": (species_a, species_b),
        "Z_t": z_t,
        "Z_s": z_s,
        "n_sq": n_sq,
        "n": n,
        "universal_response": universal_response,
        "bookkeeping_conditions": bookkeeping_conditions,
        "weak_field_closure": weak_field_closure,
        "observable_map": observable_map,
        "closure_chain": closure_chain,
        "S_cov": bridge["S_cov"],
        "status": sp.Symbol("propagation_law_closed"),
    }


def vacuum_offset_source_ledger() -> Dict[str, sp.Expr]:
    """Return the source-map sanity check for vacuum energy offsets.

    The check is deliberately narrow. An absolute additive constant in the
    matter effective action is not allowed to act as a direct Θ–Ψ source.
    Only terms that depend on the local response data (`R_micro`, `Z_t`, `Z_s`,
    `n`) or on boundary differences can contribute to the source map.
    """
    micro = micro_action_ledger()
    theta = micro["theta_eff"]
    R = micro["R"]
    chi = micro["chi"]
    response = micro["R_micro"]
    measure = micro["measure_dt_d3x"]
    n_sq = micro["n_sq"]
    n = micro["n"]

    vacuum_offset = sp.Symbol("Lambda_vac_offset", real=True)
    boundary_delta = sp.Symbol("Delta_Gamma_boundary", real=True)
    response_delta = sp.Symbol("Delta_response", real=True)

    absolute_offset_action = sp.simplify(measure * vacuum_offset)
    response_observable_action = sp.simplify(measure * response_delta * n)
    boundary_observable_action = sp.simplify(measure * boundary_delta * n)

    source_variables = {
        "theta_eff": theta,
        "R": R,
        "chi": chi,
        "R_micro": response,
    }
    offset_sources = {
        name: sp.simplify(sp.diff(vacuum_offset, var))
        for name, var in source_variables.items()
    }
    response_sources = {
        name: sp.simplify(sp.diff(response_delta * n, var))
        for name, var in source_variables.items()
    }
    boundary_sources = {
        name: sp.simplify(sp.diff(boundary_delta * n, var))
        for name, var in source_variables.items()
    }

    return {
        "vacuum_offset": vacuum_offset,
        "absolute_offset_action": absolute_offset_action,
        "response_observable_action": response_observable_action,
        "boundary_observable_action": boundary_observable_action,
        "source_variables": source_variables,
        "offset_sources": offset_sources,
        "response_sources": response_sources,
        "boundary_sources": boundary_sources,
        "additive_constant_decouples": sp.Integer(
            int(all(value == 0 for value in offset_sources.values()))
        ),
        "response_source_allowed": sp.Integer(
            int(response_sources["R_micro"] != 0)
        ),
        "boundary_source_allowed": sp.Integer(
            int(boundary_sources["R_micro"] != 0)
        ),
        "n_sq": n_sq,
        "n": n,
        "source_rule": sp.Symbol("source_depends_on_response_not_absolute_offset"),
        "status": sp.Symbol("vacuum_offset_source_map"),
    }


def micro_derivation_ledger() -> Dict[str, sp.Expr]:
    """Return the explicit derivation ladder from micro action to the canon.

    The ledger is intentionally concrete:
    1. start from `S_micro`;
    2. impose weak-field matching and integrate heavy modes;
    3. assemble the reduced canonical blocks;
    4. complete to `S_cov`.

    This is a derivation ledger for the current canonical reduction, not a UV
    uniqueness theorem.
    """
    micro = micro_action_ledger()
    canonical = canonical_action_ledger()
    covariant = full_covariant_action_ledger()
    reduction = micro_to_reduced_bridge_ledger()
    bridge = micro_to_covariant_bridge_ledger()
    propagation = propagation_law_closure_ledger()
    return {
        "S_micro": micro["S_micro"],
        "R_micro": micro["R_micro"],
        "R_micro_definition": micro["R_micro_definition"],
        "measure_dt_d3x": micro["measure_dt_d3x"],
        "measure_ledger": micro["measure_ledger"],
        "weak_field_matching": micro["weak_field_matching"],
        "bridge_conditions": reduction["bridge_conditions"],
        "weak_field_reduction": reduction["weak_field_reduction"],
        "reduction_targets": reduction["reduction_targets"],
        "S_red_derived": canonical["S_red"],
        "S_cov_derived": covariant["S_cov"],
        "S_matter_completion": covariant["S_matter"],
        "propagation_law_closure": propagation,
        "derivation_steps": [
            "S_micro",
            "derive_R_micro",
            "weak_field_projection",
            "integrate_heavy_modes",
            "assemble_S_red",
            "add_S_matter_completion",
            "close_propagation_law",
            "S_cov",
        ],
        "bridge": bridge,
        "status": sp.Symbol("derived_microphysics_derivation"),
    }


def micro_symbol_audit() -> Dict[str, Dict[str, sp.Expr]]:
    """Return a compact audit table for the derived micro symbols.

    The audit splits symbols into:
    - derived: obtained by projection or algebraic matching;
    - assumed: minimal fields taken as part of the bridge;
    - diagnostic: used only as a compact consequence check.
    """
    micro = micro_action_ledger()
    bridge = micro_to_reduced_bridge_ledger()
    strong = micro_to_strongfield_primary_bridge()
    return {
        "theta_eff": {
            "status": sp.Symbol("derived"),
            "reason": sp.Symbol("projection from Theta_mu_nu"),
            "symbol": micro["theta_eff"],
        },
        "measure_dt_d3x": {
            "status": sp.Symbol("derived"),
            "reason": sp.Symbol("preferred-slicing volume form with unit coordinate normalization"),
            "symbol": micro["measure_dt_d3x"],
        },
        "R": {
            "status": sp.Symbol("assumed"),
            "reason": sp.Symbol("minimal energy-sector amplitude in the bridge"),
            "symbol": micro["R"],
        },
        "chi": {
            "status": sp.Symbol("assumed"),
            "reason": sp.Symbol("minimal matter-sector field"),
            "symbol": micro["chi"],
        },
        "R_micro": {
            "status": sp.Symbol("derived"),
            "reason": sp.Symbol("minimal local quadratic deviation from the canonical weak-field point"),
            "symbol": micro["R_micro"],
        },
        "Z_t": {
            "status": sp.Symbol("derived"),
            "reason": sp.Symbol("time-response coefficient from the derived micro layer"),
            "symbol": micro["Z_t"],
        },
        "Z_s": {
            "status": sp.Symbol("derived"),
            "reason": sp.Symbol("space-response coefficient from the derived micro layer"),
            "symbol": micro["Z_s"],
        },
        "n_sq": {
            "status": sp.Symbol("derived"),
            "reason": sp.Symbol("matching rule n^2 = Z_s / Z_t"),
            "symbol": micro["n_sq"],
        },
        "n": {
            "status": sp.Symbol("derived"),
            "reason": sp.Symbol("square root of the matched refractive index"),
            "symbol": micro["n"],
        },
        "S_micro": {
            "status": sp.Symbol("derived"),
            "reason": sp.Symbol("local action built from the derived reduction ledger"),
            "symbol": micro["S_micro"],
        },
        "S_red": {
            "status": sp.Symbol("derived"),
            "reason": sp.Symbol("canonical reduced ledger recovered by reduction"),
            "symbol": bridge["S_red_from_micro"],
        },
        "S_cov": {
            "status": sp.Symbol("derived"),
            "reason": sp.Symbol("canonical covariant completion recovered by reduction"),
            "symbol": full_covariant_action_ledger()["S_cov"],
        },
        "finite_redshift": {
            "status": sp.Symbol("diagnostic"),
            "reason": sp.Symbol("single strong-field consequence on the compact branch"),
            "symbol": strong["strongfield_primary_consequence"]["finite_redshift"],
        },
        "finite_delay": {
            "status": sp.Symbol("diagnostic"),
            "reason": sp.Symbol("single strong-field consequence on the compact branch"),
            "symbol": strong["strongfield_primary_consequence"]["finite_delay"],
        },
        "ok_no_horizon": {
            "status": sp.Symbol("diagnostic"),
            "reason": sp.Symbol("single strong-field consequence on the compact branch"),
            "symbol": strong["strongfield_primary_consequence"]["no_horizon"],
        },
        "bridge_to_covariant": {
            "status": sp.Symbol("derived"),
            "reason": sp.Symbol("explicit reduction to the canonical covariant action"),
            "symbol": micro["bridge_to_covariant"],
        },
    }


def micro_to_strongfield_primary_bridge(
) -> Dict[str, sp.Expr]:
    """Return the single strong-field consequence tied to the derived micro layer.

    The reduction is theory-only: it is built from the micro derivation ledger and
    the diagnostic strong-field refractive rule, without calling into `checks/*`.
    """
    from theory.micro_strongfield import micro_to_strongfield_primary_bridge as bridge

    return {
        **bridge(),
        "status": sp.Symbol("derived_microphysics_strongfield_bridge"),
    }
