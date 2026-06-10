# -*- coding: utf-8 -*-
"""Anti-fit ledger for the Θ–Ψ canon.

This module records which knobs are allowed in the theory layer and which
belong only to diagnostics. It is deliberately conservative: a quantity can
enter the canon only if it is fixed before data-facing fits or derived from a
smaller set of assumptions.
"""
from __future__ import annotations

from typing import Dict

import sympy as sp


def anti_fit_ledger() -> Dict[str, object]:
    """Return the current anti-fit classification of risky model knobs."""
    canonical_fixed = {
        "primary_ontology": ("Theta_mu_nu", "Psi"),
        "theta_eff_projection": "theta_eff^2 = Theta_mu_nu u^mu u^nu",
        "phi_eff": "ln(theta_eff/theta_0)",
        "weak_field_n_slope": 2,
        "weak_field_n_profile": "asinh",
        "adm_minimal_mixings_from_action": {
            "epsK*": 0,
            "epsG*": 0,
            "alpha_k2*": 0,
        },
        "matter_response_ratio": "n^2 = Z_s/Z_t",
        "epr_correlation_rule": "E(a,b) = -a dot b",
        "particle_definition": "stable localized Theta projection of Psi sector",
    }
    diagnostic_only = {
        "weak_field_n_profile_alternatives": ("linear", "exp", "tanh2"),
        "legacy_explicit_profile_selector": "deprecated scan API; not a canon knob",
        "refractive_family_v1": "F_v1(J; J_star) = 1 - exp(-J/J_star)",
        "legacy_adm_eps_mixings": (
            "epsK_theta_R",
            "epsK_theta_phi",
            "epsK_R_phi",
            "epsG_theta_R",
            "epsG_theta_phi",
            "epsG_R_phi",
            "alpha_k2_theta",
            "alpha_k2_R",
            "alpha_k2_phi",
        ),
        "sigma_extension": "legacy/exploratory perturbation scaffold",
        "strongfield_bridge_calibration": (
            "theta_c_scale",
            "m2_scale",
            "lam_scale",
            "powerlaw2d_coefficients",
        ),
        "inverse_bao_required_multipliers": "diagnostic, not a model patch",
    }
    promotion_requirements = {
        "derive_from_action": sp.Integer(1),
        "fixed_before_fit": sp.Integer(1),
        "no_test_set_selection": sp.Integer(1),
        "no_extra_light_modes": sp.Integer(1),
        "independent_observable_prediction": sp.Integer(1),
    }
    forbidden_moves = {
        "choose_profile_after_data": sp.Integer(1),
        "retune_bridge_on_test": sp.Integer(1),
        "promote_inverse_diagnostic_to_correction": sp.Integer(1),
        "use_adm_eps_without_derivation": sp.Integer(1),
        "add_sigma_as_physical_field_without_new_canon": sp.Integer(1),
    }
    return {
        "canonical_fixed": canonical_fixed,
        "diagnostic_only": diagnostic_only,
        "promotion_requirements": promotion_requirements,
        "forbidden_moves": forbidden_moves,
        "canonical_free_profile_count": sp.Integer(0),
        "canonical_bridge_fit_parameters": sp.Integer(0),
        "data_facing_fit_allowed_after_closure": sp.Integer(1),
        "status": sp.Symbol("anti_fit_guardrails"),
    }


def canonical_promotion_check(
    *,
    derived_from_action: bool,
    fixed_before_fit: bool,
    no_test_set_selection: bool,
    no_extra_light_modes: bool,
    independent_observable_prediction: bool,
) -> Dict[str, object]:
    """Return whether a diagnostic knob may be promoted into the canon."""
    checks = {
        "derive_from_action": bool(derived_from_action),
        "fixed_before_fit": bool(fixed_before_fit),
        "no_test_set_selection": bool(no_test_set_selection),
        "no_extra_light_modes": bool(no_extra_light_modes),
        "independent_observable_prediction": bool(independent_observable_prediction),
    }
    accepted = all(checks.values())
    return {
        "checks": checks,
        "accepted_for_canon": accepted,
        "rejection_reason": sp.Symbol("missing_promotion_requirement")
        if not accepted
        else sp.Symbol("none"),
        "status": sp.Symbol("canonical_promotion_check"),
    }
