# -*- coding: utf-8 -*-
"""Derived microphysics perturbation ledger for the Θ–Ψ layer.

This module checks the microscopic action in a preferred slicing and records
the minimal conditions under which the derived micro layer reduces to the
current canonical reduced sector without introducing new propagating modes.
"""
from __future__ import annotations

from typing import Dict, Any

import sympy as sp

from theory.micro_action import micro_to_reduced_bridge_ledger


def micro_perturbation_ledger(
    *,
    a_theta: float = 1.0,
    b_theta: float = 1.0,
    a_psi: float = 1.0,
    b_psi: float = 1.0,
    z_t: float = 1.0,
    z_s: float = 1.0,
    lambda_tpsi: float = 0.0,
    lambda_tchi: float = 0.0,
) -> Dict[str, Any]:
    """Return a local perturbation ledger for the micro action.

    The ledger is intentionally minimal:
    - positive kinetic coefficients are required;
    - matter response coefficients must stay positive;
    - weak-field matching is encoded as a bridge to the reduced canon;
    - no new geometry is introduced.
    """
    bridge = micro_to_reduced_bridge_ledger()
    coeffs = {
        "A_Theta": float(a_theta),
        "B_Theta": float(b_theta),
        "A_Psi": float(a_psi),
        "B_Psi": float(b_psi),
        "Z_t": float(z_t),
        "Z_s": float(z_s),
    }

    kinetic_positive = all(v > 0.0 for v in coeffs.values())
    matter_response_positive = float(z_t) > 0.0 and float(z_s) > 0.0
    weak_field_match = (
        float(a_theta) == 1.0
        and float(b_theta) == 1.0
        and float(a_psi) == 1.0
        and float(b_psi) == 1.0
        and float(z_t) == 1.0
        and float(z_s) == 1.0
    )
    n_sq = float(z_s) / float(z_t)
    c_theta_sq = float(b_theta) / float(a_theta)
    c_psi_sq = float(b_psi) / float(a_psi)
    c_matter_sq = float(z_s) / float(z_t)

    no_new_modes = kinetic_positive and matter_response_positive
    healthy = bool(
        kinetic_positive
        and matter_response_positive
        and c_theta_sq > 0.0
        and c_psi_sq > 0.0
        and c_matter_sq > 0.0
        and n_sq > 0.0
    )

    perturbation_quadratic = {
        "K_micro_diag": [sp.Symbol("A_Theta"), sp.Symbol("A_Psi"), sp.Symbol("Z_t")],
        "G_micro_diag": [sp.Symbol("B_Theta"), sp.Symbol("B_Psi"), sp.Symbol("Z_s")],
        "c_theta_sq": c_theta_sq,
        "c_psi_sq": c_psi_sq,
        "c_matter_sq": c_matter_sq,
    }

    return {
        "bridge": bridge,
        "coeffs": coeffs,
        "lambda_tpsi": float(lambda_tpsi),
        "lambda_tchi": float(lambda_tchi),
        "weak_field_match": weak_field_match,
        "n_sq": n_sq,
        "no_new_modes": no_new_modes,
        "healthy": healthy,
        "perturbation_quadratic": perturbation_quadratic,
        "status": sp.Symbol("derived_micro_perturbation"),
    }


def casimir_boundary_ledger(
    *,
    plate_distance: float = 1.0,
    z_t: float = 1.0,
    z_s: float = 1.0,
) -> Dict[str, Any]:
    """Return the minimal Casimir sanity check for the micro matter sector.

    This is a boundary effect of the matter modes governed by
    ``L_matter = 1/2 Z_t chi_dot^2 - 1/2 Z_s |grad chi|^2 - ...``.
    For homogeneous response coefficients the mode dispersion is
    ``omega = sqrt(Z_s/Z_t) |k|`` in the preferred-slicing coordinates, so the
    parallel-plate Casimir energy scales by the same factor and returns the
    standard result at the unit weak-field point.
    """
    a = float(plate_distance)
    zt = float(z_t)
    zs = float(z_s)
    if a <= 0.0:
        raise ValueError("plate_distance must be positive")

    response_positive = zt > 0.0 and zs > 0.0
    n_sq = zs / zt if response_positive else float("nan")
    c_eff = n_sq ** 0.5 if response_positive else float("nan")

    standard_energy_per_area = -(sp.pi**2) / (720 * sp.Symbol("a") ** 3)
    standard_pressure = -(sp.pi**2) / (240 * sp.Symbol("a") ** 4)
    pi_sq = float((sp.pi**2).evalf())
    energy_per_area = -pi_sq * c_eff / (720.0 * a**3)
    pressure = -pi_sq * c_eff / (240.0 * a**4)

    return {
        "boundary_type": sp.Symbol("parallel_plates"),
        "matter_sector": sp.Symbol("chi_boundary_modes"),
        "new_ontology": 0,
        "plate_distance": a,
        "Z_t": zt,
        "Z_s": zs,
        "n_sq": n_sq,
        "c_eff": c_eff,
        "standard_energy_per_area": standard_energy_per_area,
        "standard_pressure": standard_pressure,
        "energy_per_area": energy_per_area,
        "pressure": pressure,
        "weak_field_standard_limit": bool(
            response_positive and zt == 1.0 and zs == 1.0 and c_eff == 1.0
        ),
        "healthy": bool(response_positive and n_sq > 0.0),
        "status": sp.Symbol("casimir_boundary_sanity"),
    }
