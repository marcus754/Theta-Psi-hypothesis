# -*- coding: utf-8 -*-
"""Derived particle ledger for the Θ–Ψ interpretation.

The module treats a particle as a stable localized Θ-projection of a Ψ-sector,
not as a third primitive ontology.
"""
from __future__ import annotations

from dataclasses import dataclass
from typing import Dict, Iterable, Tuple

import sympy as sp


@dataclass(frozen=True)
class QuantumLabels:
    """Minimal conserved labels identifying a repeatable particle sector."""

    mass: float
    charge: float
    spin: float


@dataclass(frozen=True)
class LocalizationProfile:
    """Minimal localization data for a Θ-event footprint."""

    norm: float = 1.0
    width: float = 1.0
    energy: float = 1.0


def localization_condition(profile: LocalizationProfile) -> Dict[str, object]:
    """Return whether a profile can represent a localized finite-energy event."""
    finite_norm = profile.norm > 0.0
    finite_width = profile.width > 0.0
    finite_energy = profile.energy > 0.0
    return {
        "finite_norm": finite_norm,
        "finite_width": finite_width,
        "finite_energy": finite_energy,
        "localized": bool(finite_norm and finite_width and finite_energy),
    }


def conserved_label_condition(labels_before: QuantumLabels, labels_after: QuantumLabels) -> Dict[str, object]:
    """Return whether a sector keeps the same particle-identifying labels."""
    return {
        "mass_conserved": labels_before.mass == labels_after.mass,
        "charge_conserved": labels_before.charge == labels_after.charge,
        "spin_conserved": labels_before.spin == labels_after.spin,
        "labels_conserved": labels_before == labels_after,
    }


def detector_event_map(
    *,
    coupling: float = 1.0,
    threshold: float = 0.5,
    local_density: float = 1.0,
) -> Dict[str, object]:
    """Return a minimal detector click map for a localized Θ-projection."""
    if threshold <= 0.0:
        raise ValueError("threshold must be positive")
    response = float(coupling) * float(local_density)
    return {
        "coupling": float(coupling),
        "threshold": float(threshold),
        "local_density": float(local_density),
        "response": response,
        "click": response >= float(threshold),
        "event_type": sp.Symbol("Theta_local_detector_event"),
    }


def multiparticle_state_ledger(*, entangled: bool = False) -> Dict[str, object]:
    """Return the state distinction between product and nonseparable Ψ sectors."""
    if entangled:
        factorization = sp.Symbol("Psi_AB != Psi_A_tensor_Psi_B")
        independent_complete_states = sp.Integer(0)
        local_events = sp.Integer(2)
    else:
        factorization = sp.Symbol("Psi_AB = Psi_A_tensor_Psi_B")
        independent_complete_states = sp.Integer(1)
        local_events = sp.Integer(2)
    return {
        "state": sp.Symbol("Psi_AB"),
        "factorization": factorization,
        "entangled": int(entangled),
        "independent_complete_states": independent_complete_states,
        "local_theta_events": local_events,
        "status": sp.Symbol("multiparticle_state_distinction"),
    }


def classical_limit_ledger(
    *,
    action_scale: float = 1000.0,
    hbar_scale: float = 1.0,
    decoherence_strength: float = 1.0,
) -> Dict[str, object]:
    """Return a minimal classical trajectory sanity check."""
    if hbar_scale <= 0.0:
        raise ValueError("hbar_scale must be positive")
    phase_ratio = float(action_scale) / float(hbar_scale)
    stationary_phase = phase_ratio > 100.0
    decohered = float(decoherence_strength) > 0.0
    return {
        "action_scale": float(action_scale),
        "hbar_scale": float(hbar_scale),
        "phase_ratio": phase_ratio,
        "stationary_phase": stationary_phase,
        "decohered": decohered,
        "classical_trajectory_limit": bool(stationary_phase and decohered),
        "status": sp.Symbol("classical_particle_limit"),
    }


def particle_ledger(
    *,
    labels: QuantumLabels = QuantumLabels(mass=1.0, charge=-1.0, spin=0.5),
    profile: LocalizationProfile = LocalizationProfile(),
) -> Dict[str, object]:
    """Return the derived Θ–Ψ definition of a particle."""
    localization = localization_condition(profile)
    labels_check = conserved_label_condition(labels, labels)
    detector = detector_event_map(local_density=profile.norm / profile.width)
    classical = classical_limit_ledger()
    return {
        "definition": sp.Symbol("stable_localized_Theta_projection_of_Psi_sector"),
        "primary_ontology": (sp.Symbol("Theta_mu_nu"), sp.Symbol("Psi")),
        "new_ontology": sp.Integer(0),
        "psi_sector": sp.Symbol("Psi_sector"),
        "theta_projection": sp.Symbol("Theta_localization_anchor"),
        "labels": labels,
        "localization": localization,
        "conserved_labels": labels_check,
        "detector_event": detector,
        "single_particle_state": sp.Symbol("Psi_A"),
        "entangled_pair_state": multiparticle_state_ledger(entangled=True),
        "product_pair_state": multiparticle_state_ledger(entangled=False),
        "classical_limit": classical,
        "acceptance": {
            "localized": localization["localized"],
            "labels_conserved": labels_check["labels_conserved"],
            "detector_click_possible": detector["click"],
            "classical_limit_exists": classical["classical_trajectory_limit"],
        },
        "status": sp.Symbol("derived_particle_definition"),
    }
