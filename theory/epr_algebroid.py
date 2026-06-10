# -*- coding: utf-8 -*-
"""EPR/CHSH algebroid sanity ledger for the Θ–Ψ interpretation.

This module is intentionally narrow. It does not introduce a tunable hidden
variable model. It pins one checkable picture:

- Ψ carries the nonseparable singlet state;
- Θ supplies the causal anchor;
- correlation generators live in ker rho_Theta and therefore do not project to
  controllable superluminal signals;
- signal generators, when present, must project inside the Θ-causal cone.
"""
from __future__ import annotations

from math import sqrt
from typing import Callable, Dict, Iterable, Tuple

import sympy as sp

Vector = Tuple[float, float, float]


def _dot(a: Vector, b: Vector) -> float:
    return sum(x * y for x, y in zip(a, b))


def _norm(a: Vector) -> float:
    return sqrt(_dot(a, a))


def normalize(v: Iterable[float]) -> Vector:
    """Return a unit 3-vector."""
    vec = tuple(float(x) for x in v)
    if len(vec) != 3:
        raise ValueError("measurement direction must be a 3-vector")
    length = _norm(vec)  # type: ignore[arg-type]
    if length <= 0.0:
        raise ValueError("measurement direction must be nonzero")
    return tuple(x / length for x in vec)  # type: ignore[return-value]


def singlet_correlation(a: Iterable[float], b: Iterable[float]) -> float:
    """Return the fixed singlet correlation E(a,b) = -a·b."""
    ahat = normalize(a)
    bhat = normalize(b)
    return -_dot(ahat, bhat)


def singlet_joint_probability(
    outcome_a: int,
    outcome_b: int,
    a: Iterable[float],
    b: Iterable[float],
) -> float:
    """Return P(A=outcome_a, B=outcome_b | a,b) for outcomes ±1."""
    if outcome_a not in (-1, 1) or outcome_b not in (-1, 1):
        raise ValueError("singlet outcomes must be ±1")
    corr = singlet_correlation(a, b)
    return 0.25 * (1.0 + outcome_a * outcome_b * corr)


def marginal_a(outcome_a: int, a: Iterable[float], b: Iterable[float]) -> float:
    """Return P(A=outcome_a | a,b)."""
    return sum(singlet_joint_probability(outcome_a, outcome_b, a, b) for outcome_b in (-1, 1))


def marginal_b(outcome_b: int, a: Iterable[float], b: Iterable[float]) -> float:
    """Return P(B=outcome_b | a,b)."""
    return sum(singlet_joint_probability(outcome_a, outcome_b, a, b) for outcome_a in (-1, 1))


def chsh_value(
    a0: Iterable[float],
    a1: Iterable[float],
    b0: Iterable[float],
    b1: Iterable[float],
) -> float:
    """Return the CHSH value for the fixed singlet rule."""
    return (
        singlet_correlation(a0, b0)
        + singlet_correlation(a0, b1)
        + singlet_correlation(a1, b0)
        - singlet_correlation(a1, b1)
    )


def chsh_from_correlations(e00: float, e01: float, e10: float, e11: float) -> float:
    """Return S = E00 + E01 + E10 - E11 for explicit correlations."""
    return float(e00) + float(e01) + float(e10) - float(e11)


def optimal_chsh_directions() -> Dict[str, Vector]:
    """Return one standard direction choice reaching |S| = 2√2."""
    inv_sqrt2 = 1.0 / sqrt(2.0)
    return {
        "a0": (1.0, 0.0, 0.0),
        "a1": (0.0, 1.0, 0.0),
        "b0": (inv_sqrt2, inv_sqrt2, 0.0),
        "b1": (inv_sqrt2, -inv_sqrt2, 0.0),
    }


def tsirelson_gate(s_chsh: float) -> Dict[str, object]:
    """Classify a CHSH value by local, quantum and super-quantum bounds."""
    value = abs(float(s_chsh))
    bell_bound = 2.0
    tsirelson = 2.0 * sqrt(2.0)
    return {
        "S_CHSH": float(s_chsh),
        "abs_S_CHSH": value,
        "bell_bound": bell_bound,
        "tsirelson_bound": tsirelson,
        "violates_bell": value > bell_bound,
        "within_quantum_bound": value <= tsirelson,
        "reject_super_tsirelson": value > tsirelson,
        "status": sp.Symbol("tsirelson_gate"),
    }


def pr_box_ledger() -> Dict[str, object]:
    """Return the PR-box diagnostic: no-signalling but super-quantum."""
    # CHSH convention matching chsh_from_correlations:
    # E00 = E01 = E10 = 1, E11 = -1 gives S = 4.
    correlations = {
        "E00": 1.0,
        "E01": 1.0,
        "E10": 1.0,
        "E11": -1.0,
    }
    s_chsh = chsh_from_correlations(
        correlations["E00"],
        correlations["E01"],
        correlations["E10"],
        correlations["E11"],
    )
    gate = tsirelson_gate(s_chsh)
    return {
        "correlations": correlations,
        "S_CHSH": s_chsh,
        "no_signalling": sp.Integer(1),
        "gate": gate,
        "accepted_by_theta_psi": sp.Integer(0),
        "rejection_reason": sp.Symbol("super_tsirelson_no_signalling_is_not_enough"),
        "status": sp.Symbol("pr_box_rejected"),
    }


def foliation_order_probabilities(
    a: Iterable[float],
    b: Iterable[float],
) -> Dict[str, Dict[Tuple[int, int], float]]:
    """Return joint probabilities for A-before-B and B-before-A descriptions."""
    joint = {
        (oa, ob): singlet_joint_probability(oa, ob, a, b)
        for oa in (-1, 1)
        for ob in (-1, 1)
    }
    return {
        "A_before_B": joint,
        "B_before_A": dict(joint),
    }


def ghz_mermin_ledger() -> Dict[str, object]:
    """Return the GHZ/Mermin contextuality sanity ledger.

    The four GHZ constraints multiply to -1 quantum mechanically, while any
    noncontextual predetermined ±1 valuation makes the left-hand product +1.
    """
    constraints = {
        "XYY": 1,
        "YXY": 1,
        "YYX": 1,
        "XXX": -1,
    }
    quantum_product = 1
    for value in constraints.values():
        quantum_product *= value
    hidden_valuation_product = 1
    return {
        "state": sp.Symbol("GHZ"),
        "constraints": constraints,
        "quantum_product": quantum_product,
        "noncontextual_hidden_valuation_product": hidden_valuation_product,
        "global_hidden_valuation_exists": sp.Integer(0),
        "contextuality_required": sp.Integer(1),
        "status": sp.Symbol("ghz_mermin_contextuality"),
    }


def experimental_profile_from_angles(
    angle_a0: float,
    angle_a1: float,
    angle_b0: float,
    angle_b1: float,
    *,
    correlation: Callable[[Iterable[float], Iterable[float]], float] = singlet_correlation,
) -> Dict[str, object]:
    """Return no-fit CHSH predictions for coplanar angle settings."""
    def direction(phi: float) -> Vector:
        import math

        return (math.cos(float(phi)), math.sin(float(phi)), 0.0)

    a0 = direction(angle_a0)
    a1 = direction(angle_a1)
    b0 = direction(angle_b0)
    b1 = direction(angle_b1)
    correlations = {
        "E00": correlation(a0, b0),
        "E01": correlation(a0, b1),
        "E10": correlation(a1, b0),
        "E11": correlation(a1, b1),
    }
    s_chsh = chsh_from_correlations(
        correlations["E00"],
        correlations["E01"],
        correlations["E10"],
        correlations["E11"],
    )
    return {
        "angles": {
            "a0": float(angle_a0),
            "a1": float(angle_a1),
            "b0": float(angle_b0),
            "b1": float(angle_b1),
        },
        "correlations": correlations,
        "S_CHSH": s_chsh,
        "gate": tsirelson_gate(s_chsh),
        "free_fit_parameters": sp.Integer(0),
        "status": sp.Symbol("no_fit_epr_angle_profile"),
    }


def quantum_layer_bridge_ledger() -> Dict[str, object]:
    """Return the bridge from the quantum state layer to the EPR ledger."""
    return {
        "quantum_layer": sp.Symbol("docs_04_quantum_state_layer"),
        "epr_layer": sp.Symbol("docs_14_epr_algebroid"),
        "shared_object": sp.Symbol("Psi_state"),
        "single_particle_reduction": sp.Symbol("Psi(x)"),
        "two_subsystem_reduction": sp.Symbol("Psi_AB"),
        "singlet_projection": sp.Symbol("Psi_AB_singlet"),
        "bridge_chain": (
            "Psi_state_layer",
            "two_subsystem_entanglement",
            "singlet_projection",
            "fixed_CHSH_correlations",
            "Theta_anchor_no_signalling",
        ),
        "status": sp.Symbol("quantum_to_epr_bridge"),
    }


def epr_algebroid_ledger() -> Dict[str, object]:
    """Return the minimal EPR algebroid ledger with no fit parameters."""
    dirs = optimal_chsh_directions()
    s_chsh = chsh_value(dirs["a0"], dirs["a1"], dirs["b0"], dirs["b1"])
    correlation_generator = sp.Symbol("C_AB")
    signal_generator = sp.Symbol("S_signal")

    anchor = {
        "rho_Theta": sp.Symbol("rho_Theta"),
        "correlation_generator": correlation_generator,
        "signal_generator": signal_generator,
        "rho_correlation": sp.Integer(0),
        "rho_signal_speed_sq": sp.Symbol("v_signal_sq", nonnegative=True, real=True),
        "causal_signal_condition": sp.Symbol("v_signal_sq <= 1"),
    }

    return {
        "state": sp.Symbol("Psi_AB_singlet"),
        "algebroid": sp.Symbol("A_to_M"),
        "state_bundle": sp.Symbol("E_to_M"),
        "anchor": anchor,
        "correlation_rule": sp.Symbol("E(a,b) = - dot(a,b)"),
        "free_correlation_function": sp.Integer(0),
        "free_fit_parameters": sp.Integer(0),
        "directions": dirs,
        "S_CHSH": s_chsh,
        "bell_bound": 2.0,
        "tsirelson_bound": 2.0 * sqrt(2.0),
        "violates_bell": abs(s_chsh) > 2.0,
        "respects_tsirelson": abs(s_chsh) <= 2.0 * sqrt(2.0),
        "tsirelson_gate": tsirelson_gate(s_chsh),
        "pr_box_rejection": pr_box_ledger(),
        "ghz_contextuality": ghz_mermin_ledger(),
        "quantum_layer_bridge": quantum_layer_bridge_ledger(),
        "no_signalling_target": 0.5,
        "interpretation": sp.Symbol("Psi_nonseparable_Theta_causal_anchor"),
        "status": sp.Symbol("epr_algebroid_sanity"),
    }
