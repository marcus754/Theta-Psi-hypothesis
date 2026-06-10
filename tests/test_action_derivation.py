import sympy as sp

from theory.action_derivation import (
    adm_minimal_quadratic_derivation,
    derivation_status_ledger,
    refractive_function_derivation_status,
    refractive_semigroup_principle,
    strongfield_bridge_derivation_status,
)


def test_adm_minimal_quadratic_derivation_sets_mixings_to_zero():
    d = adm_minimal_quadratic_derivation()
    assert d["status"] == sp.Symbol("adm_minimal_quadratic_derived_from_action")
    assert d["K_from_action"].shape == (3, 3)
    assert d["G_from_action"].shape == (3, 3)
    assert d["M_from_action"].shape == (3, 3)
    assert all(v == 0 for v in d["derived_zero_mixings"].values())
    assert d["cs2_from_action"] == (1, 1, 1)
    assert d["requires_extra_operators_for_nonzero_mixings"] == 1


def test_refractive_f_v1_is_not_derived_by_current_action():
    d = refractive_function_derivation_status()
    assert d["status"] == sp.Symbol("refractive_function_not_yet_derived")
    assert "F_v1" in d["not_derived_function_choice"]
    assert d["derived_argument"] == sp.Symbol("J_refr", nonnegative=True, real=True)


def test_refractive_semigroup_principle_yields_exponential_solution():
    d = refractive_semigroup_principle()
    assert d["status"] == sp.Symbol("refractive_semigroup_principle")
    assert "S(J) = 1 - F(J)" in d["residual"]
    assert str(d["semigroup_law"]).startswith("Eq(")
    assert "exp(-J/J_star)" in str(d["normalized_solution"])
    assert "1 - exp(-J/J_star)" in str(d["F_solution"])


def test_strongfield_bridge_scales_are_not_derived_without_sourced_solution():
    d = strongfield_bridge_derivation_status()
    assert d["status"] == sp.Symbol("strongfield_bridge_not_yet_derived")
    assert d["forbidden_as_canon"] == 1
    assert "theta_c_scale" in d["not_derived_scales"]


def test_derivation_status_ledger_collects_all_risky_knobs():
    d = derivation_status_ledger()
    assert d["status"] == sp.Symbol("action_derivation_status")
    assert set(d) >= {"adm", "refractive", "strongfield_bridge"}
