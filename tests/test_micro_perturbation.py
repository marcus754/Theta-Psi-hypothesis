import pytest

from theory.micro_perturbation import casimir_boundary_ledger, micro_perturbation_ledger


def test_micro_perturbation_ledger_is_healthy_at_unit_point():
    d = micro_perturbation_ledger()
    assert d["weak_field_match"] is True
    assert d["no_new_modes"] is True
    assert d["healthy"] is True
    assert d["n_sq"] == 1.0
    assert d["coeffs"]["Z_t"] == 1.0
    assert d["coeffs"]["Z_s"] == 1.0


def test_micro_perturbation_ledger_requires_positive_coefficients():
    d = micro_perturbation_ledger(a_theta=2.0, b_theta=1.0, a_psi=3.0, b_psi=4.0, z_t=5.0, z_s=6.0)
    assert d["healthy"] is True
    assert d["no_new_modes"] is True
    assert d["n_sq"] == 6.0 / 5.0


def test_micro_perturbation_ledger_flags_nonpositive_response():
    d = micro_perturbation_ledger(z_t=1.0, z_s=0.0)
    assert d["healthy"] is False
    assert d["no_new_modes"] is False


def test_casimir_boundary_ledger_returns_standard_weak_field_limit():
    d = casimir_boundary_ledger(plate_distance=2.0)
    assert d["healthy"] is True
    assert d["weak_field_standard_limit"] is True
    assert d["new_ontology"] == 0
    assert d["n_sq"] == 1.0
    assert d["c_eff"] == 1.0
    assert d["energy_per_area"] == pytest.approx(-(3.141592653589793**2) / (720.0 * 2.0**3))
    assert d["pressure"] == pytest.approx(-(3.141592653589793**2) / (240.0 * 2.0**4))


def test_casimir_boundary_ledger_scales_with_micro_response():
    d = casimir_boundary_ledger(plate_distance=1.0, z_t=4.0, z_s=9.0)
    assert d["healthy"] is True
    assert d["weak_field_standard_limit"] is False
    assert d["n_sq"] == 9.0 / 4.0
    assert d["c_eff"] == 1.5
    assert d["energy_per_area"] == pytest.approx(-1.5 * (3.141592653589793**2) / 720.0)


def test_casimir_boundary_ledger_requires_positive_plate_distance():
    with pytest.raises(ValueError):
        casimir_boundary_ledger(plate_distance=0.0)
