from theory.refractive_sector import f_v1, f_v1_prime, f_v1_double_prime
from theory.refractive_sector import (
    f_v1_with_scale,
    f_v1_prime_with_scale,
    f_v1_double_prime_with_scale,
)
from theory.action_derivation import refractive_function_derivation_status


def test_f_v1_is_monotone_and_saturating():
    vals = [f_v1(x) for x in (0.0, 0.2, 1.0, 5.0, 20.0)]
    assert vals[0] == 0.0
    assert all(vals[i] <= vals[i + 1] for i in range(len(vals) - 1))
    assert vals[-1] < 1.0
    assert vals[-1] > 0.9


def test_f_v1_prime_is_positive_and_decays():
    vals = [f_v1_prime(x) for x in (0.0, 0.2, 1.0, 5.0)]
    assert all(v > 0.0 for v in vals)
    assert all(vals[i] >= vals[i + 1] for i in range(len(vals) - 1))


def test_f_v1_double_prime_is_negative():
    vals = [f_v1_double_prime(x) for x in (0.0, 0.2, 1.0, 5.0)]
    assert all(v < 0.0 for v in vals)
    assert all(abs(vals[i]) >= abs(vals[i + 1]) for i in range(len(vals) - 1))


def test_scaled_family_remains_diagnostic_only():
    vals = [f_v1_with_scale(x, j_star=0.7) for x in (0.0, 0.2, 1.0, 5.0, 20.0)]
    assert vals[0] == 0.0
    assert vals[-1] < 1.0
    assert f_v1_prime_with_scale(0.0, j_star=0.7) > 0.0
    assert f_v1_double_prime_with_scale(0.0, j_star=0.7) < 0.0


def test_f_v1_is_representative_not_derived_canon():
    d = refractive_function_derivation_status()
    assert "F_v1" in d["not_derived_function_choice"]
