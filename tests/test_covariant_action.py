import sympy as sp

from theory.covariant_action import (
    canonical_action_ledger,
    covariant_action_summary,
    frw_reduced_lagrangian,
    full_covariant_action_ledger,
    phi_eff_canonical,
    stationary_refractive_sector,
    theta_eff_projection,
)


def test_theta_eff_projection_has_expected_structure():
    d = theta_eff_projection()
    assert "theta_eff" in d
    assert "theta_eff_sq" in d
    assert "Theta" in str(d["theta_eff_sq"])


def test_phi_eff_canonical_is_log_ratio():
    th, th0 = sp.symbols("th th0", positive=True)
    phi = phi_eff_canonical(th, th0)
    assert phi == sp.log(th / th0)


def test_frw_reduced_lagrangian_contains_theta_and_psi_kinetics():
    d = frw_reduced_lagrangian()
    text = str(d["L_eff"])
    assert "Derivative(theta" in text
    assert "Derivative(R" in text
    assert "V0" in text


def test_stationary_refractive_sector_builds_j_refr():
    d = stationary_refractive_sector()
    assert "J_refr" in d
    assert "log" in str(d["Phi_eff"])
    assert "Derivative(theta" in str(d["I_grad"])


def test_covariant_action_summary_contains_core_blocks():
    d = covariant_action_summary()
    s = str(d["S_red"])
    assert "EH" in s
    assert "L_Theta" in s
    assert "L_R" in s
    assert "S_cov" in d
    assert "S_m[g_hat,matter]" in str(d["S_matter"])
    assert "E_theta" in d
    assert "E_R" in d
    assert "rho_red" in d
    assert "p_red" in d


def test_canonical_action_ledger_contains_core_blocks():
    d = canonical_action_ledger()
    for key in ("theta_eff", "Phi_eff", "I_grad", "J_refr", "R", "L_Theta", "L_R", "U_red", "L_refr", "S_red"):
        assert key in d
    assert "log" in str(d["Phi_eff"])
    assert "F(" in str(d["L_refr"])


def test_full_covariant_action_ledger_contains_matter_completion():
    d = full_covariant_action_ledger()
    assert "S_cov" in d
    assert "S_matter" in d
    assert "g_hat" in d
