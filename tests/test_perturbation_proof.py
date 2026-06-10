from theory.perturbation_proof import perturbation_proof_ledger


def test_perturbation_proof_ledger_closes_in_nominal_regime():
    d = perturbation_proof_ledger(
        theta=1.2,
        R=0.9,
        a=1.1,
        k=0.7,
        gamma=0.8,
        Q=0.0,
        V0=0.2,
        H=0.5,
        theta_dot=0.03,
        R_dot=0.02,
        Hdot=-0.01,
    )
    assert "S_cov" in d
    assert "g_hat" in d
    assert d["routes_match"]
    assert d["health"]["healthy"]
    assert d["proof_closed"]
