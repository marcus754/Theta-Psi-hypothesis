import math
from src.linear_modes_adm import stability_adm, stability_finite_k_adm


def test_adm_with_elimination_small_couplings_is_stable():
    gamma = 1.0
    V0 = 0.1
    Q = 0.0
    a = 1.0
    k = 0.1
    # Background point in safe region
    x = 0.5
    theta = (math.sqrt(3.0)/gamma) * math.tanh(x)
    R = 0.2

    elim = {
        # very small couplings to constraints; C is SPD and carries k/a factors where needed
        'lambda_N_theta': 1e-5,
        'lambda_N_R': 1e-5,
        'lambda_N_phi': 1e-5,
        'lambda_B_theta': 1e-5,
        'lambda_B_R': 1e-5,
        'lambda_B_phi': 1e-5,
        'cN2': 1.0,
        'cB2': 1.0,
        'cNB': 0.0,
    }
    res = stability_adm(theta, R, a, k, gamma, Q, V0, ctheta2=0.9, cR2=0.9, cphi2=0.95, elim=elim)
    assert res['no_ghost']
    assert res['no_gradient_instability']

    # Finite-k dispersion is also non-tachyonic along same point
    st = stability_finite_k_adm(theta, R, a, k, gamma, Q, V0, ctheta2=0.9, cR2=0.9, cphi2=0.95, elim=elim)
    assert st['no_ghost'] and st['no_tachyon']
