import numpy as np

from theory.adm_symbolic import effective_matrices_sigma_strict_compare


def test_sigma_strict_schur_matches_deltaM_route_Qneq0():
    gamma = 1.0
    V0 = 0.2
    Q = 0.01
    a = 1.0
    k = 0.2
    theta = 0.3
    R = 0.12
    H = 0.05
    theta_dot = 0.02
    R_dot = -0.01
    Hdot = -0.001

    (K1, G1, M1), (K2, G2, M2) = effective_matrices_sigma_strict_compare(
        theta, R, a, k,
        gamma=gamma, Q=Q, V0=V0,
        H=H, theta_dot=theta_dot, R_dot=R_dot, Hdot=Hdot,
    )

    assert K1.equals(K2)
    assert G1.equals(G2)
    M1n = np.array([[float(M1[i,j]) for j in range(3)] for i in range(3)])
    M2n = np.array([[float(M2[i,j]) for j in range(3)] for i in range(3)])
    assert np.allclose(M1n, M2n, atol=1e-12, rtol=0.0)

