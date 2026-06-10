import math
import numpy as np

from theory.adm_symbolic import (
    effective_matrices_sigma_strict,
    effective_matrices_sigma_strict_schur,
)


def test_sigma_strict_schur_matches_deltaM_route_Q0():
    gamma = 1.0
    V0 = 0.2
    Q = 0.0
    a = 1.0
    k = 0.2
    theta = 0.25
    R = 0.15
    H = 0.05
    theta_dot = -0.01
    R_dot = 0.02
    Hdot = -0.001

    K1, G1, M1 = effective_matrices_sigma_strict(
        theta, R, a, k,
        gamma=gamma, Q=Q, V0=V0,
        H=H, theta_dot=theta_dot, R_dot=R_dot, Hdot=Hdot,
    )
    K2, G2, M2 = effective_matrices_sigma_strict_schur(
        theta, R, a, k,
        gamma=gamma, Q=Q, V0=V0,
        H=H, theta_dot=theta_dot, R_dot=R_dot, Hdot=Hdot,
    )

    # K and G should match exactly; M should match within numerical tolerance
    assert K1.equals(K2)
    assert G1.equals(G2)
    M1n = np.array([[float(M1[i,j]) for j in range(3)] for i in range(3)])
    M2n = np.array([[float(M2[i,j]) for j in range(3)] for i in range(3)])
    assert np.allclose(M1n, M2n, atol=1e-10, rtol=0.0)

