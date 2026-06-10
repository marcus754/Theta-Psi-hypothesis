import math
import numpy as np

from theory.adm_symbolic import mass_matrix_sigma_total, schur_blocks_sigma_strict_numeric


def test_schur_blocks_reproduce_deltaM_in_sigma_strict_Q0():
    # Background with Q=0 so ΔM has rank ≤2
    gamma = 1.0
    V0 = 0.2
    Q = 0.0
    a = 1.0
    theta = 0.3
    R = 0.1
    H = 0.05
    theta_dot = 0.01
    R_dot = -0.02
    Hdot = -0.001

    Mfull = np.array(mass_matrix_sigma_total(theta, R, a, gamma=gamma, Q=Q, V0=V0, H=H, theta_dot=theta_dot, R_dot=R_dot, Hdot=Hdot))
    Mpot_rr = V0
    Mpot = np.array([[0.0, 0.0, 0.0], [0.0, Mpot_rr, 0.0], [0.0, 0.0, 0.0]])
    deltaM = Mfull - Mpot

    B, C = schur_blocks_sigma_strict_numeric(theta, R, a, gamma=gamma, Q=Q, V0=V0, H=H, theta_dot=theta_dot, R_dot=R_dot, Hdot=Hdot)
    B = np.array(B)
    C = np.array(C)
    # Reconstruct -ΔM via B C^{-1} B^T
    target = -deltaM
    Cinvt = np.linalg.inv(C)
    recon = B @ Cinvt @ B.T
    assert np.allclose(recon, target, atol=1e-12, rtol=0.0)

