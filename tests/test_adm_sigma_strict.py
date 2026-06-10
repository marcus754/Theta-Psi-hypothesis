import math

from src.linear_modes_adm import dispersion_eigs_adm


def test_sigma_strict_matches_gravM_dispersion_block():
    # Background point
    gamma = 1.0
    V0 = 0.2
    Q = 0.0
    a = 1.0
    k = 0.15
    theta = 0.2
    R = 0.1
    H = 0.05
    theta_dot = 0.01
    R_dot = -0.02
    Hdot = -0.001

    # ω² with gravM corrections via explicit flag
    w2_grav = dispersion_eigs_adm(
        theta, R, a, k, gamma, Q, V0,
        gravM=True, H=H, theta_dot=theta_dot, R_dot=R_dot, Hdot=Hdot,
    )
    # ω² with sigma_strict elimination should match
    w2_strict = dispersion_eigs_adm(
        theta, R, a, k, gamma, Q, V0,
        elim='sigma_strict', H=H, theta_dot=theta_dot, R_dot=R_dot, Hdot=Hdot,
    )

    for a1, a2 in zip(w2_grav, w2_strict):
        assert abs(a1 - a2) < 1e-12

