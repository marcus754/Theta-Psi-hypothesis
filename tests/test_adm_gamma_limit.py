import math

from src.linear_modes_adm import cs2_eigs_adm


def test_adm_gamma_to_zero_limit_keeps_speeds():
    V0 = 0.2
    Q = 0.0
    a = 1.0
    k = 0.2
    # Choose a safe background point
    theta = 0.2
    R = 0.1
    # Two gamma values: O(1) and small
    g1 = 1.0
    g2 = 1e-6
    c1 = cs2_eigs_adm(theta, R, a, k, g1, Q, V0)
    c2 = cs2_eigs_adm(theta, R, a, k, g2, Q, V0)
    for a1, a2 in zip(c1, c2):
        assert abs(a1 - a2) < 1e-12

