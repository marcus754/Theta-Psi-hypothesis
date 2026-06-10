# -*- coding: utf-8 -*-
import numpy as np

from theory.adm_symbolic import effective_matrices_sym_v6, effective_matrices_v6


def test_effective_matrices_sym_v6_shapes():
    K, G, M = effective_matrices_sym_v6(
        theta=0.2, R=0.1, a=1.0, k=1.0, gamma=1e-3, Q=0.0, V0=0.1
    )
    assert K.shape == (4, 4)
    assert G.shape == (4, 4)
    assert M.shape == (4, 4)


def test_effective_matrices_v6_numeric_sigma_block():
    K, G, M = effective_matrices_v6(
        theta=0.2,
        R=0.1,
        a=1.0,
        k=1.0,
        gamma=1e-3,
        Q=0.0,
        V0=0.1,
        csigma2=0.85,
        msigma2=0.02,
    )
    Kd = np.array(K, dtype=float)
    Gd = np.array(G, dtype=float)
    Md = np.array(M, dtype=float)
    assert Kd.shape == (4, 4)
    assert abs(Kd[3, 3] - 1.0) < 1e-12
    assert abs(Gd[3, 3] - 0.85) < 1e-12
    assert abs(Md[3, 3] - 0.02) < 1e-12
