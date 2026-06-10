import sympy as sp

from theory.adm_symbolic import (
    kinetic_matrix_sym,
    gradient_matrix_sym,
    mass_matrix_sym,
    eliminate_lapse_shift_sym,
)


def test_symbolic_shapes_and_elimination_noop():
    th, R, a, k, g, V0 = 0.5, 0.1, 1.0, 0.2, 1.0, 0.3
    K = kinetic_matrix_sym(th, R, g)
    G = gradient_matrix_sym(th, R, a, k, g)
    M = mass_matrix_sym(th, R, a, k, g, V0)
    assert K.shape == (3, 3) and G.shape == (3, 3) and M.shape == (3, 3)
    Ke, Ge, Me = eliminate_lapse_shift_sym(K, G, M)
    assert Ke.equals(K) and Ge.equals(G) and Me.equals(M)

