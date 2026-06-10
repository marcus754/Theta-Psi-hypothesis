
"""
Linear modes (v0.3): kinetic/gradient with mixing
q = (δθ, δΨ). We model the quadratic action as
  L^(2) = 1/2 [ δq̇^T K δq̇ - δq^T G (k^2/a^2) δq ],
with
  K = [[Kθθ, KθΨ],
       [KθΨ, 1]],   where Kθθ = 3/θ^2 + γ^2,  KθΨ = ε_K
  G = [[Gθθ, GθΨ],
       [GθΨ, GΨΨ]], where Gθθ = cθ^2 Kθθ, GΨΨ = cΨ^2, GθΨ = ε_G

This captures leading mixing without a full metric derivation.
Stability requires:
  - No ghosts: K ≻ 0  (Sylvester: Kθθ>0 and det K > 0)
  - No gradient instabilities: in the high-k limit G ≻ 0 in the kinetic metric,
    equivalently eigenvalues of K^{-1} G are positive → we check via Cholesky of K and symmetrization.

Functions:
- kinetic_matrix(theta, gamma, epsK)
- gradient_matrix(theta, gamma, ctheta2, cpsi2, epsG)
- cs2_eigs(theta, gamma, ctheta2, cpsi2, epsK, epsG)
- stability(theta, gamma, ctheta2, cpsi2, epsK, epsG)
Notes:
- Conceptually Ψ is complex in the hypothesis; the current numerical scaffold
  treats Ψ as a real amplitude (phase unused) until the full ADM sector is derived.
"""

from __future__ import annotations
from typing import Dict, Tuple
import math

def kinetic_matrix(theta: float, gamma: float, epsK: float = 0.0):
    if theta == 0.0:
        raise ValueError("theta=0 is singular; use x variable to avoid exact zero.")
    Kth = 3.0/(theta*theta) + gamma*gamma
    return ((Kth, epsK),
            (epsK, 1.0))

def gradient_matrix(theta: float, gamma: float, ctheta2: float = 1.0, cpsi2: float = 1.0, epsG: float = 0.0):
    K = kinetic_matrix(theta, gamma, 0.0)
    Gth = ctheta2 * K[0][0]
    return ((Gth, epsG),
            (epsG, cpsi2))

def _det2(M):
    return M[0][0]*M[1][1] - M[0][1]*M[1][0]

def _is_pd_2x2(M):
    # Sylvester's criterion
    return (M[0][0] > 0.0) and (_det2(M) > 0.0)

def _matmul(A,B):
    return ((A[0][0]*B[0][0]+A[0][1]*B[1][0], A[0][0]*B[0][1]+A[0][1]*B[1][1]),
            (A[1][0]*B[0][0]+A[1][1]*B[1][0], A[1][0]*B[0][1]+A[1][1]*B[1][1]))

def _inv2(M):
    det = _det2(M)
    if det == 0.0:
        raise ValueError("Singular matrix")
    invdet = 1.0/det
    return (( M[1][1]*invdet, -M[0][1]*invdet),
            (-M[1][0]*invdet,  M[0][0]*invdet))

def cs2_eigs(theta: float, gamma: float, ctheta2: float = 1.0, cpsi2: float = 1.0, epsK: float = 0.0, epsG: float = 0.0):
    K = kinetic_matrix(theta, gamma, epsK)
    G = gradient_matrix(theta, gamma, ctheta2, cpsi2, epsG)
    # eigenvalues of K^{-1} G for 2x2
    Kinv = _inv2(K)
    A = _matmul(Kinv, G)
    # eigenvalues of 2x2 matrix
    tr = A[0][0] + A[1][1]
    det = _det2(A)
    disc = max(tr*tr - 4.0*det, 0.0)
    s1 = 0.5*(tr + math.sqrt(disc))
    s2 = 0.5*(tr - math.sqrt(disc))
    return (s1, s2)

def stability(theta: float, gamma: float, ctheta2: float = 1.0, cpsi2: float = 1.0, epsK: float = 0.0, epsG: float = 0.0) -> Dict[str,bool]:
    K = kinetic_matrix(theta, gamma, epsK)
    G = gradient_matrix(theta, gamma, ctheta2, cpsi2, epsG)
    no_ghost = _is_pd_2x2(K)
    lam = cs2_eigs(theta, gamma, ctheta2, cpsi2, epsK, epsG)
    no_grad = (lam[0] > 0.0) and (lam[1] > 0.0)
    return {"no_ghost": no_ghost, "no_gradient_instability": no_grad, "cs2_eigs": lam}

# Backwards-compatible aliases used by older tests/utilities ------------------

def cs2_eigenvalues(*args, **kwargs):
    return cs2_eigs(*args, **kwargs)


def stability_conditions(theta: float, gamma: float, ctheta2: float = 1.0, cpsi2: float = 1.0, epsK: float = 0.0, epsG: float = 0.0) -> Dict[str, bool]:
    return stability(theta, gamma, ctheta2, cpsi2, epsK, epsG)


def is_kinetically_stable(theta: float, theta_dot: float, psi: float, psi_dot: float, H: float, gamma: float) -> bool:
    """Return True if the kinetic matrix is positive-definite.

    The current placeholder ignores time derivatives and H; they are accepted
    for API compatibility and future extensions.
    """
    # Currently only depends on theta and gamma (no mixing with velocities)
    K = kinetic_matrix(theta, gamma, 0.0)
    return _is_pd_2x2(K)
