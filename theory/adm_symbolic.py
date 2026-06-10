# -*- coding: utf-8 -*-
"""
Symbolic + numeric scaffold for scalar perturbations (ADM route).

Goal: construct effective quadratic action for q = (δθ, δR, δφ) after
eliminating non-dynamical ADM variables (lapse δN and scalar shift potential B).

Design:
- Provide base symbolic K(θ,R,γ) and gradient/mass blocks G, M for the 3×3
  dynamical sector, matching the working diagonal core used across the repo.
- Provide a generic symbolic elimination of lapse/shift via Schur complement of
  their algebraic constraints, parameterized by a small set of coupling
  coefficients. This keeps the interface stable and enables systematic
  refinements once explicit coefficients from the full derivation are available.

Diagonal core (unchanged by default when all ADM couplings are 0):
  Kθθ = 3/θ² + γ²,  KRR = 1,  Kφφ = R²
  Gθθ = cθ² Kθθ,    GRR = cR², Gφφ = cφ² R²
  Mθθ = 0,          MRR = V0,  Mφφ = 0

Elimination model (constraints):
- Introduce algebraic variables c = (δN, B). At quadratic order they enter as
  L ⊃ − 1/2 cᵀ C c − qᵀ B c, with symmetric 2×2 C and 3×2 B depending on (a,k).
- Variation w.r.t c gives c* = − C^{-1} Bᵀ q. Substitution yields an effective
  potential/gradient block A_eff = A − B C^{-1} Bᵀ, where A = G (k²/a²) + M.
- No q̇–c couplings are included by default (keeps K unchanged), but the API is
  structured so that such terms can be added later if needed.

Important defaults:
- All ADM couplings are 0 by default; elimination is an identity map and the
  previous behavior is preserved exactly.
- Couplings are exposed so they can be tuned in studies; suggested usage is to
  set them ∝ γ (so the ΛCDM/γ→0 limit smoothly recovers the canonical form).
"""

from __future__ import annotations

from typing import Tuple, Optional, Dict, Any

import sympy as sp
import numpy as np


# ---------------------------------------------------------------------------
# Symbolic assembly (for future full derivation)

def kinetic_matrix_sym(
    theta: float | sp.Expr,
    R: float | sp.Expr,
    gamma: float | sp.Expr,
    *,
    epsK_theta_R: float | sp.Expr = 0.0,
    epsK_theta_phi: float | sp.Expr = 0.0,
    epsK_R_phi: float | sp.Expr = 0.0,
):
    theta_s = sp.sympify(theta)
    R_s = sp.sympify(R)
    gamma_s = sp.sympify(gamma)
    Kth = 3.0 / (theta_s * theta_s) + gamma_s * gamma_s
    KRR = sp.Integer(1)
    Kpp = R_s * R_s
    return sp.Matrix(
        [
            [Kth, epsK_theta_R, epsK_theta_phi],
            [epsK_theta_R, KRR, epsK_R_phi],
            [epsK_theta_phi, epsK_R_phi, Kpp],
        ]
    )


def gradient_matrix_sym(
    theta: float | sp.Expr,
    R: float | sp.Expr,
    a: float | sp.Expr,
    k: float | sp.Expr,
    gamma: float | sp.Expr,
    *,
    ctheta2: float | sp.Expr = 1.0,
    cR2: float | sp.Expr = 1.0,
    cphi2: float | sp.Expr = 1.0,
    epsG_theta_R: float | sp.Expr = 0.0,
    epsG_theta_phi: float | sp.Expr = 0.0,
    epsG_R_phi: float | sp.Expr = 0.0,
    alpha_k2_theta: float | sp.Expr = 0.0,
    alpha_k2_R: float | sp.Expr = 0.0,
    alpha_k2_phi: float | sp.Expr = 0.0,
):
    theta_s = sp.sympify(theta)
    R_s = sp.sympify(R)
    a_s = sp.sympify(a)
    k_s = sp.sympify(k)
    gamma_s = sp.sympify(gamma)
    Kth = 3.0 / (theta_s * theta_s) + gamma_s * gamma_s
    Gth = sp.sympify(ctheta2) * Kth + sp.sympify(alpha_k2_theta) * (k_s * k_s) / (a_s * a_s)
    GRR = sp.sympify(cR2) + sp.sympify(alpha_k2_R) * (k_s * k_s) / (a_s * a_s)
    Gpp = sp.sympify(cphi2) * (R_s * R_s) + sp.sympify(alpha_k2_phi) * (k_s * k_s) / (a_s * a_s)
    return sp.Matrix(
        [
            [Gth, epsG_theta_R, epsG_theta_phi],
            [epsG_theta_R, GRR, epsG_R_phi],
            [epsG_theta_phi, epsG_R_phi, Gpp],
        ]
    )


def mass_matrix_sym(
    theta: float | sp.Expr,
    R: float | sp.Expr,
    a: float | sp.Expr,
    k: float | sp.Expr,
    gamma: float | sp.Expr,
    V0: float | sp.Expr,
    *,
    Q: float | sp.Expr = 0.0,
    mtheta2: float | sp.Expr | None = None,
):
    """
    Mass matrix M for (δθ, δR, δφ) in the minimal sigma-model (EH + scalars).

    For Ψ = R e^{iφ} with conserved charge Q = a^3 R^2 φ̇, the effective
    potential contains 1/2 V0 R^2 + 1/2 Q^2/(a^6 R^2). Thus
      ∂^2 V_eff/∂R^2 = V0 + 3 Q^2/(a^6 R^4).
    We set Mθθ = 0 by default (no explicit V(θ)), Mφφ = 0 (shift symmetry).
    """
    Mth = sp.Integer(0) if mtheta2 is None else sp.sympify(mtheta2)
    V0_s = sp.sympify(V0)
    Q_s = sp.sympify(Q)
    a_s = sp.sympify(a)
    R_s = sp.sympify(R)
    MRR = V0_s + sp.Integer(3) * (Q_s * Q_s) / (a_s ** 6 * R_s ** 4)
    Mpp = sp.Integer(0)
    return sp.Matrix([[Mth, 0, 0], [0, MRR, 0], [0, 0, Mpp]])


def _to_tuple_matrix(M: sp.Matrix) -> Tuple[Tuple[float, ...], ...]:
    return tuple(tuple(float(M[i, j]) for j in range(M.shape[1])) for i in range(M.shape[0]))


# ---------------------------------------------------------------------------
# Numeric API (stable for downstream code)

def kinetic_matrix(
    theta: float,
    R: float,
    a: float,
    k: float,
    gamma: float,
    Q: float,
    V0: float,
    *,
    epsK_theta_R: float = 0.0,
    epsK_theta_phi: float = 0.0,
    epsK_R_phi: float = 0.0,
) -> Tuple[Tuple[float, ...], ...]:
    if theta == 0.0:
        raise ValueError("theta=0 is singular; use x-variable to avoid exact zero.")
    Msym = kinetic_matrix_sym(theta, R, gamma,
                              epsK_theta_R=epsK_theta_R, epsK_theta_phi=epsK_theta_phi, epsK_R_phi=epsK_R_phi)
    return _to_tuple_matrix(Msym)


def gradient_matrix(
    theta: float,
    R: float,
    a: float,
    k: float,
    gamma: float,
    Q: float,
    V0: float,
    *,
    ctheta2: float = 1.0,
    cR2: float = 1.0,
    cphi2: float = 1.0,
    epsG_theta_R: float = 0.0,
    epsG_theta_phi: float = 0.0,
    epsG_R_phi: float = 0.0,
    # Optional k-dependent diagonal corrections (mimic metric-induced terms)
    alpha_k2_theta: float = 0.0,
    alpha_k2_R: float = 0.0,
    alpha_k2_phi: float = 0.0,
) -> Tuple[Tuple[float, ...], ...]:
    Msym = gradient_matrix_sym(
        theta,
        R,
        a,
        k,
        gamma,
        ctheta2=ctheta2,
        cR2=cR2,
        cphi2=cphi2,
        epsG_theta_R=epsG_theta_R,
        epsG_theta_phi=epsG_theta_phi,
        epsG_R_phi=epsG_R_phi,
        alpha_k2_theta=alpha_k2_theta,
        alpha_k2_R=alpha_k2_R,
        alpha_k2_phi=alpha_k2_phi,
    )
    return _to_tuple_matrix(Msym)


def mass_matrix(
    theta: float,
    R: float,
    a: float,
    k: float,
    gamma: float,
    Q: float,
    V0: float,
    *,
    mtheta2: float | None = None,
) -> Tuple[Tuple[float, ...], ...]:
    Msym = mass_matrix_sym(theta, R, a, k, gamma, V0, Q=Q, mtheta2=mtheta2)
    return _to_tuple_matrix(Msym)


# ---------------------------------------------------------------------------
# Quadratic action scaffold and elimination of lapse/shift (generic Schur)

def build_lagrangian2_sym(
    K: sp.Matrix, G: sp.Matrix, M: sp.Matrix, a: float | sp.Expr, k: float | sp.Expr
) -> sp.Expr:
    """
    Return symbolic quadratic action L^(2) = 1/2 δq̇ᵀ K δq̇ − 1/2 δqᵀ [ G (k²/a²) + M ] δq.

    Variables are implicit; consumers are expected to substitute the field
    components or use this as a formal container for coefficients.
    """
    # Formal placeholders for derivatives and fields
    dq = sp.symbols("dth dR dphi")
    dqdot = sp.symbols("dthdot dRdot dphidot")
    v = sp.Matrix(dq)
    vdot = sp.Matrix(dqdot)
    a_s = sp.sympify(a); k_s = sp.sympify(k)
    block = G * (k_s * k_s / (a_s * a_s)) + M
    return sp.Rational(1, 2) * (vdot.T * K * vdot)[0] - sp.Rational(1, 2) * (v.T * block * v)[0]

def _constraint_blocks_sym(
    theta: float | sp.Expr,
    R: float | sp.Expr,
    a: float | sp.Expr,
    k: float | sp.Expr,
    gamma: float | sp.Expr,
    Q: float | sp.Expr,
    V0: float | sp.Expr,
    *,
    # Couplings of constraints to dynamical sector (columns: δN, B)
    lambda_N_theta: float | sp.Expr = 0.0,
    lambda_N_R: float | sp.Expr = 0.0,
    lambda_N_phi: float | sp.Expr = 0.0,
    lambda_B_theta: float | sp.Expr = 0.0,
    lambda_B_R: float | sp.Expr = 0.0,
    lambda_B_phi: float | sp.Expr = 0.0,
    # Quadratic form for constraints (positive-definite)
    cN2: float | sp.Expr = 1.0,
    cB2: float | sp.Expr = 1.0,
    cNB: float | sp.Expr = 0.0,
    ctheta2: float | sp.Expr = 1.0,
    cR2: float | sp.Expr = 1.0,
    cphi2: float | sp.Expr = 1.0,
    epsG_theta_R: float | sp.Expr = 0.0,
    epsG_theta_phi: float | sp.Expr = 0.0,
    epsG_R_phi: float | sp.Expr = 0.0,
    alpha_k2_theta: float | sp.Expr = 0.0,
    alpha_k2_R: float | sp.Expr = 0.0,
    alpha_k2_phi: float | sp.Expr = 0.0,
) -> Tuple[sp.Matrix, sp.Matrix, sp.Matrix, sp.Expr, sp.Expr]:
    """
    Build symbolic blocks (A, B, C) for constraints elimination with
    A = G (k^2/a^2) + M  (3×3), B (3×2), C (2×2) and return also a_s, k_s.

    Notes:
    - B = [b_N, b_B] with b_B carrying an explicit (k/a) to mimic gradient
      nature of the shift constraint; C_BB ∝ (k^2/a^2) ensures the Schur
      complement yields an A_eff affine in (k^2/a^2).
    - With all λ's = 0 the elimination is a no-op.
    - cNB off-diagonal is supported but defaults to 0 to keep A_eff even in k.
    """
    th_s = sp.sympify(theta); R_s = sp.sympify(R)
    a_s = sp.sympify(a); k_s = sp.sympify(k); g_s = sp.sympify(gamma)

    # Base blocks
    Gsym = gradient_matrix_sym(th_s, R_s, a_s, k_s, g_s,
                               ctheta2=ctheta2, cR2=cR2, cphi2=cphi2,
                               epsG_theta_R=epsG_theta_R, epsG_theta_phi=epsG_theta_phi, epsG_R_phi=epsG_R_phi,
                               alpha_k2_theta=alpha_k2_theta, alpha_k2_R=alpha_k2_R, alpha_k2_phi=alpha_k2_phi)
    Msym = mass_matrix_sym(th_s, R_s, a_s, k_s, g_s, V0, Q=Q)
    A = Gsym * (k_s * k_s / (a_s * a_s)) + Msym

    # Kθθ used as a convenient scaling for couplings (dimensionless proxy)
    Kth = 3.0 / (th_s * th_s) + g_s * g_s
    # Columns: δN and B (with k/a factor in B-column couplings)
    bN = sp.Matrix([
        sp.sympify(lambda_N_theta) * Kth,
        sp.sympify(lambda_N_R) * 1,
        sp.sympify(lambda_N_phi) * (R_s * R_s),
    ])
    bB_core = sp.Matrix([
        sp.sympify(lambda_B_theta) * Kth,
        sp.sympify(lambda_B_R) * 1,
        sp.sympify(lambda_B_phi) * (R_s * R_s),
    ])
    bB = bB_core * (k_s / a_s)
    B = sp.Matrix.hstack(bN, bB)

    # Constraint quadratic form C (symmetric 2×2)
    cN2_s = sp.sympify(cN2)
    cB2_s = sp.sympify(cB2)
    cNB_s = sp.sympify(cNB)
    C = sp.Matrix([
        [cN2_s,          cNB_s * (k_s / a_s)],
        [cNB_s * (k_s / a_s),  cB2_s * (k_s * k_s / (a_s * a_s))],
    ])
    return A, B, C, a_s, k_s


def _split_A_into_GM(A: sp.Matrix, a_s: sp.Expr, k_s: sp.Expr) -> Tuple[sp.Matrix, sp.Matrix]:
    """
    Split A(=G k^2/a^2 + M) back into (G, M) assuming affinity in k^2.

    To avoid k→0 singularities in C^{-1} when constraints depend on (k/a),
    use a two-point evaluation at k=1 and k=2 (non-zero) and solve for G, M.
    """
    A1 = sp.simplify(A.subs({k_s: sp.Integer(1)}))
    A2 = sp.simplify(A.subs({k_s: sp.Integer(2)}))
    denom = (2**2 - 1**2)
    G = sp.simplify((A2 - A1) * (a_s * a_s) / denom)
    M = sp.simplify(A1 - G * (1 / (a_s * a_s)))
    return G, M


def eliminate_lapse_shift_sym(
    K: sp.Matrix,
    G: sp.Matrix,
    M: sp.Matrix,
    *,
    # Background and scale parameters (required for nontrivial elimination)
    theta: Optional[float | sp.Expr] = None,
    R: Optional[float | sp.Expr] = None,
    a: Optional[float | sp.Expr] = None,
    k: Optional[float | sp.Expr] = None,
    gamma: Optional[float | sp.Expr] = None,
    Q: Optional[float | sp.Expr] = None,
    V0: Optional[float | sp.Expr] = None,
    # Couplings (see _constraint_blocks_sym)
    lambda_N_theta: float | sp.Expr = 0.0,
    lambda_N_R: float | sp.Expr = 0.0,
    lambda_N_phi: float | sp.Expr = 0.0,
    lambda_B_theta: float | sp.Expr = 0.0,
    lambda_B_R: float | sp.Expr = 0.0,
    lambda_B_phi: float | sp.Expr = 0.0,
    cN2: float | sp.Expr = 1.0,
    cB2: float | sp.Expr = 1.0,
    cNB: float | sp.Expr = 0.0,
    ctheta2: float | sp.Expr = 1.0,
    cR2: float | sp.Expr = 1.0,
    cphi2: float | sp.Expr = 1.0,
    epsG_theta_R: float | sp.Expr = 0.0,
    epsG_theta_phi: float | sp.Expr = 0.0,
    epsG_R_phi: float | sp.Expr = 0.0,
    alpha_k2_theta: float | sp.Expr = 0.0,
    alpha_k2_R: float | sp.Expr = 0.0,
    alpha_k2_phi: float | sp.Expr = 0.0,
) -> tuple[sp.Matrix, sp.Matrix, sp.Matrix]:
    """
    Eliminate algebraic ADM variables (δN, B) at quadratic order and return
    effective (K_eff, G_eff, M_eff) in the dynamical sector.

    Behavior:
    - If any of (theta, R, a, k, gamma, V0) is None, or all couplings λ_* are
      zero, return the inputs unchanged (backward-compatible no-op).
    - Otherwise compute A_eff = A − B C^{-1} Bᵀ and split back into (G_eff, M_eff).
    - K is left unchanged (no q̇–constraint couplings included by default).
    """
    # Quick guard for no-op: missing context or all-lambda-zero
    if (
        theta is None
        or R is None
        or a is None
        or k is None
        or gamma is None
        or V0 is None
    ):
        return K, G, M

    if all(
        float(v) == 0.0
        for v in [lambda_N_theta, lambda_N_R, lambda_N_phi, lambda_B_theta, lambda_B_R, lambda_B_phi]
    ):
        return K, G, M

    A, B, C, a_s, k_s = _constraint_blocks_sym(
        theta,
        R,
        a,
        k,
        gamma,
        (0 if Q is None else Q),
        V0,
        lambda_N_theta=lambda_N_theta,
        lambda_N_R=lambda_N_R,
        lambda_N_phi=lambda_N_phi,
        lambda_B_theta=lambda_B_theta,
        lambda_B_R=lambda_B_R,
        lambda_B_phi=lambda_B_phi,
        cN2=cN2,
        cB2=cB2,
        cNB=cNB,
        ctheta2=ctheta2,
        cR2=cR2,
        cphi2=cphi2,
        epsG_theta_R=epsG_theta_R,
        epsG_theta_phi=epsG_theta_phi,
        epsG_R_phi=epsG_R_phi,
        alpha_k2_theta=alpha_k2_theta,
        alpha_k2_R=alpha_k2_R,
        alpha_k2_phi=alpha_k2_phi,
    )
    # Schur complement
    Cinv = sp.simplify(C.inv())
    Aeff = sp.simplify(A - B * Cinv * B.T)
    Geff, Meff = _split_A_into_GM(Aeff, a_s, k_s)
    # Keep K unchanged in this model (no qdot–constraint couplings)
    return K, Geff, Meff


def effective_matrices_sigma_strict(
    theta: float | sp.Expr,
    R: float | sp.Expr,
    a: float | sp.Expr,
    k: float | sp.Expr,
    *,
    gamma: float | sp.Expr,
    Q: float | sp.Expr,
    V0: float | sp.Expr,
    H: float | sp.Expr | None = None,
    theta_dot: float | sp.Expr | None = None,
    R_dot: float | sp.Expr | None = None,
    Hdot: float | sp.Expr | None = None,
) -> Tuple[sp.Matrix, sp.Matrix, sp.Matrix]:
    """
    Strict σ-model completion: return (K,G,M) where
    - K = diag(γ²+3/θ², 1, R²)
    - G = K (luminal high-k)
    - M = V_pot + ΔM_grav (k⁰ gravitational corrections), using
      the same ΔM as in mass_matrix_sigma_total.

    This corresponds to eliminating lapse/shift in the minimal σ-model
    and retaining only k⁰ contributions to the potential block. It avoids
    any heuristic couplings and preserves the luminal high-k structure.
    """
    Ksym = kinetic_matrix_sym(theta, R, gamma)
    Gsym = gradient_matrix_sym(theta, R, a, k, gamma)
    # Compute M only if background derivatives are provided; otherwise return zeros (cs2-only use case)
    if H is None or theta_dot is None or R_dot is None or Hdot is None:
        Msym = sp.zeros(3, 3)
    else:
        M_num = mass_matrix_sigma_total(
            float(theta), float(R), float(a),
            gamma=float(gamma), Q=float(Q), V0=float(V0),
            H=float(H), theta_dot=float(theta_dot), R_dot=float(R_dot), Hdot=float(Hdot)
        )
        Msym = sp.Matrix([[M_num[0][0], M_num[0][1], M_num[0][2]],
                          [M_num[1][0], M_num[1][1], M_num[1][2]],
                          [M_num[2][0], M_num[2][1], M_num[2][2]]])
    return Ksym, Gsym, Msym


# ---------------------------------------------------------------------------
# Schur blocks (numeric) reproducing ΔM_grav in the strict σ-model

def schur_blocks_sigma_strict_numeric(
    theta: float,
    R: float,
    a: float,
    *,
    gamma: float,
    Q: float,
    V0: float,
    H: float,
    theta_dot: float,
    R_dot: float,
    Hdot: float,
) -> Tuple[Tuple[Tuple[float, ...], ...], Tuple[Tuple[float, ...], ...]]:
    """
    Construct numeric Schur blocks (B 3x2, C 2x2) such that
        A_eff = A - B C^{-1} B^T
    reproduces the k^0 gravitational correction ΔM_grav of the strict σ-model.

    We factor -ΔM as B C^{-1} B^T using the top-2 eigenpairs of (-ΔM) and a
    diagonal C with ±1 signs to accommodate indefinite cases.
    For Q=0 the ΔM built from the σ-model has rank ≤ 2 and the factorization
    is exact. For general Q the rank may be 3; then this returns the best
    rank-2 reconstruction in the eigen-subspace of the two largest |λ|.
    """
    # Base potential mass block (numeric)
    Mpot_rr = float(V0) + (0.0 if Q == 0.0 else 3.0 * (Q * Q) / (a ** 6 * (R ** 4)))
    Mpot = np.array([[0.0, 0.0, 0.0], [0.0, Mpot_rr, 0.0], [0.0, 0.0, 0.0]], dtype=float)
    Mfull = np.array(mass_matrix_sigma_total(theta, R, a, gamma=gamma, Q=Q, V0=V0, H=H, theta_dot=theta_dot, R_dot=R_dot, Hdot=Hdot))
    deltaM = Mfull - Mpot
    # Target is -ΔM in the Schur term
    target = -deltaM
    # Eigen-decomposition
    vals, vecs = np.linalg.eigh( target )
    # Pick two eigenpairs with largest absolute value
    order = np.argsort(np.abs(vals))[::-1]
    idx = order[:2]
    lam = vals[idx]
    V = vecs[:, idx]
    # Build B and C
    B = np.zeros((3, 2), dtype=float)
    Cinv_diag = np.zeros(2, dtype=float)
    for i in range(2):
        si = lam[i]
        vi = V[:, i]
        mag = float(np.sqrt(abs(si)))
        B[:, i] = mag * vi
        # Choose C so that (1/C_i) = sign(si)
        Cinv_diag[i] = (1.0 if si >= 0.0 else -1.0)
    # Return B and C where C is diagonal with entries = 1/Cinv_diag
    C = np.diag([ (1.0 if d == 0.0 else 1.0/d) for d in Cinv_diag ])
    return (tuple(tuple(float(x) for x in row) for row in B),
            tuple(tuple(float(x) for x in row) for row in C))


def effective_matrices_sigma_strict_schur(
    theta: float | sp.Expr,
    R: float | sp.Expr,
    a: float | sp.Expr,
    k: float | sp.Expr,
    *,
    gamma: float | sp.Expr,
    Q: float | sp.Expr,
    V0: float | sp.Expr,
    H: float | sp.Expr,
    theta_dot: float | sp.Expr,
    R_dot: float | sp.Expr,
    Hdot: float | sp.Expr,
) -> Tuple[sp.Matrix, sp.Matrix, sp.Matrix]:
    """
    Strict σ-model via Schur: build B,C numerically to reproduce ΔM_grav and
    form A_eff = A - B C^{-1} B^T, then split into (G_eff, M_eff). Returns
    (K, G_eff, M_eff) as sympy matrices. High-k structure remains luminal.
    """
    # Base K,G (luminal)
    Ksym = kinetic_matrix_sym(theta, R, gamma)
    Gsym = gradient_matrix_sym(theta, R, a, k, gamma)
    # Potential block without ΔM
    Mpot = mass_matrix_sym(theta, R, a, k, gamma, V0, Q=Q)
    # Symbolic Schur blocks to reproduce ΔM
    B, C = schur_blocks_sigma_strict_symbolic(theta, R, a, gamma=gamma, Q=Q, V0=V0, H=H, theta_dot=theta_dot, R_dot=R_dot, Hdot=Hdot)
    # Since B,C are k^0 (no k-dependence), Geff = G and Meff = Mpot - B C^{-1} B^T
    A_schur = sp.simplify(B * C.inv() * B.T)
    Meff = sp.simplify(Mpot - A_schur)
    return Ksym, Gsym, Meff


def effective_matrices_sigma_strict_compare(
    theta: float | sp.Expr,
    R: float | sp.Expr,
    a: float | sp.Expr,
    k: float | sp.Expr,
    *,
    gamma: float | sp.Expr,
    Q: float | sp.Expr,
    V0: float | sp.Expr,
    H: float | sp.Expr,
    theta_dot: float | sp.Expr,
    R_dot: float | sp.Expr,
    Hdot: float | sp.Expr,
) -> Tuple[Tuple[sp.Matrix, sp.Matrix, sp.Matrix], Tuple[sp.Matrix, sp.Matrix, sp.Matrix]]:
    """
    Convenience: return both strict paths (ΔM-route and Schur-route) for comparison/testing.
    """
    K1, G1, M1 = effective_matrices_sigma_strict(theta, R, a, k, gamma=gamma, Q=Q, V0=V0,
                                                 H=H, theta_dot=theta_dot, R_dot=R_dot, Hdot=Hdot)
    K2, G2, M2 = effective_matrices_sigma_strict_schur(theta, R, a, k, gamma=gamma, Q=Q, V0=V0,
                                                       H=H, theta_dot=theta_dot, R_dot=R_dot, Hdot=Hdot)
    return (K1, G1, M1), (K2, G2, M2)


def schur_blocks_sigma_strict_symbolic(
    theta: float | sp.Expr,
    R: float | sp.Expr,
    a: float | sp.Expr,
    *,
    gamma: float | sp.Expr,
    Q: float | sp.Expr,
    V0: float | sp.Expr,
    H: float | sp.Expr,
    theta_dot: float | sp.Expr,
    R_dot: float | sp.Expr,
    Hdot: float | sp.Expr,
) -> Tuple[sp.Matrix, sp.Matrix]:
    """
    Return symbolic Schur blocks (B,C) that reproduce ΔM_grav in the strict σ-model.

    Construction:
      v ≡ G · φ̇,  φ̇ = (θ̇, Ṙ, Q/(a^3 R^2)),  G = diag(γ²+3/θ², 1, R²)
      d ≡ (1/H) D_t v  (with D_t v components matching mass_matrix_sigma_total)
      C^{-1} = [[α, β],[β, 0]],  α = (3H − Ḣ/H), β = 1/H
      B = [v, d]
    Then:  B C^{-1} B^T = α v v^T + β (v d^T + d v^T) = −ΔM_grav.
    """
    th = sp.sympify(theta); Rs = sp.sympify(R); a_s = sp.sympify(a)
    g = sp.sympify(gamma); Qs = sp.sympify(Q); V0s = sp.sympify(V0)
    Hs = sp.sympify(H); thd = sp.sympify(theta_dot); Rd = sp.sympify(R_dot); Hds = sp.sympify(Hdot)
    # Metric factors and velocities
    Gth = g * g + sp.Integer(3) / (th * th)
    Gpp = Rs * Rs
    phid = Qs / (a_s ** 3 * Rs ** 2)
    vth_l = Gth * thd
    vr_l = Rd
    vp_l = Gpp * phid  # = Q/a^3
    # Time derivatives (match mass_matrix_sigma_total)
    dGth_dth = -sp.Integer(6) / (th ** 3)
    dGpp_dR = sp.Integer(2) * Rs
    Vth = sp.Integer(0)
    VR = V0s * Rs - (Qs * Qs) / (a_s ** 6 * (Rs ** 3))
    Vp = sp.Integer(0)
    vth = thd
    vr = Rd
    vp = phid
    term_th = dGth_dth * vth * vth
    term_R = sp.Integer(0)
    term_p = dGpp_dR * vr * vp
    Dt_vth_l = -sp.Integer(3) * Hs * vth_l - Vth + term_th
    Dt_vr_l = -sp.Integer(3) * Hs * vr_l - VR + term_R
    Dt_vp_l = -sp.Integer(3) * Hs * vp_l - Vp + term_p
    # Build v and d = D_t v (no 1/H here)
    v = sp.Matrix([vth_l, vr_l, vp_l])
    d = sp.Matrix([Dt_vth_l, Dt_vr_l, Dt_vp_l])
    # C^{-1} and C
    alpha = (sp.Integer(3) * Hs - Hds / Hs)
    beta = sp.Integer(1) / Hs
    Cinv = sp.Matrix([[alpha, beta], [beta, sp.Integer(0)]])
    C = sp.simplify(Cinv.inv())
    B = sp.Matrix.hstack(v, d)
    return B, C

def effective_matrices_sym(
    theta: float | sp.Expr,
    R: float | sp.Expr,
    a: float | sp.Expr,
    k: float | sp.Expr,
    gamma: float | sp.Expr,
    Q: float | sp.Expr,
    V0: float | sp.Expr,
    *,
    ctheta2: float | sp.Expr = 1.0,
    cR2: float | sp.Expr = 1.0,
    cphi2: float | sp.Expr = 1.0,
    epsK_theta_R: float | sp.Expr = 0.0,
    epsK_theta_phi: float | sp.Expr = 0.0,
    epsK_R_phi: float | sp.Expr = 0.0,
    epsG_theta_R: float | sp.Expr = 0.0,
    epsG_theta_phi: float | sp.Expr = 0.0,
    epsG_R_phi: float | sp.Expr = 0.0,
    alpha_k2_theta: float | sp.Expr = 0.0,
    alpha_k2_R: float | sp.Expr = 0.0,
    alpha_k2_phi: float | sp.Expr = 0.0,
    # Elimination couplings (default 0 → identity)
    lambda_N_theta: float | sp.Expr = 0.0,
    lambda_N_R: float | sp.Expr = 0.0,
    lambda_N_phi: float | sp.Expr = 0.0,
    lambda_B_theta: float | sp.Expr = 0.0,
    lambda_B_R: float | sp.Expr = 0.0,
    lambda_B_phi: float | sp.Expr = 0.0,
    cN2: float | sp.Expr = 1.0,
    cB2: float | sp.Expr = 1.0,
    cNB: float | sp.Expr = 0.0,
) -> Tuple[sp.Matrix, sp.Matrix, sp.Matrix]:
    """
    Convenience: build base K,G,M and apply lapse/shift elimination with the
    supplied couplings. Returns (K_eff, G_eff, M_eff) as sympy matrices.
    """
    Ksym = kinetic_matrix_sym(theta, R, gamma,
                              epsK_theta_R=epsK_theta_R, epsK_theta_phi=epsK_theta_phi, epsK_R_phi=epsK_R_phi)
    Gsym = gradient_matrix_sym(theta, R, a, k, gamma,
                               ctheta2=ctheta2, cR2=cR2, cphi2=cphi2,
                               epsG_theta_R=epsG_theta_R, epsG_theta_phi=epsG_theta_phi, epsG_R_phi=epsG_R_phi,
                               alpha_k2_theta=alpha_k2_theta, alpha_k2_R=alpha_k2_R, alpha_k2_phi=alpha_k2_phi)
    Msym = mass_matrix_sym(theta, R, a, k, gamma, V0, Q=Q)
    Ke, Ge, Me = eliminate_lapse_shift_sym(
        Ksym,
        Gsym,
        Msym,
        theta=theta,
        R=R,
        a=a,
        k=k,
        gamma=gamma,
        Q=Q,
        V0=V0,
        lambda_N_theta=lambda_N_theta,
        lambda_N_R=lambda_N_R,
        lambda_N_phi=lambda_N_phi,
        lambda_B_theta=lambda_B_theta,
        lambda_B_R=lambda_B_R,
        lambda_B_phi=lambda_B_phi,
        cN2=cN2,
        cB2=cB2,
        cNB=cNB,
        ctheta2=ctheta2,
        cR2=cR2,
        cphi2=cphi2,
        epsG_theta_R=epsG_theta_R,
        epsG_theta_phi=epsG_theta_phi,
        epsG_R_phi=epsG_R_phi,
        alpha_k2_theta=alpha_k2_theta,
        alpha_k2_R=alpha_k2_R,
        alpha_k2_phi=alpha_k2_phi,
    )
    return Ke, Ge, Me


def augment_with_sigma_sym(
    K: sp.Matrix,
    G: sp.Matrix,
    M: sp.Matrix,
    *,
    csigma2: float | sp.Expr = 1.0,
    msigma2: float | sp.Expr = 0.0,
    epsK_theta_sigma: float | sp.Expr = 0.0,
    epsK_R_sigma: float | sp.Expr = 0.0,
    epsK_phi_sigma: float | sp.Expr = 0.0,
    epsG_theta_sigma: float | sp.Expr = 0.0,
    epsG_R_sigma: float | sp.Expr = 0.0,
    epsG_phi_sigma: float | sp.Expr = 0.0,
    epsM_theta_sigma: float | sp.Expr = 0.0,
    epsM_R_sigma: float | sp.Expr = 0.0,
    epsM_phi_sigma: float | sp.Expr = 0.0,
) -> Tuple[sp.Matrix, sp.Matrix, sp.Matrix]:
    """
    Extend (θ,R,φ) matrices to (θ,R,φ,σ).

    Canonical V6 defaults:
      Kσσ=1, Gσσ=cσ², Mσσ=mσ².
    Optional symmetric mixings are exposed for future derivations.
    """
    if K.shape != (3, 3) or G.shape != (3, 3) or M.shape != (3, 3):
        raise ValueError("augment_with_sigma_sym expects 3x3 K,G,M")
    K4 = sp.Matrix([
        [K[0, 0], K[0, 1], K[0, 2], sp.sympify(epsK_theta_sigma)],
        [K[1, 0], K[1, 1], K[1, 2], sp.sympify(epsK_R_sigma)],
        [K[2, 0], K[2, 1], K[2, 2], sp.sympify(epsK_phi_sigma)],
        [sp.sympify(epsK_theta_sigma), sp.sympify(epsK_R_sigma), sp.sympify(epsK_phi_sigma), sp.Integer(1)],
    ])
    G4 = sp.Matrix([
        [G[0, 0], G[0, 1], G[0, 2], sp.sympify(epsG_theta_sigma)],
        [G[1, 0], G[1, 1], G[1, 2], sp.sympify(epsG_R_sigma)],
        [G[2, 0], G[2, 1], G[2, 2], sp.sympify(epsG_phi_sigma)],
        [sp.sympify(epsG_theta_sigma), sp.sympify(epsG_R_sigma), sp.sympify(epsG_phi_sigma), sp.sympify(csigma2)],
    ])
    M4 = sp.Matrix([
        [M[0, 0], M[0, 1], M[0, 2], sp.sympify(epsM_theta_sigma)],
        [M[1, 0], M[1, 1], M[1, 2], sp.sympify(epsM_R_sigma)],
        [M[2, 0], M[2, 1], M[2, 2], sp.sympify(epsM_phi_sigma)],
        [sp.sympify(epsM_theta_sigma), sp.sympify(epsM_R_sigma), sp.sympify(epsM_phi_sigma), sp.sympify(msigma2)],
    ])
    return K4, G4, M4


def effective_matrices_sym_v6(
    theta: float | sp.Expr,
    R: float | sp.Expr,
    a: float | sp.Expr,
    k: float | sp.Expr,
    gamma: float | sp.Expr,
    Q: float | sp.Expr,
    V0: float | sp.Expr,
    *,
    ctheta2: float | sp.Expr = 1.0,
    cR2: float | sp.Expr = 1.0,
    cphi2: float | sp.Expr = 1.0,
    csigma2: float | sp.Expr = 1.0,
    msigma2: float | sp.Expr = 0.0,
    epsK_theta_R: float | sp.Expr = 0.0,
    epsK_theta_phi: float | sp.Expr = 0.0,
    epsK_R_phi: float | sp.Expr = 0.0,
    epsG_theta_R: float | sp.Expr = 0.0,
    epsG_theta_phi: float | sp.Expr = 0.0,
    epsG_R_phi: float | sp.Expr = 0.0,
    alpha_k2_theta: float | sp.Expr = 0.0,
    alpha_k2_R: float | sp.Expr = 0.0,
    alpha_k2_phi: float | sp.Expr = 0.0,
    lambda_N_theta: float | sp.Expr = 0.0,
    lambda_N_R: float | sp.Expr = 0.0,
    lambda_N_phi: float | sp.Expr = 0.0,
    lambda_B_theta: float | sp.Expr = 0.0,
    lambda_B_R: float | sp.Expr = 0.0,
    lambda_B_phi: float | sp.Expr = 0.0,
    cN2: float | sp.Expr = 1.0,
    cB2: float | sp.Expr = 1.0,
    cNB: float | sp.Expr = 0.0,
    epsK_theta_sigma: float | sp.Expr = 0.0,
    epsK_R_sigma: float | sp.Expr = 0.0,
    epsK_phi_sigma: float | sp.Expr = 0.0,
    epsG_theta_sigma: float | sp.Expr = 0.0,
    epsG_R_sigma: float | sp.Expr = 0.0,
    epsG_phi_sigma: float | sp.Expr = 0.0,
    epsM_theta_sigma: float | sp.Expr = 0.0,
    epsM_R_sigma: float | sp.Expr = 0.0,
    epsM_phi_sigma: float | sp.Expr = 0.0,
) -> Tuple[sp.Matrix, sp.Matrix, sp.Matrix]:
    """
    Canonical V6 symbolic effective matrices for q=(δθ,δR,δφ,δσ).

    Base 3x3 sector is built by effective_matrices_sym and then augmented by
    augment_with_sigma_sym.
    """
    K3, G3, M3 = effective_matrices_sym(
        theta,
        R,
        a,
        k,
        gamma,
        Q,
        V0,
        ctheta2=ctheta2,
        cR2=cR2,
        cphi2=cphi2,
        epsK_theta_R=epsK_theta_R,
        epsK_theta_phi=epsK_theta_phi,
        epsK_R_phi=epsK_R_phi,
        epsG_theta_R=epsG_theta_R,
        epsG_theta_phi=epsG_theta_phi,
        epsG_R_phi=epsG_R_phi,
        alpha_k2_theta=alpha_k2_theta,
        alpha_k2_R=alpha_k2_R,
        alpha_k2_phi=alpha_k2_phi,
        lambda_N_theta=lambda_N_theta,
        lambda_N_R=lambda_N_R,
        lambda_N_phi=lambda_N_phi,
        lambda_B_theta=lambda_B_theta,
        lambda_B_R=lambda_B_R,
        lambda_B_phi=lambda_B_phi,
        cN2=cN2,
        cB2=cB2,
        cNB=cNB,
    )
    return augment_with_sigma_sym(
        K3,
        G3,
        M3,
        csigma2=csigma2,
        msigma2=msigma2,
        epsK_theta_sigma=epsK_theta_sigma,
        epsK_R_sigma=epsK_R_sigma,
        epsK_phi_sigma=epsK_phi_sigma,
        epsG_theta_sigma=epsG_theta_sigma,
        epsG_R_sigma=epsG_R_sigma,
        epsG_phi_sigma=epsG_phi_sigma,
        epsM_theta_sigma=epsM_theta_sigma,
        epsM_R_sigma=epsM_R_sigma,
        epsM_phi_sigma=epsM_phi_sigma,
    )


def effective_matrices_v6(
    theta: float,
    R: float,
    a: float,
    k: float,
    gamma: float,
    Q: float,
    V0: float,
    *,
    ctheta2: float = 1.0,
    cR2: float = 1.0,
    cphi2: float = 1.0,
    csigma2: float = 1.0,
    msigma2: float = 0.0,
    epsK_theta_R: float = 0.0,
    epsK_theta_phi: float = 0.0,
    epsK_R_phi: float = 0.0,
    epsG_theta_R: float = 0.0,
    epsG_theta_phi: float = 0.0,
    epsG_R_phi: float = 0.0,
    alpha_k2_theta: float = 0.0,
    alpha_k2_R: float = 0.0,
    alpha_k2_phi: float = 0.0,
    lambda_N_theta: float = 0.0,
    lambda_N_R: float = 0.0,
    lambda_N_phi: float = 0.0,
    lambda_B_theta: float = 0.0,
    lambda_B_R: float = 0.0,
    lambda_B_phi: float = 0.0,
    cN2: float = 1.0,
    cB2: float = 1.0,
    cNB: float = 0.0,
    epsK_theta_sigma: float = 0.0,
    epsK_R_sigma: float = 0.0,
    epsK_phi_sigma: float = 0.0,
    epsG_theta_sigma: float = 0.0,
    epsG_R_sigma: float = 0.0,
    epsG_phi_sigma: float = 0.0,
    epsM_theta_sigma: float = 0.0,
    epsM_R_sigma: float = 0.0,
    epsM_phi_sigma: float = 0.0,
) -> Tuple[Tuple[float, ...], ...]:
    """Numeric wrapper returning tuple-matrices for canonical V6 4x4 sector."""
    K4, G4, M4 = effective_matrices_sym_v6(
        theta,
        R,
        a,
        k,
        gamma,
        Q,
        V0,
        ctheta2=ctheta2,
        cR2=cR2,
        cphi2=cphi2,
        csigma2=csigma2,
        msigma2=msigma2,
        epsK_theta_R=epsK_theta_R,
        epsK_theta_phi=epsK_theta_phi,
        epsK_R_phi=epsK_R_phi,
        epsG_theta_R=epsG_theta_R,
        epsG_theta_phi=epsG_theta_phi,
        epsG_R_phi=epsG_R_phi,
        alpha_k2_theta=alpha_k2_theta,
        alpha_k2_R=alpha_k2_R,
        alpha_k2_phi=alpha_k2_phi,
        lambda_N_theta=lambda_N_theta,
        lambda_N_R=lambda_N_R,
        lambda_N_phi=lambda_N_phi,
        lambda_B_theta=lambda_B_theta,
        lambda_B_R=lambda_B_R,
        lambda_B_phi=lambda_B_phi,
        cN2=cN2,
        cB2=cB2,
        cNB=cNB,
        epsK_theta_sigma=epsK_theta_sigma,
        epsK_R_sigma=epsK_R_sigma,
        epsK_phi_sigma=epsK_phi_sigma,
        epsG_theta_sigma=epsG_theta_sigma,
        epsG_R_sigma=epsG_R_sigma,
        epsG_phi_sigma=epsG_phi_sigma,
        epsM_theta_sigma=epsM_theta_sigma,
        epsM_R_sigma=epsM_R_sigma,
        epsM_phi_sigma=epsM_phi_sigma,
    )
    return _to_tuple_matrix(K4), _to_tuple_matrix(G4), _to_tuple_matrix(M4)


# ---------------------------------------------------------------------------
# Default elimination couplings (pragmatic completion)

# Note: default elimination couplings were removed. Elimination now only
# accepts explicit coefficient dictionaries (see eliminate_lapse_shift_sym).


# ---------------------------------------------------------------------------
# Minimal sigma-model gravitational corrections to M (k^0 block)

def mass_matrix_sigma_total(
    theta: float,
    R: float,
    a: float,
    *,
    gamma: float,
    Q: float,
    V0: float,
    H: float,
    theta_dot: float,
    R_dot: float,
    Hdot: float,
) -> tuple[tuple[float, ...], ...]:
    """
    Return full M_ab = M_pot + ΔM_grav for (δθ, δR, δφ) in the minimal σ-model.

    - Base potential block: M_RR = V0 + 3 Q^2/(a^6 R^4), others 0.
    - Gravitational corrections (k^0):
      ΔM_ab = - (3H - Hdot/H) T_ab - (1/H) D_t(T_ab), with
        T_ab = φ̇_a φ̇_b,   φ̇_a = G_ab φ̇^b,
        φ̇^θ = theta_dot, φ̇^R = R_dot, φ̇^φ = Q/(a^3 R^2),
        G_ab = diag(γ^2 + 3/θ^2, 1, R^2).
    """
    if a <= 0.0:
        raise ValueError("Scale factor a must be positive")
    if R == 0.0:
        raise ValueError("R=0 invalid for Q dynamics")
    Gth = (gamma * gamma) + 3.0 / (theta * theta)
    Gpp = R * R
    phidot = Q / (a ** 3 * R * R)
    vth_l = Gth * theta_dot
    vr_l = R_dot
    vp_l = Gpp * phidot  # = Q/a^3
    Tthth = vth_l * vth_l
    TthR = vth_l * vr_l
    Tthp = vth_l * vp_l
    TRR = vr_l * vr_l
    TRp = vr_l * vp_l
    Tpp = vp_l * vp_l
    dGth_dth = -6.0 / (theta * theta * theta)
    dGpp_dR = 2.0 * R
    Vth = 0.0
    VR = V0 * R - (Q * Q) / (a ** 6 * (R ** 3))
    Vp = 0.0
    vth = theta_dot
    vr = R_dot
    vp = phidot
    term_th = dGth_dth * vth * vth
    term_R = 0.0
    term_p = dGpp_dR * vr * vp
    Dt_vth_l = -3.0 * H * vth_l - Vth + term_th
    Dt_vr_l = -3.0 * H * vr_l - VR + term_R
    Dt_vp_l = -3.0 * H * vp_l - Vp + term_p
    DtT_thth = 2.0 * vth_l * Dt_vth_l
    DtT_thR = Dt_vth_l * vr_l + vth_l * Dt_vr_l
    DtT_thp = Dt_vth_l * vp_l + vth_l * Dt_vp_l
    DtT_RR = 2.0 * vr_l * Dt_vr_l
    DtT_Rp = Dt_vr_l * vp_l + vr_l * Dt_vp_l
    DtT_pp = 2.0 * vp_l * Dt_vp_l
    if H == 0.0:
        raise ValueError("H=0 makes gravitational correction ill-defined")
    pref = (3.0 * H - (Hdot / H))
    dM_thth = - (pref * Tthth + (1.0 / H) * DtT_thth)
    dM_thR  = - (pref * TthR  + (1.0 / H) * DtT_thR)
    dM_thp  = - (pref * Tthp  + (1.0 / H) * DtT_thp)
    dM_RR   = - (pref * TRR   + (1.0 / H) * DtT_RR)
    dM_Rp   = - (pref * TRp   + (1.0 / H) * DtT_Rp)
    dM_pp   = - (pref * Tpp   + (1.0 / H) * DtT_pp)
    MRR_base = V0 + 3.0 * (Q * Q) / (a ** 6 * (R ** 4))
    M = (
        (0.0 + dM_thth, dM_thR, dM_thp),
        (dM_thR, MRR_base + dM_RR, dM_Rp),
        (dM_thp, dM_Rp, 0.0 + dM_pp),
    )
    return M
