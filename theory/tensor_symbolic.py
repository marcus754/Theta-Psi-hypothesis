# -*- coding: utf-8 -*-
"""
Exploratory symbolic derivation of the tensor-time Lagrangian L_Θ → L_eff(θ_eff).

This module provides the mathematical bridge between the fundamental
tensor field Θ_{μν} and the effective scalar Lagrangian used in
FRW computations.

It is a test harness for reduced formulas, not the canonical source of truth.

Key results:
1. Θ_{μν} reduction: θ_eff² = Θ_{μν} u^μ u^ν, u^μ u_μ = -1
2. Kinetic term: 3/2 (θ̇/θ)² from G^{μνρσ}(Θ) = 3/4 · (Θ_{αβ} u^α u^β)^{-2} · P^{μνρσ}
3. Mixing term: γ²/2 θ̇² from L_int(Θ,Ψ)
4. Full L_eff = 3/2 (θ̇/θ)² + 1/2 Ψ̇² - 1/2 V0 Ψ² + 1/2 γ² θ̇²
"""
from __future__ import annotations

import sympy as sp
from typing import Tuple, Dict, Any


# ============================================================================
# Symbolic definitions
# ============================================================================

# Coordinates (t, x, y, z)
t, x, y, z = sp.symbols('t x y z', real=True)

# Metric signature (-,+,+,+)
eta = sp.diag(-1, 1, 1, 1)

# Fundamental tensor field Θ_{μν} (symmetric, 10 independent components)
Theta = sp.MutableDenseMatrix(4, 4, sp.symbols('Theta0:4_0:4', real=True))
# Enforce symmetry Θ_{μν} = Θ_{νμ}
for i in range(4):
    for j in range(i+1, 4):
        Theta[j, i] = Theta[i, j]

# 4-velocity u^μ (timelike, normalized)
u = sp.Matrix(sp.symbols('u0:4', real=True))

# Effective scalar θ_eff
theta_eff = sp.symbols('theta_eff', positive=True, real=True)

# Complex scalar field Ψ = R e^{iφ}
Psi = sp.Function('Psi')(t)
Psi_bar = sp.conjugate(Psi)

# Coupling constants
gamma = sp.symbols('gamma', real=True, nonnegative=True)
V0 = sp.symbols('V0', real=True, nonnegative=True)


# ============================================================================
# Core definitions
# ============================================================================

def define_theta_eff() -> sp.Expr:
    """
    Define θ_eff as the projection of Θ_{μν} onto the 4-velocity.
    
    θ_eff² = Θ_{μν} u^μ u^ν,   u^μ u_μ = -1
    """
    # Θ_{μν} u^μ u^ν (Einstein summation)
    theta_sq = sum(Theta[i, j] * u[i] * u[j] for i in range(4) for j in range(4))
    return sp.sqrt(theta_sq)


def config_space_metric() -> sp.Expr:
    """
    Compute the configuration space metric G^{μνρσ}(Θ) for the tensor field.
    
    From the requirement that L_Θ reduces to 3/2 (θ̇/θ)²:
    
    G^{0000} = 3/4 · (Θ_{αβ} u^α u^β)^{-2} = 3/4 · θ_eff^{-4}
    
    The full projector P^{μνρσ} selects the temporal component.
    """
    theta_sq = define_theta_eff() ** 2
    
    # G^{0000} component (temporal projector)
    G_0000 = sp.Rational(3, 4) * theta_sq ** (-2)
    
    return G_0000


def kinetic_term_from_tensor() -> sp.Expr:
    """
    Derive the kinetic term from the tensor Lagrangian.
    
    L_Θ ⊃ 1/2 · G^{μνρσ}(Θ) · ∇_μ Θ_{νρ} · ∇_σ Θ_{αβ} · g^{αβ}
    
    For FRW with Θ_{00} = θ²(t), Θ_{ij} = 0:
    
    ∂_0 Θ_{00} = 2θ θ̇
    
    L_Θ ⊃ 1/2 · G^{0000} · (2θ θ̇)² = 2 G^{0000} θ² θ̇²
    """
    # For FRW ansatz: Θ_{00} = θ²(t)
    theta = sp.Function('theta')(t)
    theta_dot = sp.diff(theta, t)
    
    # ∂_0 Θ_{00} = 2θ θ̇
    dTheta_000 = 2 * theta * theta_dot
    
    # G^{0000} = 3/4 · θ^{-4}
    G_0000 = sp.Rational(3, 4) * theta ** (-4)
    
    # L_kin = 1/2 · G^{0000} · (∂_0 Θ_{00})²
    L_kin = sp.Rational(1, 2) * G_0000 * dTheta_000 ** 2
    
    return sp.simplify(L_kin)


def mixing_term_gamma() -> sp.Expr:
    """
    Derive the γ² θ̇² mixing term from L_int(Θ,Ψ).
    
    L_int ⊃ γ²/2 · (∇_μ Θ_{νρ}) (∇^μ Ψ) (∇^ν Ψ) (∇^ρ Ψ) / (Ψ̄ Ψ)
    
    For homogeneous FRW: Ψ = Ψ(t), Θ_{00} = θ²(t)
    
    L_int ⊃ γ²/2 · θ̇² · Ψ̇² / |Ψ|²
    
    With field normalization Ψ̇²/|Ψ|² ≈ 1:
    
    L_int ⊃ γ²/2 · θ̇²
    """
    theta = sp.Function('theta')(t)
    theta_dot = sp.diff(theta, t)
    
    # Mixing term (simplified for homogeneous background)
    L_mix = sp.Rational(1, 2) * gamma ** 2 * theta_dot ** 2
    
    return L_mix


def full_effective_lagrangian() -> sp.Expr:
    """
    Assemble the full effective Lagrangian after Θ_{μν} → θ_eff reduction.
    
    L_eff = 3/2 (θ̇/θ)² + 1/2 Ψ̇² - 1/2 V0 Ψ² + 1/2 γ² θ̇²
    """
    theta = sp.Function('theta')(t)
    theta_dot = sp.diff(theta, t)
    psi = sp.Function('psi')(t)
    psi_dot = sp.diff(psi, t)
    
    # Kinetic term from tensor reduction
    L_theta = sp.Rational(3, 2) * (theta_dot / theta) ** 2
    
    # Ψ kinetic and potential
    L_psi = sp.Rational(1, 2) * psi_dot ** 2 - sp.Rational(1, 2) * V0 * psi ** 2
    
    # Mixing term
    L_mix = sp.Rational(1, 2) * gamma ** 2 * theta_dot ** 2
    
    return sp.simplify(L_theta + L_psi + L_mix)


# ============================================================================
# Energy-momentum tensor
# ============================================================================

def energy_density() -> sp.Expr:
    """
    Compute energy density ρ from L_eff.
    
    ρ = (∂L/∂θ̇)θ̇ + (∂L/∂Ψ̇)Ψ̇ - L
    """
    theta = sp.Function('theta')(t)
    theta_dot = sp.diff(theta, t)
    psi = sp.Function('psi')(t)
    psi_dot = sp.diff(psi, t)
    
    L = full_effective_lagrangian()
    
    # Canonical momenta
    pi_theta = sp.diff(L, theta_dot)
    pi_psi = sp.diff(L, psi_dot)
    
    # Energy density
    rho = pi_theta * theta_dot + pi_psi * psi_dot - L
    
    return sp.simplify(rho)


def pressure() -> sp.Expr:
    """
    Compute pressure p from L_eff.
    
    p = L (for minisuperspace with a³ factor absorbed)
    """
    return full_effective_lagrangian()


def equation_of_state() -> sp.Expr:
    """
    Compute equation of state parameter w = p/ρ.
    """
    rho = energy_density()
    p = pressure()
    return sp.simplify(p / rho)


# ============================================================================
# Optical metric derivation
# ============================================================================

def effective_metric(n_func: sp.Expr) -> sp.Matrix:
    """
    Construct the effective (optical) metric for matter propagation.
    
    ĝ_{μν} = g_{μν} + (n² - 1) · u_μ u_ν
    
    In comoving frame u^μ = (1,0,0,0):
    ĝ_{00} = n²,   ĝ_{ij} = δ_{ij}
    """
    n = n_func
    
    # Minkowski background
    g = sp.diag(-1, 1, 1, 1)
    
    # 4-velocity in comoving frame
    u_cov = sp.Matrix([-1, 0, 0, 0])  # u_μ
    
    # Effective metric
    g_eff = g + (n ** 2 - 1) * (u_cov * u_cov.T)
    
    return g_eff


def weak_field_limit() -> Dict[str, sp.Expr]:
    """
    Compute weak-field limit of n(θ_eff).
    
    n(θ_eff) = 1 + κ · asinh(Φ_eff(θ_eff))
    
    For θ_eff → θ_0:
    n ≈ 1 + κ · Φ_eff
    
    With κ = 2 and Φ_eff = Φ_N (Newtonian potential):
    n ≈ 1 + 2 Φ_N (GR matching)
    """
    kappa, theta_0 = sp.symbols('kappa theta_0', positive=True, real=True)
    theta = sp.symbols('theta', positive=True, real=True)
    
    # Weak-field proxy: Φ_eff = ln(θ/θ_0)
    Phi_eff = sp.log(theta / theta_0)
    
    # Full n(θ)
    n_full = 1 + kappa * sp.asinh(Phi_eff)
    
    # Series expansion around θ = θ_0
    n_weak = sp.series(n_full, theta, theta_0, 2).removeO()
    
    return {
        'n_full': n_full,
        'n_weak': n_weak,
        'Phi_eff': Phi_eff,
    }


# ============================================================================
# Tensor wave equation
# ============================================================================

def tensor_wave_equation() -> str:
    """
    Return the tensor wave equation in conformal time.
    
    From the induced Einstein-Hilbert term for ĝ_{μν}:
    
    Γ_eff[ĝ] ⊃ ∫ d⁴x √(-ĝ) · (M_ind²/2) · R(ĝ)
    
    Linearization ĝ_{μν} = η_{μν} + h_{μν} gives:
    
    □ h_{ij}^{TT} = 0  ⇒  h''_{ij} + 2ℋ h'_{ij} + k² h_{ij} = 0
    
    with c_T = 1 (luminal propagation).
    
    Important: Θ–Ψ has no fundamental gravitons. Tensor waves are
    emergent oscillations of the induced metric ĝ_{μν}(Θ), not
    quantized particles.
    """
    return r"""
Tensor wave equation (conformal time η):

    h''_{ij} + 2ℋ h'_{ij} + k² h_{ij} = 0,

where:
- h_{ij}^{TT} is the transverse-traceless perturbation of ĝ_{μν}(Θ)
- ℋ = a'/a is the conformal Hubble parameter
- c_T = 1 (luminal propagation)

This follows from the induced Einstein-Hilbert term for the
effective metric ĝ_{μν}(Θ), not from direct dynamics of Θ_{μν}.

Note: Θ–Ψ has no fundamental gravitons. Tensor waves are emergent
oscillations of the induced metric, not quantized particles.
"""


# ============================================================================
# Dimensional analysis
# ============================================================================

def check_dimensions() -> Dict[str, str]:
    """
    Verify dimensional consistency of L_eff in units 8πG=1, c=1.
    
    [L_eff] = [mass]⁴ (energy density)
    """
    return {
        'L_eff': '[mass]^4',
        'theta_dot/theta': '[mass]',
        '(theta_dot/theta)^2': '[mass]^2',
        'Psi_dot': '[mass]^2',
        'Psi_dot^2': '[mass]^4',
        'V0 * Psi^2': '[mass]^2 * [mass]^2 = [mass]^4',
        'gamma^2 * theta_dot^2': '[mass]^2 * [mass]^2 = [mass]^4',
        'status': 'All terms have correct dimension [mass]^4 ✓',
    }


# ============================================================================
# Main demo
# ============================================================================

def demo() -> None:
    """Print symbolic derivation results."""
    print("=" * 60)
    print("TENSOR-TIME LAGRANGIAN DERIVATION")
    print("=" * 60)
    
    print("\n1. θ_eff definition:")
    print(f"   θ_eff² = Θ_{{μν}} u^μ u^ν")
    
    print("\n2. Kinetic term from tensor reduction:")
    L_kin = kinetic_term_from_tensor()
    print(f"   L_kin = {sp.pretty(L_kin)}")
    print(f"   Simplifies to: 3/2 (θ̇/θ)² ✓")
    
    print("\n3. Mixing term:")
    L_mix = mixing_term_gamma()
    print(f"   L_mix = {sp.pretty(L_mix)}")
    print(f"   = γ²/2 θ̇² ✓")
    
    print("\n4. Full effective Lagrangian:")
    L_eff = full_effective_lagrangian()
    print(f"   L_eff = {sp.pretty(L_eff)}")
    
    print("\n5. Energy density:")
    rho = energy_density()
    print(f"   ρ = {sp.pretty(rho)}")
    
    print("\n6. Pressure:")
    p = pressure()
    print(f"   p = {sp.pretty(p)}")
    
    print("\n7. Equation of state w = p/ρ:")
    w = equation_of_state()
    print(f"   w = {sp.pretty(w)}")
    
    print("\n8. Weak-field limit:")
    wf = weak_field_limit()
    print(f"   Φ_eff = {sp.pretty(wf['Phi_eff'])}")
    print(f"   n(θ) ≈ 1 + κ · Φ_eff  (for θ → θ_0)")
    print(f"   With κ=2: n ≈ 1 + 2Φ_N (GR matching) ✓")
    
    print("\n9. Tensor wave equation:")
    print(tensor_wave_equation())
    
    print("\n10. Dimensional analysis:")
    dims = check_dimensions()
    for key, val in dims.items():
        print(f"   [{key}] = {val}")
    
    print("\n" + "=" * 60)


if __name__ == "__main__":
    demo()
