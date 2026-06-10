# -*- coding: utf-8 -*-
"""
Ray tracing through the optical metric g̃_{μν}(θ) for EHT forward modeling.

This module computes observable quantities for compact objects:
- Photon ring radius (shadow size)
- Redshift z_surface from surface emission
- Time delay Δt (Shapiro-like)
- Deflection angle α

The optical metric is:
    g̃_{μν} = g_{μν} + (n²(θ) - 1) u_μ u_ν,

where n(θ) is the universal refractive index and u_μ is the 4-velocity
of the Θ-field rest frame.

For spherically symmetric static configurations:
    ds² = n²(r) dt² - dr² - r²(dθ² + sin²θ dφ²)

Null geodesics (ds²=0) give photon trajectories.
"""
from __future__ import annotations

import math
from dataclasses import dataclass
from typing import List, Tuple, Optional, Dict

from theory.optical_metric import diagnostic_refractive_index_from_theta


@dataclass
class RayTracingParams:
    """Parameters for ray tracing through optical metric."""
    
    # Background field parameters
    gamma: float = 0.002
    V0: float = 0.01
    Q: float = 0.0
    
    # Optical metric parameters
    theta_scale: float = 1.0
    kappa_wf: float = 2.0
    
    # Compact object profile parameters
    theta_c: float = 1.0  # Central θ value
    m2: float = 10.0     # Mass parameter in θ(r) profile
    lam: float = 0.1     # Self-interaction
    
    # Observer parameters
    r_obs: float = 1000.0  # Observer distance (in GM/c² units)
    i_obs: float = math.pi / 4  # Inclination angle
    
    # Integration parameters
    r_min: float = 1.0    # Minimum radius (surface)
    r_max: float = 1000.0  # Maximum radius
    n_steps: int = 1000   # Number of integration steps


@dataclass
class RayResult:
    """Result of single ray tracing."""
    
    # Impact parameter
    b: float
    
    # Deflection angle
    alpha: float
    
    # Redshift (if emitted from radius r_emit)
    z_redshift: float
    
    # Time delay (relative to straight line)
    delta_t: float
    
    # Closest approach radius
    r_min_approach: float
    
    # Did the ray hit the surface?
    hit_surface: bool


@dataclass
class EHTObservable:
    """EHT observables for a compact object."""
    
    # Shadow radius (in GM/c² units)
    shadow_radius: float
    
    # Shadow radius uncertainty
    shadow_radius_err: float
    
    # Surface redshift (if applicable)
    z_surface: float
    
    # Time delay proxy
    delay_proxy: float
    
    # Peak refractive index
    n_peak: float
    
    # Model parameters
    gamma: float
    V0: float
    Q: float


def theta_profile_schwarzschild_like(
    r: float,
    theta_c: float,
    m2: float,
    lam: float,
) -> float:
    """
    Approximate θ(r) profile for a compact object.
    
    This is a toy model inspired by scalar field configurations:
    
    θ(r) ≈ θ_c · exp(-m2 · r² / 2) / (1 + lam · r²)
    
    For r → 0: θ → θ_c
    For r → ∞: θ → 0
    
    Parameters:
    - r: radial coordinate (in GM/c² units)
    - theta_c: central field value
    - m2: mass parameter (controls exponential decay)
    - lam: self-interaction (controls power-law tail)
    
    Returns:
    - θ(r): field value at radius r
    """
    if r < 0:
        raise ValueError("r must be non-negative")
    
    # Exponential decay
    exp_factor = math.exp(-m2 * r * r / 2.0)
    
    # Power-law suppression
    pow_factor = 1.0 / (1.0 + lam * r * r)
    
    return theta_c * exp_factor * pow_factor


def n_of_r(
    r: float,
    theta_c: float,
    m2: float,
    lam: float,
    theta_scale: float,
    kappa_wf: float,
) -> float:
    """
    Refractive index profile n(r) for a compact object.

    Combines θ(r) profile with finite n(θ).
    """
    theta = theta_profile_schwarzschild_like(r, theta_c, m2, lam)
    return diagnostic_refractive_index_from_theta(
        theta,
        theta_scale=theta_scale,
        kappa_wf=kappa_wf,
        profile="asinh",
    )


def effective_potential(b: float, r: float, n_r: float) -> float:
    """
    Effective potential for photon radial motion.
    
    For the optical metric, the radial equation is:
    
    (dr/dλ)² + V_eff(r) = 1/b²,
    
    where b is the impact parameter and
    
    V_eff(r) = 1 / (n²(r) · r²)
    
    The photon sphere is at the maximum of V_eff.
    """
    if n_r <= 0 or r <= 0:
        return float('inf')
    return 1.0 / (n_r ** 2 * r ** 2)


def find_photon_sphere(
    theta_c: float,
    m2: float,
    lam: float,
    theta_scale: float,
    kappa_wf: float,
    r_min: float = 1.0,
    r_max: float = 50.0,
    tol: float = 1e-6,
) -> float:
    """
    Find the photon sphere radius (unstable circular orbit).
    
    The photon sphere is at the maximum of V_eff(r) = 1/(n²r²).
    
    Uses golden section search for the maximum.
    """
    phi = (1.0 + math.sqrt(5.0)) / 2.0  # Golden ratio
    
    def V_eff(r: float) -> float:
        n_r = n_of_r(r, theta_c, m2, lam, theta_scale, kappa_wf)
        return effective_potential(1.0, r, n_r)
    
    # Golden section search
    a, b = r_min, r_max
    c = b - (b - a) / phi
    d = a + (b - a) / phi
    
    for _ in range(100):
        if abs(b - a) < tol:
            break
        
        if V_eff(c) > V_eff(d):
            b = d
            d = a + (b - a) / phi
        else:
            a = c
            c = b - (b - a) / phi
    
    return (a + b) / 2.0


def trace_ray(
    b: float,
    params: RayTracingParams,
) -> RayResult:
    """
    Trace a single photon ray with impact parameter b.
    
    Integrates the null geodesic equation through the optical metric.
    
    The deflection angle is:
    
    α(b) = 2 · ∫_{r_min}^{∞} (dr / r² · √(1/b² - V_eff(r))) - π,
    
    where r_min is the closest approach (turning point).
    """
    r_min_approach = b  # Initial guess (straight line)
    
    # Find turning point: V_eff(r_min) = 1/b²
    # Use Newton-Raphson iteration
    def f(r: float) -> float:
        n_r = n_of_r(r, params.theta_c, params.m2, params.lam,
                     params.theta_scale, params.kappa_wf)
        return 1.0 / (n_r ** 2 * r ** 2) - 1.0 / (b ** 2)
    
    def df_dr(r: float, h: float = 1e-8) -> float:
        return (f(r + h) - f(r - h)) / (2 * h)
    
    # Newton-Raphson for turning point
    r_turn = max(b, params.r_min * 1.1)
    for _ in range(50):
        fr = f(r_turn)
        if abs(fr) < 1e-10:
            break
        dfr = df_dr(r_turn)
        if abs(dfr) < 1e-15:
            break
        r_turn = r_turn - fr / dfr
        r_turn = max(r_turn, params.r_min)
    
    r_min_approach = r_turn
    
    # Check if ray hits the surface
    hit_surface = r_min_approach <= params.r_min
    
    if hit_surface:
        # Ray hits the surface - compute redshift
        n_surface = n_of_r(params.r_min, params.theta_c, params.m2, params.lam,
                          params.theta_scale, params.kappa_wf)
        n_obs = n_of_r(params.r_obs, params.theta_c, params.m2, params.lam,
                      params.theta_scale, params.kappa_wf)
        z_redshift = n_surface / n_obs - 1.0
        
        # Time delay (approximate)
        delta_t = compute_time_delay(params.r_min, params.r_obs, params)
        
        return RayResult(
            b=b,
            alpha=0.0,  # Absorbed, no deflection
            z_redshift=z_redshift,
            delta_t=delta_t,
            r_min_approach=r_min_approach,
            hit_surface=True,
        )
    
    # Compute deflection angle using numerical integration
    alpha = compute_deflection_angle(b, r_min_approach, params)
    
    # Redshift for emission at r_min_approach
    n_emit = n_of_r(r_min_approach, params.theta_c, params.m2, params.lam,
                   params.theta_scale, params.kappa_wf)
    n_obs = n_of_r(params.r_obs, params.theta_c, params.m2, params.lam,
                  params.theta_scale, params.kappa_wf)
    z_redshift = n_emit / n_obs - 1.0
    
    # Time delay
    delta_t = compute_time_delay(r_min_approach, params.r_obs, params)
    
    return RayResult(
        b=b,
        alpha=alpha,
        z_redshift=z_redshift,
        delta_t=delta_t,
        r_min_approach=r_min_approach,
        hit_surface=False,
    )


def compute_deflection_angle(
    b: float,
    r_min: float,
    params: RayTracingParams,
) -> float:
    """
    Compute deflection angle using numerical integration.
    
    α = 2 · ∫_{r_min}^{r_max} (dr / (r² · √(1/b² - V_eff(r)))) - π
    """
    n_steps = params.n_steps
    dr = (params.r_max - r_min) / n_steps
    
    def integrand(r: float) -> float:
        n_r = n_of_r(r, params.theta_c, params.m2, params.lam,
                    params.theta_scale, params.kappa_wf)
        V = effective_potential(b, r, n_r)
        arg = 1.0 / (b ** 2) - V
        if arg <= 0:
            return 0.0
        return 1.0 / (r ** 2 * math.sqrt(arg))
    
    # Simpson's rule integration
    total = integrand(r_min) + integrand(params.r_max)
    
    for i in range(1, n_steps):
        r = r_min + i * dr
        weight = 4 if i % 2 == 1 else 2
        total += weight * integrand(r)
    
    integral = total * dr / 3.0
    
    alpha = 2.0 * integral - math.pi
    
    return alpha


def compute_time_delay(
    r_emit: float,
    r_obs: float,
    params: RayTracingParams,
) -> float:
    """
    Compute Shapiro-like time delay.
    
    Δt = ∫ (n(r) - 1) dr
    
    This is the extra time compared to propagation in vacuum (n=1).
    """
    n_steps = params.n_steps
    dr = (r_obs - r_emit) / n_steps
    
    total = 0.0
    for i in range(n_steps):
        r = r_emit + (i + 0.5) * dr
        n_r = n_of_r(r, params.theta_c, params.m2, params.lam,
                    params.theta_scale, params.kappa_wf)
        total += (n_r - 1.0)
    
    return total * dr


def compute_shadow_radius(params: RayTracingParams) -> float:
    """
    Compute the shadow radius (critical impact parameter).
    
    The shadow boundary is at the impact parameter b_c corresponding
    to the photon sphere.
    
    For the optical metric:
    
    b_c = r_ps · n(r_ps),
    
    where r_ps is the photon sphere radius.
    """
    r_ps = find_photon_sphere(
        params.theta_c,
        params.m2,
        params.lam,
        params.theta_scale,
        params.kappa_wf,
    )

    n_ps = n_of_r(r_ps, params.theta_c, params.m2, params.lam,
                 params.theta_scale, params.kappa_wf)
    
    return r_ps * n_ps


def compute_eht_observables(
    gamma: float,
    V0: float,
    Q: float,
    bridge_coeffs: Optional[Dict] = None,
) -> EHTObservable:
    """
    Compute EHT observables from background parameters.
    
    This is the main entry point for EHT forward modeling.
    
    Parameters:
    - gamma, V0, Q: background field parameters
    - bridge_coeffs: optional dictionary with bridge calibration
    
    Returns:
    - EHTObservable with shadow radius, redshift, delay, etc.
    """
    # Resolve bridge parameters
    if bridge_coeffs is None:
        bridge_coeffs = {}
    
    theta_c = bridge_coeffs.get('theta_c_scale', 1.0)
    m2 = bridge_coeffs.get('m2_scale', 10.0)
    lam = bridge_coeffs.get('lam_scale', 0.1)
    
    # Set up ray tracing parameters
    params = RayTracingParams(
        gamma=gamma,
        V0=V0,
        Q=Q,
        theta_c=theta_c,
        m2=m2,
        lam=lam,
    )
    
    # Compute shadow radius
    shadow_r = compute_shadow_radius(params)
    
    # Estimate uncertainty (numerical error)
    shadow_r_err = shadow_r * 1e-3  # 0.1% numerical uncertainty
    
    # Surface redshift
    n_surface = n_of_r(params.r_min, theta_c, m2, lam,
                       params.theta_scale, params.kappa_wf)
    n_obs = n_of_r(params.r_obs, theta_c, m2, lam,
                   params.theta_scale, params.kappa_wf)
    z_surface = n_surface / n_obs - 1.0
    
    # Time delay proxy
    delay = compute_time_delay(params.r_min, params.r_obs, params)
    
    # Peak refractive index
    n_peak = n_surface  # Maximum at center
    
    return EHTObservable(
        shadow_radius=shadow_r,
        shadow_radius_err=shadow_r_err,
        z_surface=z_surface,
        delay_proxy=delay,
        n_peak=n_peak,
        gamma=gamma,
        V0=V0,
        Q=Q,
    )


def eht_likelihood(
    obs: EHTObservable,
    shadow_obs: float,
    shadow_err: float,
    z_max: Optional[float] = None,
) -> float:
    """
    Compute EHT likelihood for observed shadow.
    
    Parameters:
    - obs: computed observables
    - shadow_obs: observed shadow radius (in GM/c²)
    - shadow_err: observational uncertainty
    - z_max: optional upper limit on surface redshift
    
    Returns:
    - log-likelihood value
    """
    # Shadow radius likelihood (Gaussian)
    chi2_shadow = ((obs.shadow_radius - shadow_obs) / shadow_err) ** 2
    logL_shadow = -0.5 * chi2_shadow
    
    # Redshift constraint (one-sided)
    logL_z = 0.0
    if z_max is not None:
        if obs.z_surface > z_max:
            # Penalize models with too high redshift
            logL_z = -0.5 * ((obs.z_surface - z_max) / 0.1) ** 2
    
    return logL_shadow + logL_z


def demo() -> None:
    """Demonstrate EHT forward modeling."""
    print("=" * 60)
    print("EHT FORWARD MODELING DEMO")
    print("=" * 60)
    
    # Default parameters (M87*-like)
    gamma = 0.002
    V0 = 0.01
    Q = 0.0
    
    print(f"\nInput parameters:")
    print(f"  γ = {gamma}")
    print(f"  V0 = {V0}")
    print(f"  Q = {Q}")
    
    # Compute observables
    obs = compute_eht_observables(gamma, V0, Q)
    
    print(f"\nEHT Observables:")
    print(f"  Shadow radius: {obs.shadow_radius:.3f} ± {obs.shadow_radius_err:.3f} GM/c²")
    print(f"  Surface redshift: z = {obs.z_surface:.3f}")
    print(f"  Time delay: Δt = {obs.delay_proxy:.3f} GM/c³")
    print(f"  Peak index: n_peak = {obs.n_peak:.3f}")
    
    # Compare with M87* observation
    # M87*: shadow radius ≈ 5.2 ± 0.9 GM/c² (EHT 2019)
    m87_shadow = 5.2
    m87_err = 0.9
    
    logL = eht_likelihood(obs, m87_shadow, m87_err)
    print(f"\nM87* comparison:")
    print(f"  Observed: {m87_shadow} ± {m87_err} GM/c²")
    print(f"  log-likelihood: {logL:.3f}")
    
    print("\n" + "=" * 60)


if __name__ == "__main__":
    demo()
