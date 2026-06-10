# -*- coding: utf-8 -*-
"""
Tests for EHT forward modeling (ray tracing through optical metric).
"""
import math
import pytest

from fitting.forward_eht import (
    RayTracingParams,
    EHTObservable,
    theta_profile_schwarzschild_like,
    n_of_r,
    effective_potential,
    find_photon_sphere,
    trace_ray,
    compute_shadow_radius,
    compute_eht_observables,
    eht_likelihood,
)


class TestThetaProfile:
    """Test θ(r) profile for compact objects."""
    
    def test_theta_at_center(self):
        """θ(r=0) should equal theta_c."""
        theta_c = 1.0
        theta = theta_profile_schwarzschild_like(0.0, theta_c, m2=10.0, lam=0.1)
        assert abs(theta - theta_c) < 1e-10
    
    def test_theta_decays_with_radius(self):
        """θ(r) should decrease with radius."""
        theta_c = 1.0
        theta_1 = theta_profile_schwarzschild_like(1.0, theta_c, m2=10.0, lam=0.1)
        theta_2 = theta_profile_schwarzschild_like(2.0, theta_c, m2=10.0, lam=0.1)
        assert theta_1 > theta_2
    
    def test_theta_positive(self):
        """θ(r) should be positive everywhere."""
        for r in [0.0, 0.1, 1.0, 5.0, 10.0]:
            theta = theta_profile_schwarzschild_like(r, theta_c=1.0, m2=10.0, lam=0.1)
            assert theta > 0
    
    def test_theta_zero_at_negative_r(self):
        """θ(r) should raise for negative r."""
        with pytest.raises(ValueError):
            theta_profile_schwarzschild_like(-1.0, 1.0, 10.0, 0.1)


class TestRefractiveIndex:
    """Test n(r) profile."""
    
    def test_n_at_center(self):
        """n(r=0) should grow with the central field."""
        n = n_of_r(0.0, theta_c=10.0, m2=10.0, lam=0.1,
                   theta_scale=1.0, kappa_wf=2.0)
        assert n > 1.0
    
    def test_n_at_infinity(self):
        """n(r→∞) should approach 1."""
        n = n_of_r(100.0, theta_c=1.0, m2=10.0, lam=0.1,
                   theta_scale=1.0, kappa_wf=2.0)
        assert abs(n - 1.0) < 0.1
    
    def test_n_monotonic(self):
        """n(r) should decrease with radius."""
        n_1 = n_of_r(1.0, theta_c=1.0, m2=10.0, lam=0.1,
                     theta_scale=1.0, kappa_wf=2.0)
        n_2 = n_of_r(2.0, theta_c=1.0, m2=10.0, lam=0.1,
                     theta_scale=1.0, kappa_wf=2.0)
        assert n_1 >= n_2
    
    def test_n_positive(self):
        """n(r) should stay finite and above unity."""
        for r in [0.0, 0.5, 1.0, 2.0, 5.0]:
            n = n_of_r(r, theta_c=1.0, m2=10.0, lam=0.1,
                       theta_scale=1.0, kappa_wf=2.0)
            assert n >= 1.0


class TestEffectivePotential:
    """Test effective potential for photon motion."""
    
    def test_potential_positive(self):
        """V_eff should be positive."""
        V = effective_potential(b=5.0, r=10.0, n_r=1.5)
        assert V > 0
    
    def test_potential_decreases_with_n(self):
        """V_eff should decrease with larger n."""
        V_1 = effective_potential(b=5.0, r=10.0, n_r=1.0)
        V_2 = effective_potential(b=5.0, r=10.0, n_r=2.0)
        assert V_1 > V_2
    
    def test_potential_decreases_with_r(self):
        """V_eff should decrease with larger r."""
        V_1 = effective_potential(b=5.0, r=5.0, n_r=1.5)
        V_2 = effective_potential(b=5.0, r=10.0, n_r=1.5)
        assert V_1 > V_2


class TestPhotonSphere:
    """Test photon sphere finding."""
    
    def test_photon_sphere_exists(self):
        """Photon sphere should exist at finite radius."""
        r_ps = find_photon_sphere(
            theta_c=1.0, m2=10.0, lam=0.1,
            theta_scale=1.0, kappa_wf=2.0,
        )
        assert r_ps > 1.0
        assert r_ps < 50.0
    
    def test_photon_sphere_stable(self):
        """Photon sphere radius should be stable under parameter changes."""
        r_ps_1 = find_photon_sphere(
            theta_c=1.0, m2=10.0, lam=0.1,
            theta_scale=1.0, kappa_wf=2.0,
        )
        r_ps_2 = find_photon_sphere(
            theta_c=1.1, m2=10.0, lam=0.1,
            theta_scale=1.0, kappa_wf=2.0,
        )
        # Should be close (within 50%)
        assert abs(r_ps_1 - r_ps_2) / r_ps_1 < 0.5


class TestRayTracing:
    """Test ray tracing."""
    
    def test_ray_misses_surface(self):
        """Ray with large b should miss the surface."""
        params = RayTracingParams(
            theta_c=1.0, m2=10.0, lam=0.1,
            r_min=1.0, r_max=1000.0,
        )
        result = trace_ray(b=20.0, params=params)
        assert not result.hit_surface
        assert result.r_min_approach > params.r_min
    
    def test_ray_hits_surface(self):
        """Ray with small b should hit the surface."""
        params = RayTracingParams(
            theta_c=1.0, m2=10.0, lam=0.1,
            r_min=1.0, r_max=1000.0,
        )
        result = trace_ray(b=0.5, params=params)
        assert result.hit_surface
    
    def test_deflection_angle_small_for_large_b(self):
        """Deflection should be small for large impact parameter."""
        params = RayTracingParams(
            theta_c=0.1, m2=1.0, lam=0.01,  # Weak field
            r_min=1.0, r_max=1000.0,
        )
        result = trace_ray(b=50.0, params=params)
        # Deflection should be modest (|α| < 1 rad)
        # Note: can be negative (deflection toward the mass)
        assert abs(result.alpha) < 1.0
    
    def test_redshift_positive(self):
        """Redshift should be positive (n_emit > n_obs)."""
        params = RayTracingParams(
            theta_c=1.0, m2=10.0, lam=0.1,
            r_min=1.0, r_max=1000.0,
        )
        result = trace_ray(b=5.0, params=params)
        assert result.z_redshift >= 0


class TestShadowRadius:
    """Test shadow radius computation."""
    
    def test_shadow_radius_positive(self):
        """Shadow radius should be positive."""
        params = RayTracingParams(
            theta_c=1.0, m2=10.0, lam=0.1,
        )
        r_shadow = compute_shadow_radius(params)
        assert r_shadow > 0
    
    def test_shadow_radius_scale(self):
        """Shadow radius should scale with compact object size."""
        params_1 = RayTracingParams(
            theta_c=1.0, m2=10.0, lam=0.1,
            r_min=1.0,
        )
        params_2 = RayTracingParams(
            theta_c=5.0, m2=10.0, lam=0.1,  # Much larger theta_c
            r_min=1.0,
        )
        r_1 = compute_shadow_radius(params_1)
        r_2 = compute_shadow_radius(params_2)
        # Larger theta_c should give comparable or larger shadow
        assert r_2 >= r_1 * 0.9


class TestEHTObservables:
    """Test full EHT observable computation."""
    
    def test_observables_complete(self):
        """Should compute all observables."""
        obs = compute_eht_observables(gamma=0.002, V0=0.01, Q=0.0)
        
        assert obs.shadow_radius > 0
        assert obs.shadow_radius_err > 0
        assert obs.z_surface >= 0
        assert obs.delay_proxy >= 0
        assert obs.n_peak >= 1.0
        assert obs.gamma == 0.002
        assert obs.V0 == 0.01
        assert obs.Q == 0.0
    
    def test_observables_with_bridge(self):
        """Should work with bridge coefficients."""
        bridge = {
            'theta_c_scale': 2.0,
            'm2_scale': 15.0,
            'lam_scale': 0.2,
        }
        obs = compute_eht_observables(gamma=0.002, V0=0.01, Q=0.0,
                                      bridge_coeffs=bridge)
        assert obs.shadow_radius > 0


class TestEHTLikelihood:
    """Test EHT likelihood computation."""
    
    def test_likelihood_at_best_fit(self):
        """Likelihood should be maximal at best fit."""
        obs = compute_eht_observables(gamma=0.002, V0=0.01, Q=0.0)
        
        # Best fit: observed = predicted
        logL_best = eht_likelihood(obs, obs.shadow_radius, obs.shadow_radius_err)
        
        # Off best fit: should be lower
        logL_off = eht_likelihood(obs, obs.shadow_radius * 1.5, obs.shadow_radius_err)
        
        assert logL_best > logL_off
    
    def test_likelihood_with_z_constraint(self):
        """Should penalize high redshift models."""
        obs = compute_eht_observables(gamma=0.002, V0=0.01, Q=0.0)
        
        # Without z constraint
        logL_no_z = eht_likelihood(obs, obs.shadow_radius, obs.shadow_radius_err,
                                   z_max=None)
        
        # With z constraint (satisfied)
        logL_with_z = eht_likelihood(obs, obs.shadow_radius, obs.shadow_radius_err,
                                     z_max=obs.z_surface * 2)
        
        # With z constraint (violated)
        logL_violated = eht_likelihood(obs, obs.shadow_radius, obs.shadow_radius_err,
                                       z_max=obs.z_surface / 2)
        
        assert logL_no_z >= logL_violated
        assert logL_with_z >= logL_violated


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
