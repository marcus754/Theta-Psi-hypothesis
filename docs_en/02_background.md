02. BACKGROUND: UNITS, LAGRANGIAN, EQUATIONS

TIME CANON

- The fundamental object of time in the hypothesis is the Theta_munu tensor (see docs/11_derivation_path.md).
- The formulas in this section use the reduced FRW variable theta ≡ theta_eff, where theta_eff^2 = Theta_munu u^mu u^nu, u^mu u_mu = -1.
- All notations of the form n(theta) in this file should be read as a shorthand for the full response n(theta_eff, Theta, Psi).
- This section describes the computational reduced circuit rather than the full tensor dynamics of Theta_munu.

CANONICAL FORMULATION (MANDATORY FOR ALL SECTIONS)

This file pertains only to the FRW module. It uses the physical basis of the hypothesis but does not replace either the ontology of time or the strong field layer.

S = ∫ d^4x √(-g) [ (M_pl^2/2)R + L_m + L_ThetaPsi + L_refr(Theta,Psi) ].

In the current canon, only Theta_munu and Psi are primary. In the full notation, Psi is treated as a complex order parameter Psi = R e^iφ; in the FRW reduction, the amplitude R and the phase charge Q = a^3 R^2 dot φ are used, and |Psi|^2 is naturally read as the amplitude density. The strong field response is parameterized directly by a composite invariant

J_refr(Theta,Psi) = Phi_eff(theta_eff) * I_grad(Theta,Psi),
I_grad(Theta,Psi) = l_theta^2 * |∇_⊥theta|^2 + l_Psi^2 * |∇_⊥Psi|^2.

On a homogeneous FRW background:
I_grad = 0, J_refr = 0.

Therefore, the background cosmology itself should not include the strong field refractive branch. The working form of the index remains
n(theta_eff, Theta, Psi) = 1 + 𝒩(J_refr).

The invalid background branch is excluded by parameterization:
A = 3 tanh^2( sqrt(A_raw/3) ) in [0,3),
H^2 = (ρ_bg + V/2 + ρ_phase)/(3-A).

The section defines only the base units, the homogeneous Lagrangian, and the FRW equations used in all numerical checks of the background module.

---
