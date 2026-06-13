01. Interpretation and physical basis

Time Canon

- Within the hypothesis framework, the Theta_munu tensor is adopted as the fundamental object of time (see docs/11_derivation_path.md).
- In the formulas of this section, the reduced FRW variable theta ≡ theta_eff is used, where theta_eff^2 = Theta_munu u^mu u^nu, u^mu u_mu = -1.
- All records of the form n(theta) in this file should be read as a shorthand for the full response n(theta_eff, Theta, Psi).
- This section describes the computational reduced circuit, not the full tensor dynamics of Theta_munu.

Canonical formulation (mandatory for all sections)

In the project, Theta_munu and Psi are accepted as the basis; the reduced circuit does not require an external threshold and does not use manual branch cutting.
S = ∫ d^4x √(-g) [ (M_pl^2/2)R + L_m + L_ThetaPsi + L_refr(Theta,Psi) ].

Source of the regime (without Weyl threshold in the canon):
I_grad(Theta,Psi) = l_theta^2 · |∇_⊥theta|^2 + l_Psi^2 · |∇_⊥Psi|^2, ∇_⊥ = spatial part of the gradient (orthogonal to local u^mu).

Meaning:
- In FRW (homogeneous background): I_grad=0, so the refractive sector receives no source and remains switched off.
- In strong field (large spatial gradients): I_grad>0, and the refractive sector turns on automatically as a direct consequence of the interaction between Theta and Psi.
- Time in the canon is understood as a physical object; therefore, the allowed branch must preserve the positivity of the temporal regime and not allow its physical vanishing.

Working reduced implementation of this principle:
J_refr(Theta,Psi)=Phi_eff(theta_eff)·I_grad(Theta,Psi),
n(theta_eff,Theta,Psi)=1+𝒩(J_refr).

The invalid background branch is excluded by parameterization (rather than a filter):
A = 3 tanh^2( sqrt(A_raw/3) ) in [0,3),
H^2 = (ρ_bg + V/2 + ρ_phase)/(3-A).

The regime transition is understood as a direct consequence of the Theta and Psi interaction, not as an external switch.
The section fixes the physical meaning of the hypothesis based on Theta_munu, Psi. Observational checks and comparisons with external theories belong to a separate observational layer.

---

Axioms and entities

The goal of the section is to fix a minimal set of principles and an object model, without which the Theta–Psi hypothesis loses unambiguousness.

1) Time field Theta
- In the working version of the hypothesis, it is assumed that Theta determines the local rate of causality. The speed of light c is interpreted as the limiting speed of propagation of causal perturbations in Theta.
- Restriction: |theta| < theta⋆, where theta⋆ = √3 / |gamma|. As the boundary is approached, the local rate of processes asymptotically slows down.

2) Energy substrate Psi — complex
- Psi = R e^iφ. Amplitude R ≥ 0, phase φ ∈ [0,2π).
- U(1) invariant: Q = a³ R² φ̇ — constant of motion ("phase charge").

3) Asymmetry Theta↔Psi and the arrow of time
- In the working version of the hypothesis, Theta and Psi are asymmetric: Theta sets the causality background, Psi fills the relief with energy.
- It is assumed that the 3H friction in the background equations generates a macroscopic arrow of time.

4) Unified constant gamma
- gamma is introduced as a dimensionless coupling of time and energy; physically gamma ≲ 10⁻¹⁸ (see 02_background.md).

5) Regimes as limits
- In the working reading of the hypothesis, the Newtonian limit, "dark" effects, and quantum wave nature are viewed as different regimes of a single system (Theta, Psi).

Analogue interpretation (vacuum refraction)
- Geometry is viewed as a possible language for an effective description, not as the primary ontology of the hypothesis.
- The hypothesis uses a physical picture of vacuum refraction: it is assumed that a single scalar n(theta_eff) determines the local passage of time dtau = dt / n(theta_eff).
- The only physical propagation speed is the speed of light, c = 1. No "two speeds". The inhomogeneity of ln n(theta_eff) (its gradient) provides gravitational effects.
- Classic consequences:
  a = ∇ ln n(theta_eff),; t = ∫ n(theta_eff) dl,; Δnu/nu ≈ Δ ln n(theta_eff),; α ≈ ∫ ∇_⊥ ln n(theta_eff) dl.

What exactly is the source of gravity

- The source of gravity in Theta–Psi is the spatial inhomogeneity of the temporal regime set by the Theta_munu field.
- In the reduced circuit, this is expressed through the inhomogeneity of theta_eff(x) and the derivative vacuum response n(theta_eff(x)).
- Therefore, the observed "force of gravity" is not a separate fundamental force, but the response of bodies to the gradient of the local passage of causal processes:
a = ∇ ln n(theta_eff).

- In this sense, mass/energy does not create a force directly, but through feedback with Theta creates an inhomogeneous causal-temporal relief.

What exactly is refracted

- Not only a ray of light is refracted.
- In the strict language of the hypothesis, the causal conductivity of the vacuum changes, that is, how the vacuum unfolds local processes at a given Theta.
- Light manifests this as bending and delay.
- Particles manifest this as a change in trajectory.
- Clocks manifest this as a change in the passage rate.

In other words, "vacuum refraction" means refraction not of a separate object, but of the regime of causal propagation itself.

Status of the current version of the hypothesis

At the current stage, the following have been completed:
- reduced action;
- background cosmology;
- scalar perturbations;
- stationary strong-field regime;
- basic stability checks;
- EPR/CHSH check.

Partially completed:
- observational checks;
- weak-field limit;
- tensor perturbations.

Not completed:
- full dynamics of the tensor temporal field;
- general quantum limit;
- origin of the response function from first principles.

---
