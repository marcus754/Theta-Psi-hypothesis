04. QUANTUM GAUSSIAN SUPERSTRUCTURE

TIME CANON

- The fundamental object of time in the hypothesis is the Theta_munu tensor (see docs/11_derivation_path.md).
- The formulas in this section use the reduced FRW variable theta ≡ theta_eff, where theta_eff^2 = Theta_munu u^mu u^nu, u^mu u_mu = -1.
- All notations of the form n(theta) in this file should be read as n(theta_eff, Theta, Psi).
- This section describes the computational reduced circuit rather than the full tensor dynamics of Theta_munu.

CANONICAL FORMULATION (MANDATORY FOR ALL SECTIONS)

In the project, Theta_munu and Psi are taken as the base; the reduced circuit does not require an external threshold and does not use manual branch cutting.
S = ∫ d^4x √(-g) [ (M_pl^2/2)R + L_m + L_ThetaPsi + M_*^4 F(J_refr) ],
J_refr(Theta,Psi) = Phi_eff(theta_eff) * I_grad(Theta,Psi).

Regime source (without Weyl threshold in the canon):
I_grad(Theta,Psi) = l_theta^2 * |∇_⊥theta|^2 + l_Psi^2 * |∇_⊥Psi|^2,
∇_⊥ = spatial part of the gradient (orthogonal to the local u^mu).

Meaning:
- In FRW (homogeneous background): I_grad=0, therefore J_refr=0.
- In strong field (large spatial gradients): I_grad>0, therefore the refractive sector is enabled directly via F(J_refr).
- There is no separate regime field in the canon.

The vacuum index is defined by the finite working form:
n(theta_eff, Theta, Psi) = 1 + 𝒩(J_refr),
𝒩(0)=0, 𝒩(J) >= 0, 𝒩(J) is monotonic and finite.

The invalid background branch is excluded by parameterization (rather than a filter):
A = 3 tanh^2( sqrt(A_raw/3) ) in [0,3),
H^2 = (ρ_bg + V/2 + ρ_phase)/(3-A).

In comparison with GR, it is not a "manual inclusion" that is checked, but a direct response to spatial inhomogeneities of Theta_munu, Psi.

Idea: we quantize linear modes dq_k = (d_theta_k, d_R_k, deltaφ_k). The state is Gaussian, fully defined by the covariance matrix sum(t,k).

Quadratic Hamiltonian (after Lapse/Shift elimination):
L^(2) = 1/2 · dq̇ᵀ K dq̇ − 1/2 · dqᵀ A dq,
A = G (k²/a²) + M.

Switching to canonical variables x = (dq, π) with π = K dq̇, we obtain the Hamilton matrix
H = 1/2 · xᵀ H_x x, H_x = [[ A, 0 ], [ 0, K⁻¹ ]].

Canonics and HUP
[x_i, x_j] = i Ω_ij ⇒ sum + i Ω / 2 ≥ 0,
Ω = [[0, I], [−I, 0]].

Evolution of sum (linear):
Ẋ = Ω H_x X ⇒ suṁ = Ω H_x sum + sum H_xᵀ Ωᵀ.
This is an ODE for each mode (t,k) and allows for the calculation of squeezing, coherence, and entanglement of (k,−k) pairs.

Expected effects: squeezing at k/a ≲ H, growth of entanglement, stabilization in the k/a ≫ H regime. The source of squeezing is the non-adiabaticity of the background Theta–Psi regime.

Metrics
- Entanglement entropy (for mode pairs), measures of coherence/decoherence (by sum).
- Energy density of quantum fluctuations: ⟨H⟩ = 1/2 Tr(H_x sum) (attention to ultraviolet regularization: k-cutoff or additive counterterms).

CASIMIR AS A BOUNDARY CONSISTENCY CHECK

The Casimir effect in this framework does not add a new entity to the Theta_munu, Psi basis. It belongs to the derivative microphysical layer of the matter sector: boundaries change the allowed modes of the field χ, and vacuum energy appears as the difference of spectra for different boundary data.

Minimal local notation for the matter sector:
L_matter^micro =
  1/2 Z_t χ_dot² - 1/2 Z_s |∇χ|² - V_matter(χ).

For homogeneous response coefficients, the dispersion of boundary modes in the coordinates of local time slicing has the form:
ω_k = sqrt(Z_s/Z_t) |k|.

Therefore, for parallel plates at a distance a, the Casimir contribution scales by the same factor:
E_C/A = - pi²/(720 a³) * sqrt(Z_s/Z_t),
P_C = - pi²/(240 a⁴) * sqrt(Z_s/Z_t).

In the weak-field limit:
R_micro -> 0,
Z_t -> 1,
Z_s -> 1,
n² = Z_s/Z_t -> 1,

the standard result is obtained:
E_C/A = - pi²/(720 a³),
P_C = - pi²/(240 a⁴).

Meaning of the check: Casimir is a boundary quantum contribution of matter modes in the same Theta–Psi response layer. It should not include the FRW refractive sector by itself and should not introduce independent vacuum energy as a third ontology. The mathematical boundary consistency check is established in theory.micro_perturbation.casimir_boundary_registry().
