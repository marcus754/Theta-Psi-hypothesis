03. SCALAR PERTURBATIONS (ADM)

TIME CANON

- The fundamental object of time in the hypothesis is the Theta_munu tensor (see docs/11_derivation_path.md).
- The formulas in this section use the reduced FRW variable theta ≡ theta_eff, where theta_eff^2 = Theta_munu u^mu u^nu, u^mu u_mu = -1.
- All notations of the form n(theta) in this file should be read as n(theta_eff, Theta, Psi).
- This section describes the computational reduced circuit rather than the full tensor dynamics of Theta_munu.

CANONICAL FORMULATION (MANDATORY FOR ALL SECTIONS)

In the project, Theta_munu and Psi are taken as the base; the reduced circuit does not require an external threshold and does not use manual branch cutting.
S = ∫ d^4x √(-g) [ (M_pl^2/2)R + L_m + L_ThetaPsi
                   + M_*^4 F(J_refr) ],
J_refr(Theta,Psi) = Phi_eff(theta_eff) * I_grad(Theta,Psi).

Regime source (without Weyl threshold in the canon):
I_grad(Theta,Psi) = l_theta^2 * |∇_⊥theta|^2 + l_Psi^2 * |∇_⊥Psi|^2, ∇_⊥ = spatial part of the gradient (orthogonal to the local u^mu).

Meaning:
- In FRW (homogeneous background): I_grad=0, therefore J_refr=0.
- In strong field (large spatial gradients): I_grad>0, therefore the refractive sector is enabled directly via F(J_refr).
- There is no separate regime field; the canonical scalar/ADM circuit contains only perturbations of the reduced fields d_theta, d_R, deltaφ.

The vacuum index is defined by the finite working form:
n(theta_eff, Theta, Psi) = 1 + 𝒩(J_refr),
𝒩(0)=0, 𝒩(J) >= 0, 𝒩(J) is monotonic and finite.

The invalid background branch is excluded by parameterization (rather than a filter):
A = 3 tanh^2( sqrt(A_raw/3) ) in [0,3),
H^2 = (ρ_bg + V/2 + ρ_phase)/(3-A).

In comparison with GR, it is not a "manual inclusion" that is checked, but a direct response to spatial inhomogeneities of Theta_munu, Psi.

BRIEF OVERVIEW
- Perturbations q = (d_theta, d_R, deltaφ). Quadratic action: L^(2) = ½ dq̇ᵀKdq̇ − ½ dqᵀ[ G (k²/a²) + M ]dq.
- Minimal high-k model: K = G = diag(gamma² + 3/theta², 1, R²) ⇒ scalar modes propagate at the speed of light (c_s = 1).
- Masses: M_RR = V0 + 3 Q²/(a⁶ R⁴); M_thetatheta = M_φφ = 0; k⁰-corrections ΔM_ab(H, Ḣ, thetȧ, Ṙ, Q).
- Elimination of d_N, B via Schur: A_eff = A − B C⁻¹ Bᵀ; we reconstruct G_eff, M_eff.
- Health conditions: K≻0; eig(K⁻¹G_eff)>0 (high-k); ω²(k)>0 (finite-k); no superluminality.

Goal: matrices K(a,k), G(a,k), M(a,k) for (d_theta, d_R, deltaφ) after eliminating algebraic ADM variables (Lapse/Shift), with the ability to check ghosts/gradients and finite-k dispersion ω²(k).

Implementation
- Effective action without constraints:
  L^(2) = 1/2 · dq̇ᵀ K dq̇ − 1/2 · dqᵀ [ G (k²/a²) + M ] dq, dq = (d_theta, d_R, deltaφ).

- Elimination of Lapse/Shift (algebraic variables c = (d_N, B)) is introduced as
  L^(2) ⊃ − dqᵀ B c − 1/2 · cᵀ C c.

  Stationarity with respect to c gives c* = − C⁻¹ Bᵀ dq and the effective block
  A_eff = A − B C⁻¹ Bᵀ, where A = G (k²/a²) + M.

  Then we reconstruct G_eff and M_eff as parts of A_eff proportional to (k²/a²) and constant, respectively.
  In the minimal strict regime, we use an equivalent form without explicit mixing: K = diag(gamma² + 3/theta², 1, R²), G = K (luminal),
  M = M_pot + ΔM_grav(k⁰). ΔM is implemented by default via symbolic Schur: M_eff = M_pot − B C⁻¹ Bᵀ, where
  B = [G·φ̇, D_t(G·φ̇)], and C^-1 = [[3H − Ḣ/H, 1/H], [1/H, 0]]. Old APIs named sigma_strict are considered legacy names for this block.

- The diagonal block is now considered a minimal consistent implementation of the current action, not an approximation for fitting:
  K = G = diag(gamma² + 3/theta², 1, R²),
  M = diag(0, V0 + 3Q²/(a⁶R⁴), 0).
  Therefore, in the canon ε_K = ε_G = α_k² = 0. Non-zero mixings are allowed only as a diagnostic auxiliary regime or after adding explicit operators to the action.
- Numerically, generalized eigenvalues are calculated via Cholesky reduction (symmetrized K⁻¹G), see src/linear_modes_adm.py.

Working expressions (diagonal approximation, high-k):
- Kinetics K
  - Kthetatheta = 3/theta² + gamma²
  - KRR = 1
  - Kφφ = R²

- Gradients G (coefficients of k²/a²)
  - Gthetatheta = Kthetatheta (minimal luminal model, c_s = 1)
  - GRR = 1
  - Gφφ = R²

- Masses M
  - Mthetatheta = 0 (default)
  - MRR = V0 + 3 Q²/(a⁶ R⁴)
  - Mφφ = 0

Interface
- Symbolic/analytical expressions K, G, M are defined in theory/adm_symbolic.py. The old name elim='sigma_strict' is a legacy mode for the minimal Schur block.
- Fast checks and eigenvalues are implemented in src/linear_modes_adm.py:
  - cs2_eigs_adm(...) — c_s² as eigenvalues of the generalized problem G v = λ K v.
  - dispersion_eigs_adm(...) — ω²(k) from K⁻¹ [ G (k²/a²) + M ].
  - stability_adm(...), stability_finite_k_adm(...) — flags K≻0, c_s²>0, ω²>0.
  - The parameter elim=... is available everywhere to enable Lapse/Shift elimination.

Further work
- Full derivation without diagonal reduction (controlled mixings K_ij, G_ij, M_ij will appear).

PRACTICAL NOTES (STRICT PROFILE)

- In health checks with use_adm=True, finite k (ω²>0) along the trajectory is enabled by default.
- The strict Schur regime (elim='sigma_strict', legacy name) does not use heuristic "default" eliminations; only K=G and k⁰-corrections to masses.

ANALYTICAL CONDITIONS (APPROXIMATION)

- Ghosts (K≻0): without mixing, Kthetatheta>0, KRR=1, Kφφ=R²>0 are sufficient. With mixings ε_K — Sylvester's conditions for 2×2 and 3×3 minors are met.
- Gradients: eigenvalues of K⁻¹G are positive; without mixing, these are requirements c_theta², c_R², c_φ²>0. Subluminality is implemented as c_i²≤1.
- Finite k: ω² — eigenvalues of K⁻¹ [ G (k²/a²) + M ]. With V0≥0 and positive c_i², we get ω²(k)>0 for k>0.

Practical range (demo): c_i²∈(0,1], |ε_K|, |ε_G|≪1, V0∈[1e−3,1]. See also 05_consistency_checks.md.

Note on hyperbolicity and subluminality
- Hyperbolicity is ensured by the positive definiteness of K and the positivity of the eigenvalues of K⁻¹G (high-k) and K⁻¹[ G(k²/a²) + M ] (finite-k).
- Subluminality is set by the requirement max eig(K⁻¹G) ≤ 1 (with a small tol). These conditions correspond to the checks implemented in the code.

GRAVITATIONAL k^0-CORRECTIONS

For the minimal model (G_ab = diag(gamma² + 3/theta², 1, R²)), after eliminating d_N, B, a k^0-contribution appears in the mass block:
ΔM_ab = − (3H − Ḣ/H) T_ab − (1/H) D_t T_ab,
T_ab = φ̇_a φ̇_b, φ̇_a = G_ab φ̇^b,
φ̇^theta = thetȧ,
φ̇^R = Ṙ,
φ̇^φ = Q/(a³ R²).

The resulting M = M_pot + ΔM with M_RR^pot = V0 + 3Q²/(a⁶R⁴). In the code: theory/adm_symbolic.mass_matrix_sigma_total.
