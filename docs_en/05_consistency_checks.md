05. HEALTH CRITERIA AND CUTOFFS

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
- In FRW (homogeneous background): I_grad=0, therefore J_refr=0 and the refractive sector receives no source.
- In strong field (large spatial gradients): I_grad>0, therefore the refractive sector is enabled directly as a function F(J_refr).
- There is no separate regime field σ in the canon. This removes the third entity and leaves only Theta_munu and Psi as primary.

The finiteness of the vacuum index is ensured by the model construction:
n(theta_eff, Theta, Psi) = 1 + 𝒩(J_refr).

Numerical regularization coefficients apply only to the integrator and are not recorded as a physical entity.

The invalid background branch is excluded by parameterization (rather than a filter):
A = 3 tanh^2( sqrt(A_raw/3) ) in [0,3),
H^2 = (ρ_bg + V/2 + ρ_phase)/(3-A).

In comparison with GR, it is not a "manual inclusion" that is checked, but a direct response F(J_refr) to spatial inhomogeneities of Theta_munu, Psi.

BASIC CHECKS (USED IN SCANS)
- K ≻ 0 (no ghosts)
- eig(K⁻¹G) > 0 (no gradient instabilities)
- Subluminality: c_s² ≤ 1 (with tol margin)
- Finite k: ω²(k) > 0 for all modes (if enabled)
- Early Universe (BBN): H(a) monotonically decreases for a ∈ [1e−10, 1e−6] (optional)

The default mode is the strict mode (only consequences of the Lagrangian and equations of motion). Heuristic filters (effective theory limit and PPN-proxy) are disabled by default and can be enabled explicitly for fast scans; they are not considered part of the model's physics.

The project includes a unified check aggregator checks/health.py, which takes a background trajectory and parameters and returns ok_* flags.

Note: the canonical ADM circuit uses fields q=(d_theta, d_R, deltaφ). Old sigma-extensions and strict sigma-model checks are not considered part of the canon and may remain only as outdated test frameworks until removed from the API.

---

APPENDIX: ANALYTICAL PARAMETER DOMAIN (INTEGRATED FROM FORMER SECTION 09)

ANALYTICAL HEALTH CONDITIONS AND PARAMETER DOMAIN

Simple and transparent conditions ensuring the "health" of the scalar sector under the current approximation are collected here, and the working parameter ranges are specified. They are useful as a priori constraints before comparison with data.

Notation: q=(d_theta, d_R, deltaφ), K, G, M from L^(2) = 1/2 dq̇ᵀKdq̇ − 1/2 dqᵀ[ G (k²/a²) + M ]dq.
Diagonal kernels:
Kthetatheta = 3/theta² + gamma², KRR = 1, Kφφ = R²,
Gthetatheta = c_theta² Kthetatheta, GRR = c_R², Gφφ = c_φ² R²,
MRR = V0, Mthetatheta = Mφφ = 0.

1) ABSENCE OF GHOSTS (K ≻ 0)

Without mixing (ε=0) — sufficient:
Kthetatheta > 0 (always for |theta| < √3/|gamma|), KRR = 1 > 0, Kφφ = R² > 0 ⇒ R ≠ 0.

With mixings (small ε_K): by Sylvester for 3×3

Kthetatheta > 0,
det [[Kthetatheta, εthetaR], [εthetaR, KRR]] > 0 ⇒ εthetaR² < Kthetatheta * KRR,
det K > 0 ⇒ |εthetaφ|, |εRφ| are sufficiently small (explicit formulas omitted).

In practice: keep |ε_K,ij| ≪ min(Kthetatheta, KRR, Kφφ).

2) ABSENCE OF GRADIENT INSTABILITY (EIGENVALUES OF K⁻¹G ARE POSITIVE)

Without mixing: eigenvalues c_s² = c_theta², c_R², c_φ² ⇒ c_theta² > 0, c_R² > 0, c_φ² > 0 are required.

Subluminality imposes c_i² ≤ 1 (with a small margin >1 to account for rounding in numerical calculations).

With mixings (small ε_G): similar to point 1, but for the matrix K⁻¹G. It is sufficient to require smallness of ε_G and positivity of diagonal c_i².

3) POSITIVITY OF DISPERSION AT FINITE k (ω² > 0)

ω² — eigenvalues of K⁻¹ [ G (k²/a²) + M ].

With a diagonal and V0 ≥ 0: the contribution of G (k²/a²) is positive, and M is non-negative ⇒ ω² > 0 ∀ k>0. At the point k→0, one mode (R) gains mass V0, the others remain zero (expected for theta and phase φ).

4) THETA REGULARIZATION AND LIMITS

The working domain for the background ensures |theta| < theta⋆ = √3/|gamma| through the substitution theta = (√3/gamma) tanh x. This eliminates the parameterization singularity and maintains Kthetatheta > 0.

5) EFFECTIVE THEORY LIMIT (BACKGROUND)

Rough domain of applicability of the effective theory for the background:
|thetȧ/theta| < Λ, |gamma thetȧ| < Λ, Λ ≪ 1 (in Planck units).

Typically in numerical runs, we take Λ ≈ 1e−1 as a hard cutoff.

6) PPN-PROXY (LOCAL TESTS)

Fast heuristic (not a strict limit):
ε_loc ≃ 1 / (gamma² theta² + 3) < ε_max.

The threshold ε_max is chosen as a soft filter (default 0.5 in code tests); strict PPN boundaries will require a full derivation of the weak-field effective description. This proxy should not be interpreted as a final result.

6.1) STATUS OF WEAK-FIELD CHECKS

- Reference: GPS clock correction, Hafele-Keating, Pound-Rebka, Shapiro, solar redshift, solar deflection, Mercury, Lense-Thirring, binary pulsar timing, LLR, atomic limit, twin-paradox, c_T=1.

7) FINAL WORKING DOMAIN "RECIPE"

- Background: |theta| < √3/|gamma|, |thetȧ/theta| < Λ, |gamma thetȧ| < Λ.
- Kinetics: Kthetatheta > 0, R ≠ 0; for ε_K ≠ 0 — check minors by Sylvester.
- Gradients: c_i² ∈ (0, 1] and small ε_G.
- Finite k: V0 ≥ 0 ⇒ ω²(k) > 0 for k > 0.
- Local tests: ε_loc < ε_max (proxy), to be refined in the future full PPN derivation.
- Satellite clocks and GPS: reference check with a published standard of about 38 mus/day.
- Laboratory and transport clocks: Pound-Rebka, Hafele-Keating.
- Orbital and relativistic dynamics: Shapiro, solar redshift, solar deflection, Mercury, Lense-Thirring, binary pulsar timing, lunar laser ranging (LLR).
- Atomic / chemical limit: in the weak-field limit n(theta_eff) → 1, gamma_PPN = 1, beta_PPN = 1; no additional shifts for standard relativistic chemistry.

Acceptable parameter ranges in code units for demonstration numerical runs only:
gamma ≈ O(1), V0 ∈ [1e−3, 1], Q ≈ 0 (or |Q| ≪ 1),
ctheta², cR², cphi² ∈ (0,1], |ε_K|, |ε_G| ≪ 1.

Matching with reality requires gamma_phys ≲ 1e−18; this condition is applied a priori in data analysis.
