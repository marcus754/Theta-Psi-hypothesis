08. FOUNDATIONS — CORE AND PROVABLE CONDITIONS

TIME CANON

- The fundamental object of time in the hypothesis is the Theta_munu tensor (see docs/11_derivation_path.md).
- The formulas in this section use the reduced FRW variable theta ≡ theta_eff, where theta_eff^2 = Theta_munu u^mu u^nu, u^mu u_mu = -1.
- All notations of the form n(theta) in this file should be read as n(theta_eff, Theta, Psi).
- This section describes the computational reduced circuit rather than the full tensor dynamics of Theta_munu.

CANONICAL FORMULATION (MANDATORY FOR ALL SECTIONS)

This file fixes specifically the mathematical bedrock of the hypothesis. In the current canon, only Theta_munu and Psi are primary, and the strong field response is built directly via a composite invariant

J_refr(Theta,Psi) = Phi_eff(theta_eff) * I_grad(Theta,Psi),
I_grad(Theta,Psi) = l_theta^2 * |∇_⊥theta|^2 + l_Psi^2 * |∇_⊥Psi|^2.

On a FRW background:
I_grad = 0, J_refr = 0.

The vacuum index is defined by the finite working form:
n(theta_eff, Theta, Psi) = 1 + 𝒩(J_refr),
𝒩(0)=0, 𝒩(J) >= 0, 𝒩(J) is monotonic and finite.

The invalid background branch is excluded by parameterization:
A = 3 tanh^2( sqrt(A_raw/3) ) in [0,3),
H^2 = (ρ_bg + V/2 + ρ_phase) / (3 - A).

Goal: to fix the "bedrock" — formulas and sufficient conditions upon which all phenomenology based on Theta_munu and Psi is built, with transparent code ↔ formula links.

UNITS AND SCALES
- Chosen 8πG=1, c=1 ⇒ H² = ρ_tot / 3. Conversion of H0 from km/s/Mpc: see 02_background.md.
- gamma is dimensionless in these units. Realistic regime: gamma ≲ 10⁻¹⁸.

HOMOGENEOUS LAGRANGIAN AND REGULARIZATION
Fields theta(t), Psi(t) = R e^iφ:

L_eff = 3/2 (thetȧ/theta)² + 1/2 Ṙ² − 1/2 V0 R² + 1/2 gamma² thetȧ²,
Q = a³ R² φ̇,
ρ_phase = p_phase = 1/2 · Q²/(a⁶ R²).

The equation of motion for theta has a pole at |theta| = √3/|gamma|. Regularization
theta = (√3/gamma) tanh x, thetȧ = (√3/gamma) sech²x · ẋ
eliminates the pole: |x|→∞ corresponds to the boundary in theta. See 02_background.md.

From a physics perspective, an even more transparent variable is the canonical logarithm χ ≡ √3 ln theta. Then the kinetics become simply ½χ̇², and the boundary theta→0 corresponds to χ→−∞ without special behavior in the Lagrangian. This transition fixes the normalization of the dynamics, because the coefficient 3/2 before (thetȧ/theta)² ensures a canonical metric on the logarithmic coordinate. One can freely work in this χ-parameterization or in x = artanh(gamma theta / √3) — both map the same physics, but χ clearly shows that there is no physical singularity, the coordinate simply goes to infinity. Symbolic and numerical formulas are provided in background/frw_symbolic.py, where the χ-interpretation is taken as the source, and the x-regularization remains convenient for numerical schemes.

CANONICAL FINITE PARAMETERIZATION (WITHOUT MANUAL BRANCH CUTTING)
To ensure the "undesirable" branch is not cut post-factum but is absent in the mathematical formulation itself, we use a finite working form of the effective load in the Friedmann equation:

A_raw = A_raw(x, y_N, u_N) >= 0,
χ_A = sqrt(A_raw / 3),
A = 3 tanh^2(χ_A) in [0, 3),
H^2 = (ρ_bg + V/2 + ρ_phase) / (3 - A).

Properties:
- at A_raw << 1, we get A ≈ A_raw (weak-field regime is not distorted);
- the region A >= 3 is excluded by construction, not by a "filter after calculation."

SCALAR PERTURBATIONS AND ADM ELIMINATION OF CONSTRAINTS
After linearization with respect to dq = (d_theta, d_R, deltaφ), we obtain the quadratic action

L^(2) = 1/2 · dq̇ᵀ K dq̇ − 1/2 · dqᵀ A dq,
A(a,k) = G (k²/a²) + M.

Lapse/shift (algebraic) variables c = (d_N, B) introduce
L^(2) ⊃ − dqᵀ B c − 1/2 · cᵀ C c,
and are eliminated by stationarity c* = − C⁻¹ Bᵀ dq ⇒

A_eff = A − B C⁻¹ Bᵀ,
G_eff and M_eff as parts of A_eff at (k²/a²) and the constant term.

Symbolic implementation: theory/adm_symbolic.py (effective_matrices_sym).
Numerical API: src/linear_modes_adm.py (elim=... parameter).

Base diagonal kernels (consistent with L_eff and Psi = R e^iφ):
Kthetatheta = 3/theta² + gamma², KRR = 1, Kφφ = R²;
Gthetatheta = c_theta² * Kthetatheta, GRR = c_R², Gφφ = c_φ² * R²;
MRR = V0, Mthetatheta = Mφφ = 0.

Mixings and α_k² adjustments are allowed parametrically and controlled in the interface (see theory/adm_symbolic.py).

SUFFICIENT CONDITIONS FOR HYPERBOLICITY AND SUBLUMINALITY
1) Absence of ghosts (K ≻ 0):
   - Without mixing: Kthetatheta > 0, KRR = 1 > 0, Kφφ = R² > 0 ⇒ K ≻ 0.
   - With mixing: Sylvester conditions on the principal minors of K.
2) Positivity of the "gradient operator" in the kinetic metric:
   - In high-k: positive eigenvalues of K⁻¹ G_eff.
   - At finite k: positivity of eigenvalues of K⁻¹ [ G_eff (k²/a²) + M_eff ].
3) Subluminality of scalar modes: max eig(K⁻¹ G_eff) ≤ 1 + tol.

Note on Schur complement:
If C ≻ 0, then A_eff = A − B C⁻¹ Bᵀ is symmetric. For sufficiently small ‖B‖ and positive diagonals G and M, the positivity of A_eff is preserved. This is reflected in stability tests with small mixings (see tests/test_linear_modes_adm.py and the new test with elim).

CONSISTENCY: WEAK FIELD AND TENSOR SECTOR
- The PPN layer in the current implementation is considered a weak-field diagnostic module. The interface checks/ppn_full.compute_ppn_params returns the current working configuration before the implementation of the full PN derivation.
- Tensor sector: α_T = 0 ⇒ c_T = 1. See docs/07_tensor_optical_ppn.md and tests/test_tensor_speed.py.

REALISTIC REGIME AND BACKGROUND INVARIANCE TO GAMMA
For gamma ≈ 10⁻¹⁸ and an initial state without dynamics (Psi ≈ 0, thetȧ ≈ 0), the background H(a) coincides with the native inert background limit within numerical precision. See tests/test_gamma_limit.py.

WEAK-FIELD MATCHING AS A BOUNDARY CONDITION
The optical index is introduced in a finite working form and simultaneously fixed by the weak field:

n(theta_eff) = 1 + 𝒩(Phi_eff(theta_eff)),
𝒩(0)=0, 𝒩 >= 0, 𝒩 is monotonic and finite,
where Phi_eff(theta_eff) >= 0 is a smooth weak-field working form.

At Phi_eff -> 0:
n(theta_eff) ≈ 1 + κ Phi_eff.

The canonical choice κ=2 sets the weak-field response normalization n ≈ 1 + 2Phi_eff (at c=1) as an observational boundary condition.

MINIMAL TWO-ENTITY MODEL FOR THE TIME FIELD AND SUBSTRATE
We choose the minimally-coupled action of two scalars in the reduced field-space:
- Field metric: G_thetatheta(theta) = gamma² + 3/theta², G_RR = 1, G_φφ = R².
- Gradients (high-k): c_s = 1 for all modes ⇒ G = K.
- Mass block: M_RR = V0 + 3 Q²/(a⁶ R⁴), M_thetatheta = M_φφ = 0 (for V(theta) = const, V(φ) is absent).

Eliminating d_N, B gives a standard quadratic functional with K = G_ab and G = G_ab on gradients; finite-k effects and all mass block dynamics are collected in M (see theory/adm_symbolic.py).

CAUSALITY
- The carrier of causality is the photon; c=1 fixes the light cone.
- The tensor sector is unmodified (EH) ⇒ c_T = 1.
- Scalar modes in high-k have c_s = 1 (K = G). Any case of c_s² > 1 is prohibited by veto.
- At finite k, k^0-corrections ΔM_ab are taken into account, which do not change the high-k speed.
