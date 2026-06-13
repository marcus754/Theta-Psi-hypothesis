00. Theta–Psi Overview with formulas (v2)

The file serves as a module map. It fixes only the basis of the hypothesis and the connections between layers:
- core: what is assumed about time, gravity, and refractive response;
- background: how the reduced FRW circuit is built on top of the core;
- observables: how observables are extracted from the background and strong field solutions.

Comparison with external models is not part of this file's physical basis. It belongs only to a separate observational layer.

Time Canon

- Within the hypothesis framework, the Theta_munu tensor is adopted as the fundamental object of time (see docs/11_derivation_path.md).
- In the formulas of this section, the reduced FRW variable theta ≡ theta_eff is used, where theta_eff^2 = Theta_munu u^mu u^nu, u^mu u_mu = -1.
- All records of the form n(theta) in this file should be read as n(theta_eff, Theta, Psi).
- This section describes the computational reduced circuit, not the full tensor dynamics of Theta_munu.

Canonical formulation (mandatory for all sections)

In the project, Theta_munu and Psi are accepted as the basis; the reduced circuit does not require an external threshold and does not use manual branch cutting.
S = ∫ d^4x √(-g) [ (M_pl^2/2)R + L_m + L_ThetaPsi + L_refr(Theta,Psi) ].

Source of the regime (without Weyl threshold in the canon):
I_grad(Theta,Psi) = l_theta² |∇_⊥theta|² + l_Psi² |∇_⊥Psi|²
∇_⊥ = spatial part of the gradient (orthogonal to local u^mu).

Meaning:
- In FRW (homogeneous background): I_grad=0, so the refractive sector receives no source and remains switched off.
- In strong field (large spatial gradients): I_grad>0, and the refractive sector turns on automatically as a direct consequence of the interaction between Theta and Psi.
- Time in the hypothesis is understood as a physical object; therefore, in the working canon, the allowed dynamics should not transition the temporal regime into zero or a negative branch.

Working reduced implementation of this principle:
J_refr(Theta,Psi)=Phi_eff(theta_eff)·I_grad(Theta,Psi),
n(theta_eff,Theta,Psi)=1+𝒩(J_refr).

The invalid background branch is excluded by parameterization (rather than a filter):
A = 3 tanh^2( sqrt(A_raw/3) ) in [0,3),
H^2 = (ρ_bg + V/2 + ρ_phase)/(3-A).

The regime transition is understood as a direct consequence of the Theta and Psi interaction, not as an external switch.

Core (seed)
- Units: 8πG=1, c=1. Space is flat; "geometry" is a language, not an ontology.
- Causal structure is determined by a common propagation cone; the only physical speed is c=1 (dx/dtau=1).
- Vacuum refraction: dtau = dt / n(theta_eff). This is an extension of coordinate time and optical path length, not a "second speed".
- Physical branch of the model: the time index remains finite on the realized branch; therefore, a true event horizon does not form.
- Effects of a single n(theta_eff): a = −∇ ln n; t = ∫ n(theta_eff) dl; Δnu/nu ≈ Δ ln n(theta_eff); α ≈ ∫ ∇_⊥ ln n(theta_eff) dl.
- Homogeneous L_eff: 3/2(thetȧ/theta)² + 1/2Psi̇² − 1/2V0Psi² + 1/2gamma²thetȧ²; Psi=R e^iφ, Q=a³R²φ̇.
- ADM high‑k: K=G=diag(gamma²+3/theta², 1, R²); M_RR = V0 + 3Q²/(a⁶R⁴); ΔM_ab — k⁰‑corrections (H, Ḣ, thetȧ, Ṙ, Q).
- Scalars are luminal in high‑k; superluminal is prohibited. Tensors: c_T=1.
- Compatibility: universality of n, isotropy, almost no dispersion; weak field alignment remains a condition on observables, not part of the ontology.
- Code is a verification tool; reference sources for defining the model are this document and background/frw_symbolic.py.

This file is a brief human summary of the hypothesis based on Theta_munu, Psi with reference formulas. The code serves as a verification tool, but the canonical sources for defining the model are considered to be the consistent records here and in background/frw_symbolic.py.

Modularity of the hypothesis

1. Hypothesis basis:
   Theta_munu, Psi, positivity of time, J_refr, n, g_hat_munu.
2. Gravitational layer: gravity as a response to the inhomogeneity of the temporal regime and the causal conductivity of the vacuum.
3. Background cosmology: a separate FRW reduction from the basis, but not the entire basis.
4. Strong field: a separate stationary and wave layer using the same J_refr.
5. Observables / fitting: a separate layer for comparison with data and external baselines.

Units and parameters
- Units: chosen such that the Friedmann equation takes the form H^2 = ρ_tot / 3 (i.e., 8πG=1, c=1).
- Model parameters: gamma — coupling of time and energy (in reality gamma_phys ≲ 10⁻¹⁸), V0 — quadratic potential of the substrate, Q — phase charge of complex Psi.
- H0 conversion from km/s/Mpc into these units: see background/components.py (omegas_from_km_s_Mpc).

Causality
- Causal structure is determined by a common propagation cone. The only physical speed is the speed of light, c=1.
- Gravity as vacuum refraction: the vacuum acquires an effective refractive index with the index n(theta_eff), which extends the coordinate time: dtau = dt / n(theta_eff). This is an extension of the path time, not a change in the maximum speed.
- Classic effects follow from this: t = ∫ n(theta_eff) dl (Shapiro delay), Δnu/nu ≈ Δ ln n(theta_eff) (redshift), α ≈ ∫ ∇_⊥ ln n(theta_eff) dl (lensing).
- The tensor sector is unmodified: c_T = 1. Scalar modes in high-k propagate at the speed of light (K=G ⇒ c_s = 1). Any superluminality (c_s²>1) is prohibited.

Interpretation (analogue picture)
- Geometry is a language. Physics is here: vacuum with a refractive index n(theta_eff) extends coordinate time (dtau = dt / n(theta_eff)).
- The only physical speed is c=1; refraction is the optical path length t = ∫ n(theta_eff) dl, not "another speed of light".
- All gravitational effects are reduced to the gradient of ln n(theta_eff): free fall, redshift, Shapiro delay, lensing.

Minimal working class of the refractive sector

To ensure that the transition
Theta, Psi -> J_refr(Theta,Psi) -> n(theta_eff,Theta,Psi) -> g_hat_munu

is fixed not only by words, the canon takes a minimal class
L_refr = M_*^4 · F(J_refr),
J_refr = Phi_eff(theta_eff) · I_grad(Theta,Psi),
Phi_eff(theta_eff) = ln(theta_eff / theta_0).

Meaning:
- Phi_eff — canonical weak field indicator of the temporal regime;
- I_grad measures the strength of the spatial inhomogeneity of the Theta and Psi interaction;
- F defines a universal vacuum response, the same for light, clocks, and massive bodies through a common effective metric g_hat_munu.

In working record, this is written as
n(theta_eff,Theta,Psi)=1+𝒩(J_refr).

Thus, n acts not as an arbitrary knob, but as a compressed record of the universal refractive response of the vacuum to the direct invariant of the Theta and Psi interaction.

Canonical normalized representative-form:
F_v1(J) = 1 - exp(-J).

It was chosen because it immediately gives:
- F(0)=0 and the absence of artificial background activation;
- F'(J)>0 and monotonic enhancement of the regime;
- F''(J)<0 and smooth saturation without a threshold;
- activation of the strong field branch directly through the growth of I_grad, without a third field.

But this is not yet a derivation from the action: the action fixes J_refr and the requirements for F, not a single function F(J). An explicit J_* is not needed here and is considered a normalization of J_refr.

Minimal lemmas of the canon:
- at I_grad=0, the refractive sector on FRW does not receive a source;
- at 𝒩'(J)>0 and finite saturation, always 1 <= n(J) < ∞;
- at finite energy and theta(t0)>0, the state theta=0 is not reached in finite time;
- universality of g_hat_munu makes free fall a consequence, not a separate postulate.

Where gravity comes from
- The source of gravity in this frame is not a separate force of attraction, but the spatial inhomogeneity of the temporal regime set by Theta_munu.
- In the reduced description, this is read as the inhomogeneity of theta_eff(x) and, consequently, n(theta_eff(x)).
- If n depends on position, then the local course of processes, frequencies, and optical path length depends on position as well. This very gradient gives the observed acceleration:
a = -∇ ln n(theta_eff).
- Therefore, gravity here is not a primary force entity, but a response of bodies and signals to the inhomogeneous causal-temporal structure of the vacuum.

What "vacuum refraction" means
- This is not a claim that "only light breaks". First and foremost, the vacuum conductivity itself changes with respect to causal processes.
- Light manifests this as ray bending and increased flight time.
- Particles manifest this as a change in trajectory.
- Clocks manifest this as a slowing of the local passage of time.
- Therefore, the word "refraction" in the hypothesis should be read broadly: it is not a single photon that is refracted, but the vacuum regime of causality transfer and process unfolding.

Positivity of time and compact objects
Physical thesis of the hypothesis: time is not an empty coordinate and therefore cannot physically vanish along a permissible trajectory. In the current reduced circuit, this is coded, in particular, by the condition:
1 ≤ n(theta_eff) < ∞.

Consequences:
1 + z = n_emit / n_obs < ∞,
Δt = ∫_Γ n(theta_eff) dl < ∞.

That is, a signal from the surface of a supercompact object can always exit, but at a very large n_emit, it will be monstrously redshifted and practically unobservable. The regularized form here is a working record of a more general requirement: the temporal regime remains physically permissible and positive.

Background module

This section refers only to the FRW reduction. It does not exhaust the hypothesis and should not be read as its complete canon.

Homogeneous Lagrangian and main quantities
Fields: theta(t) — the norm of the time field Theta; R(t) — the amplitude of the complex energy substrate Psi = R e^iφ. Effective homogeneous Lagrangian:
L_eff = 3/2 (thetȧ/theta)^2 + 1/2 Ṙ^2 − 1/2 V0 R^2 + 1/2 gamma^2 thetȧ^2.

From it follow the energy density and pressure (accurate to factor a^3):
ρ = (∂L/∂thetȧ) thetȧ + (∂L/∂Psi̇) Psi̇ − L,
p = L.

Useful identities: ρ − p = V0 R^2; ρ + p = gamma^2 thetȧ^2 + Ṙ^2 + 3 (thetȧ/theta)^2.
The expansion is given by
H^2 = (ρ_fields + ρ_m + ρ_r + ρ_Λ) / 3,
where ρ_fields is the Theta–Psi contribution, and background components are set by Omegas (matter, radiation, Λ).

Regularizing variable x
The equation of motion for theta contains a singular coordinate point at 1 − gamma²theta²/3 = 0. We introduce a safe variable
theta = (√3/gamma) tanh x,
thetȧ = (√3/gamma) sech²x · ẋ,

then it carries over to the limit |x| → ∞. This is a coordinate (numerical) regularization of the equations, not a separate physical postulate about the horizon. The equation in x (at an arbitrary H(t)):
ẍ + 3H ẋ = ẋ² ( 2 tanh x + coth x ).

Complex Psi = R e^iφ
We introduce the charge Q = a³ R² φ̇ (integral of motion). Then
R̈ + 3H Ṙ − Q²/(a⁶ R³) + V0 R = 0,
φ̇ = Q / (a³ R²),
ρ_phase = p_phase = 1/2 · Q²/(a⁶ R²).

In the limit Q → 0, a real prototype remains with a purely "spring" term V0 R.

Integration by N = ln a
For numerical stability, it is convenient to move to N. Let y = dx/dN, u = dR/dN. Then
ρ_fields = A(x,y,u) · H² + 1/2 V0 R², where
A(x,y,u) = 3/2 · (sech⁴x / tanh²x) · y² + 1/2 · u² + 3/2 · sech⁴x · y².
H² = [ ρ_m(a) + ρ_r(a) + ρ_Λ(a) + 1/2 V0 R² + 1/2 Q²/(a⁶ R²) ] / (3 − A).
Ḣ = −1/2 [ (gamma² thetȧ² + Ṙ² + 3 (thetȧ/theta)²) + ρ_m + 4/3 ρ_r + Q²/(a⁶ R²) ].

From the equations for ẍ and R̈, we obtain equations in N for y′(N), u′(N) and t′(N)=1/H (see background/frw_background.py).

Scalar perturbations (ADM)
Quadratic action for scalars in basic form:
L^(2) = 1/2 · dq̇ᵀ K dq̇ − 1/2 · dqᵀ [ G (k²/a²) + M ] dq,
dq = (d_theta, d_R, deltaφ).

Diagonal two-entity high-k model (consistent with the background):
Kthetatheta = 3/theta² + gamma²,
KRR = 1,
Kφφ = R²,

Gthetatheta = Kthetatheta,
GRR = 1,
Gφφ = R²,

MRR = V0 + 3 Q²/(a⁶ R⁴), Mthetatheta = Mφφ = 0.

Mixings (small ε) and diagonal corrections of the form + α_i (k²/a²) are allowed and explicitly given in the interface (theory/adm_symbolic.py).

Health conditions:
- No ghosts: K ≻ 0 (positive definite).
- No gradient instabilities: eigenvalues of K⁻¹ G are positive.
- Subluminality: c_s² ≤ 1 (taking into account numerical error).
- For finite k: ω²(k) — eigenvalues of K⁻¹ [ G (k²/a²) + M ] are positive.

Note: the code implements a numerically stable estimate of eigenvalues via Cholesky reduction (see src/linear_modes_adm.py). The symbolic framework for future full elimination of Lapse/Shift is in theory/adm_symbolic.py.

Lapse/Shift elimination (c = (d_N, B)):
L^(2) ⊃ − dqᵀ B c − 1/2 · cᵀ C c ⇒ c* = − C⁻¹ Bᵀ dq,
A_eff = A − B C⁻¹ Bᵀ, A = G (k²/a²) + M.

From here, G_eff and M_eff are recovered (terms at k²/a² and constants). In the numerical API, passing elim=... includes this step (see src/linear_modes_adm.py).

Equivalence principle: "why" it works here

In this frame, equivalence is not postulated but follows from a single "time index" of the vacuum n(theta_eff):
dtau = dt / n(theta_eff),
L ≈ (m/2n) v² − m/n ⇒ a = −∇ ln n(theta_eff),

the mass cancels out, and all bodies fall equally (m_grav = m_inert). Photons see the same n(theta_eff) — redshift, Shapiro delay, and lensing follow from the same function. See docs/01_interpretation.md for details.

Health checks and cutoffs
- Finite k: positivity of ω²(k) for all modes (if enabled).
- (Optional) heuristics for fast scans: EFT limit and PPN proxy; disabled by default, only conditions from the Lagrangian and equations of motion are considered primary.
- BBN: monotonic decrease of H(a) on a ∈ [1e−10, 1e−6] (radiation era).
- Superluminality: max eig(K⁻¹ G) ≤ 1 + tol (2×2 or ADM 3×3).

All checks are collected in a single aggregator checks/health.py.

Observables and fit
- Background: H(z), distances (SN, BAO). For example — fitting/ with demo CC-data and light Metropolis.
- For real H0, use omegas_from_km_s_Mpc conversion.
- Parameter scan (gamma, V0[, Q]) with a "hard stop" on health: scripts/run_grid_scan.py (--health, --adm, --finite_k, --lnN).

Future development
1) Full ADM-derivation of K, G, M with Lapse/Shift elimination; analytical proof of hyperbolicity.
2) Connecting real SN/BAO/CC and separate observational evaluation of the background without mixing it with the physical basis.
3) Expanding checks (more accurate PPN, EFT boundaries) and automation (CI-stop upon violations).
