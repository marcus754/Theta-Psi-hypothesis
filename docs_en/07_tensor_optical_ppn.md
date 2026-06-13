07. TENSOR SECTOR, OPTICAL METRIC, AND PPN

TIME CANON

- The fundamental object of time in the hypothesis is the Theta_mu_nu tensor (see docs/11_derivation_path.md).
- The formulas in this section use the reduced FRW variable theta ≡ theta_eff, where theta_eff^2 = Theta_mu_nu u^mu u^nu, u^mu u_mu = -1.
- All notations of the form n(theta) in this file should be read as n(theta_eff, Theta, Psi).
- This section describes the computational reduced circuit rather than the full tensor dynamics of Theta_mu_nu.

CANONICAL FORMULATION (MANDATORY FOR ALL SECTIONS)

In the project, Theta_mu_nu and Psi are taken as the base; the reduced circuit does not require an external threshold and does not use manual branch cutting.

S = ∫ d^4x √(-g) [ (M_pl^2/2)R + L_m + L_ThetaPsi + M_*^4 F(J_refr) ],
J_refr(Theta,Psi) = Phi_eff(theta_eff) * I_grad(Theta,Psi).

Regime source (without Weyl threshold in the canon):
I_grad(Theta,Psi) = l_theta^2 * |∇_⊥theta|^2 + l_Psi^2 * |∇_⊥Psi|^2,
∇_⊥ = spatial part of the gradient (orthogonal to the local u^mu).

Meaning:
- In FRW (homogeneous background): I_grad=0, therefore J_refr=0.
- In strong field (large spatial gradients): I_grad>0, therefore the refractive sector is enabled directly via F(J_refr).
- There is no separate regime field in the canon.

The vacuum index is limited by construction:
n(theta_eff, Theta, Psi) = 1 + 𝒩(J_refr).

The invalid background branch is excluded by parameterization (rather than a filter):
A = 3 tanh^2( sqrt(A_raw/3) ) in [0,3),
H^2 = (ρ_bg + V/2 + ρ_phase)/(3-A).

In observational checks, it is not a "manual inclusion" that is evaluated, but a direct response to spatial inhomogeneities of Theta_mu_nu, Psi.
The section fixes the link: c_T = 1, universal optical metric n(theta_eff), and weak-field PPN limit gamma=beta=1.

---

1) SPEED OF TENSOR WAVES (c_T)

Key question: does the model change the speed of gravitational waves (c_T), strongly constrained by LIGO/Virgo observations (|c_T/c − 1| ≲ 10⁻¹⁵)?

In the current version of the hypothesis based on Theta_mu_nu, Psi, the gravitational sector is unmodified; the fields (theta, Psi) belong to "matter" (source for background and scalar perturbations). In the effective parameterization of modified gravity, this corresponds to α_T = 0, which strictly gives c_T = 1 for the tensor sector on a FRW background.

Practical conclusion:
- The properties of the tensor sector with c_T^2 = 1 are defined in the code (theory/tensor.py).
- A test (tests/test_tensor_speed.py) checks that c_T^2 = 1 is returned.
- If in the future operators changing the gravity kinetics are added, they must be explicitly included and immediately checked for compatibility with LIGO/Virgo constraints, i.e., α_T ≈ 0 in late epochs.

---

2) OPTICAL METRIC FROM n(theta_eff) AND EMERGENT TENSOR WAVES

Goal: to fix "from first principles" the link:
- one vacuum time index n(theta_eff) (refraction), principle of least time;
- a single optical metric for all matter;
- the emergence of the tensor sector (two TT polarizations, c_T = 1) as a universal low-energy response of the medium.

2.1) One n(theta_eff) and Fermat's Principle

Postulate: the local passage of proper time for matter and light is stretched by the same index
dtau = dt / n(theta_eff).

The physical formulation takes a bounded index:
1 ≤ n(theta_eff) < ∞.

Hence, in a strong field, we obtain not a "true horizon", but a horizonless limit with a very large but finite redshift.
Rays and trajectories minimize the optical length (Fermat's principle)

L_opt[Γ] = ∫_Γ n(theta_eff) dl,

 which already gives three classic weak-field effects from a single n(theta_eff): redshift, Shapiro delay, lensing (see 01_interpretation.md).

2.2) Optical Metric g_hat(n, u)

The variational principle is conveniently written covariantly via a Gordon-type "optical metric". Let
- u_mu ≡ ∂_mu theta / √(∂theta·∂theta) — a preferred timelike direction ("clocks");
- the base background metric in locally inertial coordinates — eta.

Then the effective metric, along whose geodesics light and matter travel, can be written as
g_hat_mu_nu(theta) = eta_mu_nu + (n(theta_eff)^2 − 1) u_mu u_nu

(equivalent notations are also found with a coefficient for u_mu u_nu depending on the chosen normalization; the fact of existence of a single g_hat encoding n(theta_eff) is significant). In a locally comoving frame u^mu = (1,0,0,0) and g_hat_00 = n^2, so dtau^2 = dt^2 / n^2.

Consequence: "vacuum refraction" reduces to geodesics of one and the same g_hat(theta) for all matter fields. This is the strict form of the principle of least time in terms of action.

For compact objects, the observability criterion:
1+z = n_emit/n_obs < ∞, Δt = ∫ n dl < ∞.

2.3) Emergent Tensors from a Single Propagation Law

Further, general EFT logic applies. If all matter and Psi interact via a single g_hat(theta):
S_total = S_ThetaPsi[theta, Psi] + sum_i S_matter[Phi_i; g_hat(theta)],

then the quantum effective functional at low energies inevitably contains geometric invariants of g_hat (Sakharov mechanism):
Γ_eff[g_hat] ⊃ ∫ d^4x √−g_hat [ Λ_ind^4 + (M_ind^2/2) R(g_hat) + O(∂^4) ].

Minimality argument: the only relevant (by derivatives) second-order tensor invariant is the curvature scalar R(g_hat). Its coefficient M_ind^2 is determined by integrating the fast modes of matter/Psi. This does not add a new ontology: g_hat remains a compact notation for the fact that all processes use the same clock-sector response.

From the linear variation of the EH term with respect to the TT part deltag_hat_ij^TT, the standard gravitational wave equation on an FRW background follows (in conformal time):
h''_ij + 2 H_conf h'_ij + k^2 h_ij = 0,

i.e., two transverse tensor polarizations with c_T = 1.

Result: "observable tensority" is the TT perturbations of the composite metric g_hat(theta) generated by the same n(theta_eff). Theta remains a scalar-clock; it is not Theta that is tensor, but the collective object g_hat.

WEAK-FIELD OBSERVABLE MATCHING (CONFORMAL FORM)

Let the conformal coefficient A(theta_eff) = n(theta_eff)^2 in statics and isotropy be representable as
A(theta_eff) = e^2Phi(x) = 1 + 2Phi + 2Phi^2 + O(Phi^3),

where Phi is the weak potential. Then the optical metric

g_hat_00 = A = 1 + 2Phi + 2Phi^2 + ...,
g_hat_ij = -A delta_ij = -(1 - 2Phi) delta_ij + O(Phi^2).

Comparing with the observational weak-field parameterization in isotropic gauge,
g_00 = 1 + 2Phi + 2beta Phi^2 + ...,
g_ij = −(1 − 2gamma Phi) delta_ij + ...,

we get the required weak-field coefficients gamma = 1, beta = 1. In Theta–Psi, this is read not as an ontological link to geometry, but as a sanity condition: the same n(theta_eff) must give the correct coefficients for redshift, delay, lensing, and orbital dynamics.

2.4) Compatibility with the Theta_mu_nu, Psi Base

- Homogeneous limit: the minisuperspace Lagrangian (see background/frw_symbolic.py) describes the zero modes of Theta and Psi amplitudes on an FRW background; meanwhile, the leading tensor dynamics come from the induced EH term for g_hat and give c_T=1 (see docs/07_tensor_optical_ppn.md and tests).
- Scalar perturbations: the current K, G, M matrices of the minimal two-entity reduced circuit are the sector of scalar phonons of the same medium. Hyperbolicity and subluminality conditions (K≻0, eig(K⁻¹G)≤1, ω²>0) are already automated in checks/.

2.5) What to check and where the risks are

- PPN/dipole: derive gamma_PPN, beta_PPN and the absence of dipole radiation as a consequence of universal n and single g_hat.
- c_T and polarizations: c_T=1 (already fixed), absence of dominant scalar/vector impurities in GW.
- Superluminality/ghosts: controlled by conditions on K, G, M (in our case — checks/health.py).
- Λ_ind and UV cutoff: phenomenological control (EFT regime |thetȧ/theta|, |gamma thetȧ| < Λ).

2.6) Roadmap to "the Stone"

1) Formally link the current n(theta_eff) with an explicit g_hat(theta) and write out the TT equation on FRW (short calculation, c_T=1).
2) Bring the PPN module to formulas (gamma_PPN, beta_PPN), add tests.
3) Develop the ADM branch: scalar sector on the background of g_hat(theta), checks of ω²(k)>0 and prohibition of superluminality during Lapse/Shift elimination.

This note fixes the foundation: "one n(theta_eff) → one g_hat(theta) → tensors as emergent TT modes with c_T=1". All classic effects and observed tensor gravity follow from the same object — the vacuum time index.

---

3) PPN: SANITY WEAK-FIELD TEST gamma_PPN = beta_PPN = 1

Goal: to show that with a universal coupling of matter via the same time index n(theta_eff), weak-field observables do not break: PPN parameters in a weak, static field take values gamma=1, beta=1, and dipole radiation is absent. This is a check of compatibility with measurements, not an ontological basis for the hypothesis.

3.1) Universal Effective Metric

In the weak field, statics, and isotropy, it is convenient to take a universal conformal form (consistent with the equivalence principle):
g_munu(x) = A(theta_eff(x)) eta_munu, A(theta_eff) = n(theta_eff)^2,

where eta_mu_nu is Minkowski, and A(theta_eff)>0 encodes the same n(theta_eff) for all matter fields. This form:
- fixes proper time stretching dtau^2 = dt^2/n^2;
- gives matching coefficients for g_00 and g_ij perturbations at the first order, which leads to gamma=1.
- in the full (not just PN) formulation, a bounded index 1 ≤ n(theta_eff) < ∞ is used, which excludes a true horizon and leaves only a finite extreme redshift for supercompact objects.

Note. In Theta–Psi, the primary entity is not the metric, but the single time index n(theta_eff). The metric notation is needed only for a compact accounting of the fact that light, clocks, and massive matter use the same local propagation law.

3.2) Linearization and Comparison with PPN

We write A(theta_eff) = 1 + 2Phi(x) + O(Phi^2), where Phi is the weak potential (proportional to ln n). Then

g_00 = 1 + 2Phi + O(Phi^2),
g_ij = -(1 - 2Phi)delta_ij + O(Phi^2).

Observational PPN notation in isotropic gauge:
g_00 = 1 + 2Phi + 2betaPhi^2 + O(Phi^3),
g_ij = -(1 - 2gammaPhi)delta_ij + O(Phi^2).

Comparing the first order, we get gamma = 1.
At the second order, in the absence of additional non-local/non-linear sources in A(theta_eff), beta = 1 holds (the conformal form scales g_00 and g_ij equally and does not introduce special quadratic corrections).

Result: gamma_PPN = 1, beta_PPN = 1.

3.3) Absence of Dipole Radiation

The scalar "charging" of matter is identical (universal coupling), so there is no separation of bodies by "scalar charge" — the dipole moment of the system does not radiate at the leading order. Radiative losses begin with the quadrupole channel, which is compatible with pulsar tests.

3.4) Link with n(theta_eff) and Classic Tests

From section 01_interpretation.md: classic weak-field effects follow from a single n(theta_eff) via Deltanu/nu approx Deltaln n, t = int n(theta_eff)dl, alpha approx int nabla_perp ln n dl. Choosing the normalization ln n approx -2Phi aligns the numerical coefficients of redshift, Shapiro delay, and lensing with observables and is compatible with gamma=1.

3.5) Limitations and Conditions of Applicability

- Statics/weak field/isotropy: |Phi| ll 1, strong anisotropic theta gradients are absent.
- Universality of coupling: all matter uses the same n(theta_eff) and the same effective propagation record g(theta); otherwise, the equivalence principle is violated and PPN parameters deviate from observed values.
- High-order conformal/disformal corrections are possible outside the weak field — they must preserve local tests (parameter constraints).
- The horizonless strong-field regime is formulated as n(theta_eff) being bounded from above; then 1+z and ∫ n dl are finite.

Thus, with a universal n(theta_eff) and a single effective propagation record g = n(theta_eff)^2 eta, weak-field checks pass: gamma=beta=1, scalar dipole is absent.
