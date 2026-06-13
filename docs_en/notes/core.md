CORE (Θ–Ψ)

Unified working document: physical core + operational plan + observational benchmarks.
meaning.md — author's manifesto and is not edited.

1. PHYSICAL CORE

- The Theta_munu time tensor field and the Psi energy substrate are accepted as the two fundamental entities in the hypothesis.

- Both Theta_munu and Psi are primary but functionally asymmetric: Theta_munu defines the structure of causality and dynamics, while in the reduced layer, the amplitude R, which is part of the complex Psi = R e^{iφ}, is dynamized.

- In the full version of the canon, Psi is treated as a complex order parameter Psi = R e^{iφ}; in the FRW reduction, the amplitude R and a possible phase charge Q = a^3 R^2 \dot φ remain. In this sense, |Psi|^2 is naturally read as the amplitude density rather than as a separate new entity.

- The scalar theta in the calculations is the reduced projection theta_eff (rather than an independent ontology).

- Operational definition of reduction: theta_eff^2 = Theta_munu u^mu u^nu, u^mu u_mu = -1.

- Canonical definition of Phi_eff: Phi_eff^can = ln(theta_eff/theta_0).

- The smooth Phi_eff used in some solvers is considered only as numerical regularization of the positive branch and does not introduce new physics.

- Canonical regime: FRW/weak-field reproduces verified local observational tests; the Theta–Psi difference in the refractive sector is enabled through spatial inhomogeneities, i.e., via I_grad. Weak-field reference checks: GPS, Hafele-Keating, Pound-Rebka, Shapiro, solar redshift, solar deflection, Mercury, Lense-Thirring, binary pulsar timing, LLR, atomic limit, twin-paradox, c_T=1.

- Strong-field reference check: public event-level report on M87*, Sgr A* and GWTC-3 echo-scale. Conservative compact-bound scan remains diagnostic and currently has no allowed points.

- On a homogeneous background I_grad=0, so FRW by itself does not include the refractive sector.

- Time in the hypothesis is understood not as a scale but as a positive physical object; the allowed dynamics should not transition the time regime into zero or the negative branch.

- The working branch of the reduced circuit implements this principle, in particular, via the finite index 1 <= n(theta_eff) < ∞.

- Weak-field matching is fixed as a sanity condition on observables: n ≈ 1 + 2 Phi_eff in the region of weak inhomogeneities.

- Strong-field is enabled directly via the composite invariant of the regime J_refr = Phi_eff(theta_eff) I_grad(Theta, Psi).

- The minimal stationary strong-field sector is defined by the functional L_stat = 1/2 I_grad - U(theta, Psi) + αF(J_refr) and leads to a quasi-linear matrix system for theta'' and Psi''.

- A unified covariant completion covering both FRW and strong-field:
S_cov = S_red + S_m[g_hat, matter]

- g and g_hat in this completion are not primary entities of the hypothesis. They are used as a technical covariant notation for the unified law of process propagation in the Theta–Psi medium. The ontological basis remains solely Theta_munu and Psi.

- The horizon is not introduced as a physical object of Theta–Psi. In the strong regime, the canonical time variable remains on the positive branch, and the observed effect is read as a large but finite slowing down of processes: finite n, redshift, and delay.

- Derived microphysics layer, derived from the Theta_munu, Psi base, sets the local action S_micro with response-scalar R_micro, derived as a minimal local quadratic deviation scalar from locality / positivity / weak-field normalization, and with the preferred-slicing measure measure_dt_d3x fixed as derived volume normalization. Matching rule:
n^2 = Z_s / Z_t

- The working bridge is read as S_micro -> S_red -> S_cov; S_red remains the reduced/covariant completion, while S_micro is the derived layer from the Theta_munu, Psi basis.

- Algebraic reduction is fixed by R_micro -> 0 and the canonical sector weights, not by free coefficient functions or free potentials.

- Reduction-check: A_Theta, B_Theta, A_Psi, B_Psi, Z_t, Z_s -> 1, n -> 1, local quadratic blocks are reconstructed, heavy modes are integrated out.

- In code the derived micro layer is pinned to the canonical ledgers: S_red_from_micro = S_red, S_cov_from_micro = S_cov.

- Ordered derivation ladder is fixed as S_micro -> weak_field_projection -> integrate_heavy_modes -> assemble_S_red -> add_S_matter_completion -> S_cov.

- Micro perturbation check lives in theory/micro_perturbation.py and requires positive kinetic/response coefficients, positive sound-speed squares, and n² = Z_s / Z_t > 0.

- Single strong-field consequence of the derived micro layer is the compact branch summary with finite redshift, finite delay, and ok_no_horizon == True. The explicit derived strong-field ledger is micro_strongfield_primary_ledger() and micro_to_strongfield_primary_bridge().

- Stationary strong-field branch is closed by strongfield_branch_closure_ledger() and test_strongfield_branch_closure_pins_horizonless_criteria. Closure means: positive theta_eff, finite n, finite redshift/delay, outgoing characteristic, finite energy, regular center, elliptic radial operator, forward-observable map, and explicit falsification criteria.

- Dynamic collapse is closed as a canonical evolution contract by dynamic_collapse_closure_ledger() and test_dynamic_collapse_closure_contract_pins_evolution_criteria. Closure means finite-energy initial data, evolution of Theta–Psi, tracked min theta_eff and max n, preserved outgoing characteristic, no blow-up before branch decision, hyperbolic evolution operator, and only two allowed outcomes: settle to stationary branch or disperse to weak-field. A numerical PDE solver is an implementation of this contract, not a new canon.

- Symbol audit classifies the micro layer as: derived = projection/matching/reduction outputs, including measure_dt_d3x and R_micro; assumed = R, chi; diagnostic = the compact-branch consequence.

- Closure of points 1-3 is fixed as follows:

S_micro
  -> universal Z_t(R_micro), Z_s(R_micro)
  -> n^2 = Z_s/Z_t
  -> g_hat as technical propagation record
  -> weak-field observable coefficients

- Point 1 (g/g_hat): g and g_hat do not add a third entity if all matter-sectors receive them only as a record of the same local propagation law derived from Z_t, Z_s.

- Point 2 (universality): in the current derived micro layer, universality is ensured by the fact that one common set of response coefficients Z_t(R_micro), Z_s(R_micro) is introduced for the matter-sector rather than separate coefficients for different substances or fields.

- Point 3 (weak-field): the weak-field test is considered closed at the canonical level if Z_t -> 1, Z_s -> 1, n -> 1 when R_micro -> 0, and the first non-zero correction to n is normalized as n ≈ 1 + 2 Phi_eff and thereby gives the correct observed coefficients for redshift, delay, lensing, and orbital dynamics.

- Closure status: points 1-3 are closed in the canon via propagation_law_closure_ledger() and the test test_propagation_law_closure_pins_points_1_2_3. This ledger prohibits an independent g_hat ontology, fixes one common Z_t, Z_s for matter-sectors, and establishes the weak-field normalization before fit.

- Vacuum-offset sanity: the absolute additive constant in the effective matter action is not considered a direct Theta–Psi source. The source map depends on response/boundary differences (R_micro, Z_t, Z_s, n, boundary data) rather than an arbitrary Lambda_vac_offset. Casimir is preserved as a boundary spectral difference. Code fixation: vacuum_offset_source_ledger() and test_vacuum_offset_source_ledger_*.

- EPR/CHSH sanity: the algebroid layer reads Theta as a causal anchor/foliation and Psi_AB as an indecomposable singlet-state. The first test fixes without fit-parameters: E(a,b) = -a·b, |S_CHSH| = 2√2, no-signalling marginals are equal to 1/2, and the correlation generator lies in ker ρ_Theta and is not a signal. The layer also explicitly rejects PR-box/super-Tsirelson correlations, fixes GHZ/Mermin contextuality, foliation-order independence, and the no-fit angle profile. Code fixation: theory/epr_algebroid.py, tests/test_epr_algebroid.py, docs/13_epr_algebroid.md.

- Particle sanity: a particle is defined as a derived object rather than a third ontology: particle = stable localized Theta-projection of a Psi-sector. Minimal checks: finite localized footprint, conserved labels (mass, charge, spin), detector event map, distinction between product and entangled multiparticle sectors, classical trajectory limit. Code fixation: theory/particle.py, tests/test_particle.py, docs/14_particle.md.

- Anti-fit guardrails: the weak-field n(theta_eff) profile is closed as a canonical asinh with a fixed slope of 2; linear/exp/tanh2 alternatives remain diagnostic-only. Saturating family F_v1(J) with removable J_* normalization, ADM eps*/alpha* mixings, sigma-extension, bridge calibration coefficients, and inverse-BAO required multipliers are also considered diagnostic-only until derived from the action and fixed before fit. In the canon now: canonical_free_profile_count = 0, canonical_bridge_fit_parameters = 0. Code fixation: theory/anti_fit.py, tests/test_anti_fit.py.

- Action-derivation status: from the current minimal action, the diagonal ADM high-k block K = G = diag(3/theta² + gamma², 1, R²) and M = diag(0, V0 + 3Q²/(a⁶R⁴), 0) already follow. Therefore, canonical epsK*, epsG*, alpha_k2* are zero; non-zero values require new operators in the action. The F_v1(J) form and strong-field bridge scales are not yet derived: the action sets the argument J_refr and requirements for F but does not select the F(J) function itself, and the bridge requires solving sourced static Euler-Lagrange equations. Code fixation: theory/action_derivation.py, tests/test_action_derivation.py.

- Refractive semigroup principle: if we set S(J) = 1 - F(J) and require S(J1 + J2) = S(J1)S(J2), then with continuity we obtain S(J) = exp(-J/J_*), and thus F(J) = 1 - exp(-J/J_*). J_* here is not physics but the normalization of the remaining response. This is no longer a choice of form but a law for the remaining response. Code fixation: theory/action_derivation.py, tests/test_action_derivation.py.

- Base homogeneous Lagrangian in reduced amplitude form:
L_eff = 1.5 * (thetȧ/theta)^2 + 0.5 * Ṙ^2 - 0.5 * V0 * R^2 + 0.5 * gamma^2 * thetȧ^2,
where theta = theta_eff — FRW reduction of the Theta_munu tensor, and R — reduced amplitude representative of the complex Psi = R e^{iφ}.

- Energy and pressure:
rho = Σ_i (∂L/∂q̇_i) q̇_i - L
p = L

- Mathematical criterion "without true horizon" (at finite energy):
L_eff ⊃ (3/2) (thetȧ/theta)^2 => rho ⊃ (3/2) (thetȧ/theta)^2.

If along the physical trajectory rho(t) <= rho_max < ∞, then
|d ln theta / dt| = |thetȧ/theta| <= sqrt(2rho_max/3) = K.

Integrating:
ln theta(t) >= ln theta(t0) - K|t-t0|,
theta(t) >= theta(t0) exp(-K|t-t0|) > 0 for any finite t.

Consequence: the state theta=0 is not reached in finite time at finite energy. This is the substantive form of the thesis about the positivity of time as an object. Within the n=n(theta) interpretation, this gives a dynamic barrier to the formation of a true horizon: an asymptotic regime of strong slowing down is obtained rather than a physically reachable "no return" surface.

- SSOT according to the formulas of the reduced circuit: background/frw_symbolic.py.

- Canonical physical layer (tensor time): docs/11_derivation_path.md.

2. NUMERICAL SCHEME

- background/frw_integrator.py — adaptive RK45 in the variable x = artanh(gamma theta / sqrt(3)).

- background/frw_background.py — FRW circuit with scale factor and Omega-components; the invalid branch is konstrukтивно removed via smooth A saturation.

- scripts/run_grid_scan.py — base scan + health-filters.

3. OBSERVATIONAL BENCHMARKS FOR gamma

- Working frame by scale:
  - laboratory: sensitivity roughly down to gamma ~ 10^-6;
  - astrophysics: gamma ~ 10^-10 ... 10^-8;
  - cosmology/early Universe: gamma <= 10^-18 (key limit);
  - supercompact objects: check finite n(theta) and finite redshift/delay.

- Strict calibration requires real tables and explicit unit normalization in the fit-pipeline.
