09. SUPERCOMPACT OBJECTS WITHOUT EVENT HORIZONS

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
- In FRW (homogeneous background): I_grad=0, therefore J_refr=0 and the refractive sector is disabled without an external switch.
- In strong field (large spatial gradients): I_grad>0, therefore the refractive sector is enabled directly via F(J_refr).
- There is no separate regime field σ in the canon; only Theta_munu and Psi remain primary.

The vacuum index in the solver implementation is written via the regularized response:

n(theta_eff, Theta, Psi) = 1 + 𝒩(J_refr),
𝒩(0)=0; 𝒩'(J)>0; 𝒩(J) < ∞ on the physical branch.

In the solver implementation, a numerical integrator limiter may be used, but it is not part of the physical formula.

The invalid background branch is excluded by parameterization (rather than a filter):
A = 3 tanh^2( sqrt(A_raw/3) ) in [0,3),
H^2 = (ρ_bg + V/2 + ρ_phase)/(3-A).

In comparison with GR, it is not a "manual inclusion" that is checked, but the dynamic response F(J_refr) arising from spatial inhomogeneities of Theta–Psi.

Goal: to describe the strong-field regime such that it follows from the equations themselves: the model does not produce a true event horizon (as in a standard black hole), but a supercompact object with a very large but finite redshift. This section describes the steady strong-field regime; the question of object formation is considered separately.

1) INITIAL ASSUMPTION

In the optical notation
dtau = dt / n(theta_eff),
physically allowed states satisfy
1 ≤ n(theta_eff) < ∞.

Canonical finite form (working in the code):
n(theta_eff, Theta, Psi) = 1 + 𝒩(J_refr).

This means:
- an infinite redshift does not arise on the physical branch;
- all local clocks slow down in a finite manner.

2) WHY THERE IS NO EVENT HORIZON IN THIS MODEL

For an effective metric of type g_hat with a positive and finite n:
- the coefficient of the time part does not vanish;
- the light characteristics do not lose the outgoing branch due to n alone.

Result: an event horizon does not form in this branch. Light from the surface can escape outside, but at large n_surface, it will be extremely redshifted.

Key point: the absence of a horizon here is obtained not from an external physical ceiling, but because the physical branch does not reach n → ∞.
n_guard is needed only for calculation stability.

2.0) STATUS OF THE STATIONARY BRANCH CLOSURE

The stationary compact branch is considered closed in the canon only as a set of verifiable criteria, not as a verbal thesis. Code source:
theory.micro_strongfield.strongfield_branch_closure_registry()

Criteria:
theta_eff > 0,
n finite,
z_surface finite,
delay finite,
outgoing characteristic exists,
finite energy,
regular center,
elliptic radial operator.

This registry also fixes the direct map:
n_profile -> z_surface -> delay -> ray/image map -> echo transfer.

And explicit branch rejection conditions:
theta_eff -> 0,
n -> infinity,
energy -> infinity,
loss of ellipticity,
loss of outgoing characteristic,
incompatibility with direct observables.

Important: this closes the stationary strong field branch. Dynamical collapse is closed by a separate evolution contract:
theory.micro_strongfield.dynamic_collapse_closure_registry()

Its criteria:
finite-energy initial data,
theta_eff(t,r) > 0,
n(t,r) finite,
E_total(t) finite,
outgoing characteristic exists,
hyperbolic evolution operator,
no blow-up before branch decision.

Allowed outcomes:
settle to stationary branch,
disperse to weak field.

Forbidden outcomes:
theta_eff -> 0,
n -> infinity,
energy blow-up,
loss of hyperbolicity,
loss of outgoing characteristic.

The numerical PDE solver is considered an implementation of this contract, not a separate canonical assumption.

2.1) DYNAMICAL CRITERION OF PRE-LIMIT REGIME UNREACHABILITY

The homogeneous sector contains the term
L_eff ⊃ (3/2) (thetȧ/theta)^2,
meaning in the energy density
ρ ⊃ (3/2) (thetȧ/theta)^2.

For any trajectory with finite energy ρ(t) <= ρ_max < ∞:
|d ln theta / dt| = |thetȧ/theta| <= sqrt(2ρ_max/3) = K.

Integrating with respect to time:
ln theta(t) >= ln theta(t0) - K|t-t0|,
theta(t) >= theta(t0) exp(-K|t-t0|) > 0 for any finite t.

Consequently, theta=0 is not reached in finite time with finite energy. If the "true horizon" requires a limit state theta->0 (or equivalently an unreachable limit of causal transmission), then such a horizon is not realized as a physically reachable surface. An asymptotic regime remains in the model: very large redshift and large delay, but without a strict boundary of no return.

3) OBSERVABLES IN STRONG FIELD (APPROXIMATE INDICATORS)

For the surface of a supercompact object:
1 + z_surface = n_surface / n_obs < ∞,
Δt = ∫_Γ n(theta_eff) dl < ∞.

The following indicators are used in the code:
- n_surface, n_peak;
- z_surface = n_surface / n_obs - 1 (in the demo n_obs=1);
- delay = sum (n_i - n_obs) along the discretized trajectory.

3.1) WHAT IS CONSIDERED PRIMARY HERE

For the theory, there is one primary strong-field effect: finite redshift in the absence of a true horizon. In practice, this is read as the combination
z_surface + delay + ok_no_horizon.

The meaning of this combination:
- z_surface grows with compactness but remains finite;
- delay grows along with the same regime;
- ok_no_horizon=True fixes that the branch does not transition into a normal black hole structure.

ring_uas and echo_delay_ms remain downstream checks, not a separate physical entity.

MINIMAL STATIONARY SPHERICAL SECTOR
For the first step in the strong field, we fix the radial profile theta(r):
theta''(r) + (2/r) theta'(r) + ∂U/∂theta = 0,
U(theta) = 1/2 m² theta² + 1/4 λ theta⁴, m²>0, λ≥0.

Boundary conditions:
theta(0) = theta_c, theta'(0) = 0, theta(r→∞) → 0.

After solving the ODE:
n(r) = n(theta(r)),
z_surface = n(0) - 1,
Δt = ∫_0^r_max (n(r)-1) dr.

4) COMPUTATIONAL IMPLEMENTATION

- Bounded index profile: theory/optical_metric.py (refractive_index_from_theta, n_bounds_ok).
- No-horizon check: checks/no_horizon.py (no_horizon_from_n_values, no_horizon_from_theta).
- General cutoff block: checks/health.py (check_n_bounds, check_no_horizon).
- Scan and report: scripts/run_grid_scan.py, checks/compact_report.py.
- Stationary ODE solver: checks/compact_star.py.
- Stationary sector scan: scripts/run_compact_star_scan.py.
- Exclusion map (target-based): scripts/run_compact_star_exclusion.py.
- compact-likelihood supports two modes:
  - mode="gaussian" (deprecated: symmetric point fitting);
  - mode="constraints" (primary: soft one-sided constraints).
- Constraint informativeness scan mode="constraints":
  scripts/run_compact_constraint_scan.py (fraction of allowed area and exclusion map).

4.1) STATUS OF STRONG FIELD CHECKS

- REFERENCE: public forward event report for M87*, Sgr A* and GWTC-3 echo-scale.
- DIAGNOSTIC ARCHIVE: conservative compact-bound scan from the EHT/LVK envelope.
- Public event-level comparator currently passes within the chosen tolerance.
- Conservative compact-bound scan currently has no allowed points and remains archival.

5) WHAT FALSIFIES THE BRANCH

In this formulation, the branch is rejected if simultaneously required:
- observationally, z/delays incompatible with a finite physical branch;
- or dynamics require exceeding the finite n(theta_eff) profile;
- or the full likelihood (CC+SN+BAO+compact) does not allow an area with ok_no_horizon=True;
- or the constraint scan shows an empty allowed region at physical thresholds for z_surface/delays.

This turns the "horizonless black star" into a testable hypothesis rather than an interpretational comment.

---

APPENDIX: REGISTRY OF STRONG FIELD OBSERVATIONS AND THEIR TRANSLATION INTO INDICATORS

REGISTRY OF STRONG FIELD OBSERVATIONS AND THEIR TRANSLATION INTO INDICATORS

Goal: to fix how observational strong field results (EHT/LVK) are matched with the working indicators z_surface, delay, n_peak in the diagnostic Theta–Psi circuit.

1) SOURCES
- EHT M87* (2019-04-10):
  - article: EHT Collaboration, ApJL 875 L1.
  - URL: https://iopscience.iop.org/article/10.3847/2041-8213/ab0ec7
  - used: ring/shadow size and consistency of strong-field geometry.
- EHT Sgr A* (2022-05-12):
  - article: EHT Collaboration, ApJL 930 L12/L15.
  - URL: https://iopscience.iop.org/article/10.3847/2041-8213/ac6674
  - used: ring size and consistency of morphology with a compact object.
- LVK GWTC-3 tests of GR (2022-04-28, LIGO-P2100275):
  - URL: https://dcc.ligo.org/LIGO-P2100275/public
  - used: absence of statistically significant echo-like post-merger signals.

2) TRANSLATION RULE (CURRENT VERSION)
Since there are no direct publication intervals specifically for z_surface, delay, n_peak, a conservative translation is applied:
- EHT data -> lower bounds on z_surface;
- LVK data (no-echo) -> lower bound on delay;
- EHT+LVK consistency envelope -> soft upper limit for n_peak.

This is a diagnostic matching, not a direct observable-to-observable fit. Consequently, these constraints should not be used as an automatic PASS/FAIL falsifiability criterion; their place is in diagnostic scans and triage.

3) CURRENT AGGREGATED CONSTRAINTS
Aggregation is performed by the script scripts/run_build_compact_constraints_from_bounds.py according to the rule:
- for lower-bounds, the maximum value is taken,
- for upper-bounds, the minimum value is taken.

Current output file:
- fitting/data/compact_constraints_observational_proxy.json (archival diagnostic output).

4) WHAT TO IMPROVE NEXT
- replace the inferred transformation with an explicit forward-model model -> image/ring/GW-observable;
- add VLBI/GRAVITY sources and new LVK catalogs;
- add a machine-readable mapping-id and translation formula to each entry.
