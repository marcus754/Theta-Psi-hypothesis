06. OBSERVABLES AND DATA COMPARISON

TIME CANON

- The fundamental object of time in the hypothesis is the Theta_munu tensor (see docs/11_derivation_path.md).
- The formulas in this section use the reduced FRW variable theta ≡ theta_eff, where theta_eff^2 = Theta_munu u^mu u^nu, u^mu u_mu = -1.
- All notations of the form n(theta) in this file should be read as n(theta_eff, Theta, Psi).
- This section describes the computational reduced circuit rather than the full tensor dynamics of Theta_munu.

CANONICAL FORMULATION (MANDATORY FOR ALL SECTIONS)

In the current reduced computational layer, the strong field response is built directly from the Theta_munu, Psi base:
J_refr(Theta,Psi) = Phi_eff(theta_eff) * I_grad(Theta,Psi),
L_refr = M_*^4 F(J_refr),
n = n(theta_eff, Theta, Psi).

For data in the repository, it is important to distinguish:
- EMPIRICAL observables: direct data-facing values used in falsifiability;
- DIAGNOSTIC observables: internal diagnostics and scans for development.

BRIEF OVERVIEW
- Background observables: H(z), distances (SN, BAO). Units: 8πG=1, c=1; convert H0 via omegas_from_km_s_Mpc.
- Growth of structures: D(a) from D'' + (2 + d ln H / d ln a) D' − (3/2) Ω_m mu D = 0; in the basic version mu=1 (tensor sector unmodified).
- Causality: c_T=1; scalars are luminal in high-k; superluminality is forbidden.
- Comparison: first, a "health" slice (K≻0, c_s²≤1, ω²>0, PPN/BBN), then fit CC/BAO/SN.
- Practice: construction of H(a) and D(a) from Omegas; reports on units and H0 conversion are mandatory.

Comparison goals:
- Background: H(z), BAO, SN — check expansion and effective energy.
- Local tests: PPN restrictions, gravitational waves (c_T=1), Solar System tests.
- Satellite synchronization: GPS clock correction as a reference check; the sign and magnitude of the correction coincide with the standard relativistic estimation.
  In the data-facing report, we use the published standard of about 38 mus/day for GPS.
- Laboratory and transport synchronization: Pound-Rebka, Hafele-Keating, and the twin paradox. All three are formatted as reference checks.
- Orbital weak-field dynamics: Lense-Thirring, binary pulsar timing, LLR.
- Atomic / chemical limit: Theta–Psi should not change standard relativistic chemistry; in this regime, n(theta_eff) -> 1, gamma_PPN = 1, beta_PPN = 1 is required.
- Cosmological perturbations: spectra, c_s², stability.

Reading status of this block:
- REFERENCE: GPS, Hafele-Keating, Pound-Rebka, Shapiro, solar redshift, solar limb deflection, Mercury, Lense-Thirring, binary pulsar timing, LLR;
- REFERENCE: atomic limit, twin-paradox, c_T=1;
- WEAK FIELD consistency only.

Strong field status:
- REFERENCE: public event-level report for M87*, Sgr A*, and GWTC-3 echo-scale;
- DIAGNOSTIC: conservative compact-bound scan from the EHT/LVK envelope;
- The compact-bound scan currently has no allowed points.

Practice:
- Use scripts/run_grid_scan.py for a preliminary slice of parameters (gamma, V0, Q) with automatic checks.
- Next — choice of MCMC/scanner for fitting H(z), BAO, SN, with fixed H0 (via omegas_from_km_s_Mpc).

DEMONSTRATION PIPELINES (IN CODE UNITS)

- CC (H(z)): fitting/model.py::Hz_dimless and fitting/data_io.py::load_cc_dimless.
- SN (mu(z)): fitting/model.py::mu_dimless and fitting/data_io.py::load_sn_dimless.
- BAO (D_V(z)): fitting/model.py::DV_dimless and fitting/data_io.py::load_bao_dimless.
- Combined log-likelihood and MCMC: python -m fitting.run_mcmc --steps 200 --output results/fit_demo.txt.
- For real H0, pass --H0_km_s_Mpc 70.0 (omegas_from_km_s_Mpc will be used internally), or manually --H0_dimless.
- To save the trace, specify --save_csv results/chain.csv (format: step, gamma, V0, Q, logpost, accept).
- Mercury perihelion check:
  - analytical PPN block: python scripts/run_mercury_precession_check.py;
  - numerical orbit integration: python scripts/run_mercury_numeric_check.py.
- Joint parameter region map:
  python scripts/run_joint_region_scan.py (CC+SN+BAO+compact_star + diagnostic compact scan).

HOW TO PREPARE REAL DATA (OFFLINE)

1) Download public tables (example):
   - CC (H(z)): cosmic chronometers compilation (z, H, σ_H) in km/s/Mpc.
   - SN: Pantheon/Pantheon+ (z, mu, σ_mu).
   - BAO: D_V(z) table (or convert to D_V) with (z, DV, σ_DV).
2) Convert formats to simple CSV (without network access):
   - python -m fitting.convert_public --cc CC_public.csv --H0_km_s_Mpc 70.0 --outdir fitting/data/real_dimless
   - python -m fitting.convert_public --sn Pantheon.csv --outdir fitting/data/real_dimless
   - python -m fitting.convert_public --bao BAO_DV.csv --outdir fitting/data/real_dimless
   - strong field empirical-level events (EHT ring catalog -> direct-model at the event level):
     python -m fitting.convert_public --forward_events strong_field_events.csv --outdir fitting/data/real_dimless
     where strong_field_events.csv has columns:
     name, observable, target, sigma, mass_msun, distance_m, source. In the empirical circuit, only observable=ring_uas is currently supported.
3) Run combined fit on converted data (replace paths --data/--sn/--bao).

ONE-STEP STRONG FIELD PIPELINE

If a direct run from the event table to reports is needed:
python scripts/run_forward_public_pipeline.py --forward_events_csv strong_field_events.csv --outdir results/public_forward

The following will be automatically collected:
- forward_event_catalog.json
- forward_event_scan.csv/.md
- forward_event_report.json/.md
- empirical_checks.md
- pipeline_summary.md

Note:
- echo-like values and conservative compact-bound scans are not part of the weak field empirical circuit; for strong field, they remain diagnostic artifacts.

Note. If BAO DV is given in Mpc, multiply by a dimensionless scale (e.g., H0/c) via --bao_scale_dimless during conversion.

LINEAR GROWTH OF STRUCTURES (QUASI-STATIONARY SUB-HORIZON LIMIT)

In the variable N=ln a, the growth equation for the factor D(a) takes the form
D'' + (2 + d ln H / d ln a) D' − (3/2) Ω_m(a) mu(a,k) D = 0,

where primes are derivatives with respect to N, and mu(a,k) is the gravity modification (in GR mu≡1).
For background cosmology with E(a)^2 = Ω_r a^-4 + Ω_m a^-3 + Ω_Λ, we obtain Ω_m(a) = Ω_m a^-3 / E^2(a),
d ln H / d ln a = (1/2) d ln E^2 / d ln a = (a / (2 E^2)) dE^2/da,
dE^2/da = −4 Ω_r a^-5 − 3 Ω_m a^-4.

In the current implementation, the base version is GR (tensor sector unmodified), so mu(a,k)=1.
Implementation: see src/growth.py, function growth_factor(a_min, a_max, omegas, nsteps, mu_of_a=None)
returns (a_grid, D_norm) with normalization D(a=1)=1.

RELATION TO P(k)

In the first approximation, the spectrum today can be assembled as P(k, z=0) ∝ D^2(a=1) * T^2(k) * P_prim(k),

where T(k) is the transfer function (can be borrowed from ΛCDM at early stages), and P_prim(k) ∝ k^n_s.
Differences from ΛCDM at the growth level are reflected in D(a) (and, if desired, in mu(a,k)).

Note: for a strict calculation of P(k) and CMB observables, a full Boltzmann code is needed; the current goal is to ensure a correct background and stable modes for initial phenomenological comparisons.
