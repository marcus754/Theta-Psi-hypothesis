10. FALSIFIABLE CRITERIA

TIME CANON

- The fundamental object of time in the hypothesis is the Theta_munu tensor (see docs/11_derivation_path.md).
- The formulas in this section use the reduced FRW variable theta ≡ theta_eff, where theta_eff^2 = Theta_munu u^mu u^nu, u^mu u_mu = -1.
- All notations of the form n(theta) in this file should be read as n(theta_eff, Theta, Psi).
- This section describes the computational reduced circuit rather than the full tensor dynamics of Theta_munu.

CANONICAL FORMULATION (MANDATORY FOR ALL SECTIONS)

In the current canon, only two entities are primary: Theta_munu and Psi. The strong field response is set directly by the composite invariant of the regime

J_refr(Theta,Psi) = Phi_eff(theta_eff) * I_grad(Theta,Psi),
I_grad(Theta,Psi) = l_theta^2 * |∇_⊥theta|^2 + l_Psi^2 * |∇_⊥Psi|^2.

The refractive sector is defined via

L_refr = M_*^4 F(J_refr),
n = n(theta_eff, Theta, Psi),
g_hat_munu = g_munu + (n^2 - 1) u_mu u_nu.

Key separation for this file:
- EMPIRICAL checks rely only on direct data-facing observables and are used in scripts/run_falsifiability_check.py;
- DIAGNOSTIC checks may use internal diagnostics, scans, and calibration artifacts but are not considered a strict falsification criterion.

0) WORKING ORDER

In the current version, the following are already fixed:
- a single reduced canon;
- covariant form of the action;
- link between background and strong-field regimes;
- scalar perturbation sector;
- base set of stability checks.

Remaining unfinished are:
- full dynamics of the tensor time field;
- general non-stationary strong-field regime;
- strict derivation of all quantum effects from the base action.

1) WEAK-FIELD TEST PACKAGE

Criteria are checked via scripts/run_prediction_suite.py and scripts/run_falsifiability_check.py.

- Mercury (anomalous precession):
  - target value 42.98"/century,
  - allowed deviation: |pull| <= 3σ.
- Redshift at the Sun's surface:
  - range 1.8e-6 .. 2.4e-6.
- Shapiro (roundtrip, near conjunction):
  - range 1.0e-4 .. 4.0e-4 seconds.
- Light deflection at the Sun's limb:
  - range 1.70 .. 1.80 arcseconds.

2) JOINT REGION CRITERION

According to scripts/run_joint_region_scan.py:

- a non-zero 3σ region must exist (N_3sigma >= 1),
- meanwhile, the 3σ region must not occupy almost the entire grid (f_3sigma <= 0.20) — otherwise the model is non-falsifiable on this set of tests.

3) DIAGNOSTIC STRONG FIELD SCANS

Scans like scripts/run_compact_constraint_scan.py and inferred transitions from strong_field_bounds_registry.json remain in the repository only as archival diagnostics. They may use inferred transitions of type observation -> inferred indicator, so they are not included in the PASS/FAIL for falsifiability.

4) MACHINE-READABLE SPECIFICATION

Threshold values for empirical falsifiability are moved to:
- fitting/data/prediction_targets.json.

Automatic check:
- python scripts/run_falsifiability_check.py.

For the current version of the canon, two classes of direct distinguishing tests are prioritized:
1. strong field saturation tests on the primary compact-channel (z_surface, delay, ok_no_horizon), using ring_uas and echo_delay_ms as downstream forward checks;
2. background compatibility tests on CC+SN+BAO, but only as a second layer after closing points 1-3 from the section above.

If the result is FAIL, the version is considered unfinished and requires refinement on the empirical circuit. Diagnostic reports may fall or diverge without an automatic falsification status.
