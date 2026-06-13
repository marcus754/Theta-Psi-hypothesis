STATUS (V2.4)

The file answers only the question: which artifacts are to be considered the source of truth.

Legend:
- Canonical — can be used in derivations.
- Provisional — working layer, not canon.
- Exploratory — research/legacy layer.

| Status | Path | Purpose |
|---|---|---|
| Canonical | docs/notes/meaning.md | Author's manifesto of the hypothesis (do not edit automatically) |
| Canonical | docs/notes/core.md | Unified summary: core + plan + observational benchmarks |
| Canonical | docs/11_derivation_path.md | Canonical physical layer: tensor time Theta_munu |
| Canonical | background/frw_symbolic.py | SSOT via formulas L/EOM/rho/p |
| Canonical | docs/notes/frw_reduced_eom.md | Readable form of the background equations |
| Canonical | background/frw_integrator.py | Adaptive integration in variable x |
| Canonical | background/frw_background.py | FRW integrator with Omega-components |
| Canonical | scripts/run_grid_scan.py | Base parametric scan |
| Canonical | theory/adm_symbolic.py | ADM-scalar sector of the reduced theta_eff circuit |
| Canonical | theory/covariant_action.py | Symbolic ledger of the full covariant Theta–Psi action |
| Canonical | theory/perturbation_proof.py | Full covariant completion and perturbation ledger |
| Canonical | theory/micro_action.py | Derived microphysics ledger and local matching rule |
| Canonical | theory/micro_strongfield.py | Single strong-field consequence of the derived micro layer |
| Canonical | theory/micro_perturbation.py | Derived microphysics perturbation ledger |
| Provisional | theory/twofield_stability.py | Minimal FRW/stationary stability diagnostics |
| Exploratory | theory/tensor_symbolic.py | Symbolic check of the reduced formulas |
| Exploratory | theory/sigma_gate.py | Legacy gate for the old sigma-layer |
| Canonical | theory/optical_metric.py | finite n(theta_eff), strong-field branch |
| Provisional | docs/00_overview.md ... docs/10_falsifiability.md | Reduced computational theta_eff layer (Theta_munu canon defined in docs/11_derivation_path.md) |
| Canonical | scripts/run_validation_suite.py | End-to-end validation |
| Canonical | scripts/run_falsifiability_check.py | Version PASS/FAIL criteria |
| Provisional | checks/ | Veto/health checks (part phenomenological) |
| Provisional | fitting/forward_transfer.py | explicit/event-level forward |
| Provisional | fitting/strong_field_link.py | background->compact link-model |
| Provisional | fitting/data/strong_field_bounds_registry.json | strong-field registry (some entries exploratory) |
| Provisional | docs/13_publication_draft.md | External paper-level draft for publication |
| Exploratory | — | no dedicated exploratory modules currently |

USAGE RULES

1. Take formulas from background/frw_symbolic.py, not from paraphrases.
2. If Canonical and Provisional diverge, Canonical takes priority.
3. Keep the docs/notes documents in a minimal set: meaning.md, core.md, frw_reduced_eom.md, STATUS.md.
4. Time ontology: Theta_munu and Psi are primary entities; theta in scripts and FRW docs is treated as theta_eff.
