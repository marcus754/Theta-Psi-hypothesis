# -*- coding: utf-8 -*-
"""
Run full validation pipeline for Θ–Ψ and collect a single summary report.
"""
from __future__ import annotations

import argparse
import subprocess
import sys
from pathlib import Path
from typing import Sequence


EMPIRICAL_STEPS = {
    "prediction_suite",
    "forward_event_scan",
    "forward_event_report",
    "joint_region_scan",
    "falsifiability_check",
}

DIAGNOSTIC_STEPS = {
    "forward_observables_scan",
    "likelihood_comparison",
    "joint_strongfield_profile_scan",
    "inverse_bao_fit",
    "parameter_envelope",
    "distinctive_predictions",
    "head2head_comparison",
    "killer_tests",
}


def parse_args(argv: Sequence[str] | None = None) -> argparse.Namespace:
    p = argparse.ArgumentParser(description="Run full validation suite")
    p.add_argument("--output_md", type=str, default="results/validation_suite.md")
    p.add_argument("--forward_event_catalog_json", type=str, default="fitting/data/forward_event_catalog.json")
    return p.parse_args(argv)


def _run(cmd: list[str]) -> tuple[int, str]:
    proc = subprocess.run(cmd, capture_output=True, text=True)
    text = (proc.stdout or "") + (proc.stderr or "")
    return proc.returncode, text.strip()


def build_steps(py: str, forward_event_catalog_json: str) -> list[tuple[str, list[str]]]:
    return [
        ("prediction_suite", [py, "-m", "scripts.run_prediction_suite", "--output_json", "results/prediction_suite.json", "--output_md", "results/prediction_suite.md"]),
        ("forward_observables_scan", [py, "-m", "scripts.run_forward_observables_scan", "--n_theta_c", "16", "--n_m2", "16", "--output_csv", "results/forward_observables_scan.csv", "--output_md", "results/forward_observables_scan.md"]),
        ("forward_event_scan", [py, "-m", "scripts.run_forward_event_scan", "--event_catalog_json", forward_event_catalog_json, "--n_theta_c", "16", "--n_m2", "16", "--output_csv", "results/forward_event_scan.csv", "--output_md", "results/forward_event_scan.md"]),
        ("forward_event_report", [py, "-m", "scripts.run_forward_event_report", "--scan_csv", "results/forward_event_scan.csv", "--event_catalog_json", forward_event_catalog_json, "--output_json", "results/forward_event_report.json", "--output_md", "results/forward_event_report.md"]),
        ("likelihood_comparison", [py, "-m", "scripts.run_likelihood_comparison", "--steps", "180", "--compact_sanity_weight", "1.0", "--output_csv", "results/likelihood_comparison.csv", "--output_md", "results/likelihood_comparison.md"]),
        ("head2head_comparison", [py, "-m", "scripts.run_head2head_comparison", "--theta_samples", "120", "--output_md", "results/head2head_comparison.md", "--output_json", "results/head2head_comparison.json"]),
        ("inverse_bao_fit", [py, "-m", "scripts.run_inverse_bao_fit", "--head2head_json", "results/head2head_comparison.json", "--output_csv", "results/inverse_bao_fit.csv", "--output_json", "results/inverse_bao_fit.json", "--output_md", "results/inverse_bao_fit.md"]),
        ("joint_strongfield_profile_scan", [py, "-m", "scripts.run_joint_strongfield_profile_scan", "--ng", "6", "--nv", "6", "--n_profile", "4", "--output_csv", "results/joint_strongfield_profile_scan.csv", "--output_md", "results/joint_strongfield_profile_scan.md"]),
        ("joint_region_scan", [py, "-m", "scripts.run_joint_region_scan", "--ng", "25", "--nv", "25", "--output_csv", "results/joint_region_scan.csv", "--output_md", "results/joint_region_scan.md"]),
        ("parameter_envelope", [py, "-m", "scripts.run_parameter_envelope_report", "--joint_csv", "results/joint_region_scan.csv", "--profile_csv", "results/joint_strongfield_profile_scan.csv", "--output_json", "results/parameter_envelope.json", "--output_md", "results/parameter_envelope.md"]),
        ("distinctive_predictions", [py, "-m", "scripts.run_distinctive_predictions_report", "--prediction_json", "results/prediction_suite.json", "--forward_event_report_json", "results/forward_event_report.json", "--joint_profile_csv", "results/joint_strongfield_profile_scan.csv", "--output_md", "results/distinctive_predictions.md"]),
        ("falsifiability_check", [py, "-m", "scripts.run_falsifiability_check", "--targets_json", "fitting/data/prediction_targets.json", "--prediction_json", "results/prediction_suite.json", "--joint_csv", "results/joint_region_scan.csv", "--forward_event_report_json", "results/forward_event_report.json", "--output_json", "results/falsifiability_check.json", "--output_md", "results/falsifiability_check.md"]),
        ("killer_tests", [py, "-m", "scripts.run_killer_tests_report", "--falsifiability_json", "results/falsifiability_check.json", "--head2head_json", "results/head2head_comparison.json", "--forward_event_report_json", "results/forward_event_report.json", "--output_md", "results/killer_tests.md"]),
    ]


def main(argv: Sequence[str] | None = None) -> int:
    args = parse_args(argv)
    py = sys.executable
    out_md = Path(args.output_md)
    out_tmp = out_md.with_suffix(out_md.suffix + ".tmp")
    steps = build_steps(py, args.forward_event_catalog_json)

    if out_md.exists():
        out_md.unlink()
    if out_tmp.exists():
        out_tmp.unlink()

    results: list[dict] = []
    total = len(steps)
    for idx, (name, cmd) in enumerate(steps, start=1):
        print(f"[validation] Step {idx}/{total} START {name}", flush=True)
        rc, out = _run(cmd)
        results.append({"name": name, "rc": rc, "output": out})
        status = "PASS" if rc == 0 else "FAIL"
        tail = out.splitlines()[-1] if out else ""
        print(f"[validation] Step {idx}/{total} {status} {name} (rc={rc})", flush=True)
        if tail:
            print(f"[validation]   {tail}", flush=True)

    out_md.parent.mkdir(parents=True, exist_ok=True)
    lines = ["# Validation Suite\n\n"]
    overall = all(r["rc"] == 0 for r in results)
    lines.append(f"- Overall command status: {'PASS' if overall else 'FAIL'}\n")
    lines.append("\n## Classes\n")
    lines.append("- `empirical`: direct data-facing checks used in falsifiability.\n")
    lines.append("- `diagnostic`: proxy / scan / calibration reports kept for triage only.\n")
    lines.append("\n## Steps\n")
    for r in results:
        if r["name"] in EMPIRICAL_STEPS:
            cls = "empirical"
        elif r["name"] in DIAGNOSTIC_STEPS:
            cls = "diagnostic"
        else:
            cls = "unclassified"
        lines.append(f"- {r['name']} [{cls}]: {'PASS' if r['rc']==0 else 'FAIL'} (rc={r['rc']})\n")
    lines.append("\n## Notes\n")
    for r in results:
        snippet = r["output"].splitlines()[-1] if r["output"] else ""
        lines.append(f"- {r['name']}: {snippet}\n")
    out_tmp.write_text("".join(lines), encoding="utf-8")
    out_tmp.replace(out_md)

    print(f"Wrote validation suite report to {args.output_md}")
    return 0 if overall else 1


if __name__ == "__main__":
    raise SystemExit(main())
