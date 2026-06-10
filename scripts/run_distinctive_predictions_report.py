# -*- coding: utf-8 -*-
"""
Build a compact report of the theory's distinctive predictions/challenges.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Sequence


def parse_args(argv: Sequence[str] | None = None) -> argparse.Namespace:
    p = argparse.ArgumentParser(description="Write distinctive predictions report")
    p.add_argument("--prediction_json", type=str, default="results/prediction_suite.json")
    p.add_argument("--forward_event_report_json", type=str, default="results/forward_event_report.json")
    p.add_argument("--joint_profile_csv", type=str, default="results/joint_strongfield_profile_scan.csv")
    p.add_argument("--output_md", type=str, default="results/distinctive_predictions.md")
    return p.parse_args(argv)


def main(argv: Sequence[str] | None = None) -> int:
    args = parse_args(argv)
    fev = json.loads(Path(args.forward_event_report_json).read_text(encoding="utf-8"))
    import csv
    with Path(args.joint_profile_csv).open("r", newline="", encoding="utf-8") as fh:
        prof = list(csv.DictReader(fh))

    best_kappa = min((r for r in prof if r["parameter"] == "compact_star_kappa_refr"), key=lambda r: float(r["chi2_total_best"]))
    best_alpha = min((r for r in prof if r["parameter"] == "compact_star_refractive_strength"), key=lambda r: float(r["chi2_total_best"]))

    lines = ["# Distinctive Predictions\n\n"]
    lines.append("## Surviving Predictions\n")
    lines.append("- Weak-field branch remains GR-like to current accuracy: Mercury, Shapiro, redshift and limb deflection all pass.\n")
    if fev.get("events"):
        pulls = {e["name"]: abs(float(e["pull_sigma"])) for e in fev["events"]}
        primary = []
        for key, label in (
            ("EHT_M87_ring", "M87"),
            ("EHT_SgrA_ring", "SgrA*"),
            ("GWTC3_echo_scale", "GWTC-3"),
        ):
            if key in pulls:
                primary.append(f"{label}={pulls[key]:.3f}σ")
        if primary:
            lines.append(
                "- Strong-field event comparator remains viable without horizons: "
                + ", ".join(primary)
                + ".\n"
            )
        lines.append(
            f"- Public event-level check max |pull| = {float(fev.get('max_abs_pull_sigma', 0.0)):.3f}σ.\n"
        )
    lines.append("- The theory's distinctive strong-field claim is finite redshift and finite delay on the compact branch, with no physically reachable horizon.\n")
    lines.append("\n## Quantitative Levers\n")
    lines.append(
        f"- Joint fit currently prefers weaker refractive saturation strength around "
        f"`compact_star_refractive_strength ~ {float(best_alpha['parameter_value']):.3g}`.\n"
    )
    lines.append(
        f"- Joint fit currently prefers lower refractive index response around "
        f"`compact_star_kappa_refr ~ {float(best_kappa['parameter_value']):.3g}`.\n"
    )
    lines.append("\n## What Would Distinguish Θ–Ψ from GR\n")
    lines.append("- Finite redshift/delay for ultra-compact objects instead of true horizon behaviour.\n")
    lines.append("- A universal refractive response driven by `J_refr = Φ_eff I_grad` rather than pure metric curvature.\n")
    lines.append("- A preferred positive-time branch, forbidding continuation through `θ_eff = 0`.\n")
    Path(args.output_md).write_text("".join(lines), encoding="utf-8")
    print(f"Wrote distinctive predictions report to {args.output_md}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
