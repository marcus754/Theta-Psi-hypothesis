# -*- coding: utf-8 -*-
"""
Summarize the current empirically allowed / preferred parameter envelope.
"""
from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path
from typing import Sequence


def parse_args(argv: Sequence[str] | None = None) -> argparse.Namespace:
    p = argparse.ArgumentParser(description="Write parameter envelope report")
    p.add_argument("--joint_csv", type=str, default="results/joint_region_scan.csv")
    p.add_argument("--profile_csv", type=str, default="results/joint_strongfield_profile_scan.csv")
    p.add_argument("--output_json", type=str, default="results/parameter_envelope.json")
    p.add_argument("--output_md", type=str, default="results/parameter_envelope.md")
    return p.parse_args(argv)


def main(argv: Sequence[str] | None = None) -> int:
    args = parse_args(argv)
    with Path(args.joint_csv).open("r", newline="", encoding="utf-8") as fh:
        joint = list(csv.DictReader(fh))
    with Path(args.profile_csv).open("r", newline="", encoding="utf-8") as fh:
        prof = list(csv.DictReader(fh))

    allowed = [r for r in joint if int(r.get("allowed_3sigma", "0")) == 1]
    gamma_vals = [float(r["gamma"]) for r in allowed]
    v0_vals = [float(r["V0"]) for r in allowed]
    envelope = {
        "joint_allowed_3sigma": {
            "count": len(allowed),
            "gamma_min": min(gamma_vals) if gamma_vals else None,
            "gamma_max": max(gamma_vals) if gamma_vals else None,
            "V0_min": min(v0_vals) if v0_vals else None,
            "V0_max": max(v0_vals) if v0_vals else None,
        },
        "profile_best_by_parameter": {},
    }
    for par_name in sorted({r["parameter"] for r in prof}):
        block = [r for r in prof if r["parameter"] == par_name]
        best = min(block, key=lambda r: float(r["chi2_total_best"]))
        envelope["profile_best_by_parameter"][par_name] = {
            "value": float(best["parameter_value"]),
            "gamma_best": float(best["gamma_best"]),
            "V0_best": float(best["V0_best"]),
            "chi2_total_best": float(best["chi2_total_best"]),
        }

    Path(args.output_json).write_text(json.dumps(envelope, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    lines = ["# Parameter Envelope\n\n"]
    lines.append("## Joint 3σ Envelope\n")
    lines.append(f"- Count: {envelope['joint_allowed_3sigma']['count']}\n")
    lines.append(f"- gamma: [{envelope['joint_allowed_3sigma']['gamma_min']}, {envelope['joint_allowed_3sigma']['gamma_max']}]\n")
    lines.append(f"- V0: [{envelope['joint_allowed_3sigma']['V0_min']}, {envelope['joint_allowed_3sigma']['V0_max']}]\n")
    lines.append("\n## Strong-Field Profile Best Values\n")
    for name, d in envelope["profile_best_by_parameter"].items():
        lines.append(f"- {name}: value={d['value']}, gamma={d['gamma_best']}, V0={d['V0_best']}, chi2={d['chi2_total_best']}\n")
    Path(args.output_md).write_text("".join(lines), encoding="utf-8")
    print(f"Wrote parameter envelope JSON to {args.output_json}")
    print(f"Wrote parameter envelope report to {args.output_md}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
