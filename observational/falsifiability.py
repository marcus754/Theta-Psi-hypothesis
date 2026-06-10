# -*- coding: utf-8 -*-
"""
Evaluate falsifiability criteria for Θ–Ψ from generated artifacts.
"""
from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path
from typing import Sequence


def parse_args(argv: Sequence[str] | None = None) -> argparse.Namespace:
    p = argparse.ArgumentParser(description="Run falsifiability checks")
    p.add_argument("--targets_json", type=str, default="fitting/data/prediction_targets.json")
    p.add_argument("--prediction_json", type=str, default="results/prediction_suite.json")
    p.add_argument("--joint_csv", type=str, default="results/joint_region_scan.csv")
    p.add_argument("--forward_event_report_json", type=str, default=None)
    p.add_argument("--output_json", type=str, default="results/falsifiability_check.json")
    p.add_argument("--output_md", type=str, default="results/falsifiability_check.md")
    return p.parse_args(argv)


def _load_csv(path: str) -> list[dict]:
    with Path(path).open("r", newline="", encoding="utf-8") as fh:
        return list(csv.DictReader(fh))


def main(argv: Sequence[str] | None = None) -> int:
    args = parse_args(argv)
    targets = json.loads(Path(args.targets_json).read_text(encoding="utf-8"))
    pred = json.loads(Path(args.prediction_json).read_text(encoding="utf-8"))
    joint = _load_csv(args.joint_csv)
    fev = None
    if args.forward_event_report_json:
        p = Path(args.forward_event_report_json)
        if p.exists():
            fev = json.loads(p.read_text(encoding="utf-8"))

    checks = []

    mt = targets["mercury_arcsec_per_century"]
    m_target = float(mt["target"])
    m_sigma = float(mt["sigma"])
    max_pull = float(mt["max_abs_pull_sigma"])
    m_ana = float(pred["mercury_analytic_arcsec_per_century"])
    m_num = float(pred["mercury_numeric_arcsec_per_century"])
    p_ana = (m_ana - m_target) / m_sigma
    p_num = (m_num - m_target) / m_sigma
    checks.append({"name": "mercury_analytic", "value": p_ana, "pass": abs(p_ana) <= max_pull})
    checks.append({"name": "mercury_numeric", "value": p_num, "pass": abs(p_num) <= max_pull})

    for key in ("solar_redshift_z", "solar_limb_deflection_arcsec", "shapiro_roundtrip_s"):
        v = float(pred[key])
        lo = float(targets[key]["min"])
        hi = float(targets[key]["max"])
        checks.append({"name": key, "value": v, "pass": lo <= v <= hi, "min": lo, "max": hi})

    jt = targets["joint_region"]
    n = len(joint)
    n3 = sum(1 for r in joint if int(r.get("allowed_3sigma", "0")) == 1)
    frac3 = n3 / max(n, 1)
    checks.append(
        {
            "name": "joint_allowed_3sigma_fraction",
            "value": frac3,
            "pass": (n3 >= int(jt["min_allowed_3sigma_points"])) and (frac3 <= float(jt["max_allowed_3sigma_fraction"])),
            "count": n3,
            "total": n,
        }
    )

    if fev is not None and "forward_events" in targets:
        ft = targets["forward_events"]
        event_pull = float(fev.get("max_abs_pull_sigma", 0.0))
        n_events = int(fev.get("n_events", 0))
        checks.append(
            {
                "name": "forward_events_max_abs_pull",
                "value": event_pull,
                "pass": (n_events >= int(ft.get("min_events", 1))) and (event_pull <= float(ft["max_abs_pull_sigma"])),
                "count": n_events,
                "total": n_events,
            }
        )

    overall = all(bool(c["pass"]) for c in checks)
    out = {
        "overall_pass": overall,
        "checks": checks,
        "inputs": {
            "targets_json": args.targets_json,
            "prediction_json": args.prediction_json,
            "joint_csv": args.joint_csv,
            "forward_event_report_json": args.forward_event_report_json,
        },
    }

    out_json = Path(args.output_json)
    out_json.parent.mkdir(parents=True, exist_ok=True)
    out_json.write_text(json.dumps(out, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")

    lines = ["# Falsifiability Check\n\n", f"- Overall: {'PASS' if overall else 'FAIL'}\n", "\n## Checks\n"]
    for c in checks:
        line = f"- {c['name']}: {'PASS' if c['pass'] else 'FAIL'} (value={c['value']})"
        if "min" in c and "max" in c:
            line += f", range=[{c['min']}, {c['max']}]"
        if "count" in c and "total" in c:
            line += f", count={c['count']}/{c['total']}"
        lines.append(line + "\n")

    out_md = Path(args.output_md)
    out_md.parent.mkdir(parents=True, exist_ok=True)
    out_md.write_text("".join(lines), encoding="utf-8")
    print(f"Wrote falsifiability JSON to {args.output_json}")
    print(f"Wrote falsifiability report to {args.output_md}")
    return 0 if overall else 2


if __name__ == "__main__":
    raise SystemExit(main())
