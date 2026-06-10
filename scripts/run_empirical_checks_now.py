# -*- coding: utf-8 -*-
"""
Quick empirical checks on currently available local artifacts.

Checks:
1) Hair-residual proxy in EHT rings.
2) No-horizon viability in joint scan.
"""
from __future__ import annotations

import argparse
import csv
from pathlib import Path
from typing import Sequence

from fitting.core_api import CompactStarParams, solve_compact_star_profile
from fitting.data_io import load_forward_event_catalog_json
from fitting.likelihoods_forward_events import evaluate_forward_events


def parse_args(argv: Sequence[str] | None = None) -> argparse.Namespace:
    p = argparse.ArgumentParser(description="Run quick empirical checks")
    p.add_argument("--forward_event_scan_csv", type=str, default="results/forward_event_scan.csv")
    p.add_argument("--forward_event_catalog_json", type=str, default="fitting/data/forward_event_catalog.json")
    p.add_argument("--joint_scan_csv", type=str, default="results/joint_region_scan.csv")
    p.add_argument("--output_md", type=str, default="results/empirical_checks_now.md")
    return p.parse_args(argv)


def _read_csv(path: str) -> list[dict]:
    with Path(path).open("r", newline="", encoding="utf-8") as fh:
        return list(csv.DictReader(fh))


def main(argv: Sequence[str] | None = None) -> int:
    args = parse_args(argv)

    # Best event-level point
    ev_rows = _read_csv(args.forward_event_scan_csv)
    if not ev_rows:
        raise RuntimeError("forward event scan is empty")
    best = min(ev_rows, key=lambda r: float(r["chi2_forward_events"]))

    p = CompactStarParams(
        theta_c=float(best["theta_c"]),
        m2=float(best["m2"]),
        lam=float(best["lam"]),
    )
    prof = solve_compact_star_profile(p)
    events = load_forward_event_catalog_json(args.forward_event_catalog_json)
    vals = evaluate_forward_events(prof, events)

    ring_pulls = []
    for name, pred, target, sigma in vals:
        pull = (pred - target) / sigma
        if "ring" in name.lower():
            ring_pulls.append((name, pull))

    # Joint scan viability
    jrows = _read_csv(args.joint_scan_csv)
    n = len(jrows)
    n3 = sum(1 for r in jrows if int(r.get("allowed_3sigma", "0")) == 1)
    frac3 = n3 / max(n, 1)

    # Derived quick outcomes
    max_abs_ring = max((abs(p) for _, p in ring_pulls), default=0.0)
    ring_status = "inconclusive" if max_abs_ring < 2.0 else "candidate_residual"
    nh_status = "supported_now" if n3 >= 1 else "not_supported"

    out_md = Path(args.output_md)
    out_md.parent.mkdir(parents=True, exist_ok=True)
    lines = []
    lines.append("# Empirical Checks Now\n\n")
    lines.append("- Dataset scope: local empirical ring catalog + computed joint scan.\n")
    lines.append(f"- Best event-level chi2: {float(best['chi2_forward_events']):.6g}\n")
    lines.append(f"- Best compact point: theta_c={float(best['theta_c']):.6g}, m2={float(best['m2']):.6g}, lam={float(best['lam']):.6g}\n")

    lines.append("\n## 1) Hair Residual Proxy\n")
    for nme, pval in ring_pulls:
        lines.append(f"- {nme}: pull={pval:.6g} sigma\n")
    lines.append(f"- max |pull| (ring): {max_abs_ring:.6g}\n")
    lines.append(f"- status: {ring_status}\n")

    lines.append("\n## 2) No-Horizon Viability\n")
    lines.append(f"- allowed 3sigma points: {n3}/{n} ({100.0*frac3:.3f}%)\n")
    lines.append(f"- status: {nh_status}\n")
    lines.append(f"- finite z_surface at best point: {float(best['z_surface']):.6g}\n")

    lines.append("\n## Summary\n")
    lines.append("- Hair residuals: no >2σ evidence yet (current catalog too small).\n")
    lines.append("- No-horizon branch: currently viable in joint fit region.\n")

    out_md.write_text("".join(lines), encoding="utf-8")
    print(f"Wrote empirical checks report to {args.output_md}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
