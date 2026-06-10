# -*- coding: utf-8 -*-
"""
Build empirical event-by-event residual report at best point from forward_event_scan CSV.
"""
from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path
from typing import Sequence

from fitting.core_api import CompactStarParams, solve_compact_star_profile
from fitting.data_io import load_forward_event_catalog_json
from fitting.likelihoods_forward_events import evaluate_forward_events


def parse_args(argv: Sequence[str] | None = None) -> argparse.Namespace:
    p = argparse.ArgumentParser(description="Forward event-by-event report")
    p.add_argument("--scan_csv", type=str, default="results/forward_event_scan_v3.csv")
    p.add_argument("--event_catalog_json", type=str, default="fitting/data/forward_event_catalog.json")
    p.add_argument("--r_max", type=float, default=20.0)
    p.add_argument("--dr", type=float, default=0.02)
    p.add_argument("--theta_scale", type=float, default=1.0)
    p.add_argument("--output_json", type=str, default="results/forward_event_report_v3.json")
    p.add_argument("--output_md", type=str, default="results/forward_event_report_v3.md")
    return p.parse_args(argv)


def main(argv: Sequence[str] | None = None) -> int:
    args = parse_args(argv)
    rows = list(csv.DictReader(Path(args.scan_csv).open("r", newline="", encoding="utf-8")))
    if not rows:
        raise RuntimeError("scan CSV is empty")
    best = min(rows, key=lambda r: float(r["chi2_forward_events"]))

    p = CompactStarParams(
        theta_c=float(best["theta_c"]),
        m2=float(best["m2"]),
        lam=float(best["lam"]),
        r_max=args.r_max,
        dr=args.dr,
        theta_scale=args.theta_scale,
    )
    prof = solve_compact_star_profile(p)
    events = load_forward_event_catalog_json(args.event_catalog_json)
    vals = evaluate_forward_events(prof, events)

    evt = []
    max_abs_pull = 0.0
    for name, pred, target, sigma in vals:
        pull = (pred - target) / sigma
        max_abs_pull = max(max_abs_pull, abs(pull))
        evt.append({
            "name": name,
            "predicted": pred,
            "target": target,
            "sigma": sigma,
            "pull_sigma": pull,
        })

    out = {
        "best_point": {
            "theta_c": float(best["theta_c"]),
            "m2": float(best["m2"]),
            "lam": float(best["lam"]),
            "chi2_forward_events": float(best["chi2_forward_events"]),
            "z_surface": float(best["z_surface"]),
            "delay_proxy": float(best["delay_proxy"]),
        },
        "events": evt,
        "max_abs_pull_sigma": max_abs_pull,
        "n_events": len(evt),
    }

    out_json = Path(args.output_json)
    out_json.parent.mkdir(parents=True, exist_ok=True)
    out_json.write_text(json.dumps(out, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")

    lines = []
    lines.append("# Forward Event Report\n\n")
    lines.append(f"- scan_csv: {args.scan_csv}\n")
    lines.append(f"- catalog: {args.event_catalog_json}\n")
    lines.append(f"- best chi2_forward_events: {out['best_point']['chi2_forward_events']:.6g}\n")
    lines.append(f"- max |pull|: {max_abs_pull:.6g}\n")
    lines.append("\n## Events\n")
    for e in evt:
        lines.append(
            f"- {e['name']}: pred={e['predicted']:.6g}, target={e['target']:.6g}, "
            f"sigma={e['sigma']:.6g}, pull={e['pull_sigma']:.6g}\n"
        )

    out_md = Path(args.output_md)
    out_md.parent.mkdir(parents=True, exist_ok=True)
    out_md.write_text("".join(lines), encoding="utf-8")
    print(f"Wrote forward-event report JSON to {args.output_json}")
    print(f"Wrote forward-event report to {args.output_md}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
