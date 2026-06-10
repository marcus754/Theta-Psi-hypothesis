# -*- coding: utf-8 -*-
"""
Scan empirical event-level forward likelihood over compact-star grid.
"""
from __future__ import annotations

import argparse
import csv
from pathlib import Path
from typing import Sequence

from fitting.core_api import CompactStarParams, solve_compact_star_profile
from fitting.data_io import load_forward_event_catalog_json
from fitting.likelihoods_forward_events import loglike_forward_events


def _linspace(a: float, b: float, n: int) -> list[float]:
    if n <= 1:
        return [a]
    h = (b - a) / (n - 1)
    return [a + i * h for i in range(n)]


def parse_args(argv: Sequence[str] | None = None) -> argparse.Namespace:
    p = argparse.ArgumentParser(description="Forward event-level scan")
    p.add_argument("--theta_c_min", type=float, default=0.1)
    p.add_argument("--theta_c_max", type=float, default=3.0)
    p.add_argument("--n_theta_c", type=int, default=20)
    p.add_argument("--m2_min", type=float, default=0.05)
    p.add_argument("--m2_max", type=float, default=2.0)
    p.add_argument("--n_m2", type=int, default=20)
    p.add_argument("--lam", type=float, default=0.001)
    p.add_argument("--event_catalog_json", type=str, default="fitting/data/forward_event_catalog.json")
    p.add_argument("--r_max", type=float, default=20.0)
    p.add_argument("--dr", type=float, default=0.02)
    p.add_argument("--theta_scale", type=float, default=1.0)
    p.add_argument("--psi_grad_weight", type=float, default=1.0)
    p.add_argument("--theta_psi_coupling", type=float, default=1.0)
    p.add_argument("--refractive_strength", type=float, default=1.0)
    p.add_argument("--refractive_j_star", type=float, default=1.0)
    p.add_argument("--kappa_refr", type=float, default=1.0)
    p.add_argument("--output_csv", type=str, default="results/forward_event_scan.csv")
    p.add_argument("--output_md", type=str, default="results/forward_event_scan.md")
    return p.parse_args(argv)


def main(argv: Sequence[str] | None = None) -> int:
    args = parse_args(argv)
    events = load_forward_event_catalog_json(args.event_catalog_json)

    rows = []
    for tc in _linspace(args.theta_c_min, args.theta_c_max, args.n_theta_c):
        for m2 in _linspace(args.m2_min, args.m2_max, args.n_m2):
            p = CompactStarParams(
                theta_c=tc,
                m2=m2,
                lam=args.lam,
                r_max=args.r_max,
                dr=args.dr,
                theta_scale=args.theta_scale,
                psi_grad_weight=args.psi_grad_weight,
                theta_psi_coupling=args.theta_psi_coupling,
                refractive_strength=args.refractive_strength,
                refractive_j_star=args.refractive_j_star,
                kappa_refr=args.kappa_refr,
            )
            prof = solve_compact_star_profile(p)
            ll = loglike_forward_events(prof, events)
            rows.append(
                {
                    "theta_c": tc,
                    "m2": m2,
                    "lam": args.lam,
                    "psi_grad_weight": args.psi_grad_weight,
                    "theta_psi_coupling": args.theta_psi_coupling,
                    "refractive_strength": args.refractive_strength,
                    "refractive_j_star": args.refractive_j_star,
                    "kappa_refr": args.kappa_refr,
                    "z_surface": float(prof["z_surface"]),
                    "delay_proxy": float(prof["delay_proxy"]),
                    "chi2_forward_events": -2.0 * float(ll),
                }
            )

    best = min(rows, key=lambda r: float(r["chi2_forward_events"]))

    out_csv = Path(args.output_csv)
    out_csv.parent.mkdir(parents=True, exist_ok=True)
    with out_csv.open("w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()))
        w.writeheader()
        w.writerows(rows)

    lines = []
    lines.append("# Forward Event Scan\n\n")
    lines.append("- Components: compact_star -> explicit forward ring transfer -> empirical event-level chi2\n")
    lines.append(f"- Grid: theta_c x m2 = {args.n_theta_c} x {args.n_m2}, lam={args.lam}\n")
    lines.append(
        f"- Strong-field closure: psi_grad_weight={args.psi_grad_weight}, "
        f"theta_psi_coupling={args.theta_psi_coupling}, "
        f"refractive_strength={args.refractive_strength}, "
        f"refractive_j_star={args.refractive_j_star}, "
        f"kappa_refr={args.kappa_refr}\n"
    )
    lines.append(f"- Catalog: {args.event_catalog_json}\n")
    lines.append("\n## Best point\n")
    for k in (
        "theta_c",
        "m2",
        "lam",
        "psi_grad_weight",
        "theta_psi_coupling",
        "refractive_strength",
        "refractive_j_star",
        "kappa_refr",
        "z_surface",
        "delay_proxy",
        "chi2_forward_events",
    ):
        lines.append(f"- {k}: {best[k]:.6g}\n")

    out_md = Path(args.output_md)
    out_md.parent.mkdir(parents=True, exist_ok=True)
    out_md.write_text("".join(lines), encoding="utf-8")
    print(f"Wrote forward-event scan CSV to {args.output_csv}")
    print(f"Wrote forward-event scan report to {args.output_md}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
