# -*- coding: utf-8 -*-
"""
Scan compact constraints over (gamma, V0) grid and report exclusion power.
"""
from __future__ import annotations

import argparse
import csv
from pathlib import Path
from typing import Sequence

from fitting.core_api import Omegas
from fitting.data_io import load_compact_targets_json
from fitting.likelihoods_compact import CompactTargets, compact_metrics_model, loglike_compact


def parse_args(argv: Sequence[str] | None = None) -> argparse.Namespace:
    p = argparse.ArgumentParser(description="Compact constraints scan")
    p.add_argument("--compact_targets_json", type=str, default="fitting/data/compact_constraints_observational_proxy.json")
    p.add_argument("--gmin", type=float, default=0.1)
    p.add_argument("--gmax", type=float, default=10.0)
    p.add_argument("--ng", type=int, default=40)
    p.add_argument("--vmin", type=float, default=0.01)
    p.add_argument("--vmax", type=float, default=1.0)
    p.add_argument("--nv", type=int, default=40)
    p.add_argument("--theta_scale", type=float, default=1.0)
    p.add_argument("--output_csv", type=str, default="results/compact_constraint_scan.csv")
    p.add_argument("--output_md", type=str, default="results/compact_constraint_scan.md")
    return p.parse_args(argv)


def _linspace(a: float, b: float, n: int) -> list[float]:
    if n <= 1:
        return [a]
    h = (b - a) / (n - 1)
    return [a + i * h for i in range(n)]


def main(argv: Sequence[str] | None = None) -> int:
    args = parse_args(argv)
    om = Omegas(H0=1.0)
    cj = load_compact_targets_json(args.compact_targets_json)
    target = CompactTargets(
        mode=str(cj.get("mode", "gaussian")),
        z_surface=cj.get("z_surface"),
        sigma_z_surface=cj.get("sigma_z_surface"),
        delay_proxy=cj.get("delay_proxy"),
        sigma_delay_proxy=cj.get("sigma_delay_proxy"),
        n_peak=cj.get("n_peak"),
        sigma_n_peak=cj.get("sigma_n_peak"),
        z_surface_min=cj.get("z_surface_min"),
        sigma_z_surface_min=cj.get("sigma_z_surface_min"),
        delay_proxy_min=cj.get("delay_proxy_min"),
        sigma_delay_proxy_min=cj.get("sigma_delay_proxy_min"),
        n_peak_max=cj.get("n_peak_max"),
        sigma_n_peak_max=cj.get("sigma_n_peak_max"),
    )

    rows = []
    for g in _linspace(args.gmin, args.gmax, args.ng):
        for v in _linspace(args.vmin, args.vmax, args.nv):
            m = compact_metrics_model(
                gamma=g,
                V0=v,
                Q=0.0,
                omegas=om,
                theta_scale=args.theta_scale,
            )
            ll = loglike_compact(m, target)
            chi2 = -2.0 * float(ll)
            rows.append(
                {
                    "gamma": g,
                    "V0": v,
                    "z_surface": float(m["z_surface"]),
                    "delay_proxy": float(m["delay_proxy"]),
                    "n_peak": float(m["n_peak"]),
                    "chi2_compact": chi2,
                    "allowed": int(chi2 <= 1e-12),
                }
            )

    allowed = [r for r in rows if r["allowed"] == 1]
    excluded = [r for r in rows if r["allowed"] == 0]
    best = min(rows, key=lambda r: r["chi2_compact"])
    worst = max(rows, key=lambda r: r["chi2_compact"])

    out_csv = Path(args.output_csv)
    out_csv.parent.mkdir(parents=True, exist_ok=True)
    with out_csv.open("w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()))
        w.writeheader()
        w.writerows(rows)

    out_md = Path(args.output_md)
    lines = []
    lines.append("# Compact Constraints Scan\n\n")
    lines.append(f"- Grid: gamma in [{args.gmin}, {args.gmax}] x{args.ng}, V0 in [{args.vmin}, {args.vmax}] x{args.nv}\n")
    lines.append(f"- Total points: {len(rows)}\n")
    lines.append(f"- Allowed (chi2=0): {len(allowed)} ({100.0*len(allowed)/len(rows):.1f}%)\n")
    lines.append(f"- Excluded (chi2>0): {len(excluded)} ({100.0*len(excluded)/len(rows):.1f}%)\n")
    lines.append("\n## Extremes\n")
    lines.append(f"- Best: gamma={best['gamma']:.6g}, V0={best['V0']:.6g}, chi2={best['chi2_compact']:.6g}\n")
    lines.append(f"- Worst: gamma={worst['gamma']:.6g}, V0={worst['V0']:.6g}, chi2={worst['chi2_compact']:.6g}\n")
    out_md.write_text("".join(lines), encoding="utf-8")

    print(f"Wrote compact constraint scan CSV to {args.output_csv}")
    print(f"Wrote compact constraint scan report to {args.output_md}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
