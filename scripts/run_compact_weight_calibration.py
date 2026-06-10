# -*- coding: utf-8 -*-
"""
Calibrate compact sanity weight so it does not dominate strong-field fit.
"""
from __future__ import annotations

import argparse
import csv
from pathlib import Path
from typing import Sequence

from fitting.core_api import Omegas
from fitting.data_io import (
    load_cc_dimless,
    load_sn_dimless,
    load_bao_dimless,
    load_compact_targets_json,
    load_compact_star_targets_json,
)
from fitting.likelihoods_compact import CompactTargets
from fitting.likelihoods_compact_star import CompactStarTargets
from fitting.run_mcmc import Params, loglikes_components


def parse_args(argv: Sequence[str] | None = None) -> argparse.Namespace:
    p = argparse.ArgumentParser(description="Calibrate compact sanity weight")
    p.add_argument("--gamma", type=float, default=0.4855549418483553)
    p.add_argument("--v0", type=float, default=0.06648278115869659)
    p.add_argument("--q", type=float, default=0.0)
    p.add_argument("--compact_targets_json", type=str, default="fitting/data/compact_constraints_observational_proxy.json")
    p.add_argument("--compact_star_targets_json", type=str, default="fitting/data/compact_star_observational_proxy.json")
    p.add_argument("--output_csv", type=str, default="results/compact_weight_calibration.csv")
    p.add_argument("--output_md", type=str, default="results/compact_weight_calibration.md")
    return p.parse_args(argv)


def main(argv: Sequence[str] | None = None) -> int:
    args = parse_args(argv)
    om = Omegas(H0=1.0)
    cc = load_cc_dimless("fitting/data/cc_demo_dimless.csv")
    sn = load_sn_dimless("fitting/data/sn_demo_dimless.csv")
    bao = load_bao_dimless("fitting/data/bao_demo_dimless.csv")
    cj = load_compact_targets_json(args.compact_targets_json)
    compact = CompactTargets(
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
    csj = load_compact_star_targets_json(args.compact_star_targets_json)
    compact_star = CompactStarTargets(
        z_surface=csj["z_surface"],
        sigma_z_surface=csj["sigma_z_surface"],
        delay_proxy=csj["delay_proxy"],
        sigma_delay_proxy=csj["sigma_delay_proxy"],
    )
    p = Params(gamma=args.gamma, V0=args.v0, Q=args.q)

    weight_grid = [1.0, 0.1, 0.05, 0.01, 0.005, 0.001, 0.0005, 0.0001]
    rows = []
    for w in weight_grid:
        comps = loglikes_components(
            p,
            om=om,
            cc=cc,
            sn=sn,
            bao=bao,
            compact=compact,
            compact_weight=w,
            compact_star=compact_star,
            compact_star_link_model=True,
            compact_star_link_background=True,
        )
        chi2_compact = -2.0 * float(comps.get("compact", 0.0))
        chi2_cstar = -2.0 * float(comps.get("compact_star", 0.0))
        rows.append(
            {
                "compact_weight": w,
                "chi2_compact": chi2_compact,
                "chi2_compact_star": chi2_cstar,
                "ratio_compact_to_star": chi2_compact / max(chi2_cstar, 1e-12),
            }
        )

    # Criterion: compact sanity <= 5x compact_star contribution
    recommended = None
    for r in rows:
        if r["ratio_compact_to_star"] <= 5.0:
            recommended = r
            break
    if recommended is None:
        recommended = rows[-1]

    out_csv = Path(args.output_csv)
    out_csv.parent.mkdir(parents=True, exist_ok=True)
    with out_csv.open("w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()))
        w.writeheader()
        w.writerows(rows)

    out_md = Path(args.output_md)
    lines = []
    lines.append("# Compact Weight Calibration\n\n")
    lines.append(f"- Reference params: gamma={args.gamma:.6g}, V0={args.v0:.6g}, Q={args.q:.6g}\n")
    lines.append("- Criterion: chi2(compact) / chi2(compact_star) <= 5\n")
    lines.append(f"- Recommended compact_weight: {recommended['compact_weight']:.6g}\n")
    lines.append(f"- chi2_compact: {recommended['chi2_compact']:.6g}\n")
    lines.append(f"- chi2_compact_star: {recommended['chi2_compact_star']:.6g}\n")
    lines.append(f"- ratio: {recommended['ratio_compact_to_star']:.6g}\n")
    out_md.write_text("".join(lines), encoding="utf-8")
    print(f"Wrote calibration CSV to {args.output_csv}")
    print(f"Wrote calibration report to {args.output_md}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
