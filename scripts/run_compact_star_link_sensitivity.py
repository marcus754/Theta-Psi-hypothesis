# -*- coding: utf-8 -*-
"""
Sensitivity scan for compact-star link-model coefficients.
"""
from __future__ import annotations

import argparse
import csv
from pathlib import Path
from typing import Sequence

from fitting.run_mcmc import Params
from fitting.likelihoods_compact_star import CompactStarTargets
from fitting.run_mcmc import loglikes_components
from fitting.core_api import Omegas


def parse_args(argv: Sequence[str] | None = None) -> argparse.Namespace:
    p = argparse.ArgumentParser(description="Sensitivity of compact-star link coefficients")
    p.add_argument("--output_csv", type=str, default="results/compact_star_link_sensitivity.csv")
    p.add_argument("--output_md", type=str, default="results/compact_star_link_sensitivity.md")
    return p.parse_args(argv)


def main(argv: Sequence[str] | None = None) -> int:
    args = parse_args(argv)
    om = Omegas(H0=1.0)
    p0 = Params(gamma=0.5, V0=0.08, Q=0.0)
    target = CompactStarTargets(z_surface=350.0, sigma_z_surface=175.0, delay_proxy=700.0, sigma_delay_proxy=350.0)

    coeff_grid = [0.5, 1.0, 1.5, 2.0]
    rows = []
    for cg in coeff_grid:
        for cv in coeff_grid:
            comps = loglikes_components(
                p0,
                om=om,
                compact_star=target,
                compact_star_link_model=True,
                compact_star_theta_c_from_gamma=1.5 * cg,
                compact_star_m2_from_v0=10.0 * cv,
                compact_star_lam_from_q2=0.1,
            )
            ll = float(comps["compact_star"])
            rows.append(
                {
                    "theta_c_from_gamma": 1.5 * cg,
                    "m2_from_v0": 10.0 * cv,
                    "compact_star_loglike": ll,
                    "compact_star_chi2": -2.0 * ll,
                }
            )

    out_csv = Path(args.output_csv)
    out_csv.parent.mkdir(parents=True, exist_ok=True)
    with out_csv.open("w", newline="", encoding="utf-8") as fh:
        wr = csv.DictWriter(fh, fieldnames=list(rows[0].keys()))
        wr.writeheader()
        wr.writerows(rows)

    best = min(rows, key=lambda r: r["compact_star_chi2"])
    lines = []
    lines.append("# Compact-Star Link Sensitivity\n\n")
    lines.append(f"- Rows scanned: {len(rows)}\n")
    lines.append(f"- Reference params: gamma={p0.gamma}, V0={p0.V0}, Q={p0.Q}\n")
    lines.append("\n## Best coefficients\n")
    lines.append(f"- theta_c_from_gamma: {best['theta_c_from_gamma']:.6g}\n")
    lines.append(f"- m2_from_v0: {best['m2_from_v0']:.6g}\n")
    lines.append(f"- compact_star_chi2: {best['compact_star_chi2']:.6g}\n")
    out_md = Path(args.output_md)
    out_md.parent.mkdir(parents=True, exist_ok=True)
    out_md.write_text("".join(lines), encoding="utf-8")
    print(f"Wrote sensitivity CSV to {args.output_csv}")
    print(f"Wrote sensitivity report to {args.output_md}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
