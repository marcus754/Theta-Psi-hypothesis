# -*- coding: utf-8 -*-
"""
Build exclusion map for stationary compact-star scan.
"""
from __future__ import annotations

import argparse
import csv
from pathlib import Path
from typing import Sequence

from fitting.data_io import load_compact_star_targets_json


def parse_args(argv: Sequence[str] | None = None) -> argparse.Namespace:
    p = argparse.ArgumentParser(description="Compact-star exclusion map")
    p.add_argument("--scan_csv", type=str, default="results/compact_star_scan.csv")
    p.add_argument("--targets_json", type=str, default="fitting/data/compact_star_demo_targets.json")
    p.add_argument("--chi2_max", type=float, default=9.0, help="Allowed if chi2 <= chi2_max")
    p.add_argument("--output_csv", type=str, default="results/compact_star_exclusion.csv")
    p.add_argument("--output_md", type=str, default="results/compact_star_exclusion.md")
    return p.parse_args(argv)


def main(argv: Sequence[str] | None = None) -> int:
    args = parse_args(argv)
    tgt = load_compact_star_targets_json(args.targets_json)
    rows = []
    with Path(args.scan_csv).open("r", newline="", encoding="utf-8") as fh:
        rd = csv.DictReader(fh)
        for r in rd:
            z = float(r["z_surface"])
            d = float(r["delay_proxy"])
            chi2 = (
                ((z - tgt["z_surface"]) / tgt["sigma_z_surface"]) ** 2
                + ((d - tgt["delay_proxy"]) / tgt["sigma_delay_proxy"]) ** 2
            )
            r2 = dict(r)
            r2["chi2_compact_star"] = chi2
            r2["allowed"] = chi2 <= args.chi2_max
            rows.append(r2)

    if rows:
        out_csv = Path(args.output_csv)
        out_csv.parent.mkdir(parents=True, exist_ok=True)
        flds = list(rows[0].keys())
        with out_csv.open("w", newline="", encoding="utf-8") as fh:
            wr = csv.DictWriter(fh, fieldnames=flds)
            wr.writeheader()
            wr.writerows(rows)

    n_all = len(rows)
    n_ok = sum(1 for r in rows if bool(r["allowed"]))
    best = min(rows, key=lambda r: float(r["chi2_compact_star"])) if rows else None

    lines = []
    lines.append("# Compact-Star Exclusion Map\n\n")
    lines.append(f"- Input scan: `{args.scan_csv}`\n")
    lines.append(f"- Targets: `{args.targets_json}`\n")
    lines.append(f"- Threshold: chi2 <= {args.chi2_max}\n")
    lines.append(f"- Total rows: {n_all}\n")
    lines.append(f"- Allowed rows: {n_ok}\n")
    lines.append(f"- Excluded rows: {n_all - n_ok}\n")
    if best is not None:
        lines.append("\n## Best Row\n")
        lines.append(f"- theta_c: {best['theta_c']}\n")
        lines.append(f"- m2: {best['m2']}\n")
        lines.append(f"- lam: {best['lam']}\n")
        lines.append(f"- z_surface: {float(best['z_surface']):.6g}\n")
        lines.append(f"- delay_proxy: {float(best['delay_proxy']):.6g}\n")
        lines.append(f"- chi2: {float(best['chi2_compact_star']):.6g}\n")

    out_md = Path(args.output_md)
    out_md.parent.mkdir(parents=True, exist_ok=True)
    out_md.write_text("".join(lines), encoding="utf-8")
    print(f"Wrote exclusion CSV to {args.output_csv}")
    print(f"Wrote exclusion report to {args.output_md}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
