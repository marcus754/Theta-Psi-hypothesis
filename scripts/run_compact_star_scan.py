# -*- coding: utf-8 -*-
"""Stationary compact-star grid scan for the strong-field closure branch."""
from __future__ import annotations

import argparse
import csv
from pathlib import Path
from typing import Sequence

from fitting.core_api import scan_compact_star_grid


def parse_args(argv: Sequence[str] | None = None) -> argparse.Namespace:
    p = argparse.ArgumentParser(description="Run stationary compact-star scan")
    p.add_argument("--theta_c", type=float, nargs="+", required=True)
    p.add_argument("--m2", type=float, nargs="+", required=True)
    p.add_argument("--lam", type=float, nargs="+", required=True)
    p.add_argument("--r_max", type=float, default=20.0)
    p.add_argument("--dr", type=float, default=0.02)
    p.add_argument("--theta_scale", type=float, default=1.0)
    p.add_argument("--output", type=str, default="results/compact_star_scan.csv")
    return p.parse_args(argv)


def write_rows(rows: list[dict], path: str) -> None:
    if not rows:
        return
    out = Path(path)
    out.parent.mkdir(parents=True, exist_ok=True)
    fields = sorted(rows[0].keys())
    with out.open("w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=fields)
        w.writeheader()
        w.writerows(rows)


def main(argv: Sequence[str] | None = None) -> int:
    args = parse_args(argv)
    rows = scan_compact_star_grid(
        args.theta_c,
        args.m2,
        args.lam,
        r_max=args.r_max,
        dr=args.dr,
        theta_scale=args.theta_scale,
    )
    write_rows(rows, args.output)
    print(f"Wrote {len(rows)} rows to {args.output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
