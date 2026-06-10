# -*- coding: utf-8 -*-
"""
Collect the highest-risk tests for Θ–Ψ in one place.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Sequence


def parse_args(argv: Sequence[str] | None = None) -> argparse.Namespace:
    p = argparse.ArgumentParser(description="Write killer-tests report")
    p.add_argument("--falsifiability_json", type=str, default="results/falsifiability_check.json")
    p.add_argument("--head2head_json", type=str, default="results/head2head_comparison.json")
    p.add_argument("--forward_event_report_json", type=str, default="results/forward_event_report.json")
    p.add_argument("--output_md", type=str, default="results/killer_tests.md")
    return p.parse_args(argv)


def main(argv: Sequence[str] | None = None) -> int:
    args = parse_args(argv)
    fals = json.loads(Path(args.falsifiability_json).read_text(encoding="utf-8"))
    h2h = json.loads(Path(args.head2head_json).read_text(encoding="utf-8")) if Path(args.head2head_json).exists() else None
    fev = json.loads(Path(args.forward_event_report_json).read_text(encoding="utf-8"))

    lines = ["# Killer Tests\n\n"]
    lines.append("Это не список того, что уже красиво выглядит, а того, где теория может быть убита.\n\n")
    lines.append("## Current Status\n")
    lines.append(f"- empirical falsifiability: {'PASS' if fals['overall_pass'] else 'FAIL'}\n")
    lines.append(f"- forward-event max pull: {fev.get('max_abs_pull_sigma')}\n")
    if h2h is not None:
        dchi2 = float(h2h["delta_theta_minus_lcdm"]["chi2"])
        lines.append(f"- head-to-head Δχ²(Θ–Ψ - LCDM): {dchi2}\n")
    lines.append("\n## Real Killer Tests\n")
    lines.append("- Full covariant action consistency: one action must generate FRW, weak-field and compact-star limits without manual branch switching.\n")
    lines.append("- Strong-field horizon replacement: finite redshift/delay must survive event-level compact-object transfer, not only stationary closure.\n")
    lines.append("- Head-to-head against ΛCDM/GR on unified data: surviving is not enough; theory must either win somewhere cleanly or at least not need hidden flexibility.\n")
    lines.append("- Future cluster / lensing tests: if Θ–Ψ is extended toward dark-sector claims, Bullet-like systems become killer tests.\n")
    lines.append("\n## Interpretation\n")
    if h2h is not None and float(h2h["delta_theta_minus_lcdm"]["chi2"]) < 0.0:
        lines.append("- Right now the theory survives empirical contour and beats the current flat ΛCDM/GR comparator in this repo's head-to-head setup.\n")
    else:
        lines.append("- Right now the theory survives empirical contour but does not yet establish a robust global advantage over the baseline.\n")
    Path(args.output_md).write_text("".join(lines), encoding="utf-8")
    print(f"Wrote killer-tests report to {args.output_md}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
