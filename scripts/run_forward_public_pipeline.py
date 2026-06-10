# -*- coding: utf-8 -*-
"""
End-to-end strong-field forward pipeline from public-style CSV input.

Flow:
1) Convert CSV -> forward_event_catalog.json
2) Run forward event scan
3) Build event-level report
4) Run quick empirical checks
"""
from __future__ import annotations

import argparse
import subprocess
import sys
from pathlib import Path
from typing import Sequence


def parse_args(argv: Sequence[str] | None = None) -> argparse.Namespace:
    p = argparse.ArgumentParser(description="Run forward pipeline from public strong-field CSV")
    p.add_argument("--forward_events_csv", type=str, required=True)
    p.add_argument("--outdir", type=str, default="results/public_forward")
    p.add_argument("--joint_scan_csv", type=str, default="results/joint_region_scan.csv")
    p.add_argument("--n_theta_c", type=int, default=16)
    p.add_argument("--n_m2", type=int, default=16)
    return p.parse_args(argv)


def _run(cmd: list[str]) -> tuple[int, str]:
    proc = subprocess.run(cmd, capture_output=True, text=True)
    text = ((proc.stdout or "") + (proc.stderr or "")).strip()
    return proc.returncode, text


def main(argv: Sequence[str] | None = None) -> int:
    args = parse_args(argv)
    py = sys.executable
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    catalog_json = outdir / "forward_event_catalog.json"
    scan_csv = outdir / "forward_event_scan.csv"
    scan_md = outdir / "forward_event_scan.md"
    report_json = outdir / "forward_event_report.json"
    report_md = outdir / "forward_event_report.md"
    checks_md = outdir / "empirical_checks.md"
    summary_md = outdir / "pipeline_summary.md"

    steps: list[tuple[str, list[str]]] = [
        (
            "convert_forward_events",
            [
                py,
                "-m",
                "fitting.convert_public",
                "--forward_events",
                args.forward_events_csv,
                "--outdir",
                str(outdir),
            ],
        ),
        (
            "forward_event_scan",
            [
                py,
                "-m",
                "scripts.run_forward_event_scan",
                "--event_catalog_json",
                str(catalog_json),
                "--n_theta_c",
                str(args.n_theta_c),
                "--n_m2",
                str(args.n_m2),
                "--output_csv",
                str(scan_csv),
                "--output_md",
                str(scan_md),
            ],
        ),
        (
            "forward_event_report",
            [
                py,
                "-m",
                "scripts.run_forward_event_report",
                "--scan_csv",
                str(scan_csv),
                "--event_catalog_json",
                str(catalog_json),
                "--output_json",
                str(report_json),
                "--output_md",
                str(report_md),
            ],
        ),
    ]
    if Path(args.joint_scan_csv).exists():
        steps.append(
            (
                "empirical_checks",
                [
                    py,
                    "-m",
                    "scripts.run_empirical_checks_now",
                    "--forward_event_scan_csv",
                    str(scan_csv),
                    "--forward_event_catalog_json",
                    str(catalog_json),
                    "--joint_scan_csv",
                    args.joint_scan_csv,
                    "--output_md",
                    str(checks_md),
                ],
            )
        )

    results: list[dict] = []
    for name, cmd in steps:
        rc, out = _run(cmd)
        results.append({"name": name, "rc": rc, "output": out})
        if rc != 0:
            break

    lines = []
    lines.append("# Forward Public Pipeline\n\n")
    ok = all(r["rc"] == 0 for r in results) and len(results) == len(steps)
    lines.append(f"- Status: {'PASS' if ok else 'FAIL'}\n")
    lines.append(f"- Input CSV: {args.forward_events_csv}\n")
    lines.append(f"- Output dir: {outdir}\n")
    lines.append("\n## Steps\n")
    for r in results:
        lines.append(f"- {r['name']}: {'PASS' if r['rc']==0 else 'FAIL'} (rc={r['rc']})\n")
    lines.append("\n## Notes\n")
    for r in results:
        snippet = r["output"].splitlines()[-1] if r["output"] else ""
        lines.append(f"- {r['name']}: {snippet}\n")
    summary_md.write_text("".join(lines), encoding="utf-8")

    print(f"Wrote summary to {summary_md}")
    return 0 if ok else 1


if __name__ == "__main__":
    raise SystemExit(main())
