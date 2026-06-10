# -*- coding: utf-8 -*-
"""Single entrypoint for strong-field workflows."""
from __future__ import annotations

import argparse
import subprocess
import sys
from typing import Sequence


def parse_args(argv: Sequence[str] | None = None) -> argparse.Namespace:
    p = argparse.ArgumentParser(description="Run strong-field workflow")
    p.add_argument("--mode", choices=("quick", "public"), default="quick")
    p.add_argument("--forward_events_csv", type=str, default=None)
    p.add_argument("--outdir", type=str, default="results/strong_field")
    return p.parse_args(argv)


def _run(cmd: list[str]) -> int:
    return int(subprocess.run(cmd).returncode)


def main(argv: Sequence[str] | None = None) -> int:
    args = parse_args(argv)
    py = sys.executable

    if args.mode == "public":
        if not args.forward_events_csv:
            raise SystemExit("For --mode public pass --forward_events_csv <path.csv>")
        return _run(
            [
                py,
                "-m",
                "scripts.run_forward_public_pipeline",
                "--forward_events_csv",
                args.forward_events_csv,
                "--outdir",
                args.outdir,
            ]
        )

    # quick mode: uses in-repo defaults and produces strong-field event artifacts
    rc = _run(
        [
            py,
            "-m",
            "scripts.run_forward_event_scan",
            "--event_catalog_json",
            "fitting/data/forward_event_catalog.json",
            "--output_csv",
            "results/forward_event_scan.csv",
            "--output_md",
            "results/forward_event_scan.md",
        ]
    )
    if rc != 0:
        return rc
    return _run(
        [
            py,
            "-m",
            "scripts.run_forward_event_report",
            "--scan_csv",
            "results/forward_event_scan.csv",
            "--event_catalog_json",
            "fitting/data/forward_event_catalog.json",
            "--output_json",
            "results/forward_event_report.json",
            "--output_md",
            "results/forward_event_report.md",
        ]
    )


if __name__ == "__main__":
    raise SystemExit(main())
