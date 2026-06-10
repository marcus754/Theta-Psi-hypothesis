# -*- coding: utf-8 -*-
"""
Build compact constraints JSON from a strong-field bounds registry.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Sequence

from fitting.data_io import load_strong_field_bounds_json


def parse_args(argv: Sequence[str] | None = None) -> argparse.Namespace:
    p = argparse.ArgumentParser(description="Build compact constraints from bounds registry")
    p.add_argument("--bounds_json", type=str, default="fitting/data/strong_field_bounds_registry.json")
    p.add_argument("--output_json", type=str, default="fitting/data/compact_constraints_observational_proxy.json")
    p.add_argument("--output_md", type=str, default="results/compact_constraints_build.md")
    return p.parse_args(argv)


def main(argv: Sequence[str] | None = None) -> int:
    args = parse_args(argv)
    reg = load_strong_field_bounds_json(args.bounds_json)

    z_lowers = [e for e in reg["entries"] if e["type"] == "lower" and e["observable"] == "z_surface"]
    d_lowers = [e for e in reg["entries"] if e["type"] == "lower" and e["observable"] == "delay_proxy"]
    n_uppers = [e for e in reg["entries"] if e["type"] == "upper" and e["observable"] == "n_peak"]

    out = {"mode": "constraints"}
    if z_lowers:
        # strongest lower bound
        e = max(z_lowers, key=lambda x: x["value"])
        out["z_surface_min"] = float(e["value"])
        out["sigma_z_surface_min"] = float(e["sigma"])
    if d_lowers:
        e = max(d_lowers, key=lambda x: x["value"])
        out["delay_proxy_min"] = float(e["value"])
        out["sigma_delay_proxy_min"] = float(e["sigma"])
    if n_uppers:
        # strongest upper bound
        e = min(n_uppers, key=lambda x: x["value"])
        out["n_peak_max"] = float(e["value"])
        out["sigma_n_peak_max"] = float(e["sigma"])

    out["provenance"] = {
        "bounds_json": args.bounds_json,
        "n_entries": len(reg["entries"]),
        "method": "max lower / min upper aggregation",
    }

    out_json = Path(args.output_json)
    out_json.parent.mkdir(parents=True, exist_ok=True)
    out_json.write_text(json.dumps(out, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")

    lines = []
    lines.append("# Compact Constraints Build\n\n")
    lines.append(f"- Input registry: {args.bounds_json}\n")
    lines.append(f"- Entries: {len(reg['entries'])}\n")
    lines.append(f"- Output: {args.output_json}\n")
    lines.append("\n## Sources\n")
    for e in reg["entries"]:
        lines.append(
            f"- {e['name']}: {e['source']} ({e['observed_on']}) {e['source_url']}\n"
        )
    lines.append("\n## Aggregated constraints\n")
    for k in (
        "z_surface_min",
        "sigma_z_surface_min",
        "delay_proxy_min",
        "sigma_delay_proxy_min",
        "n_peak_max",
        "sigma_n_peak_max",
    ):
        if k in out:
            lines.append(f"- {k}: {out[k]}\n")

    out_md = Path(args.output_md)
    out_md.parent.mkdir(parents=True, exist_ok=True)
    out_md.write_text("".join(lines), encoding="utf-8")
    print(f"Wrote compact constraints JSON to {args.output_json}")
    print(f"Wrote compact constraints report to {args.output_md}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
