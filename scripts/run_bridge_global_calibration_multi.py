# -*- coding: utf-8 -*-
"""
Multi-point global bridge calibration.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Sequence

from fitting.data_io import load_forward_event_catalog_json
from fitting.strong_field_link import calibrate_global_bridge_coeffs_multi


def parse_args(argv: Sequence[str] | None = None) -> argparse.Namespace:
    p = argparse.ArgumentParser(description="Calibrate global bridge coefficients over multiple background points")
    p.add_argument("--forward_event_catalog_json", type=str, default="fitting/data/forward_event_catalog.json")
    p.add_argument("--gamma", type=float, default=0.0020339622)
    p.add_argument("--V0", type=float, default=0.012364227)
    p.add_argument("--Q", type=float, default=0.0)
    p.add_argument("--points_json", type=str, default=None, help="Optional JSON with explicit points list")
    p.add_argument("--grid_gamma_factors", type=str, default="0.7,1.0,1.4")
    p.add_argument("--grid_v0_factors", type=str, default="0.7,1.0,1.4")
    p.add_argument("--weight_center", type=float, default=2.0)
    p.add_argument("--rounds", type=int, default=40)
    p.add_argument("--reg_strength", type=float, default=0.2)
    p.add_argument("--theta_c_scale0", type=float, default=1.0)
    p.add_argument("--m2_scale0", type=float, default=10.0)
    p.add_argument("--lam_scale0", type=float, default=0.1)
    p.add_argument("--output_json", type=str, default="fitting/data/bridge_coeffs_global_multi.json")
    p.add_argument("--output_md", type=str, default="results/bridge_global_calibration_multi.md")
    return p.parse_args(argv)


def _parse_floats(raw: str) -> list[float]:
    out: list[float] = []
    for t in str(raw).split(","):
        s = t.strip()
        if s:
            out.append(float(s))
    if not out:
        raise ValueError("empty factors list")
    return out


def _build_points_from_grid(args: argparse.Namespace) -> list[dict]:
    gf = _parse_floats(args.grid_gamma_factors)
    vf = _parse_floats(args.grid_v0_factors)
    base_g = float(args.gamma)
    base_v = float(args.V0)
    points: list[dict] = []
    for gk in gf:
        for vk in vf:
            g = max(1e-6, base_g * float(gk))
            v = max(1e-6, base_v * float(vk))
            w = float(args.weight_center) if (abs(gk - 1.0) < 1e-12 and abs(vk - 1.0) < 1e-12) else 1.0
            points.append(
                {
                    "gamma": g,
                    "V0": v,
                    "Q": float(args.Q),
                    "weight": float(w),
                }
            )
    return points


def main(argv: Sequence[str] | None = None) -> int:
    args = parse_args(argv)
    events = load_forward_event_catalog_json(args.forward_event_catalog_json)

    points: list[dict]
    if args.points_json:
        points = json.loads(Path(args.points_json).read_text(encoding="utf-8"))
    else:
        points = _build_points_from_grid(args)

    coeffs = calibrate_global_bridge_coeffs_multi(
        points=points,
        events=events,
        rounds=int(args.rounds),
        reg_strength=float(args.reg_strength),
        theta_c_scale0=float(args.theta_c_scale0),
        m2_scale0=float(args.m2_scale0),
        lam_scale0=float(args.lam_scale0),
    )

    payload = {
        **coeffs,
        "points": points,
        "forward_event_catalog_json": args.forward_event_catalog_json,
    }
    out_json = Path(args.output_json)
    out_json.parent.mkdir(parents=True, exist_ok=True)
    out_json.write_text(json.dumps(payload, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")

    out_md = Path(args.output_md)
    out_md.parent.mkdir(parents=True, exist_ok=True)
    lines = []
    lines.append("# Global Bridge Calibration (Multi-point)\n\n")
    lines.append("- Mode: fixed global coefficients from multiple background points\n")
    lines.append(f"- Events: {len(events)}\n")
    lines.append(f"- Points: {len(points)}\n")
    lines.append(f"- rounds: {int(args.rounds)}\n")
    lines.append(f"- reg_strength: {float(args.reg_strength):.8g}\n\n")
    lines.append("## Coefficients\n")
    lines.append(f"- theta_c_scale: {float(coeffs['theta_c_scale']):.8g}\n")
    lines.append(f"- m2_scale: {float(coeffs['m2_scale']):.8g}\n")
    lines.append(f"- lam_scale: {float(coeffs['lam_scale']):.8g}\n")
    lines.append(f"- calibrated_chi2_proxy: {float(coeffs['calibrated_chi2_proxy']):.8g}\n")
    out_md.write_text("".join(lines), encoding="utf-8")

    print(f"Wrote {out_json}")
    print(f"Wrote {out_md}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
