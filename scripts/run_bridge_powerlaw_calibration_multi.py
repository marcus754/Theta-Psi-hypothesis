# -*- coding: utf-8 -*-
"""
Calibrate functional bridge model (powerlaw2d) over multiple background points.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Sequence

from fitting.data_io import load_forward_event_catalog_json
from fitting.strong_field_link import calibrate_bridge_powerlaw2d_multi


def parse_args(argv: Sequence[str] | None = None) -> argparse.Namespace:
    p = argparse.ArgumentParser(description="Calibrate powerlaw2d bridge over multiple points")
    p.add_argument("--forward_event_catalog_json", type=str, default="fitting/data/forward_event_catalog.json")
    p.add_argument("--gamma", type=float, default=0.0020339622)
    p.add_argument("--V0", type=float, default=0.012364227)
    p.add_argument("--Q", type=float, default=0.0)
    p.add_argument("--grid_gamma_factors", type=str, default="0.7,1.0,1.4")
    p.add_argument("--grid_v0_factors", type=str, default="0.7,1.0,1.4")
    p.add_argument("--weight_center", type=float, default=2.0)
    p.add_argument("--rounds", type=int, default=36)
    p.add_argument("--reg_strength", type=float, default=0.05)
    p.add_argument("--gamma_ref", type=float, default=0.0020339622)
    p.add_argument("--v0_ref", type=float, default=0.012364227)
    p.add_argument("--theta0", type=float, default=1.0)
    p.add_argument("--m20", type=float, default=10.0)
    p.add_argument("--lam0", type=float, default=0.1)
    p.add_argument("--output_json", type=str, default="fitting/data/bridge_coeffs_powerlaw2d.json")
    p.add_argument("--output_md", type=str, default="results/bridge_powerlaw2d_calibration.md")
    return p.parse_args(argv)


def _parse_floats(raw: str) -> list[float]:
    out: list[float] = []
    for t in str(raw).split(","):
        s = t.strip()
        if s:
            out.append(float(s))
    if not out:
        raise ValueError("empty factor list")
    return out


def main(argv: Sequence[str] | None = None) -> int:
    args = parse_args(argv)
    events = load_forward_event_catalog_json(args.forward_event_catalog_json)
    gf = _parse_floats(args.grid_gamma_factors)
    vf = _parse_floats(args.grid_v0_factors)
    points: list[dict] = []
    for gk in gf:
        for vk in vf:
            g = max(1e-6, float(args.gamma) * float(gk))
            v = max(1e-6, float(args.V0) * float(vk))
            w = float(args.weight_center) if (abs(gk - 1.0) < 1e-12 and abs(vk - 1.0) < 1e-12) else 1.0
            points.append(
                {
                    "gamma": g,
                    "V0": v,
                    "Q": float(args.Q),
                    "weight": float(w),
                }
            )

    coeffs = calibrate_bridge_powerlaw2d_multi(
        points=points,
        events=events,
        rounds=int(args.rounds),
        reg_strength=float(args.reg_strength),
        gamma_ref=float(args.gamma_ref),
        V0_ref=float(args.v0_ref),
        theta0=float(args.theta0),
        m20=float(args.m20),
        lam0=float(args.lam0),
    )
    payload = {**coeffs, "points": points, "forward_event_catalog_json": args.forward_event_catalog_json}
    out_json = Path(args.output_json)
    out_json.parent.mkdir(parents=True, exist_ok=True)
    out_json.write_text(json.dumps(payload, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")

    out_md = Path(args.output_md)
    out_md.parent.mkdir(parents=True, exist_ok=True)
    lines = []
    lines.append("# Bridge Calibration (powerlaw2d)\n\n")
    lines.append(f"- points: {len(points)}\n")
    lines.append(f"- rounds: {int(args.rounds)}\n")
    lines.append(f"- reg_strength: {float(args.reg_strength):.6g}\n")
    lines.append(f"- calibrated_chi2_proxy: {float(coeffs['calibrated_chi2_proxy']):.6f}\n\n")
    lines.append("## theta_c\n")
    lines.append(f"- a0: {float(coeffs['theta_c']['a0']):.8g}\n")
    lines.append(f"- ag: {float(coeffs['theta_c']['ag']):.8g}\n")
    lines.append(f"- av: {float(coeffs['theta_c']['av']):.8g}\n")
    lines.append("## m2\n")
    lines.append(f"- a0: {float(coeffs['m2']['a0']):.8g}\n")
    lines.append(f"- ag: {float(coeffs['m2']['ag']):.8g}\n")
    lines.append(f"- av: {float(coeffs['m2']['av']):.8g}\n")
    lines.append("## lam\n")
    lines.append(f"- a0: {float(coeffs['lam']['a0']):.8g}\n")
    lines.append(f"- ag: {float(coeffs['lam']['ag']):.8g}\n")
    lines.append(f"- av: {float(coeffs['lam']['av']):.8g}\n")
    out_md.write_text("".join(lines), encoding="utf-8")
    print(f"Wrote {out_json}")
    print(f"Wrote {out_md}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
