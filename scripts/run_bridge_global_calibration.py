# -*- coding: utf-8 -*-
"""
One-shot global bridge calibration for Theta-Psi strong-field link.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Sequence

from fitting.data_io import load_forward_event_catalog_json
from fitting.strong_field_link import calibrate_global_bridge_coeffs


def parse_args(argv: Sequence[str] | None = None) -> argparse.Namespace:
    p = argparse.ArgumentParser(description="Calibrate global bridge coefficients once (no event-level inner tuning)")
    p.add_argument("--forward_event_catalog_json", type=str, default="fitting/data/forward_event_catalog.json")
    p.add_argument("--gamma", type=float, default=0.0020)
    p.add_argument("--V0", type=float, default=0.0124)
    p.add_argument("--Q", type=float, default=0.0)
    p.add_argument("--rounds", type=int, default=36)
    p.add_argument("--reg_strength", type=float, default=0.2)
    p.add_argument("--theta_c_scale0", type=float, default=1.0)
    p.add_argument("--m2_scale0", type=float, default=10.0)
    p.add_argument("--lam_scale0", type=float, default=0.1)
    p.add_argument("--output_json", type=str, default="fitting/data/bridge_coeffs_global.json")
    p.add_argument("--output_md", type=str, default="results/bridge_global_calibration.md")
    return p.parse_args(argv)


def main(argv: Sequence[str] | None = None) -> int:
    args = parse_args(argv)
    events = load_forward_event_catalog_json(args.forward_event_catalog_json)
    coeffs = calibrate_global_bridge_coeffs(
        gamma=float(args.gamma),
        V0=float(args.V0),
        Q=float(args.Q),
        events=events,
        rounds=int(args.rounds),
        reg_strength=float(args.reg_strength),
        theta_c_scale0=float(args.theta_c_scale0),
        m2_scale0=float(args.m2_scale0),
        lam_scale0=float(args.lam_scale0),
    )

    out_json = Path(args.output_json)
    out_json.parent.mkdir(parents=True, exist_ok=True)
    out_json.write_text(json.dumps(coeffs, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")

    out_md = Path(args.output_md)
    out_md.parent.mkdir(parents=True, exist_ok=True)
    lines = []
    lines.append("# Global Bridge Calibration\n\n")
    lines.append("- Mode: fixed global coefficients (no event-level inner tuning)\n")
    lines.append(f"- Events: {len(events)}\n")
    lines.append(f"- Calibration point: gamma={float(args.gamma):.8g}, V0={float(args.V0):.8g}, Q={float(args.Q):.8g}\n")
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
