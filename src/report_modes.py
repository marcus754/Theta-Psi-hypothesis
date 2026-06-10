# -*- coding: utf-8 -*-
"""
Generate a background + modes report for the Θ–Ψ model.

Integrates the background with RK45 and evaluates linear-mode sound speeds
along the trajectory using the current 2×2 scaffold (δθ, δΨ) from src.linear_modes.

CLI example:
  python -m src.report_modes --gamma 1.0 --V0 0.1 --H0_km_s_Mpc 70 --Q 0.0 \
      --t1 2.0 --h0 1e-2 --output results/modes_report.csv
"""
from __future__ import annotations

import argparse
import csv
from pathlib import Path
from typing import Sequence

from background.components import omegas_from_km_s_Mpc, Omegas
from background.frw_background import (
    Params as BGParams,
    integrate_background_rk45,
    theta_from_x,
    theta_dot_from_x,
)
from src.linear_modes import cs2_eigs


def parse_args(argv: Sequence[str] | None = None) -> argparse.Namespace:
    p = argparse.ArgumentParser(description="Background + linear modes report")
    p.add_argument("--gamma", type=float, required=True)
    p.add_argument("--V0", type=float, required=True)
    p.add_argument("--Q", type=float, default=0.0)
    p.add_argument("--H0_km_s_Mpc", type=float, default=None)
    p.add_argument("--a0", type=float, default=1e-3)
    p.add_argument("--x0", type=float, default=0.1)
    p.add_argument("--xdot0", type=float, default=0.0)
    p.add_argument("--psi0", type=float, default=0.05)
    p.add_argument("--psidot0", type=float, default=0.0)
    p.add_argument("--t1", type=float, default=2.0)
    p.add_argument("--h0", type=float, default=1e-2)
    p.add_argument("--ctheta2", type=float, default=1.0)
    p.add_argument("--cpsi2", type=float, default=1.0)
    p.add_argument("--output", type=str, default="results/modes_report.csv")
    return p.parse_args(argv)


def main(argv: Sequence[str] | None = None) -> int:
    args = parse_args(argv)

    if args.H0_km_s_Mpc is not None:
        om = omegas_from_km_s_Mpc(args.H0_km_s_Mpc)
    else:
        om = Omegas()

    bgp = BGParams(gamma=args.gamma, V0=args.V0, Q=args.Q, omegas=om)
    y0 = (args.a0, args.x0, args.xdot0, args.psi0, args.psidot0)
    out = integrate_background_rk45(0.0, y0, t1=args.t1, p=bgp, h0=args.h0)

    rows = []
    for t, a, H, x, theta, xdot, psi in zip(
        out["t"], out["a"], out["H"], out["x"], out["theta"], out["theta_dot"], out["psi"],
    ):
        # 2×2 scaffold uses only theta and gamma
        s1, s2 = cs2_eigs(theta, args.gamma, args.ctheta2, args.cpsi2)
        rows.append({
            "t": t,
            "a": a,
            "H": H,
            "theta": theta,
            "psi": psi,
            "cs2_1": s1,
            "cs2_2": s2,
        })

    outpath = Path(args.output)
    outpath.parent.mkdir(parents=True, exist_ok=True)
    with outpath.open("w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        w.writeheader()
        w.writerows(rows)

    min_c = min(min(r["cs2_1"], r["cs2_2"]) for r in rows)
    max_c = max(max(r["cs2_1"], r["cs2_2"]) for r in rows)
    print(f"Saved {len(rows)} rows to {outpath}")
    print(f"min c_s^2 = {min_c:.6g}, max c_s^2 = {max_c:.6g}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

