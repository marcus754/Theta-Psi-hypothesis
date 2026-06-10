# -*- coding: utf-8 -*-
"""
Report ADM mode diagnostics along a background trajectory.

Outputs CSV with background (t,a,H,theta,R) and, for a user-specified k-grid,
the finite-k dispersions ω^2 (eigenvalues of K^{-1}[G k^2/a^2 + M]). Also
prints min/max c_s^2 (high-k) over the trajectory.

Example:
  python -m src.report_adm --gamma 1.0 --V0 0.1 --H0_km_s_Mpc 70 \
      --t1 2.0 --h0 1e-2 --kmin 1e-2 --kmax 1.0 --nk 5 \
      --output results/adm_report.csv
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
)
from src.linear_modes_adm import cs2_eigs_adm, dispersion_eigs_adm


def parse_args(argv: Sequence[str] | None = None) -> argparse.Namespace:
    p = argparse.ArgumentParser(description="ADM modes report (3x3)")
    p.add_argument("--gamma", type=float, required=True)
    p.add_argument("--V0", type=float, required=True)
    p.add_argument("--Q", type=float, default=0.0)
    p.add_argument("--H0_km_s_Mpc", type=float, default=None)
    p.add_argument("--a0", type=float, default=1e-3)
    p.add_argument("--x0", type=float, default=0.1)
    p.add_argument("--xdot0", type=float, default=0.0)
    p.add_argument("--R0", type=float, default=0.05)
    p.add_argument("--Rdot0", type=float, default=0.0)
    p.add_argument("--t1", type=float, default=2.0)
    p.add_argument("--h0", type=float, default=1e-2)
    p.add_argument("--ctheta2", type=float, default=1.0)
    p.add_argument("--cR2", type=float, default=1.0)
    p.add_argument("--cphi2", type=float, default=1.0)
    p.add_argument("--epsKtr", type=float, default=0.0)
    p.add_argument("--epsKtphi", type=float, default=0.0)
    p.add_argument("--epsKRphi", type=float, default=0.0)
    p.add_argument("--epsGtr", type=float, default=0.0)
    p.add_argument("--epsGtphi", type=float, default=0.0)
    p.add_argument("--epsGRphi", type=float, default=0.0)
    p.add_argument("--alpha_k2_theta", type=float, default=0.0)
    p.add_argument("--alpha_k2_R", type=float, default=0.0)
    p.add_argument("--alpha_k2_phi", type=float, default=0.0)
    p.add_argument("--kmin", type=float, default=1e-2)
    p.add_argument("--kmax", type=float, default=1.0)
    p.add_argument("--nk", type=int, default=5)
    p.add_argument("--output", type=str, default="results/adm_report.csv")
    return p.parse_args(argv)


def main(argv: Sequence[str] | None = None) -> int:
    args = parse_args(argv)

    om = omegas_from_km_s_Mpc(args.H0_km_s_Mpc) if args.H0_km_s_Mpc is not None else Omegas()
    p = BGParams(gamma=args.gamma, V0=args.V0, Q=args.Q, omegas=om)
    y0 = (args.a0, args.x0, args.xdot0, args.R0, args.Rdot0)
    out = integrate_background_rk45(0.0, y0, t1=args.t1, p=p, h0=args.h0)

    # High-k sound speeds along the traj (independent of k in the current scaffold)
    cs2_min = float("inf"); cs2_max = 0.0

    # k-grid for finite-k dispersions
    import numpy as np
    ks = np.geomspace(args.kmin, args.kmax, max(args.nk, 2))

    rows = []
    for i in range(len(out["t"])):
        t = out["t"][i]; a = out["a"][i]; H = out["H"][i]
        th = out["theta"][i]; R = out["psi"][i]
        c1, c2, c3 = cs2_eigs_adm(
            th, R, a, 1.0, args.gamma, args.Q, args.V0,
            ctheta2=args.ctheta2, cR2=args.cR2, cphi2=args.cphi2,
            epsK_theta_R=args.epsKtr, epsK_theta_phi=args.epsKtphi, epsK_R_phi=args.epsKRphi,
            epsG_theta_R=args.epsGtr, epsG_theta_phi=args.epsGtphi, epsG_R_phi=args.epsGRphi,
            alpha_k2_theta=args.alpha_k2_theta, alpha_k2_R=args.alpha_k2_R, alpha_k2_phi=args.alpha_k2_phi,
        )
        cs2_min = min(cs2_min, c1, c2, c3)
        cs2_max = max(cs2_max, c1, c2, c3)
        for k in ks:
            w2 = dispersion_eigs_adm(
                th, R, a, float(k), args.gamma, args.Q, args.V0,
                ctheta2=args.ctheta2, cR2=args.cR2, cphi2=args.cphi2,
                epsK_theta_R=args.epsKtr, epsK_theta_phi=args.epsKtphi, epsK_R_phi=args.epsKRphi,
                epsG_theta_R=args.epsGtr, epsG_theta_phi=args.epsGtphi, epsG_R_phi=args.epsGRphi,
                alpha_k2_theta=args.alpha_k2_theta, alpha_k2_R=args.alpha_k2_R, alpha_k2_phi=args.alpha_k2_phi,
            )
            rows.append({
                "t": t, "a": a, "H": H, "theta": th, "R": R, "k": float(k),
                "w2_1": w2[0], "w2_2": w2[1], "w2_3": w2[2],
            })

    outpath = Path(args.output)
    outpath.parent.mkdir(parents=True, exist_ok=True)
    with outpath.open("w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        w.writeheader(); w.writerows(rows)
    print(f"Saved {len(rows)} rows to {outpath}")
    print(f"min c_s^2 = {cs2_min:.6g}, max c_s^2 = {cs2_max:.6g}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

