# -*- coding: utf-8 -*-
"""Run background continuity equation check."""
from __future__ import annotations

import argparse
from typing import Sequence

from fitting.core_api import BackgroundParams as Params, integrate_background_lnN, omegas_from_km_s_Mpc
from checks.background_conservation import continuity_report_from_background


def parse_args(argv: Sequence[str] | None = None) -> argparse.Namespace:
    p = argparse.ArgumentParser(description="Run FRW continuity check")
    p.add_argument("--H0_km_s_Mpc", type=float, default=70.0)
    p.add_argument("--gamma", type=float, default=1e-3)
    p.add_argument("--V0", type=float, default=0.1)
    p.add_argument("--Q", type=float, default=0.0)
    p.add_argument("--a0", type=float, default=1e-3)
    p.add_argument("--N1", type=float, default=0.0)
    p.add_argument("--hN", type=float, default=0.02)
    p.add_argument("--x0", type=float, default=0.1)
    p.add_argument("--xdot0", type=float, default=0.0)
    p.add_argument("--psi0", type=float, default=0.05)
    p.add_argument("--psidot0", type=float, default=0.0)
    return p.parse_args(argv)


def main(argv: Sequence[str] | None = None) -> int:
    args = parse_args(argv)
    om = omegas_from_km_s_Mpc(args.H0_km_s_Mpc)
    p = Params(gamma=args.gamma, V0=args.V0, Q=args.Q, omegas=om)
    tr = integrate_background_lnN(
        a0=args.a0,
        x0=args.x0,
        xdot0=args.xdot0,
        psi0=args.psi0,
        psidot0=args.psidot0,
        N1=args.N1,
        p=p,
        hN=args.hN,
    )
    rep = continuity_report_from_background(tr)
    print(f"samples={rep.n_samples}")
    print(f"max_abs_residual={rep.max_abs_residual:.6e}")
    print(f"max_rel_residual={rep.max_rel_residual:.6e}")
    print(f"mean_rel_residual={rep.mean_rel_residual:.6e}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
