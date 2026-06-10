# -*- coding: utf-8 -*-
"""
Compute Mercury anomalous perihelion precession from Θ–Ψ PPN baseline.
"""
from __future__ import annotations

import argparse
from typing import Sequence

from checks.ppn_full import compute_ppn_params
from checks.mercury_precession import (
    MercuryObservation,
    mercury_precession_result,
)


def parse_args(argv: Sequence[str] | None = None) -> argparse.Namespace:
    p = argparse.ArgumentParser(description="Mercury perihelion check (PPN)")
    p.add_argument("--gamma_model", type=float, default=1.0)
    p.add_argument("--theta0", type=float, default=0.0)
    p.add_argument("--obs", type=float, default=42.98, help="Observed anomalous precession [arcsec/century]")
    p.add_argument("--sigma", type=float, default=0.04, help="Observed sigma [arcsec/century]")
    return p.parse_args(argv)


def main(argv: Sequence[str] | None = None) -> int:
    args = parse_args(argv)
    ppn = compute_ppn_params(gamma_model=args.gamma_model, theta0=args.theta0)
    obs = MercuryObservation(
        arcsec_per_century=float(args.obs),
        sigma_arcsec_per_century=float(args.sigma),
    )
    res = mercury_precession_result(gamma_ppn=ppn.gamma, beta_ppn=ppn.beta, obs=obs)

    print(f"PPN gamma={ppn.gamma:.6g}, beta={ppn.beta:.6g}")
    print(f"Predicted anomalous precession: {res.arcsec_per_century:.6f} arcsec/century")
    print(f"Observed anomalous precession:  {obs.arcsec_per_century:.6f} ± {obs.sigma_arcsec_per_century:.6f} arcsec/century")
    if res.pull_sigma is not None:
        print(f"Pull: {res.pull_sigma:.3f} sigma")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
