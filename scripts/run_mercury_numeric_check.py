# -*- coding: utf-8 -*-
"""
Numerical Mercury perihelion check for Θ–Ψ weak-field sector.
"""
from __future__ import annotations

import argparse
from typing import Sequence

from checks.ppn_full import compute_ppn_params
from checks.mercury_precession import MercuryObservation
from checks.mercury_orbit_numeric import mercury_precession_numeric


def parse_args(argv: Sequence[str] | None = None) -> argparse.Namespace:
    p = argparse.ArgumentParser(description="Numerical Mercury precession from orbit integration")
    p.add_argument("--gamma_model", type=float, default=1.0)
    p.add_argument("--theta0", type=float, default=0.0)
    p.add_argument("--obs", type=float, default=42.98)
    p.add_argument("--sigma", type=float, default=0.04)
    p.add_argument("--n_orbits", type=int, default=14)
    p.add_argument("--steps_per_orbit", type=int, default=30000)
    return p.parse_args(argv)


def main(argv: Sequence[str] | None = None) -> int:
    args = parse_args(argv)
    ppn = compute_ppn_params(gamma_model=args.gamma_model, theta0=args.theta0)
    obs = MercuryObservation(args.obs, args.sigma)
    res = mercury_precession_numeric(
        gamma_ppn=ppn.gamma,
        beta_ppn=ppn.beta,
        n_orbits=args.n_orbits,
        steps_per_orbit=args.steps_per_orbit,
    )
    pull = (res.arcsec_per_century - obs.arcsec_per_century) / obs.sigma_arcsec_per_century

    print(f"PPN gamma={ppn.gamma:.6g}, beta={ppn.beta:.6g}")
    print(f"Numerical anomalous precession: {res.arcsec_per_century:.6f} arcsec/century")
    print(f"Observed anomalous precession:  {obs.arcsec_per_century:.6f} ± {obs.sigma_arcsec_per_century:.6f} arcsec/century")
    print(f"Perihelia used: {res.perihelion_count}")
    print(f"Pull: {pull:.3f} sigma")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
