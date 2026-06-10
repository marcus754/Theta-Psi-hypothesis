# -*- coding: utf-8 -*-
"""
Scan forward EHT/GW observables over compact-star parameter grid.
"""
from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path
from typing import Sequence

from fitting.core_api import CompactStarParams, solve_compact_star_profile
from fitting.forward_observables import (
    forward_observables_from_compact_metrics,
    calibrate_v1_coefficients,
)
from fitting.forward_transfer import explicit_forward_observables
from fitting.likelihoods_forward import ForwardTargets, loglike_forward
from fitting.data_io import load_forward_targets_json


def _linspace(a: float, b: float, n: int) -> list[float]:
    if n <= 1:
        return [a]
    h = (b - a) / (n - 1)
    return [a + i * h for i in range(n)]


def parse_args(argv: Sequence[str] | None = None) -> argparse.Namespace:
    p = argparse.ArgumentParser(description="Forward observables scan (EHT/GW observables)")
    p.add_argument("--theta_c_min", type=float, default=0.1)
    p.add_argument("--theta_c_max", type=float, default=3.0)
    p.add_argument("--n_theta_c", type=int, default=20)
    p.add_argument("--m2_min", type=float, default=0.05)
    p.add_argument("--m2_max", type=float, default=2.0)
    p.add_argument("--n_m2", type=int, default=20)
    p.add_argument("--lam", type=float, default=0.001)
    p.add_argument("--targets_json", type=str, default="fitting/data/forward_observational_proxy.json")
    p.add_argument("--mode", type=str, default="explicit", choices=["explicit", "v1_calibrated", "v1_fixed"])
    p.add_argument("--r_max", type=float, default=20.0)
    p.add_argument("--dr", type=float, default=0.02)
    p.add_argument("--theta_scale", type=float, default=1.0)
    p.add_argument("--psi_grad_weight", type=float, default=1.0)
    p.add_argument("--theta_psi_coupling", type=float, default=1.0)
    p.add_argument("--refractive_strength", type=float, default=1.0)
    p.add_argument("--refractive_j_star", type=float, default=1.0)
    p.add_argument("--kappa_refr", type=float, default=1.0)
    p.add_argument("--k_ring_m87", type=float, default=3.0)
    p.add_argument("--k_ring_sgra", type=float, default=3.6)
    p.add_argument("--k_echo_delay", type=float, default=1.0)
    p.add_argument("--ref_theta_c", type=float, default=1.5)
    p.add_argument("--ref_m2", type=float, default=1.0)
    p.add_argument("--ref_lam", type=float, default=0.1)
    p.add_argument("--ring_barrier_r_code", type=float, default=3.0)
    p.add_argument("--ring_band_index", type=int, default=2)
    p.add_argument("--ring_n_vac_tol", type=float, default=0.01)
    p.add_argument("--ring_r_min_code", type=float, default=2.0)
    p.add_argument("--echo_mass_msun", type=float, default=60.0)
    p.add_argument("--echo_barrier_r_code", type=float, default=3.0)
    p.add_argument("--echo_log_enhancement_coeff", type=float, default=0.6366197723675814)
    p.add_argument("--output_csv", type=str, default="results/forward_observables_scan.csv")
    p.add_argument("--output_md", type=str, default="results/forward_observables_scan.md")
    return p.parse_args(argv)


def main(argv: Sequence[str] | None = None) -> int:
    args = parse_args(argv)
    tj = load_forward_targets_json(args.targets_json)
    tgt = ForwardTargets(
        ring_m87_uas=tj["ring_m87_uas"],
        sigma_ring_m87_uas=tj["sigma_ring_m87_uas"],
        ring_sgra_uas=tj["ring_sgra_uas"],
        sigma_ring_sgra_uas=tj["sigma_ring_sgra_uas"],
        echo_delay_ms=tj.get("echo_delay_ms"),
        sigma_echo_delay_ms=tj.get("sigma_echo_delay_ms"),
    )
    k_cal = None
    if args.mode == "v1_calibrated":
        rp = CompactStarParams(
            theta_c=args.ref_theta_c,
            m2=args.ref_m2,
            lam=args.ref_lam,
            r_max=args.r_max,
            dr=args.dr,
            theta_scale=args.theta_scale,
            psi_grad_weight=args.psi_grad_weight,
            theta_psi_coupling=args.theta_psi_coupling,
            refractive_strength=args.refractive_strength,
            refractive_j_star=args.refractive_j_star,
            kappa_refr=args.kappa_refr,
        )
        ro = solve_compact_star_profile(rp)
        ref_met = {
            "z_surface": float(ro["z_surface"]),
            "delay_proxy": float(ro["delay_proxy"]),
            "n_peak": float(ro["n_peak"]),
        }
        k_cal = calibrate_v1_coefficients(
            reference_metrics=ref_met,
            target_ring_m87_uas=tgt.ring_m87_uas,
            target_ring_sgra_uas=tgt.ring_sgra_uas,
            target_echo_delay_ms=tgt.echo_delay_ms,
        )

    rows = []
    for tc in _linspace(args.theta_c_min, args.theta_c_max, args.n_theta_c):
        for m2 in _linspace(args.m2_min, args.m2_max, args.n_m2):
            p = CompactStarParams(
                theta_c=tc,
                m2=m2,
                lam=args.lam,
                r_max=args.r_max,
                dr=args.dr,
                theta_scale=args.theta_scale,
                psi_grad_weight=args.psi_grad_weight,
                theta_psi_coupling=args.theta_psi_coupling,
                refractive_strength=args.refractive_strength,
                refractive_j_star=args.refractive_j_star,
                kappa_refr=args.kappa_refr,
            )
            cm = solve_compact_star_profile(p)
            met = {
                "z_surface": float(cm["z_surface"]),
                "delay_proxy": float(cm["delay_proxy"]),
                "n_peak": float(cm["n_peak"]),
            }
            if args.mode == "explicit":
                obs = explicit_forward_observables(
                    cm,
                    ring_barrier_r_code=args.ring_barrier_r_code,
                    ring_band_index=args.ring_band_index,
                    ring_n_vac_tol=args.ring_n_vac_tol,
                    ring_r_min_code=args.ring_r_min_code,
                    echo_mass_msun=args.echo_mass_msun,
                    echo_barrier_r_code=args.echo_barrier_r_code,
                    echo_log_enhancement_coeff=args.echo_log_enhancement_coeff,
                )
            elif args.mode == "v1_calibrated":
                obs = forward_observables_from_compact_metrics(met, **k_cal)
            else:
                obs = forward_observables_from_compact_metrics(
                    met,
                    k_ring_m87=args.k_ring_m87,
                    k_ring_sgra=args.k_ring_sgra,
                    k_echo_delay=args.k_echo_delay,
                )
            ll = loglike_forward(obs, tgt)
            rows.append(
                {
                    "theta_c": tc,
                    "m2": m2,
                    "lam": args.lam,
                    "psi_grad_weight": args.psi_grad_weight,
                    "theta_psi_coupling": args.theta_psi_coupling,
                    "refractive_strength": args.refractive_strength,
                    "refractive_j_star": args.refractive_j_star,
                    "kappa_refr": args.kappa_refr,
                    "z_surface": float(cm["z_surface"]),
                    "delay_proxy": float(cm["delay_proxy"]),
                    "ring_m87_uas": obs["ring_m87_uas"],
                    "ring_sgra_uas": obs["ring_sgra_uas"],
                    "chi2_forward": -2.0 * float(ll),
                }
            )

    best = min(rows, key=lambda r: float(r["chi2_forward"]))

    out_csv = Path(args.output_csv)
    out_csv.parent.mkdir(parents=True, exist_ok=True)
    with out_csv.open("w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()))
        w.writeheader()
        w.writerows(rows)

    lines = []
    lines.append("# Forward Observables Scan\n\n")
    lines.append("- Components: compact_star -> forward(EHT/GW observables) -> chi2\n")
    lines.append(f"- Forward mode: {args.mode}\n")
    lines.append(f"- Grid: theta_c x m2 = {args.n_theta_c} x {args.n_m2}, lam={args.lam}\n")
    lines.append(
        f"- Strong-field closure: psi_grad_weight={args.psi_grad_weight}, "
        f"theta_psi_coupling={args.theta_psi_coupling}, "
        f"refractive_strength={args.refractive_strength}, "
        f"refractive_j_star={args.refractive_j_star}, "
        f"kappa_refr={args.kappa_refr}\n"
    )
    lines.append(f"- Targets file: {args.targets_json}\n")
    lines.append("\n## Best point\n")
    for k in (
        "theta_c",
        "m2",
        "lam",
        "psi_grad_weight",
        "theta_psi_coupling",
        "refractive_strength",
        "refractive_j_star",
        "kappa_refr",
        "z_surface",
        "delay_proxy",
        "ring_m87_uas",
        "ring_sgra_uas",
        "chi2_forward",
    ):
        lines.append(f"- {k}: {best[k]:.6g}\n")

    out_md = Path(args.output_md)
    out_md.parent.mkdir(parents=True, exist_ok=True)
    out_md.write_text("".join(lines), encoding="utf-8")
    print(f"Wrote forward scan CSV to {args.output_csv}")
    print(f"Wrote forward scan report to {args.output_md}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
