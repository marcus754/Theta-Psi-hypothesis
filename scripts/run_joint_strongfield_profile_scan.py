# -*- coding: utf-8 -*-
"""
Profile scan of new strong-field parameters atop the joint-fit contour.
"""
from __future__ import annotations

import argparse
import csv
from pathlib import Path
from typing import Sequence

from fitting.core_api import Omegas
from fitting.data_io import (
    load_bao_dimless,
    load_cc_dimless,
    load_compact_star_targets_json,
    load_compact_targets_json,
    load_forward_event_catalog_json,
    load_forward_targets_json,
    load_sn_dimless,
)
from fitting.likelihoods_compact import CompactTargets
from fitting.likelihoods_compact_star import CompactStarTargets
from fitting.likelihoods_forward import ForwardTargets
from fitting.run_mcmc import Params, loglikes_components


def _linspace(a: float, b: float, n: int) -> list[float]:
    if n <= 1:
        return [float(a)]
    h = (float(b) - float(a)) / float(n - 1)
    return [float(a) + i * h for i in range(n)]


def _to_compact_targets(d: dict) -> CompactTargets:
    return CompactTargets(
        mode=str(d.get("mode", "gaussian")),
        z_surface=d.get("z_surface"),
        sigma_z_surface=d.get("sigma_z_surface"),
        delay_proxy=d.get("delay_proxy"),
        sigma_delay_proxy=d.get("sigma_delay_proxy"),
        n_peak=d.get("n_peak"),
        sigma_n_peak=d.get("sigma_n_peak"),
        z_surface_min=d.get("z_surface_min"),
        sigma_z_surface_min=d.get("sigma_z_surface_min"),
        delay_proxy_min=d.get("delay_proxy_min"),
        sigma_delay_proxy_min=d.get("sigma_delay_proxy_min"),
        n_peak_max=d.get("n_peak_max"),
        sigma_n_peak_max=d.get("sigma_n_peak_max"),
    )


def parse_args(argv: Sequence[str] | None = None) -> argparse.Namespace:
    p = argparse.ArgumentParser(description="Profile scan of strong-field parameters over joint fit")
    p.add_argument("--gmin", type=float, default=0.1)
    p.add_argument("--gmax", type=float, default=10.0)
    p.add_argument("--ng", type=int, default=10)
    p.add_argument("--vmin", type=float, default=0.01)
    p.add_argument("--vmax", type=float, default=1.0)
    p.add_argument("--nv", type=int, default=10)
    p.add_argument("--q", type=float, default=0.0)
    p.add_argument("--n_profile", type=int, default=5)
    p.add_argument("--compact_targets_json", type=str, default="fitting/data/compact_constraints_observational_proxy.json")
    p.add_argument("--compact_star_targets_json", type=str, default="fitting/data/compact_star_observational_proxy.json")
    p.add_argument("--forward_targets_json", type=str, default="fitting/data/forward_observational_proxy.json")
    p.add_argument("--forward_event_catalog_json", type=str, default="fitting/data/forward_event_catalog.json")
    p.add_argument("--output_csv", type=str, default="results/joint_strongfield_profile_scan.csv")
    p.add_argument("--output_md", type=str, default="results/joint_strongfield_profile_scan.md")
    return p.parse_args(argv)


def main(argv: Sequence[str] | None = None) -> int:
    args = parse_args(argv)
    om = Omegas(H0=1.0)
    cc = load_cc_dimless("fitting/data/cc_demo_dimless.csv")
    sn = load_sn_dimless("fitting/data/sn_demo_dimless.csv")
    bao = load_bao_dimless("fitting/data/bao_demo_dimless.csv")
    compact = _to_compact_targets(load_compact_targets_json(args.compact_targets_json))
    csj = load_compact_star_targets_json(args.compact_star_targets_json)
    compact_star = CompactStarTargets(
        z_surface=csj["z_surface"],
        sigma_z_surface=csj["sigma_z_surface"],
        delay_proxy=csj["delay_proxy"],
        sigma_delay_proxy=csj["sigma_delay_proxy"],
    )
    fj = load_forward_targets_json(args.forward_targets_json)
    forward = ForwardTargets(
        ring_m87_uas=fj["ring_m87_uas"],
        sigma_ring_m87_uas=fj["sigma_ring_m87_uas"],
        ring_sgra_uas=fj["ring_sgra_uas"],
        sigma_ring_sgra_uas=fj["sigma_ring_sgra_uas"],
        echo_delay_ms=fj.get("echo_delay_ms"),
        sigma_echo_delay_ms=fj.get("sigma_echo_delay_ms"),
    )
    forward_events = load_forward_event_catalog_json(args.forward_event_catalog_json)

    base_values = {
        "compact_star_psi_grad_weight": 1.0,
        "compact_star_theta_psi_coupling": 1.0,
        "compact_star_refractive_strength": 1.0,
        "compact_star_refractive_j_star": 1.0,
        "compact_star_kappa_refr": 1.0,
    }
    scan_ranges = {
        "compact_star_psi_grad_weight": (0.25, 2.0),
        "compact_star_theta_psi_coupling": (0.25, 2.0),
        "compact_star_refractive_strength": (0.25, 2.0),
        "compact_star_refractive_j_star": (0.25, 2.0),
        "compact_star_kappa_refr": (0.25, 2.0),
    }

    rows: list[dict] = []
    for par_name, (vmin, vmax) in scan_ranges.items():
        for par_value in _linspace(vmin, vmax, args.n_profile):
            best_row = None
            best_chi2 = None
            sf_kwargs = dict(base_values)
            sf_kwargs[par_name] = float(par_value)
            for g in _linspace(args.gmin, args.gmax, args.ng):
                for v in _linspace(args.vmin, args.vmax, args.nv):
                    p = Params(gamma=g, V0=v, Q=args.q)
                    comps = loglikes_components(
                        p,
                        om=om,
                        cc=cc,
                        sn=sn,
                        bao=bao,
                        compact=compact,
                        compact_weight=1.0,
                        compact_star=compact_star,
                        compact_star_weight=1.0,
                        compact_star_link_model=True,
                        compact_star_link_background=True,
                        compact_star_theta_c_from_gamma=1.5,
                        compact_star_m2_from_v0=10.0,
                        compact_star_lam_from_q2=0.1,
                        forward=forward,
                        forward_weight=1.0,
                        forward_events=forward_events,
                        forward_events_weight=1.0,
                        **sf_kwargs,
                    )
                    chi2_total = -2.0 * sum(float(x) for x in comps.values())
                    if best_chi2 is None or chi2_total < best_chi2:
                        best_chi2 = chi2_total
                        best_row = {
                            "parameter": par_name,
                            "parameter_value": float(par_value),
                            "gamma_best": float(g),
                            "V0_best": float(v),
                            "chi2_total_best": float(chi2_total),
                        }
                        for k, ll in comps.items():
                            best_row[f"chi2_{k}"] = -2.0 * float(ll)
            assert best_row is not None
            rows.append(best_row)

    out_csv = Path(args.output_csv)
    out_csv.parent.mkdir(parents=True, exist_ok=True)
    fields = sorted({k for row in rows for k in row.keys()})
    with out_csv.open("w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=fields)
        w.writeheader()
        w.writerows(rows)

    lines = ["# Joint Strong-Field Profile Scan\n\n"]
    lines.append(f"- gamma grid: [{args.gmin}, {args.gmax}] x{args.ng}\n")
    lines.append(f"- V0 grid: [{args.vmin}, {args.vmax}] x{args.nv}\n")
    lines.append(f"- Profile samples per parameter: {args.n_profile}\n\n")
    lines.append("## Best By Parameter\n")
    for par_name in scan_ranges:
        block = [r for r in rows if r["parameter"] == par_name]
        best = min(block, key=lambda r: float(r["chi2_total_best"]))
        lines.append(
            f"- {par_name}: best value={best['parameter_value']:.6g}, "
            f"gamma={best['gamma_best']:.6g}, V0={best['V0_best']:.6g}, "
            f"chi2={best['chi2_total_best']:.6g}\n"
        )

    out_md = Path(args.output_md)
    out_md.write_text("".join(lines), encoding="utf-8")
    print(f"Wrote profile scan CSV to {args.output_csv}")
    print(f"Wrote profile scan report to {args.output_md}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
