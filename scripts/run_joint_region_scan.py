# -*- coding: utf-8 -*-
"""
Joint region scan for gamma,V0 with CC+SN+BAO+compact_star+compact constraints.
"""
from __future__ import annotations

import argparse
import csv
from pathlib import Path
from typing import Sequence

from fitting.core_api import Omegas
from fitting.data_io import (
    load_cc_dimless,
    load_sn_dimless,
    load_bao_dimless,
    load_compact_targets_json,
    load_compact_star_targets_json,
    load_forward_targets_json,
    load_forward_event_catalog_json,
)
from fitting.likelihoods_compact import CompactTargets
from fitting.likelihoods_compact_star import CompactStarTargets
from fitting.likelihoods_forward import ForwardTargets
from fitting.run_mcmc import Params, loglikes_components


def _linspace(a: float, b: float, n: int) -> list[float]:
    if n <= 1:
        return [a]
    h = (b - a) / (n - 1)
    return [a + i * h for i in range(n)]


def parse_args(argv: Sequence[str] | None = None) -> argparse.Namespace:
    p = argparse.ArgumentParser(description="Joint parameter-region scan")
    p.add_argument("--gmin", type=float, default=0.1)
    p.add_argument("--gmax", type=float, default=10.0)
    p.add_argument("--ng", type=int, default=30)
    p.add_argument("--vmin", type=float, default=0.01)
    p.add_argument("--vmax", type=float, default=1.0)
    p.add_argument("--nv", type=int, default=30)
    p.add_argument("--q", type=float, default=0.0)
    p.add_argument("--compact_targets_json", type=str, default="fitting/data/compact_constraints_observational_proxy.json")
    p.add_argument("--compact_star_targets_json", type=str, default="fitting/data/compact_star_observational_proxy.json")
    p.add_argument("--forward_targets_json", type=str, default="fitting/data/forward_observational_proxy.json")
    p.add_argument("--compact_weight", type=float, default=1.0)
    p.add_argument("--compact_star_weight", type=float, default=1.0)
    p.add_argument("--forward_weight", type=float, default=1.0)
    p.add_argument("--forward_event_catalog_json", type=str, default="fitting/data/forward_event_catalog.json")
    p.add_argument("--forward_events_weight", type=float, default=1.0)
    p.add_argument("--forward_k_ring_m87", type=float, default=3.0)
    p.add_argument("--forward_k_ring_sgra", type=float, default=3.6)
    p.add_argument("--forward_k_echo_delay", type=float, default=1.0)
    p.add_argument("--output_csv", type=str, default="results/joint_region_scan.csv")
    p.add_argument("--output_md", type=str, default="results/joint_region_scan.md")
    return p.parse_args(argv)


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

    rows = []
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
                compact_weight=args.compact_weight,
                compact_star=compact_star,
                compact_star_weight=args.compact_star_weight,
                compact_star_link_model=True,
                compact_star_link_background=True,
                compact_star_theta_c_from_gamma=1.5,
                compact_star_m2_from_v0=10.0,
                compact_star_lam_from_q2=0.1,
                forward=forward,
                forward_weight=args.forward_weight,
                forward_k_ring_m87=args.forward_k_ring_m87,
                forward_k_ring_sgra=args.forward_k_ring_sgra,
                forward_k_echo_delay=args.forward_k_echo_delay,
                forward_events=forward_events,
                forward_events_weight=args.forward_events_weight,
            )
            loglike_total = sum(float(x) for x in comps.values())
            chi2_total = -2.0 * loglike_total
            row = {
                "gamma": g,
                "V0": v,
                "Q": args.q,
                "chi2_total": chi2_total,
            }
            for k, ll in comps.items():
                row[f"chi2_{k}"] = -2.0 * float(ll)
            rows.append(row)

    best = min(rows, key=lambda r: float(r["chi2_total"]))
    chi2_best = float(best["chi2_total"])
    for r in rows:
        dchi2 = float(r["chi2_total"]) - chi2_best
        r["dchi2"] = dchi2
        r["allowed_1sigma"] = int(dchi2 <= 2.30)
        r["allowed_2sigma"] = int(dchi2 <= 6.18)
        r["allowed_3sigma"] = int(dchi2 <= 11.83)

    out_csv = Path(args.output_csv)
    out_csv.parent.mkdir(parents=True, exist_ok=True)
    fields = sorted({k for row in rows for k in row.keys()})
    with out_csv.open("w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=fields)
        w.writeheader()
        w.writerows(rows)

    n = len(rows)
    n1 = sum(int(r["allowed_1sigma"]) for r in rows)
    n2 = sum(int(r["allowed_2sigma"]) for r in rows)
    n3 = sum(int(r["allowed_3sigma"]) for r in rows)

    out_md = Path(args.output_md)
    lines = []
    lines.append("# Joint Region Scan\n\n")
    lines.append("- Components: CC + SN + BAO + compact_star(linked) + compact(constraints) + forward(EHT/GW)\n")
    lines.append(f"- Grid: gamma in [{args.gmin}, {args.gmax}] x{args.ng}, V0 in [{args.vmin}, {args.vmax}] x{args.nv}\n")
    lines.append(f"- Best point: gamma={best['gamma']:.6g}, V0={best['V0']:.6g}, chi2={chi2_best:.6g}\n")
    lines.append(f"- Allowed 1σ: {n1}/{n} ({100.0*n1/n:.2f}%)\n")
    lines.append(f"- Allowed 2σ: {n2}/{n} ({100.0*n2/n:.2f}%)\n")
    lines.append(f"- Allowed 3σ: {n3}/{n} ({100.0*n3/n:.2f}%)\n")
    out_md.write_text("".join(lines), encoding="utf-8")

    print(f"Wrote joint scan CSV to {args.output_csv}")
    print(f"Wrote joint scan report to {args.output_md}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
