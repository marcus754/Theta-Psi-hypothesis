# -*- coding: utf-8 -*-
"""
BAO-only inverse target scan for the Θ–Ψ background.

This scan answers a concrete question:
where would BAO alone place the (gamma, V0) background, and how far is that
from the current head-to-head / joint best-fit regions?
"""
from __future__ import annotations

import argparse
import csv
import json
import math
from pathlib import Path
from typing import Sequence

import numpy as np

from fitting.core_api import Omegas
from fitting.data_io import load_cc_dimless, load_sn_dimless
from fitting.head2head import load_bao_vector_from_desi
from fitting.model import Hz_dimless


def _linspace(a: float, b: float, n: int) -> list[float]:
    if n <= 1:
        return [float(a)]
    h = (float(b) - float(a)) / float(n - 1)
    return [float(a) + i * h for i in range(n)]


def _interp(x: list[float], y: list[float], xq: float) -> float:
    import bisect

    if xq <= x[0]:
        return float(y[0])
    if xq >= x[-1]:
        return float(y[-1])
    j = bisect.bisect_left(x, xq)
    x0, x1 = x[j - 1], x[j]
    y0, y1 = y[j - 1], y[j]
    t = (xq - x0) / max(x1 - x0, 1e-18)
    return float(y0 + t * (y1 - y0))


def _theta_bao_vector(gamma: float, V0: float, *, bao, cc, sn, alpha: float, omegas: Omegas) -> tuple[np.ndarray, dict[float, float]]:
    zmax = max(max(cc.z), max(sn.z), max(bao.z))
    nint = 1200
    zs = [zmax * i / nint for i in range(nint + 1)]
    Es = Hz_dimless(
        zs,
        gamma=float(gamma),
        V0=float(V0),
        Q=0.0,
        omegas=omegas,
        a_min=1e-3,
        dN=0.02,
    )
    invE = [1.0 / max(float(e), 1e-18) for e in Es]
    dc = [0.0] * len(zs)
    s = 0.0
    for i in range(1, len(zs)):
        dz = zs[i] - zs[i - 1]
        s += 0.5 * (invE[i] + invE[i - 1]) * dz
        dc[i] = s
    out = []
    e_by_z: dict[float, float] = {}
    for z, obs in zip(bao.z, bao.obs):
        Ez = _interp(zs, Es, z)
        Dm = _interp(zs, dc, z)
        e_by_z[float(z)] = float(Ez)
        if obs == "DM_over_rs":
            out.append(alpha * Dm)
        elif obs == "DH_over_rs":
            out.append(alpha / max(Ez, 1e-18))
        elif obs == "DV_over_rs":
            out.append(alpha * ((z / max(Ez, 1e-18)) * (Dm * Dm)) ** (1.0 / 3.0))
        else:
            raise ValueError(f"Unsupported BAO observable: {obs}")
    return np.asarray(out, dtype=float), e_by_z


def parse_args(argv: Sequence[str] | None = None) -> argparse.Namespace:
    p = argparse.ArgumentParser(description="Inverse BAO target scan for Θ–Ψ")
    p.add_argument("--head2head_json", type=str, default="results/head2head_comparison.json")
    p.add_argument("--bao_mean", type=str, default="fitting/data/raw_public/desi_2024_gaussian_bao_ALL_GCcomb_mean.txt")
    p.add_argument("--bao_cov", type=str, default="fitting/data/raw_public/desi_2024_gaussian_bao_ALL_GCcomb_cov.txt")
    p.add_argument("--cc", type=str, default="fitting/data/real_dimless/cc_Ez_zlt1.csv")
    p.add_argument("--sn", type=str, default="fitting/data/real_dimless/sn_dimless.csv")
    p.add_argument("--bao_rd_mpc", type=float, default=147.09)
    p.add_argument("--bao_h0_km_s_mpc", type=float, default=70.0)
    p.add_argument("--gmin", type=float, default=0.001)
    p.add_argument("--gmax", type=float, default=0.05)
    p.add_argument("--ng", type=int, default=24)
    p.add_argument("--vmin", type=float, default=0.001)
    p.add_argument("--vmax", type=float, default=0.03)
    p.add_argument("--nv", type=int, default=24)
    p.add_argument("--output_csv", type=str, default="results/inverse_bao_target_scan.csv")
    p.add_argument("--output_json", type=str, default="results/inverse_bao_target_scan.json")
    p.add_argument("--output_md", type=str, default="results/inverse_bao_target_scan.md")
    return p.parse_args(argv)


def main(argv: Sequence[str] | None = None) -> int:
    args = parse_args(argv)
    h2h = json.loads(Path(args.head2head_json).read_text(encoding="utf-8"))
    bao = load_bao_vector_from_desi(args.bao_mean, args.bao_cov)
    cc = load_cc_dimless(args.cc)
    sn = load_sn_dimless(args.sn)
    alpha = 299792.458 / (float(args.bao_h0_km_s_mpc) * float(args.bao_rd_mpc))
    omegas = Omegas(H0=1.0, Omega_r=8.4e-5, Omega_b=0.049, Omega_c=0.251, Omega_L=0.7)

    rows: list[dict] = []
    best = None
    for gamma in _linspace(args.gmin, args.gmax, args.ng):
        for V0 in _linspace(args.vmin, args.vmax, args.nv):
            pred, e_by_z = _theta_bao_vector(gamma, V0, bao=bao, cc=cc, sn=sn, alpha=alpha, omegas=omegas)
            resid = pred - bao.y
            chi2 = float(resid @ bao.cov_inv @ resid)
            row = {"gamma": float(gamma), "V0": float(V0), "chi2_bao": chi2}
            # Store a few shape diagnostics.
            for z in sorted(set(float(z) for z in bao.z)):
                row[f"E_z{z:.3f}"] = float(e_by_z[z])
            rows.append(row)
            if best is None or chi2 < best["chi2_bao"]:
                best = dict(row)

    assert best is not None
    theta_gamma = float(h2h["theta_psi"]["gamma"])
    theta_V0 = float(h2h["theta_psi"]["V0"])
    theta_bao = float(h2h["theta_psi"]["chi2_bao"])
    delta_to_theta = theta_bao - float(best["chi2_bao"])

    # nearest grid point to current head2head theta solution
    nearest = min(rows, key=lambda r: (float(r["gamma"]) - theta_gamma) ** 2 + (float(r["V0"]) - theta_V0) ** 2)

    out_csv = Path(args.output_csv)
    out_csv.parent.mkdir(parents=True, exist_ok=True)
    fields = sorted({k for row in rows for k in row.keys()})
    with out_csv.open("w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=fields)
        w.writeheader()
        w.writerows(rows)

    out = {
        "best_bao_only": best,
        "theta_head2head": {
            "gamma": theta_gamma,
            "V0": theta_V0,
            "chi2_bao": theta_bao,
        },
        "nearest_grid_to_theta_head2head": nearest,
        "improvement_vs_theta_head2head": delta_to_theta,
        "scan_window": {
            "gamma": [float(args.gmin), float(args.gmax), int(args.ng)],
            "V0": [float(args.vmin), float(args.vmax), int(args.nv)],
        },
    }
    out_json = Path(args.output_json)
    out_json.write_text(json.dumps(out, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")

    lines = ["# Inverse BAO Target Scan\n\n"]
    lines.append(f"- scan gamma: [{args.gmin}, {args.gmax}] x{args.ng}\n")
    lines.append(f"- scan V0: [{args.vmin}, {args.vmax}] x{args.nv}\n")
    lines.append(
        f"- BAO-only best: gamma={best['gamma']:.8g}, V0={best['V0']:.8g}, chi2_bao={best['chi2_bao']:.6f}\n"
    )
    lines.append(
        f"- current head-to-head Θ–Ψ: gamma={theta_gamma:.8g}, V0={theta_V0:.8g}, chi2_bao={theta_bao:.6f}\n"
    )
    lines.append(f"- BAO-only improvement over current Θ–Ψ point: {delta_to_theta:.6f}\n")
    lines.append("\n## Reading\n")
    lines.append("- If BAO-only best lies far from the current Θ–Ψ head-to-head point, the background tension is structural, not just under-optimized.\n")
    lines.append("- If BAO-only best is close in parameter space but much better in chi2, the current global fit is likely paying for SN/CC compatibility.\n")
    lines.append("- Compare `E_z...` columns in the CSV to see which parts of E(z) the BAO-only optimum shifts most strongly.\n")
    out_md = Path(args.output_md)
    out_md.write_text("".join(lines), encoding="utf-8")

    print(f"Wrote inverse BAO target scan CSV to {args.output_csv}")
    print(f"Wrote inverse BAO target scan JSON to {args.output_json}")
    print(f"Wrote inverse BAO target scan report to {args.output_md}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
