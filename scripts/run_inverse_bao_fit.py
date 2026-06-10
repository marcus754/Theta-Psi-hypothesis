# -*- coding: utf-8 -*-
"""
Inverse BAO diagnostic for the Θ–Ψ hypothesis.

Goal:
- locate where Θ–Ψ loses to the flat LCDM baseline in BAO;
- compute the required correction factors for DM/DH/DV observables;
- summarize what this implies for the background E(z)=H(z)/H0.
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


def parse_args(argv: Sequence[str] | None = None) -> argparse.Namespace:
    p = argparse.ArgumentParser(description="Inverse BAO diagnostic for Θ–Ψ")
    p.add_argument("--head2head_json", type=str, default="results/head2head_comparison.json")
    p.add_argument("--bao_mean", type=str, default="fitting/data/raw_public/desi_2024_gaussian_bao_ALL_GCcomb_mean.txt")
    p.add_argument("--bao_cov", type=str, default="fitting/data/raw_public/desi_2024_gaussian_bao_ALL_GCcomb_cov.txt")
    p.add_argument("--cc", type=str, default="fitting/data/real_dimless/cc_Ez_zlt1.csv")
    p.add_argument("--sn", type=str, default="fitting/data/real_dimless/sn_dimless.csv")
    p.add_argument("--bao_rd_mpc", type=float, default=147.09)
    p.add_argument("--bao_h0_km_s_mpc", type=float, default=70.0)
    p.add_argument("--output_csv", type=str, default="results/inverse_bao_fit.csv")
    p.add_argument("--output_json", type=str, default="results/inverse_bao_fit.json")
    p.add_argument("--output_md", type=str, default="results/inverse_bao_fit.md")
    return p.parse_args(argv)


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


def _theta_bao_predictions(
    *,
    gamma: float,
    V0: float,
    bao,
    cc,
    sn,
    alpha: float,
    omegas: Omegas,
) -> tuple[np.ndarray, dict[float, float]]:
    zmax = max(max(cc.z), max(sn.z), max(bao.z))
    nint = 1200
    zs = [zmax * i / nint for i in range(nint + 1)]
    Es = Hz_dimless(
        zs,
        gamma=gamma,
        V0=V0,
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


def _lcdm_bao_predictions(
    *,
    omega_m: float,
    bao,
    cc,
    sn,
    alpha: float,
    omegas: Omegas,
) -> np.ndarray:
    zmax = max(max(cc.z), max(sn.z), max(bao.z))
    nint = 1200
    zs = [zmax * i / nint for i in range(nint + 1)]
    Es = [math.sqrt(max(float(omegas.Omega_r) * (1.0 + z) ** 4 + float(omega_m) * (1.0 + z) ** 3 + float(omegas.Omega_L), 1e-18)) for z in zs]
    invE = [1.0 / max(float(e), 1e-18) for e in Es]
    dc = [0.0] * len(zs)
    s = 0.0
    for i in range(1, len(zs)):
        dz = zs[i] - zs[i - 1]
        s += 0.5 * (invE[i] + invE[i - 1]) * dz
        dc[i] = s

    out = []
    for z, obs in zip(bao.z, bao.obs):
        Ez = _interp(zs, Es, z)
        Dm = _interp(zs, dc, z)
        if obs == "DM_over_rs":
            out.append(alpha * Dm)
        elif obs == "DH_over_rs":
            out.append(alpha / max(Ez, 1e-18))
        elif obs == "DV_over_rs":
            out.append(alpha * ((z / max(Ez, 1e-18)) * (Dm * Dm)) ** (1.0 / 3.0))
        else:
            raise ValueError(f"Unsupported BAO observable: {obs}")
    return np.asarray(out, dtype=float)


def main(argv: Sequence[str] | None = None) -> int:
    args = parse_args(argv)
    h2h = json.loads(Path(args.head2head_json).read_text(encoding="utf-8"))
    bao = load_bao_vector_from_desi(args.bao_mean, args.bao_cov)
    cc = load_cc_dimless(args.cc)
    sn = load_sn_dimless(args.sn)
    alpha = 299792.458 / (float(args.bao_h0_km_s_mpc) * float(args.bao_rd_mpc))
    omegas = Omegas(H0=1.0, Omega_r=8.4e-5, Omega_b=0.049, Omega_c=0.251, Omega_L=0.7)

    theta_pred, theta_E = _theta_bao_predictions(
        gamma=float(h2h["theta_psi"]["gamma"]),
        V0=float(h2h["theta_psi"]["V0"]),
        bao=bao,
        cc=cc,
        sn=sn,
        alpha=alpha,
        omegas=omegas,
    )
    lcdm_pred = _lcdm_bao_predictions(
        omega_m=float(h2h["lcdm_gr"]["omega_m"]),
        bao=bao,
        cc=cc,
        sn=sn,
        alpha=alpha,
        omegas=omegas,
    )

    sig_diag = np.sqrt(np.diag(np.linalg.inv(bao.cov_inv)))
    rows: list[dict] = []
    obs_groups: dict[str, list[float]] = {"DV_over_rs": [], "DM_over_rs": [], "DH_over_rs": []}
    for i, (z, obs, data_y, sig, th, lc) in enumerate(zip(bao.z, bao.obs, bao.y, sig_diag, theta_pred, lcdm_pred)):
        th_pull = (float(th) - float(data_y)) / max(float(sig), 1e-18)
        lc_pull = (float(lc) - float(data_y)) / max(float(sig), 1e-18)
        required_ratio = float(data_y) / max(float(th), 1e-18)
        frac_error = (float(th) - float(data_y)) / max(abs(float(data_y)), 1e-18)
        row = {
            "index": i,
            "z": float(z),
            "observable": str(obs),
            "data": float(data_y),
            "theta_pred": float(th),
            "lcdm_pred": float(lc),
            "theta_pull_diag": float(th_pull),
            "lcdm_pull_diag": float(lc_pull),
            "frac_error_theta": float(frac_error),
            "required_multiplier": float(required_ratio),
            "required_frac_shift": float(required_ratio - 1.0),
        }
        if str(obs) == "DH_over_rs":
            # DH ∝ 1/E, so to match data we need E_target/E_theta = theta_pred/data
            row["required_E_multiplier"] = float(th) / max(float(data_y), 1e-18)
            row["theta_Ez"] = float(theta_E[float(z)])
        rows.append(row)
        obs_groups[str(obs)].append(float(frac_error))

    rth = theta_pred - bao.y
    rlc = lcdm_pred - bao.y
    chi2_theta = float(rth @ bao.cov_inv @ rth)
    chi2_lcdm = float(rlc @ bao.cov_inv @ rlc)
    delta = chi2_theta - chi2_lcdm

    summary = {
        "chi2_theta": chi2_theta,
        "chi2_lcdm": chi2_lcdm,
        "delta_chi2_bao": delta,
        "mean_frac_error": {
            key: (sum(vals) / len(vals) if vals else None) for key, vals in obs_groups.items()
        },
        "head2head_total_delta_chi2": float(h2h["delta_theta_minus_lcdm"]["chi2"]),
        "needed_total_improvement_to_parity": float(h2h["delta_theta_minus_lcdm"]["chi2"]),
        "needed_bao_improvement_to_parity": float(delta),
        "interpretation": {
            "DM_over_rs": "positive required_multiplier means D_M must increase when current prediction is low",
            "DH_over_rs": "required_E_multiplier > 1 means H(z) must increase locally; < 1 means H(z) must decrease",
            "DV_over_rs": "tracks combined transverse/radial correction and is mostly low by about two percent",
        },
    }

    out_csv = Path(args.output_csv)
    out_csv.parent.mkdir(parents=True, exist_ok=True)
    with out_csv.open("w", newline="", encoding="utf-8") as fh:
        fields = sorted({k for row in rows for k in row.keys()}) if rows else []
        w = csv.DictWriter(fh, fieldnames=fields)
        w.writeheader()
        w.writerows(rows)

    out_json = Path(args.output_json)
    out_json.write_text(json.dumps({"summary": summary, "rows": rows}, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")

    worst = sorted(rows, key=lambda r: abs(float(r["theta_pull_diag"])), reverse=True)[:5]
    lines = ["# Inverse BAO Fit\n\n"]
    lines.append(f"- BAO chi2 Θ–Ψ: {chi2_theta:.6f}\n")
    lines.append(f"- BAO chi2 flat ΛCDM: {chi2_lcdm:.6f}\n")
    lines.append(f"- Δchi2_BAO (Θ–Ψ - ΛCDM): {delta:.6f}\n")
    lines.append(f"- Needed total improvement to head-to-head parity: {float(h2h['delta_theta_minus_lcdm']['chi2']):.6f}\n")
    lines.append("\n## Mean Fractional Error By Observable\n")
    for key, val in summary["mean_frac_error"].items():
        lines.append(f"- {key}: {val:.6%}\n")
    lines.append("\n## Worst BAO Points\n")
    for r in worst:
        line = (
            f"- z={r['z']:.3f} {r['observable']}: "
            f"data={r['data']:.4f}, theta={r['theta_pred']:.4f}, lcdm={r['lcdm_pred']:.4f}, "
            f"theta_pull={r['theta_pull_diag']:.3f}, required_shift={r['required_frac_shift']:.3%}"
        )
        if "required_E_multiplier" in r:
            line += f", required_E_multiplier={r['required_E_multiplier']:.3f}"
        lines.append(line + "\n")
    lines.append("\n## Reading\n")
    lines.append("- `DM_over_rs` and `DV_over_rs` are mostly low, so the background needs larger transverse/integrated distances in the BAO redshift range.\n")
    lines.append("- `DH_over_rs` changes sign across redshift, so the problem is not a single overall rescaling of H(z), but the shape of E(z).\n")
    lines.append("- This is an inverse target for the background sector: adjust Θ–Ψ so BAO distances rise by about two percent while preserving SN and CC.\n")

    out_md = Path(args.output_md)
    out_md.write_text("".join(lines), encoding="utf-8")
    print(f"Wrote inverse BAO CSV to {args.output_csv}")
    print(f"Wrote inverse BAO JSON to {args.output_json}")
    print(f"Wrote inverse BAO report to {args.output_md}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
