# -*- coding: utf-8 -*-
"""
Stability scan for head-to-head Theta-Psi vs flat LCDM/GR.

Runs scripts.run_head2head_comparison over a grid of:
- strong_field_weight
- bridge_forward_calibration_rounds
- bridge_forward_calibration_reg

Collects delta chi2/AIC/BIC and component chi2 contributions.
"""
from __future__ import annotations

import argparse
import csv
import json
import subprocess
import sys
from pathlib import Path
from typing import Sequence


def _parse_grid_floats(raw: str) -> list[float]:
    out: list[float] = []
    for t in str(raw).split(","):
        s = t.strip()
        if not s:
            continue
        out.append(float(s))
    if not out:
        raise ValueError("empty float grid")
    return out


def _parse_grid_ints(raw: str) -> list[int]:
    out: list[int] = []
    for t in str(raw).split(","):
        s = t.strip()
        if not s:
            continue
        out.append(int(s))
    if not out:
        raise ValueError("empty int grid")
    return out


def parse_args(argv: Sequence[str] | None = None) -> argparse.Namespace:
    p = argparse.ArgumentParser(description="Stability scan for head2head comparison")
    p.add_argument("--cc", type=str, default="fitting/data/real_dimless/cc_Ez_zlt1.csv")
    p.add_argument("--sn", type=str, default="fitting/data/real_dimless/sn_dimless.csv")
    p.add_argument("--bao_mean", type=str, default="fitting/data/raw_public/desi_2024_gaussian_bao_ALL_GCcomb_mean.txt")
    p.add_argument("--bao_cov", type=str, default="fitting/data/raw_public/desi_2024_gaussian_bao_ALL_GCcomb_cov.txt")
    p.add_argument("--theta_samples", type=int, default=260)
    p.add_argument("--strong_field_weights", type=str, default="0.5,1.0,2.0")
    p.add_argument("--bridge_rounds", type=str, default="0,1,2")
    p.add_argument("--bridge_regs", type=str, default="0.03,0.1,0.3")
    p.add_argument("--workdir", type=str, default="results/head2head_stability_runs")
    p.add_argument("--output_csv", type=str, default="results/head2head_stability_scan.csv")
    p.add_argument("--output_md", type=str, default="results/head2head_stability_scan.md")
    return p.parse_args(argv)


def _run_one(
    *,
    cc: str,
    sn: str,
    bao_mean: str,
    bao_cov: str,
    theta_samples: int,
    strong_field_weight: float,
    bridge_rounds: int,
    bridge_reg: float,
    out_json: Path,
    out_md: Path,
) -> dict:
    cmd = [
        sys.executable,
        "-m",
        "scripts.run_head2head_comparison",
        "--cc",
        cc,
        "--sn",
        sn,
        "--bao_mean",
        bao_mean,
        "--bao_cov",
        bao_cov,
        "--include_strong_field",
        "--strong_field_weight",
        str(float(strong_field_weight)),
        "--bridge_forward_calibration_rounds",
        str(int(bridge_rounds)),
        "--bridge_forward_calibration_reg",
        str(float(bridge_reg)),
        "--penalize_bridge_latent",
        "--theta_samples",
        str(int(theta_samples)),
        "--output_md",
        str(out_md),
        "--output_json",
        str(out_json),
    ]
    proc = subprocess.run(cmd, check=False, capture_output=True, text=True)
    if proc.returncode != 0:
        raise RuntimeError(
            f"head2head failed rc={proc.returncode}\nSTDOUT:\n{proc.stdout}\nSTDERR:\n{proc.stderr}"
        )
    data = json.loads(out_json.read_text(encoding="utf-8"))
    delta = data["delta_theta_minus_lcdm"]
    th = data["theta_psi"]
    lc = data["lcdm_gr"]
    row = {
        "strong_field_weight": float(strong_field_weight),
        "bridge_rounds": int(bridge_rounds),
        "bridge_reg": float(bridge_reg),
        "delta_chi2": float(delta["chi2"]),
        "delta_AIC": float(delta["AIC"]),
        "delta_BIC": float(delta["BIC"]),
        "theta_chi2_total": float(th["chi2_total"]),
        "lcdm_chi2_total": float(lc["chi2_total"]),
        "theta_chi2_cc": float(th["chi2_cc"]),
        "theta_chi2_sn": float(th["chi2_sn"]),
        "theta_chi2_bao": float(th["chi2_bao"]),
        "theta_chi2_sf": float(th["chi2_sf"]),
        "lcdm_chi2_cc": float(lc["chi2_cc"]),
        "lcdm_chi2_sn": float(lc["chi2_sn"]),
        "lcdm_chi2_bao": float(lc["chi2_bao"]),
        "lcdm_chi2_sf": float(lc["chi2_sf"]),
        "theta_gamma": float(th["gamma"]),
        "theta_V0": float(th["V0"]),
    }
    return row


def main(argv: Sequence[str] | None = None) -> int:
    args = parse_args(argv)
    weights = _parse_grid_floats(args.strong_field_weights)
    rounds_grid = _parse_grid_ints(args.bridge_rounds)
    regs = _parse_grid_floats(args.bridge_regs)

    workdir = Path(args.workdir)
    workdir.mkdir(parents=True, exist_ok=True)

    rows: list[dict] = []
    for w in weights:
        for br in rounds_grid:
            for reg in regs:
                tag = f"w{w:g}_br{br:d}_reg{reg:g}"
                out_json = workdir / f"{tag}.json"
                out_md = workdir / f"{tag}.md"
                row = _run_one(
                    cc=args.cc,
                    sn=args.sn,
                    bao_mean=args.bao_mean,
                    bao_cov=args.bao_cov,
                    theta_samples=int(args.theta_samples),
                    strong_field_weight=float(w),
                    bridge_rounds=int(br),
                    bridge_reg=float(reg),
                    out_json=out_json,
                    out_md=out_md,
                )
                rows.append(row)
                print(
                    f"[scan] w={w:g} br={br} reg={reg:g} "
                    f"delta_chi2={row['delta_chi2']:.3f} delta_BIC={row['delta_BIC']:.3f}"
                )

    rows.sort(key=lambda r: (r["delta_BIC"], r["delta_AIC"], r["delta_chi2"]))
    out_csv = Path(args.output_csv)
    out_csv.parent.mkdir(parents=True, exist_ok=True)
    fields = list(rows[0].keys()) if rows else []
    with out_csv.open("w", newline="", encoding="utf-8") as fh:
        wr = csv.DictWriter(fh, fieldnames=fields)
        wr.writeheader()
        wr.writerows(rows)

    best_bic = rows[0]
    best_chi2 = min(rows, key=lambda r: r["delta_chi2"])
    no_bridge = [r for r in rows if int(r["bridge_rounds"]) == 0]
    with_bridge = [r for r in rows if int(r["bridge_rounds"]) > 0]

    lines: list[str] = []
    lines.append("# Head2Head Stability Scan\n\n")
    lines.append(f"- runs: {len(rows)}\n")
    lines.append(f"- strong_field_weights: {weights}\n")
    lines.append(f"- bridge_rounds: {rounds_grid}\n")
    lines.append(f"- bridge_regs: {regs}\n")
    lines.append(f"- theta_samples(each): {int(args.theta_samples)}\n\n")
    lines.append("## Best by BIC\n")
    lines.append(f"- delta_chi2: {best_bic['delta_chi2']:.6f}\n")
    lines.append(f"- delta_AIC: {best_bic['delta_AIC']:.6f}\n")
    lines.append(f"- delta_BIC: {best_bic['delta_BIC']:.6f}\n")
    lines.append(
        f"- setup: w={best_bic['strong_field_weight']:.6g}, "
        f"bridge_rounds={int(best_bic['bridge_rounds'])}, "
        f"bridge_reg={best_bic['bridge_reg']:.6g}\n\n"
    )
    lines.append("## Best by chi2\n")
    lines.append(f"- delta_chi2: {best_chi2['delta_chi2']:.6f}\n")
    lines.append(f"- delta_AIC: {best_chi2['delta_AIC']:.6f}\n")
    lines.append(f"- delta_BIC: {best_chi2['delta_BIC']:.6f}\n")
    lines.append(
        f"- setup: w={best_chi2['strong_field_weight']:.6g}, "
        f"bridge_rounds={int(best_chi2['bridge_rounds'])}, "
        f"bridge_reg={best_chi2['bridge_reg']:.6g}\n\n"
    )

    if no_bridge and with_bridge:
        mean_nb = sum(r["delta_chi2"] for r in no_bridge) / len(no_bridge)
        mean_wb = sum(r["delta_chi2"] for r in with_bridge) / len(with_bridge)
        lines.append("## Bridge impact\n")
        lines.append(f"- avg delta_chi2 (bridge_rounds=0): {mean_nb:.6f}\n")
        lines.append(f"- avg delta_chi2 (bridge_rounds>0): {mean_wb:.6f}\n")
        lines.append(f"- improvement from bridge fit: {mean_wb - mean_nb:.6f}\n\n")

    lines.append("## Top-5 (by BIC)\n")
    for i, r in enumerate(rows[:5], start=1):
        lines.append(
            f"- {i}. w={r['strong_field_weight']:.6g}, br={int(r['bridge_rounds'])}, "
            f"reg={r['bridge_reg']:.6g}, dchi2={r['delta_chi2']:.3f}, "
            f"dAIC={r['delta_AIC']:.3f}, dBIC={r['delta_BIC']:.3f}\n"
        )

    out_md = Path(args.output_md)
    out_md.parent.mkdir(parents=True, exist_ok=True)
    out_md.write_text("".join(lines), encoding="utf-8")
    print(f"Wrote {out_csv}")
    print(f"Wrote {out_md}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
