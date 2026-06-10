# -*- coding: utf-8 -*-
"""
Train/test falsification protocol for Theta-Psi vs GR baseline.

Pipeline:
1) Split CC/SN/BAO/events into train/test by index modulo.
2) Calibrate fixed global bridge coefficients on train events.
3) Fit Theta-Psi and GR baseline on train split.
4) Evaluate both models on test split with fixed train parameters (no retune).
5) Emit pass/fail verdict by predefined falsification criteria.
"""
from __future__ import annotations

import argparse
import csv
import json
import math
import subprocess
import sys
from pathlib import Path
from typing import Sequence

import numpy as np

from fitting.core_api import Omegas
from fitting.data_io import load_cc_dimless, load_forward_event_catalog_json, load_sn_dimless
from fitting.head2head import load_bao_vector_from_desi
from scripts.run_head2head_comparison import _lcdm_total_chi2, _theta_total_chi2


def parse_args(argv: Sequence[str] | None = None) -> argparse.Namespace:
    p = argparse.ArgumentParser(description="Run train/test falsification protocol")
    p.add_argument("--cc", type=str, default="fitting/data/real_dimless/cc_Ez_zlt1.csv")
    p.add_argument("--sn", type=str, default="fitting/data/real_dimless/sn_dimless.csv")
    p.add_argument("--bao_mean", type=str, default="fitting/data/raw_public/desi_2024_gaussian_bao_ALL_GCcomb_mean.txt")
    p.add_argument("--bao_cov", type=str, default="fitting/data/raw_public/desi_2024_gaussian_bao_ALL_GCcomb_cov.txt")
    p.add_argument("--forward_event_catalog_json", type=str, default="fitting/data/forward_event_catalog.json")
    p.add_argument("--bridge_base_gamma", type=float, default=0.0020339622)
    p.add_argument("--bridge_base_v0", type=float, default=0.012364227)
    p.add_argument("--bridge_model", type=str, choices=["constant", "powerlaw2d"], default="constant")
    p.add_argument("--bridge_recenter_on_train", action=argparse.BooleanOptionalAction, default=False)
    p.add_argument("--bridge_rounds", type=int, default=36)
    p.add_argument("--bridge_reg_strength", type=float, default=0.2)
    p.add_argument("--split_mod", type=int, default=3)
    p.add_argument("--test_fold", type=int, default=1)
    p.add_argument("--comparison_mode", type=str, choices=["lcdm_fit", "same_content"], default="same_content")
    p.add_argument("--omega_r", type=float, default=8.4e-5)
    p.add_argument("--omega_b", type=float, default=0.049)
    p.add_argument("--omega_c", type=float, default=0.251)
    p.add_argument("--omega_l", type=float, default=0.7)
    p.add_argument("--theta_samples_train", type=int, default=220)
    p.add_argument("--strong_field_weight", type=float, default=1.0)
    p.add_argument("--bic_fail_threshold", type=float, default=10.0)
    p.add_argument("--c_sn_tol", type=float, default=1.0)
    p.add_argument("--c_bao_tol", type=float, default=1.0)
    p.add_argument("--c_min_sf_improvement", type=float, default=0.5)
    p.add_argument("--workdir", type=str, default="results/falsification_protocol")
    p.add_argument("--output_md", type=str, default="results/falsification_protocol/report.md")
    p.add_argument("--output_json", type=str, default="results/falsification_protocol/report.json")
    return p.parse_args(argv)


def _split_idx(n: int, mod: int, test_fold: int) -> tuple[list[int], list[int]]:
    tr, te = [], []
    for i in range(n):
        if (i % mod) == test_fold:
            te.append(i)
        else:
            tr.append(i)
    return tr, te


def _split_csv_rows(path: Path, mod: int, test_fold: int) -> tuple[list[dict], list[dict], list[str]]:
    with path.open("r", newline="", encoding="utf-8") as fh:
        rd = csv.DictReader(fh)
        rows = list(rd)
        fields = rd.fieldnames or []
    tr_idx, te_idx = _split_idx(len(rows), mod, test_fold)
    tr = [rows[i] for i in tr_idx]
    te = [rows[i] for i in te_idx]
    return tr, te, fields


def _write_csv(path: Path, rows: list[dict], fields: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as fh:
        wr = csv.DictWriter(fh, fieldnames=fields)
        wr.writeheader()
        wr.writerows(rows)


def _split_bao(mean_path: Path, cov_path: Path, mod: int, test_fold: int, outdir: Path) -> dict:
    recs: list[tuple[float, float, str]] = []
    with mean_path.open("r", encoding="utf-8") as fh:
        for ln in fh:
            t = ln.strip()
            if (not t) or t.startswith("#"):
                continue
            p = t.split()
            if len(p) < 3:
                continue
            recs.append((float(p[0]), float(p[1]), str(p[2])))
    cov = np.loadtxt(str(cov_path), dtype=float)
    if cov.ndim == 1:
        n = int(round(math.sqrt(float(cov.size))))
        cov = cov.reshape((n, n))
    if cov.shape[0] != len(recs):
        raise ValueError(f"BAO mean/cov mismatch: {len(recs)} vs {cov.shape}")
    tr_idx, te_idx = _split_idx(len(recs), mod, test_fold)
    out = {}
    for name, idx in (("train", tr_idx), ("test", te_idx)):
        mpath = outdir / f"bao_{name}_mean.txt"
        cpath = outdir / f"bao_{name}_cov.txt"
        mpath.parent.mkdir(parents=True, exist_ok=True)
        with mpath.open("w", encoding="utf-8") as fh:
            for i in idx:
                z, y, obs = recs[i]
                fh.write(f"{z:.12g} {y:.12g} {obs}\n")
        subcov = cov[np.ix_(idx, idx)]
        np.savetxt(cpath, subcov)
        out[name] = {"mean": str(mpath), "cov": str(cpath), "n": len(idx)}
    return out


def _split_events(path: Path, mod: int, test_fold: int, outdir: Path) -> dict:
    events = load_forward_event_catalog_json(path)
    tr_idx, te_idx = _split_idx(len(events), mod, test_fold)
    out = {}
    for name, idx in (("train", tr_idx), ("test", te_idx)):
        events_sub = [events[i] for i in idx]
        p = outdir / f"events_{name}.json"
        p.parent.mkdir(parents=True, exist_ok=True)
        p.write_text(json.dumps({"events": events_sub}, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
        out[name] = {"path": str(p), "n": len(events_sub)}
    return out


def _run(cmd: list[str]) -> None:
    proc = subprocess.run(cmd, check=False, capture_output=True, text=True)
    if proc.returncode != 0:
        raise RuntimeError(f"Command failed: {' '.join(cmd)}\nSTDOUT:\n{proc.stdout}\nSTDERR:\n{proc.stderr}")


def _run_bridge_calibration(
    *,
    bridge_model: str,
    events_train_path: str,
    gamma: float,
    V0: float,
    rounds: int,
    reg_strength: float,
    bridge_json: Path,
    bridge_md: Path,
) -> None:
    if str(bridge_model) == "powerlaw2d":
        _run(
            [
                sys.executable,
                "-m",
                "scripts.run_bridge_powerlaw_calibration_multi",
                "--forward_event_catalog_json",
                str(events_train_path),
                "--gamma",
                str(float(gamma)),
                "--V0",
                str(float(V0)),
                "--grid_gamma_factors",
                "0.7,1.0,1.4",
                "--grid_v0_factors",
                "0.7,1.0,1.4",
                "--weight_center",
                "2.0",
                "--rounds",
                str(int(rounds)),
                "--reg_strength",
                str(float(reg_strength)),
                "--gamma_ref",
                str(float(gamma)),
                "--v0_ref",
                str(float(V0)),
                "--output_json",
                str(bridge_json),
                "--output_md",
                str(bridge_md),
            ]
        )
    else:
        _run(
            [
                sys.executable,
                "-m",
                "scripts.run_bridge_global_calibration",
                "--forward_event_catalog_json",
                str(events_train_path),
                "--gamma",
                str(float(gamma)),
                "--V0",
                str(float(V0)),
                "--rounds",
                str(int(rounds)),
                "--reg_strength",
                str(float(reg_strength)),
                "--output_json",
                str(bridge_json),
                "--output_md",
                str(bridge_md),
            ]
        )


def _run_train_head2head(
    *,
    comparison_mode: str,
    omega_r: float,
    omega_b: float,
    omega_c: float,
    omega_l: float,
    cc_train_path: Path,
    sn_train_path: Path,
    bao_train_mean: str,
    bao_train_cov: str,
    events_train_path: str,
    bridge_json: Path,
    strong_field_weight: float,
    theta_samples_train: int,
    train_md: Path,
    train_json: Path,
) -> None:
    _run(
        [
            sys.executable,
            "-m",
            "scripts.run_head2head_comparison",
            "--comparison_mode",
            str(comparison_mode),
            "--omega_r",
            str(float(omega_r)),
            "--omega_b",
            str(float(omega_b)),
            "--omega_c",
            str(float(omega_c)),
            "--omega_l",
            str(float(omega_l)),
            "--cc",
            str(cc_train_path),
            "--sn",
            str(sn_train_path),
            "--bao_mean",
            str(bao_train_mean),
            "--bao_cov",
            str(bao_train_cov),
            "--forward_event_catalog_json",
            str(events_train_path),
            "--bridge_coeffs_json",
            str(bridge_json),
            "--include_strong_field",
            "--strong_field_weight",
            str(float(strong_field_weight)),
            "--theta_samples",
            str(int(theta_samples_train)),
            "--output_md",
            str(train_md),
            "--output_json",
            str(train_json),
        ]
    )


def main(argv: Sequence[str] | None = None) -> int:
    args = parse_args(argv)
    workdir = Path(args.workdir)
    split_dir = workdir / "splits"
    split_dir.mkdir(parents=True, exist_ok=True)

    cc_tr, cc_te, cc_fields = _split_csv_rows(Path(args.cc), int(args.split_mod), int(args.test_fold))
    sn_tr, sn_te, sn_fields = _split_csv_rows(Path(args.sn), int(args.split_mod), int(args.test_fold))
    cc_tr_path = split_dir / "cc_train.csv"
    cc_te_path = split_dir / "cc_test.csv"
    sn_tr_path = split_dir / "sn_train.csv"
    sn_te_path = split_dir / "sn_test.csv"
    _write_csv(cc_tr_path, cc_tr, cc_fields)
    _write_csv(cc_te_path, cc_te, cc_fields)
    _write_csv(sn_tr_path, sn_tr, sn_fields)
    _write_csv(sn_te_path, sn_te, sn_fields)
    bao_paths = _split_bao(Path(args.bao_mean), Path(args.bao_cov), int(args.split_mod), int(args.test_fold), split_dir)
    events_paths = _split_events(Path(args.forward_event_catalog_json), int(args.split_mod), int(args.test_fold), split_dir)

    bridge_json = workdir / "bridge_coeffs_train.json"
    bridge_md = workdir / "bridge_calibration_train.md"
    _run_bridge_calibration(
        bridge_model=str(args.bridge_model),
        events_train_path=events_paths["train"]["path"],
        gamma=float(args.bridge_base_gamma),
        V0=float(args.bridge_base_v0),
        rounds=int(args.bridge_rounds),
        reg_strength=float(args.bridge_reg_strength),
        bridge_json=bridge_json,
        bridge_md=bridge_md,
    )

    train_json = workdir / "head2head_train.json"
    train_md = workdir / "head2head_train.md"
    _run_train_head2head(
        comparison_mode=str(args.comparison_mode),
        omega_r=float(args.omega_r),
        omega_b=float(args.omega_b),
        omega_c=float(args.omega_c),
        omega_l=float(args.omega_l),
        cc_train_path=cc_tr_path,
        sn_train_path=sn_tr_path,
        bao_train_mean=bao_paths["train"]["mean"],
        bao_train_cov=bao_paths["train"]["cov"],
        events_train_path=events_paths["train"]["path"],
        bridge_json=bridge_json,
        strong_field_weight=float(args.strong_field_weight),
        theta_samples_train=int(args.theta_samples_train),
        train_md=train_md,
        train_json=train_json,
    )

    if bool(args.bridge_recenter_on_train):
        tr0 = json.loads(train_json.read_text(encoding="utf-8"))
        g0 = float(tr0["theta_psi"]["gamma"])
        v0 = float(tr0["theta_psi"]["V0"])
        _run_bridge_calibration(
            bridge_model=str(args.bridge_model),
            events_train_path=events_paths["train"]["path"],
            gamma=g0,
            V0=v0,
            rounds=int(args.bridge_rounds),
            reg_strength=float(args.bridge_reg_strength),
            bridge_json=bridge_json,
            bridge_md=bridge_md,
        )
        _run_train_head2head(
            comparison_mode=str(args.comparison_mode),
            omega_r=float(args.omega_r),
            omega_b=float(args.omega_b),
            omega_c=float(args.omega_c),
            omega_l=float(args.omega_l),
            cc_train_path=cc_tr_path,
            sn_train_path=sn_tr_path,
            bao_train_mean=bao_paths["train"]["mean"],
            bao_train_cov=bao_paths["train"]["cov"],
            events_train_path=events_paths["train"]["path"],
            bridge_json=bridge_json,
            strong_field_weight=float(args.strong_field_weight),
            theta_samples_train=int(args.theta_samples_train),
            train_md=train_md,
            train_json=train_json,
        )

    tr = json.loads(train_json.read_text(encoding="utf-8"))
    theta_train = tr["theta_psi"]
    lcdm_train = tr["lcdm_gr"]
    bridge_coeffs = json.loads(bridge_json.read_text(encoding="utf-8"))

    cc_test = load_cc_dimless(cc_te_path)
    sn_test = load_sn_dimless(sn_te_path)
    bao_test = load_bao_vector_from_desi(bao_paths["test"]["mean"], bao_paths["test"]["cov"])
    sf_test = load_forward_event_catalog_json(events_paths["test"]["path"])
    c_km_s = 299792.458
    alpha = c_km_s / (70.0 * 147.09)
    omegas = Omegas(H0=1.0, Omega_r=float(args.omega_r), Omega_b=float(args.omega_b), Omega_c=float(args.omega_c), Omega_L=float(args.omega_l))

    chi2_theta_test, theta_parts_test = _theta_total_chi2(
        float(theta_train["gamma"]),
        float(theta_train["V0"]),
        cc=cc_test,
        sn=sn_test,
        bao=bao_test,
        sn_cov=None,
        alpha=alpha,
        omegas=omegas,
        strong_field_events=sf_test,
        strong_field_weight=float(args.strong_field_weight),
        bridge_coeffs=bridge_coeffs,
        bridge_forward_calibration_rounds=0,
        bridge_forward_calibration_reg=0.0,
    )

    if str(args.comparison_mode) == "same_content":
        omega_m_test = float(omegas.Omega_m)
    else:
        omega_m_test = float(lcdm_train["omega_m"])

    chi2_lcdm_test, lcdm_parts_test = _lcdm_total_chi2(
        float(omega_m_test),
        cc=cc_test,
        sn=sn_test,
        bao=bao_test,
        sn_cov=None,
        alpha=alpha,
        omegas=omegas,
        strong_field_events=sf_test,
        strong_field_weight=float(args.strong_field_weight),
    )

    n_test = len(cc_test.z) + len(sn_test.z) + len(bao_test.z) + len(sf_test)
    # Fixed-parameter evaluation on test: k=0 for both.
    delta_chi2 = float(chi2_theta_test - chi2_lcdm_test)
    delta_aic = delta_chi2
    delta_bic = delta_chi2

    fail_a = bool(delta_bic > float(args.bic_fail_threshold))
    finite_ok = all(
        math.isfinite(float(x))
        for x in (
            chi2_theta_test,
            chi2_lcdm_test,
            delta_chi2,
            theta_parts_test["chi2_cc"],
            theta_parts_test["chi2_sn"],
            theta_parts_test["chi2_bao"],
            theta_parts_test["chi2_sf"],
            lcdm_parts_test["chi2_cc"],
            lcdm_parts_test["chi2_sn"],
            lcdm_parts_test["chi2_bao"],
            lcdm_parts_test["chi2_sf"],
        )
    )
    # Criterion B: hidden tuning on test is forbidden; additionally, non-finite test
    # outputs invalidate the model run.
    fail_b = (not finite_ok)
    d_cc = float(theta_parts_test["chi2_cc"] - lcdm_parts_test["chi2_cc"])
    d_sn = float(theta_parts_test["chi2_sn"] - lcdm_parts_test["chi2_sn"])
    d_bao = float(theta_parts_test["chi2_bao"] - lcdm_parts_test["chi2_bao"])
    d_sf = float(theta_parts_test["chi2_sf"] - lcdm_parts_test["chi2_sf"])
    # Criterion C (quantitative): fail if SN and BAO both degrade beyond tolerances
    # and strong-field does not compensate with sufficient improvement.
    fail_c = (
        (d_sn > float(args.c_sn_tol))
        and (d_bao > float(args.c_bao_tol))
        and (d_sf > -float(args.c_min_sf_improvement))
    )
    verdict = "FAIL" if (fail_a or fail_b or fail_c) else "PASS"

    report = {
        "verdict": verdict,
        "criteria": {
            "A_delta_BIC_gt_threshold": {"failed": fail_a, "delta_BIC_test": delta_bic, "threshold": float(args.bic_fail_threshold)},
            "B_hidden_tuning_required": {
                "failed": fail_b,
                "note": ("non-finite test outputs" if (not finite_ok) else "fixed bridge on test; no inner event calibration"),
            },
            "C_no_joint_consistency_bao_sn_sf": {
                "failed": fail_c,
                "delta_components": {"cc": d_cc, "sn": d_sn, "bao": d_bao, "sf": d_sf},
                "thresholds": {
                    "sn_tol": float(args.c_sn_tol),
                    "bao_tol": float(args.c_bao_tol),
                    "min_sf_improvement": float(args.c_min_sf_improvement),
                },
            },
        },
        "split": {
            "mod": int(args.split_mod),
            "test_fold": int(args.test_fold),
            "n_cc_train": len(cc_tr),
            "n_cc_test": len(cc_te),
            "n_sn_train": len(sn_tr),
            "n_sn_test": len(sn_te),
            "n_bao_train": int(bao_paths["train"]["n"]),
            "n_bao_test": int(bao_paths["test"]["n"]),
            "n_events_train": int(events_paths["train"]["n"]),
            "n_events_test": int(events_paths["test"]["n"]),
        },
        "train_fit": tr,
        "test_eval_fixed": {
            "n_obs_test": n_test,
            "theta": {
                "gamma": float(theta_train["gamma"]),
                "V0": float(theta_train["V0"]),
                "chi2_total": float(chi2_theta_test),
                **{k: float(v) for k, v in theta_parts_test.items()},
            },
            "lcdm_gr": {
                "omega_m_eval": float(omega_m_test),
                "chi2_total": float(chi2_lcdm_test),
                **{k: float(v) for k, v in lcdm_parts_test.items()},
            },
            "delta_theta_minus_lcdm": {"chi2": delta_chi2, "AIC": delta_aic, "BIC": delta_bic},
        },
        "artifacts": {
            "bridge_model": str(args.bridge_model),
            "bridge_coeffs_train_json": str(bridge_json),
            "head2head_train_json": str(train_json),
        },
    }

    out_json = Path(args.output_json)
    out_json.parent.mkdir(parents=True, exist_ok=True)
    out_json.write_text(json.dumps(report, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")

    lines = []
    lines.append("# Falsification Protocol Report\n\n")
    lines.append(f"- Verdict: **{verdict}**\n")
    lines.append(f"- Test delta chi2: {delta_chi2:.6f}\n")
    lines.append(f"- Test delta BIC: {delta_bic:.6f}\n")
    lines.append(f"- Criterion A (delta_BIC > {float(args.bic_fail_threshold):.3g}): {'FAIL' if fail_a else 'PASS'}\n")
    lines.append(f"- Criterion B (hidden tuning required): {'FAIL' if fail_b else 'PASS'}\n")
    lines.append(
        f"- Criterion C (dSN>{float(args.c_sn_tol):.3g} and dBAO>{float(args.c_bao_tol):.3g} "
        f"and dSF>-{float(args.c_min_sf_improvement):.3g}): {'FAIL' if fail_c else 'PASS'}\n\n"
    )
    lines.append("## Test Components\n")
    lines.append(
        f"- Theta chi2: total={chi2_theta_test:.6f}, cc={theta_parts_test['chi2_cc']:.6f}, "
        f"sn={theta_parts_test['chi2_sn']:.6f}, bao={theta_parts_test['chi2_bao']:.6f}, sf={theta_parts_test['chi2_sf']:.6f}\n"
    )
    lines.append(
        f"- LCDM/GR chi2: total={chi2_lcdm_test:.6f}, cc={lcdm_parts_test['chi2_cc']:.6f}, "
        f"sn={lcdm_parts_test['chi2_sn']:.6f}, bao={lcdm_parts_test['chi2_bao']:.6f}, sf={lcdm_parts_test['chi2_sf']:.6f}\n"
    )
    lines.append(
        f"- Delta components (Theta-LCDM): cc={d_cc:.6f}, sn={d_sn:.6f}, bao={d_bao:.6f}, sf={d_sf:.6f}\n"
    )
    lines.append("\n## Split Sizes\n")
    lines.append(
        f"- CC train/test: {len(cc_tr)}/{len(cc_te)}; SN train/test: {len(sn_tr)}/{len(sn_te)}; "
        f"BAO train/test: {bao_paths['train']['n']}/{bao_paths['test']['n']}; "
        f"Events train/test: {events_paths['train']['n']}/{events_paths['test']['n']}\n"
    )

    out_md = Path(args.output_md)
    out_md.parent.mkdir(parents=True, exist_ok=True)
    out_md.write_text("".join(lines), encoding="utf-8")
    print(f"Wrote {out_md}")
    print(f"Wrote {out_json}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
