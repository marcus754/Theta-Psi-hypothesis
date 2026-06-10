# -*- coding: utf-8 -*-
"""
K-fold wrapper over train/test falsification protocol.
"""
from __future__ import annotations

import argparse
import json
import math
import subprocess
import sys
from pathlib import Path
from typing import Sequence


def parse_args(argv: Sequence[str] | None = None) -> argparse.Namespace:
    p = argparse.ArgumentParser(description="Run k-fold falsification protocol")
    p.add_argument("--comparison_mode", type=str, choices=["same_content", "lcdm_fit"], default="same_content")
    p.add_argument("--split_mod", type=int, default=3)
    p.add_argument("--theta_samples_train", type=int, default=180)
    p.add_argument("--bridge_rounds", type=int, default=32)
    p.add_argument("--bridge_reg_strength", type=float, default=0.2)
    p.add_argument("--bridge_model", type=str, choices=["constant", "powerlaw2d"], default="constant")
    p.add_argument("--bridge_recenter_on_train", action=argparse.BooleanOptionalAction, default=False)
    p.add_argument("--workdir", type=str, default="results/falsification_kfold")
    p.add_argument("--output_md", type=str, default="results/falsification_kfold/report.md")
    p.add_argument("--output_json", type=str, default="results/falsification_kfold/report.json")
    return p.parse_args(argv)


def _run(cmd: list[str]) -> None:
    proc = subprocess.run(cmd, check=False, capture_output=True, text=True)
    if proc.returncode != 0:
        raise RuntimeError(f"Command failed: {' '.join(cmd)}\nSTDOUT:\n{proc.stdout}\nSTDERR:\n{proc.stderr}")


def main(argv: Sequence[str] | None = None) -> int:
    args = parse_args(argv)
    workdir = Path(args.workdir)
    workdir.mkdir(parents=True, exist_ok=True)

    folds: list[dict] = []
    fail_a = 0
    fail_b = 0
    fail_c = 0
    fail_verdict = 0
    for fold in range(int(args.split_mod)):
        fold_dir = workdir / f"fold_{fold}"
        md = fold_dir / "report.md"
        js = fold_dir / "report.json"
        cmd = [
            sys.executable,
            "-m",
            "scripts.run_falsification_protocol",
            "--comparison_mode",
            str(args.comparison_mode),
            "--split_mod",
            str(int(args.split_mod)),
            "--test_fold",
            str(int(fold)),
            "--theta_samples_train",
            str(int(args.theta_samples_train)),
            "--bridge_rounds",
            str(int(args.bridge_rounds)),
            "--bridge_reg_strength",
            str(float(args.bridge_reg_strength)),
            "--bridge_model",
            str(args.bridge_model),
            "--workdir",
            str(fold_dir),
            "--output_md",
            str(md),
            "--output_json",
            str(js),
        ]
        cmd.append("--bridge_recenter_on_train" if bool(args.bridge_recenter_on_train) else "--no-bridge_recenter_on_train")
        _run(cmd)
        rep = json.loads(js.read_text(encoding="utf-8"))
        d_bic = float(rep["criteria"]["A_delta_BIC_gt_threshold"]["delta_BIC_test"])
        fa = bool(rep["criteria"]["A_delta_BIC_gt_threshold"]["failed"])
        fb = bool(rep["criteria"]["B_hidden_tuning_required"]["failed"])
        fc = bool(rep["criteria"]["C_no_joint_consistency_bao_sn_sf"]["failed"])
        if fa:
            fail_a += 1
        if fb:
            fail_b += 1
        if fc:
            fail_c += 1
        if str(rep["verdict"]).upper() == "FAIL":
            fail_verdict += 1
        folds.append(
            {
                "fold": fold,
                "verdict": rep["verdict"],
                "delta_BIC_test": d_bic,
                "delta_chi2_test": float(rep["test_eval_fixed"]["delta_theta_minus_lcdm"]["chi2"]),
                "fail_A": fa,
                "fail_B": fb,
                "fail_C": fc,
            }
        )

    db = [float(f["delta_BIC_test"]) for f in folds]
    db_fin = [x for x in db if math.isfinite(x)]
    n_nan = len(db) - len(db_fin)
    mean_db = (sum(db_fin) / max(len(db_fin), 1)) if db_fin else float("nan")
    max_db = max(db_fin) if db_fin else float("nan")
    min_db = min(db_fin) if db_fin else float("nan")
    fail6 = sum(1 for x in db_fin if x > 6.0)
    fail10 = sum(1 for x in db_fin if x > 10.0)
    all_pass10 = (fail10 == 0)

    payload = {
        "comparison_mode": str(args.comparison_mode),
        "split_mod": int(args.split_mod),
        "theta_samples_train": int(args.theta_samples_train),
        "bridge_rounds": int(args.bridge_rounds),
        "bridge_reg_strength": float(args.bridge_reg_strength),
        "bridge_model": str(args.bridge_model),
        "bridge_recenter_on_train": bool(args.bridge_recenter_on_train),
        "folds": folds,
        "summary": {
            "mean_delta_BIC_test": mean_db,
            "min_delta_BIC_test": min_db,
            "max_delta_BIC_test": max_db,
            "n_folds_delta_BIC_gt_6": fail6,
            "n_folds_delta_BIC_gt_10": fail10,
            "n_folds_nonfinite_delta_BIC": n_nan,
            "all_folds_pass_threshold_10": all_pass10,
            "n_folds_verdict_FAIL": fail_verdict,
            "n_folds_fail_A": fail_a,
            "n_folds_fail_B": fail_b,
            "n_folds_fail_C": fail_c,
        },
    }

    out_json = Path(args.output_json)
    out_json.parent.mkdir(parents=True, exist_ok=True)
    out_json.write_text(json.dumps(payload, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")

    lines = []
    lines.append("# Falsification K-Fold Report\n\n")
    lines.append(f"- comparison_mode: {args.comparison_mode}\n")
    lines.append(f"- split_mod (k): {int(args.split_mod)}\n")
    lines.append(f"- bridge_model: {args.bridge_model}\n")
    lines.append(f"- bridge_recenter_on_train: {'on' if bool(args.bridge_recenter_on_train) else 'off'}\n")
    lines.append(f"- mean delta_BIC_test: {mean_db:.6f}\n")
    lines.append(f"- min/max delta_BIC_test: {min_db:.6f} / {max_db:.6f}\n")
    lines.append(f"- folds with delta_BIC_test > 6: {fail6}\n")
    lines.append(f"- folds with delta_BIC_test > 10: {fail10}\n")
    lines.append(f"- folds with non-finite delta_BIC_test: {n_nan}\n")
    lines.append(f"- all folds pass threshold 10: {'yes' if all_pass10 else 'no'}\n\n")
    lines.append(f"- folds with verdict FAIL: {fail_verdict}\n")
    lines.append(f"- folds fail criterion A/B/C: {fail_a}/{fail_b}/{fail_c}\n\n")
    lines.append("## Folds\n")
    for f in folds:
        lines.append(
            f"- fold={int(f['fold'])}: verdict={f['verdict']}, "
            f"delta_BIC_test={float(f['delta_BIC_test']):.6f}, "
            f"delta_chi2_test={float(f['delta_chi2_test']):.6f}, "
            f"fails(A/B/C)={int(bool(f['fail_A']))}/{int(bool(f['fail_B']))}/{int(bool(f['fail_C']))}\n"
        )

    out_md = Path(args.output_md)
    out_md.parent.mkdir(parents=True, exist_ok=True)
    out_md.write_text("".join(lines), encoding="utf-8")
    print(f"Wrote {out_md}")
    print(f"Wrote {out_json}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
