# -*- coding: utf-8 -*-
"""Single entrypoint for parameter fitting."""
from __future__ import annotations

import argparse
import concurrent.futures
import re
import subprocess
import sys
from pathlib import Path
from typing import Sequence


def parse_args(argv: Sequence[str] | None = None) -> argparse.Namespace:
    p = argparse.ArgumentParser(description="Run fitting pipeline (MCMC)")
    p.add_argument("--steps", type=int, default=400)
    p.add_argument("--jobs", type=int, default=1, help="Number of parallel independent chains")
    p.add_argument("--seed", type=int, default=7)
    p.add_argument("--data", type=str, default=None)
    p.add_argument("--sn", type=str, default=None)
    p.add_argument("--bao", type=str, default=None)
    p.add_argument("--H0_km_s_Mpc", type=float, default=None)
    p.add_argument("--cc_stride", type=int, default=1)
    p.add_argument("--sn_stride", type=int, default=1)
    p.add_argument("--bao_stride", type=int, default=1)
    p.add_argument("--sn_fit_mu0", action=argparse.BooleanOptionalAction, default=True)
    p.add_argument("--joint_strong_field", action=argparse.BooleanOptionalAction, default=True)
    p.add_argument("--compact_star_targets_json", type=str, default="fitting/data/compact_star_observational_proxy.json")
    p.add_argument("--forward_targets_json", type=str, default="fitting/data/forward_observational_proxy.json")
    p.add_argument("--forward_event_catalog_json", type=str, default="fitting/data/forward_event_catalog.json")
    p.add_argument("--bridge_forward_calibration_rounds", type=int, default=1)
    p.add_argument("--bridge_forward_calibration_reg", type=float, default=0.1)
    p.add_argument("--compact_star_weight", type=float, default=1.0)
    p.add_argument("--forward_weight", type=float, default=1.0)
    p.add_argument("--forward_events_weight", type=float, default=1.0)
    p.add_argument("--prop_scale_gamma", type=float, default=0.1)
    p.add_argument("--prop_scale_v0", type=float, default=0.1)
    p.add_argument("--prop_scale_q", type=float, default=0.0)
    p.add_argument("--prior_gamma_min", type=float, default=1e-3)
    p.add_argument("--prior_gamma_max", type=float, default=10.0)
    p.add_argument("--prior_V0_min", type=float, default=1e-4)
    p.add_argument("--prior_V0_max", type=float, default=10.0)
    p.add_argument("--prior_Q_abs_max", type=float, default=10.0)
    p.add_argument("--priors_json", type=str, default=None)
    p.add_argument("--target_accept", type=float, default=0.30)
    p.add_argument("--adapt_interval", type=int, default=50)
    p.add_argument("--adapt_until_frac", type=float, default=0.5)
    p.add_argument("--output", type=str, default="results/fit_summary.txt")
    return p.parse_args(argv)


def _build_base_cmd(args: argparse.Namespace) -> list[str]:
    cmd = [
        sys.executable,
        "-m",
        "fitting.run_mcmc",
        "--cc_stride",
        str(args.cc_stride),
        "--sn_stride",
        str(args.sn_stride),
        "--bao_stride",
        str(args.bao_stride),
        "--prop_scale_gamma",
        str(args.prop_scale_gamma),
        "--prop_scale_v0",
        str(args.prop_scale_v0),
        "--prop_scale_q",
        str(args.prop_scale_q),
        "--prior_gamma_min",
        str(args.prior_gamma_min),
        "--prior_gamma_max",
        str(args.prior_gamma_max),
        "--prior_V0_min",
        str(args.prior_V0_min),
        "--prior_V0_max",
        str(args.prior_V0_max),
        "--prior_Q_abs_max",
        str(args.prior_Q_abs_max),
        "--target_accept",
        str(args.target_accept),
        "--adapt_interval",
        str(args.adapt_interval),
        "--adapt_until_frac",
        str(args.adapt_until_frac),
    ]
    if args.priors_json:
        cmd.extend(["--priors_json", str(args.priors_json)])
    cmd.extend(["--sn_fit_mu0"] if args.sn_fit_mu0 else ["--no-sn_fit_mu0"])
    if args.data:
        cmd.extend(["--data", args.data])
    if args.sn:
        cmd.extend(["--sn", args.sn])
    if args.bao:
        cmd.extend(["--bao", args.bao])
    if args.H0_km_s_Mpc is not None:
        cmd.extend(["--H0_km_s_Mpc", str(args.H0_km_s_Mpc)])
    if args.joint_strong_field:
        cmd.extend(
            [
                "--use_compact_star",
                "--compact_star_link_background",
                "--compact_star_targets_json",
                str(args.compact_star_targets_json),
                "--compact_star_weight",
                str(args.compact_star_weight),
                "--use_forward",
                "--forward_targets_json",
                str(args.forward_targets_json),
                "--forward_weight",
                str(args.forward_weight),
                "--use_forward_events",
                "--forward_event_catalog_json",
                str(args.forward_event_catalog_json),
                "--forward_events_weight",
                str(args.forward_events_weight),
                "--bridge_forward_calibration_rounds",
                str(int(args.bridge_forward_calibration_rounds)),
                "--bridge_forward_calibration_reg",
                str(float(args.bridge_forward_calibration_reg)),
            ]
        )
    return cmd


def _run_job(cmd: list[str]) -> int:
    return int(subprocess.run(cmd).returncode)


def _best_logpost(summary_text: str) -> float:
    m = re.search(r"best=\{logp:([+-]?(?:inf|[0-9]*\.?[0-9]+))", summary_text)
    if not m:
        return float("-inf")
    token = m.group(1).lower()
    if token == "inf" or token == "+inf":
        return float("inf")
    if token == "-inf":
        return float("-inf")
    return float(token)


def main(argv: Sequence[str] | None = None) -> int:
    args = parse_args(argv)
    base = _build_base_cmd(args)
    jobs = max(1, int(args.jobs))

    if jobs == 1:
        cmd = base + ["--steps", str(args.steps), "--seed", str(args.seed), "--output", args.output]
        return _run_job(cmd)

    outp = Path(args.output)
    outp.parent.mkdir(parents=True, exist_ok=True)
    steps_per_job = max(1, args.steps // jobs)
    job_files = [outp.with_suffix(outp.suffix + f".job{i}") for i in range(jobs)]
    cmds = [
        base + [
            "--steps",
            str(steps_per_job),
            "--seed",
            str(args.seed + i),
            "--output",
            str(job_files[i]),
        ]
        for i in range(jobs)
    ]

    with concurrent.futures.ThreadPoolExecutor(max_workers=jobs) as ex:
        rcs = list(ex.map(_run_job, cmds))
    if any(rc != 0 for rc in rcs):
        return 1

    best_text = ""
    best_lp = float("-inf")
    for p in job_files:
        txt = p.read_text(encoding="utf-8").strip()
        lp = _best_logpost(txt)
        if best_text == "" or lp > best_lp:
            best_lp = lp
            best_text = txt
    outp.write_text(best_text + "\n", encoding="utf-8")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
