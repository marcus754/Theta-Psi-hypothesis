# -*- coding: utf-8 -*-
"""
Observational comparison layer for short comparative fits.

Runs short Metropolis fits for:
  A: CC+SN+BAO
  B: A + compact_star (linked, primary strong-field block)
  C: B + compact (sanity, low weight)
"""
from __future__ import annotations

import argparse
import csv
import random
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Sequence

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
from fitting.run_mcmc import Params, step_proposal, logposterior_combined, loglikes_components


@dataclass
class FitResult:
    name: str
    best_logpost: float
    best_gamma: float
    best_v0: float
    best_q: float
    accept_rate: float
    comp_ll: Dict[str, float]


def _scale_compact_sigma(c: CompactTargets, scale: float) -> CompactTargets:
    if c.mode == "constraints":
        return CompactTargets(
            mode="constraints",
            z_surface_min=c.z_surface_min,
            sigma_z_surface_min=None if c.sigma_z_surface_min is None else max(1e-9, c.sigma_z_surface_min * scale),
            delay_proxy_min=c.delay_proxy_min,
            sigma_delay_proxy_min=None if c.sigma_delay_proxy_min is None else max(1e-9, c.sigma_delay_proxy_min * scale),
            n_peak_max=c.n_peak_max,
            sigma_n_peak_max=None if c.sigma_n_peak_max is None else max(1e-9, c.sigma_n_peak_max * scale),
        )
    return CompactTargets(
        mode="gaussian",
        z_surface=c.z_surface,
        sigma_z_surface=max(1e-9, c.sigma_z_surface * scale),
        delay_proxy=c.delay_proxy,
        sigma_delay_proxy=max(1e-9, c.sigma_delay_proxy * scale),
        n_peak=c.n_peak,
        sigma_n_peak=max(1e-9, c.sigma_n_peak * scale),
    )


def _run_one(
    *,
    name: str,
    steps: int,
    seed: int,
    om: Omegas,
    cc,
    sn,
    bao,
    compact: CompactTargets | None,
    compact_star: CompactStarTargets | None,
    forward: ForwardTargets | None,
    forward_events: list[dict] | None,
    compact_theta_scale: float,
    compact_weight: float,
    compact_star_theta_scale: float,
    compact_star_weight: float,
    compact_star_link_background: bool,
    forward_weight: float,
    forward_events_weight: float,
) -> FitResult:
    rng = random.Random(seed)
    p = Params(gamma=1.0, V0=0.1, Q=0.0)
    logp = logposterior_combined(
        p,
        om=om,
        cc=cc,
        sn=sn,
        bao=bao,
        compact=compact,
        compact_theta_scale=compact_theta_scale,
        compact_weight=compact_weight,
        compact_star=compact_star,
        compact_star_link_model=True,
        compact_star_theta_c_from_gamma=1.5,
        compact_star_m2_from_v0=10.0,
        compact_star_lam_from_q2=0.1,
        compact_star_theta_scale=compact_star_theta_scale,
        compact_star_weight=compact_star_weight,
        compact_star_link_background=compact_star_link_background,
        forward=forward,
        forward_weight=forward_weight,
        forward_events=forward_events,
        forward_events_weight=forward_events_weight,
    )
    best_logp = logp
    best_p = p
    accepted = 0
    for _ in range(steps):
        cand = step_proposal(p, rng)
        logc = logposterior_combined(
            cand,
            om=om,
            cc=cc,
            sn=sn,
            bao=bao,
            compact=compact,
            compact_theta_scale=compact_theta_scale,
            compact_weight=compact_weight,
            compact_star=compact_star,
            compact_star_link_model=True,
            compact_star_theta_c_from_gamma=1.5,
            compact_star_m2_from_v0=10.0,
            compact_star_lam_from_q2=0.1,
            compact_star_theta_scale=compact_star_theta_scale,
            compact_star_weight=compact_star_weight,
            compact_star_link_background=compact_star_link_background,
            forward=forward,
            forward_weight=forward_weight,
            forward_events=forward_events,
            forward_events_weight=forward_events_weight,
        )
        if (logc - logp) >= __import__("math").log(rng.random()):
            p, logp = cand, logc
            accepted += 1
            if logp > best_logp:
                best_logp = logp
                best_p = p
    comps = loglikes_components(
        best_p,
        om=om,
        cc=cc,
        sn=sn,
        bao=bao,
        compact=compact,
        compact_theta_scale=compact_theta_scale,
        compact_weight=compact_weight,
        compact_star=compact_star,
        compact_star_link_model=True,
        compact_star_theta_c_from_gamma=1.5,
        compact_star_m2_from_v0=10.0,
        compact_star_lam_from_q2=0.1,
        compact_star_theta_scale=compact_star_theta_scale,
        compact_star_weight=compact_star_weight,
        compact_star_link_background=compact_star_link_background,
        forward=forward,
        forward_weight=forward_weight,
        forward_events=forward_events,
        forward_events_weight=forward_events_weight,
    )
    return FitResult(
        name=name,
        best_logpost=best_logp,
        best_gamma=best_p.gamma,
        best_v0=best_p.V0,
        best_q=best_p.Q,
        accept_rate=accepted / max(steps, 1),
        comp_ll=comps,
    )


def parse_args(argv: Sequence[str] | None = None) -> argparse.Namespace:
    p = argparse.ArgumentParser(description="Run likelihood comparison A/B/C/D")
    p.add_argument("--steps", type=int, default=120)
    p.add_argument("--seed", type=int, default=11)
    p.add_argument("--compact_targets_json", type=str, default="fitting/data/compact_constraints_observational_proxy.json")
    p.add_argument("--compact_star_targets_json", type=str, default="fitting/data/compact_star_observational_proxy.json")
    p.add_argument("--forward_targets_json", type=str, default="fitting/data/forward_observational_proxy.json")
    p.add_argument("--forward_event_catalog_json", type=str, default="fitting/data/forward_event_catalog.json")
    p.add_argument("--compact_sanity_weight", type=float, default=0.05)
    p.add_argument("--forward_weight", type=float, default=1.0)
    p.add_argument("--forward_events_weight", type=float, default=1.0)
    p.add_argument("--output_csv", type=str, default="results/likelihood_comparison.csv")
    p.add_argument("--output_md", type=str, default="results/likelihood_comparison.md")
    return p.parse_args(argv)


def main(argv: Sequence[str] | None = None) -> int:
    args = parse_args(argv)
    om = Omegas(H0=1.0)
    cc = load_cc_dimless("fitting/data/cc_demo_dimless.csv")
    sn = load_sn_dimless("fitting/data/sn_demo_dimless.csv")
    bao = load_bao_dimless("fitting/data/bao_demo_dimless.csv")
    c = load_compact_targets_json(args.compact_targets_json)
    compact = CompactTargets(
        mode=str(c.get("mode", "gaussian")),
        z_surface=c.get("z_surface"),
        sigma_z_surface=c.get("sigma_z_surface"),
        delay_proxy=c.get("delay_proxy"),
        sigma_delay_proxy=c.get("sigma_delay_proxy"),
        n_peak=c.get("n_peak"),
        sigma_n_peak=c.get("sigma_n_peak"),
        z_surface_min=c.get("z_surface_min"),
        sigma_z_surface_min=c.get("sigma_z_surface_min"),
        delay_proxy_min=c.get("delay_proxy_min"),
        sigma_delay_proxy_min=c.get("sigma_delay_proxy_min"),
        n_peak_max=c.get("n_peak_max"),
        sigma_n_peak_max=c.get("sigma_n_peak_max"),
    )
    cs = load_compact_star_targets_json(args.compact_star_targets_json)
    compact_star = CompactStarTargets(
        z_surface=cs["z_surface"],
        sigma_z_surface=cs["sigma_z_surface"],
        delay_proxy=cs["delay_proxy"],
        sigma_delay_proxy=cs["sigma_delay_proxy"],
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

    results = []
    results.append(_run_one(
        name="A_base", steps=args.steps, seed=args.seed, om=om, cc=cc, sn=sn, bao=bao,
        compact=None, compact_star=None, forward=None, forward_events=None, compact_theta_scale=1.0, compact_weight=1.0,
        compact_star_theta_scale=1.0, compact_star_weight=1.0, compact_star_link_background=False,
        forward_weight=args.forward_weight, forward_events_weight=args.forward_events_weight
    ))
    results.append(_run_one(
        name="B_plus_compact_star", steps=args.steps, seed=args.seed + 1, om=om, cc=cc, sn=sn, bao=bao,
        compact=None, compact_star=compact_star, forward=forward, forward_events=forward_events, compact_theta_scale=1.0, compact_weight=1.0,
        compact_star_theta_scale=1.0, compact_star_weight=1.0, compact_star_link_background=True,
        forward_weight=args.forward_weight, forward_events_weight=args.forward_events_weight
    ))
    results.append(_run_one(
        name="C_plus_compact_star_plus_compact_sanity", steps=args.steps, seed=args.seed + 2,
        om=om, cc=cc, sn=sn, bao=bao,
        compact=compact, compact_star=compact_star, forward=forward, forward_events=forward_events, compact_theta_scale=1.0, compact_weight=args.compact_sanity_weight,
        compact_star_theta_scale=1.0, compact_star_weight=1.0, compact_star_link_background=True,
        forward_weight=args.forward_weight, forward_events_weight=args.forward_events_weight
    ))

    out_csv = Path(args.output_csv)
    out_csv.parent.mkdir(parents=True, exist_ok=True)
    rows = []
    for r in results:
        row = {
            "name": r.name,
            "best_logpost": r.best_logpost,
            "best_gamma": r.best_gamma,
            "best_v0": r.best_v0,
            "best_q": r.best_q,
            "accept_rate": r.accept_rate,
        }
        for k, v in r.comp_ll.items():
            row[f"loglike_{k}"] = v
            row[f"chi2_{k}"] = -2.0 * v
        if r.comp_ll:
            worst = min(r.comp_ll.items(), key=lambda kv: kv[1])[0]
            row["dominant_tension"] = worst
        rows.append(row)
    fields = sorted({k for row in rows for k in row.keys()})
    with out_csv.open("w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=fields)
        w.writeheader()
        w.writerows(rows)

    out_md = Path(args.output_md)
    out_md.parent.mkdir(parents=True, exist_ok=True)
    lines = []
    lines.append("# Likelihood Comparison v2 (A/B/C)\n\n")
    lines.append(f"- Steps per config: {args.steps}\n")
    lines.append(f"- Compact sanity weight (C): {args.compact_sanity_weight}\n")
    lines.append(f"- Forward weight (B/C): {args.forward_weight}\n")
    lines.append(f"- Forward-events weight (B/C): {args.forward_events_weight}\n")
    best = max(results, key=lambda rr: rr.best_logpost)
    lines.append(f"- Best configuration by logpost: {best.name} ({best.best_logpost:.3f})\n")
    lines.append("\n## Config Results\n")
    for r in results:
        lines.append(f"- {r.name}: logpost={r.best_logpost:.3f}, gamma={r.best_gamma:.4g}, V0={r.best_v0:.4g}, accept={r.accept_rate:.2f}\n")
        if r.comp_ll:
            worst = min(r.comp_ll.items(), key=lambda kv: kv[1])
            lines.append(f"  dominant_tension={worst[0]}, chi2={-2.0*worst[1]:.3f}\n")
    out_md.write_text("".join(lines), encoding="utf-8")
    print(f"Wrote likelihood comparison CSV to {args.output_csv}")
    print(f"Wrote likelihood comparison report to {args.output_md}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
