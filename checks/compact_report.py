# -*- coding: utf-8 -*-
"""
Compact-object report helpers for the finite optical-index regime.

This module computes simple, finite diagnostic observables for the
"black-star" (no true horizon) branch:
- n_surface, n_peak
- z_surface = n_surface / n_obs - 1
- delay_proxy = sum_i (n_i - n_obs) along a sampled trajectory

The diagnostics are intentionally lightweight and model-internal. They are
useful for scans and triage before full ray-tracing / strong-field fits.
"""
from __future__ import annotations

from pathlib import Path
from typing import Iterable, Mapping, Sequence

from theory.optical_metric import canonical_refractive_index_from_theta


def compute_compact_metrics(
    theta_values: Iterable[float],
    *,
    theta_scale: float = 1.0,
    n_obs: float = 1.0,
) -> dict:
    theta = list(theta_values)
    if not theta:
        raise ValueError("theta_values must be non-empty")
    if n_obs <= 0.0:
        raise ValueError("n_obs must be > 0")

    n_vals = [
        canonical_refractive_index_from_theta(th, theta_scale=theta_scale)
        for th in theta
    ]
    n_surface = n_vals[-1]
    n_peak = max(n_vals)
    z_surface = n_surface / n_obs - 1.0
    delay_proxy = sum((nv - n_obs) for nv in n_vals)
    return {
        "n_surface": n_surface,
        "n_peak": n_peak,
        "z_surface": z_surface,
        "z_surface_proxy": z_surface,
        "delay": delay_proxy,
        "delay_proxy": delay_proxy,
        "has_horizon": False,
    }


def write_compact_markdown(rows: Sequence[Mapping[str, object]], path: str) -> None:
    p = Path(path)
    p.parent.mkdir(parents=True, exist_ok=True)
    rows = list(rows)
    ok_rows = [r for r in rows if bool(r.get("ok", False))]
    ok_n_rows = [r for r in rows if bool(r.get("ok_n_bounds", True))]
    best = max(rows, key=lambda r: float(r.get("z_surface", r.get("z_surface_proxy", 0.0)))) if rows else None

    lines = []
    lines.append("# Compact Regime Report\n\n")
    lines.append(f"- Total rows: {len(rows)}\n")
    lines.append(f"- Health ok rows: {len(ok_rows)}\n")
    lines.append(f"- Bounded-n ok rows: {len(ok_n_rows)}\n")
    if best is not None:
        lines.append("\n## Peak Redshift\n")
        lines.append(f"- gamma: {best.get('gamma')}\n")
        lines.append(f"- V0: {best.get('V0')}\n")
        zval = float(best.get("z_surface", best.get("z_surface_proxy", 0.0)))
        dval = float(best.get("delay", best.get("delay_proxy", 0.0)))
        lines.append(f"- z_surface: {zval:.6g}\n")
        lines.append(f"- n_surface: {float(best.get('n_surface', 0.0)):.6g}\n")
        lines.append(f"- n_peak: {float(best.get('n_peak', 0.0)):.6g}\n")
        lines.append(f"- delay: {dval:.6g}\n")

    p.write_text("".join(lines), encoding="utf-8")
