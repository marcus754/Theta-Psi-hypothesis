# -*- coding: utf-8 -*-
"""Observational comparison layer for head-to-head model evaluation."""
from __future__ import annotations

import argparse
import json
import math
import random
from bisect import bisect_left
from pathlib import Path
from typing import Sequence

from fitting.core_api import CompactStarParams, Omegas, solve_compact_star_profile
from fitting.data_io import (
    load_cc_dimless,
    load_sn_dimless,
    load_forward_event_catalog_json,
)
from fitting.strong_field_link import derive_compact_star_params_from_background, resolve_bridge_scales
from fitting.likelihoods_forward_events import loglike_forward_events
from fitting.head2head import (
    BAOVectorData,
    bao_vector_chi2,
    bao_vector_model,
    _build_dc_from_E,
    gr_E_fn,
    gr_predictions_for_cc,
    gr_predictions_for_sn,
    lcdm_predictions_for_cc,
    lcdm_predictions_for_sn,
    load_bao_vector_from_desi,
    sn_chi2,
    theta_predictions_for_cc,
    theta_predictions_for_sn,
    theta_model_E_fn,
)


def parse_args(argv: Sequence[str] | None = None) -> argparse.Namespace:
    p = argparse.ArgumentParser(description="Compare Theta-Psi and flat LCDM/GR on same data")
    p.add_argument("--cc", type=str, default="fitting/data/real_dimless/cc_Ez_zlt1.csv")
    p.add_argument("--sn", type=str, default="fitting/data/real_dimless/sn_dimless.csv")
    p.add_argument("--sn_cov", type=str, default=None, help="Optional SN covariance matrix (NxN text)")
    p.add_argument("--bao_mean", type=str, default="fitting/data/raw_public/desi_2024_gaussian_bao_ALL_GCcomb_mean.txt")
    p.add_argument("--bao_cov", type=str, default="fitting/data/raw_public/desi_2024_gaussian_bao_ALL_GCcomb_cov.txt")
    p.add_argument("--bao_rd_mpc", type=float, default=147.09)
    p.add_argument("--bao_h0_km_s_mpc", type=float, default=70.0)
    p.add_argument(
        "--comparison_mode",
        type=str,
        choices=["lcdm_fit", "same_content", "gr_only"],
        default="lcdm_fit",
        help="lcdm_fit: fit Omega_m for flat LCDM baseline; same_content: use fixed Omegas for both GR and Theta-Psi; gr_only: pure GR baseline with Omega_L=0",
    )
    p.add_argument("--omega_r", type=float, default=8.4e-5)
    p.add_argument("--omega_b", type=float, default=0.049)
    p.add_argument("--omega_c", type=float, default=0.251)
    p.add_argument("--omega_l", type=float, default=0.7)
    p.add_argument("--seed", type=int, default=123)
    p.add_argument("--theta_samples", type=int, default=240)
    p.add_argument(
        "--theta_coarse_init",
        action=argparse.BooleanOptionalAction,
        default=True,
        help="Run coarse multi-start scan for (gamma,V0) before stochastic refinement",
    )
    p.add_argument("--include_strong_field", action=argparse.BooleanOptionalAction, default=True)
    p.add_argument("--forward_event_catalog_json", type=str, default="fitting/data/forward_event_catalog.json")
    p.add_argument("--strong_field_weight", type=float, default=1.0)
    p.add_argument("--bridge_forward_calibration_rounds", type=int, default=0, help="Legacy (ignored in fixed bridge mode)")
    p.add_argument("--bridge_forward_calibration_reg", type=float, default=0.0, help="Legacy (ignored in fixed bridge mode)")
    p.add_argument(
        "--bridge_coeffs_json",
        type=str,
        default="fitting/data/bridge_coeffs_global.json",
        help="Fixed global bridge coefficients JSON (theta_c_scale,m2_scale,lam_scale)",
    )
    p.add_argument(
        "--penalize_bridge_latent",
        action=argparse.BooleanOptionalAction,
        default=True,
        help="Add AIC/BIC penalty for inner bridge calibration (theta_c,m2,lam)",
    )
    p.add_argument(
        "--diagnose_bridge_dependency",
        action=argparse.BooleanOptionalAction,
        default=True,
        help="Evaluate Theta-Psi degradation at same best-fit background point when bridge rounds are forced to 0",
    )
    p.add_argument("--output_md", type=str, default="results/head2head_comparison.md")
    p.add_argument("--output_json", type=str, default="results/head2head_comparison.json")
    return p.parse_args(argv)


def _chi2_gauss_diag(model: list[float], data: list[float], sigma: list[float]) -> float:
    return sum(((m - d) / max(s, 1e-18)) ** 2 for m, d, s in zip(model, data, sigma))


def _gr_ring_uas(mass_msun: float, distance_m: float) -> float:
    # Photon-ring diameter for Schwarzschild-like GR branch: d = 2 b_ph = 6*sqrt(3) * GM/c^2.
    G = 6.67430e-11
    c = 299792458.0
    msun = 1.98847e30
    rad_to_uas = 206264806247.09636
    diam_m = 6.0 * math.sqrt(3.0) * G * mass_msun * msun / (c * c)
    return diam_m / max(distance_m, 1e-30) * rad_to_uas


def _gr_forward_events_chi2(events: list[dict]) -> float:
    chi2 = 0.0
    for e in events:
        obs = str(e["observable"])
        if obs == "ring_uas":
            pred = _gr_ring_uas(float(e["mass_msun"]), float(e["distance_m"]))
        elif obs == "echo_delay_ms":
            # Classical GR horizon branch: no echoes in this simplified comparator.
            pred = 0.0
        else:
            continue
        sig = max(float(e["sigma"]), 1e-18)
        chi2 += ((pred - float(e["target"])) / sig) ** 2
    return chi2


def _theta_forward_events_chi2(
    gamma: float,
    V0: float,
    *,
    events: list[dict],
    bridge_coeffs: dict | None = None,
    bridge_forward_calibration_rounds: int = 1,
    bridge_forward_calibration_reg: float = 0.1,
) -> float:
    scales = resolve_bridge_scales(gamma=gamma, V0=V0, bridge_coeffs=bridge_coeffs)
    d = derive_compact_star_params_from_background(
        gamma=gamma,
        V0=V0,
        Q=0.0,
        theta_c_scale=float(scales["theta_c_scale"]),
        m2_scale=float(scales["m2_scale"]),
        lam_scale=float(scales["lam_scale"]),
        # In head2head production mode we do NOT tune compact parameters per likelihood call.
        forward_events=None,
        calibrate_forward_rounds=0,
        calibrate_forward_reg_strength=0.0,
    )
    prof = solve_compact_star_profile(
        CompactStarParams(
            theta_c=float(d["theta_c"]),
            m2=float(d["m2"]),
            lam=float(d["lam"]),
            r_max=20.0,
            dr=0.02,
            theta_scale=1.0,
        )
    )
    ll = loglike_forward_events(prof, events)
    return -2.0 * float(ll)


def _theta_total_chi2(
    gamma: float,
    V0: float,
    *,
    cc,
    sn,
    bao: BAOVectorData,
    sn_cov,
    alpha,
    omegas: Omegas | None = None,
    strong_field_events: list[dict] | None = None,
    strong_field_weight: float = 1.0,
    bridge_coeffs: dict | None = None,
    bridge_forward_calibration_rounds: int = 1,
    bridge_forward_calibration_reg: float = 0.1,
) -> tuple[float, dict]:
    # Single trajectory build for CC+SN+BAO.
    from fitting.model import Hz_dimless
    import numpy as np

    zmax = max(
        max(cc.z) if cc.z else 0.0,
        max(sn.z) if sn.z else 0.0,
        max(bao.z) if bao.z else 0.0,
    )
    nint = 1200
    zs = [zmax * i / nint for i in range(nint + 1)]
    Es = Hz_dimless(
        zs, gamma=gamma, V0=V0, Q=0.0, omegas=omegas, a_min=1e-3, dN=0.02,
    )
    invE = [1.0 / max(e, 1e-18) for e in Es]
    dc = [0.0] * len(zs)
    s = 0.0
    for i in range(1, len(zs)):
        dz = zs[i] - zs[i - 1]
        s += 0.5 * (invE[i] + invE[i - 1]) * dz
        dc[i] = s

    def interp(y: list[float], xq: float) -> float:
        if xq <= zs[0]:
            return y[0]
        if xq >= zs[-1]:
            return y[-1]
        j = bisect_left(zs, xq)
        x0, x1 = zs[j - 1], zs[j]
        y0, y1 = y[j - 1], y[j]
        t = (xq - x0) / max(x1 - x0, 1e-18)
        return y0 + t * (y1 - y0)

    cc_m = [interp(Es, z) for z in cc.z]
    chi2_cc = _chi2_gauss_diag(cc_m, cc.H, cc.sigma)

    sn_m = [5.0 * math.log10(max((1.0 + z) * interp(dc, z), 1e-18)) for z in sn.z]
    chi2_sn, mu0 = sn_chi2(sn_m, sn.mu, sn.sigma, cov=sn_cov, fit_mu0=True)

    bao_out = []
    for z, o in zip(bao.z, bao.obs):
        Ez = interp(Es, z)
        Dm = interp(dc, z)
        if o == "DM_over_rs":
            bao_out.append(alpha * Dm)
        elif o == "DH_over_rs":
            bao_out.append(alpha / max(Ez, 1e-18))
        else:
            bao_out.append(alpha * ((z / max(Ez, 1e-18)) * (Dm * Dm)) ** (1.0 / 3.0))
    bao_m = np.asarray(bao_out, dtype=float)
    chi2_bao = bao_vector_chi2(bao_m, bao)
    chi2_sf = 0.0
    if strong_field_events is not None:
        chi2_sf = float(strong_field_weight) * _theta_forward_events_chi2(
            gamma,
            V0,
            events=strong_field_events,
            bridge_coeffs=bridge_coeffs,
            bridge_forward_calibration_rounds=bridge_forward_calibration_rounds,
            bridge_forward_calibration_reg=bridge_forward_calibration_reg,
        )
    total = chi2_cc + chi2_sn + chi2_bao + chi2_sf
    return {
        "total": total,
        "parts": {"chi2_cc": chi2_cc, "chi2_sn": chi2_sn, "chi2_bao": chi2_bao, "chi2_sf": chi2_sf, "mu0": mu0},
    }["total"], {
        "chi2_cc": chi2_cc,
        "chi2_sn": chi2_sn,
        "chi2_bao": chi2_bao,
        "chi2_sf": chi2_sf,
        "mu0": mu0,
    }


def _lcdm_total_chi2(
    omega_m: float,
    *,
    cc,
    sn,
    bao: BAOVectorData,
    sn_cov,
    alpha,
    omegas: Omegas | None = None,
    strong_field_events: list[dict] | None = None,
    strong_field_weight: float = 1.0,
    bridge_forward_calibration_rounds: int = 1,
    bridge_forward_calibration_reg: float = 0.1,
) -> tuple[float, dict]:
    # Shared distance grid for CC + SN + BAO.
    om = float(omega_m)
    if omegas is not None:
        om_r = max(0.0, float(omegas.Omega_r))
        ol = max(0.0, float(omegas.Omega_L))
    else:
        om_r = 0.0
        ol = 1.0 - om
    zmax = max(max(sn.z) if sn.z else 0.0, max(bao.z) if bao.z else 0.0)
    zmax = max(zmax, max(cc.z) if cc.z else 0.0)
    nint = 1200
    zs = [zmax * i / nint for i in range(nint + 1)]
    Es = [math.sqrt(max(om_r * (1.0 + z) ** 4 + om * (1.0 + z) ** 3 + ol, 1e-18)) for z in zs]
    invE = [1.0 / max(e, 1e-18) for e in Es]
    dc = [0.0] * len(zs)
    s = 0.0
    for i in range(1, len(zs)):
        dz = zs[i] - zs[i - 1]
        s += 0.5 * (invE[i] + invE[i - 1]) * dz
        dc[i] = s

    def interp(xq: float) -> float:
        if xq <= zs[0]:
            return dc[0]
        if xq >= zs[-1]:
            return dc[-1]
        j = bisect_left(zs, xq)
        x0, x1 = zs[j - 1], zs[j]
        y0, y1 = dc[j - 1], dc[j]
        t = (xq - x0) / max(x1 - x0, 1e-18)
        return y0 + t * (y1 - y0)

    def interp_E(xq: float) -> float:
        if xq <= zs[0]:
            return Es[0]
        if xq >= zs[-1]:
            return Es[-1]
        j = bisect_left(zs, xq)
        x0, x1 = zs[j - 1], zs[j]
        y0, y1 = Es[j - 1], Es[j]
        t = (xq - x0) / max(x1 - x0, 1e-18)
        return y0 + t * (y1 - y0)

    cc_m = [interp_E(z) for z in cc.z]
    chi2_cc = _chi2_gauss_diag(cc_m, cc.H, cc.sigma)
    sn_m = [5.0 * math.log10(max((1.0 + z) * interp(z), 1e-18)) for z in sn.z]
    chi2_sn, mu0 = sn_chi2(sn_m, sn.mu, sn.sigma, cov=sn_cov, fit_mu0=True)

    bao_out = []
    for z, o in zip(bao.z, bao.obs):
        Ez = interp_E(z)
        Dm = interp(z)
        if o == "DM_over_rs":
            bao_out.append(alpha * Dm)
        elif o == "DH_over_rs":
            bao_out.append(alpha / max(Ez, 1e-18))
        else:
            bao_out.append(alpha * ((z / max(Ez, 1e-18)) * (Dm * Dm)) ** (1.0 / 3.0))
    import numpy as np
    bao_m = np.asarray(bao_out, dtype=float)
    chi2_bao = bao_vector_chi2(bao_m, bao)
    chi2_sf = 0.0
    if strong_field_events is not None:
        chi2_sf = float(strong_field_weight) * _gr_forward_events_chi2(strong_field_events)
    total = chi2_cc + chi2_sn + chi2_bao + chi2_sf
    return total, {"chi2_cc": chi2_cc, "chi2_sn": chi2_sn, "chi2_bao": chi2_bao, "chi2_sf": chi2_sf, "mu0": mu0}


def _fit_lcdm(
    *,
    cc,
    sn,
    bao,
    sn_cov,
    alpha,
    strong_field_events: list[dict] | None = None,
    strong_field_weight: float = 1.0,
    omegas: Omegas | None = None,
    fit_omega_m: bool = True,
) -> tuple[float, dict]:
    if not bool(fit_omega_m):
        if omegas is None:
            raise ValueError("omegas must be provided when fit_omega_m=False")
        om_fix = max(1e-6, float(omegas.Omega_m))
        c2, parts = _lcdm_total_chi2(
            om_fix,
            cc=cc,
            sn=sn,
            bao=bao,
            sn_cov=sn_cov,
            alpha=alpha,
            omegas=omegas,
            strong_field_events=strong_field_events,
            strong_field_weight=strong_field_weight,
        )
        return c2, {"omega_m": om_fix, **parts}
    # 1D global scan + local refinement
    best = None
    for i in range(240):
        om = 0.01 + i * (0.89 / 239.0)
        c2, parts = _lcdm_total_chi2(
            om, cc=cc, sn=sn, bao=bao, sn_cov=sn_cov, alpha=alpha,
            omegas=omegas,
            strong_field_events=strong_field_events, strong_field_weight=strong_field_weight,
        )
        if best is None or c2 < best[0]:
            best = (c2, om, parts)
    assert best is not None
    c2b, omb, partsb = best
    lo = max(0.001, omb - 0.06)
    hi = min(0.99, omb + 0.06)
    for _ in range(4):
        for i in range(80):
            om = lo + i * (hi - lo) / 79.0
            c2, parts = _lcdm_total_chi2(
                om, cc=cc, sn=sn, bao=bao, sn_cov=sn_cov, alpha=alpha,
                omegas=omegas,
                strong_field_events=strong_field_events, strong_field_weight=strong_field_weight,
            )
            if c2 < c2b:
                c2b, omb, partsb = c2, om, parts
        span = (hi - lo) / 8.0
        lo = max(0.001, omb - span)
        hi = min(0.99, omb + span)
    return c2b, {"omega_m": omb, **partsb}


def _fit_gr(
    *,
    cc,
    sn,
    bao,
    sn_cov,
    alpha,
    omegas: Omegas | None = None,
    strong_field_events: list[dict] | None = None,
    strong_field_weight: float = 1.0,
) -> tuple[float, dict]:
    if omegas is None:
        raise ValueError("omegas must be provided for GR baseline")
    import numpy as np

    om_fix = max(1e-6, float(omegas.Omega_m))
    E = gr_E_fn(om_fix, omega_r=float(omegas.Omega_r))
    zmax = max(max(cc.z) if cc.z else 0.0, max(sn.z) if sn.z else 0.0, max(bao.z) if bao.z else 0.0)
    zs, dc = _build_dc_from_E(E, zmax=zmax, nint=3500)

    def interp(xq: float) -> float:
        if xq <= zs[0]:
            return dc[0]
        if xq >= zs[-1]:
            return dc[-1]
        j = bisect_left(zs, xq)
        x0, x1 = zs[j - 1], zs[j]
        y0, y1 = dc[j - 1], dc[j]
        t = (xq - x0) / max(x1 - x0, 1e-18)
        return y0 + t * (y1 - y0)

    cc_m = [E(z) for z in cc.z]
    chi2_cc = _chi2_gauss_diag(cc_m, cc.H, cc.sigma)
    sn_m = gr_predictions_for_sn(sn.z, omega_m=om_fix, omega_r=float(omegas.Omega_r))
    chi2_sn, mu0 = sn_chi2(sn_m, sn.mu, sn.sigma, cov=sn_cov, fit_mu0=True)
    bao_out = []
    for z, o in zip(bao.z, bao.obs):
        Ez = E(z)
        Dm = interp(z)
        if o == "DM_over_rs":
            bao_out.append(alpha * Dm)
        elif o == "DH_over_rs":
            bao_out.append(alpha / max(Ez, 1e-18))
        else:
            bao_out.append(alpha * ((z / max(Ez, 1e-18)) * (Dm * Dm)) ** (1.0 / 3.0))
    bao_m = np.asarray(bao_out, dtype=float)
    chi2_bao = bao_vector_chi2(bao_m, bao)
    chi2_sf = 0.0
    if strong_field_events is not None:
        chi2_sf = float(strong_field_weight) * _gr_forward_events_chi2(strong_field_events)
    total = chi2_cc + chi2_sn + chi2_bao + chi2_sf
    return total, {"omega_m": om_fix, "chi2_cc": chi2_cc, "chi2_sn": chi2_sn, "chi2_bao": chi2_bao, "chi2_sf": chi2_sf, "mu0": mu0}


def _fit_theta(
    *,
    cc,
    sn,
    bao,
    sn_cov,
    alpha,
    omegas: Omegas | None = None,
    seed: int,
    nsamples: int,
    strong_field_events: list[dict] | None = None,
    strong_field_weight: float = 1.0,
    bridge_coeffs: dict | None = None,
    bridge_forward_calibration_rounds: int = 1,
    bridge_forward_calibration_reg: float = 0.1,
    coarse_init: bool = True,
) -> tuple[float, dict]:
    rng = random.Random(seed)
    # Start around previous best known point + optional coarse multi-start.
    g = 0.03447
    v = 0.00647
    best_c2, best_parts = _theta_total_chi2(
        g, v, cc=cc, sn=sn, bao=bao, sn_cov=sn_cov, alpha=alpha,
        omegas=omegas,
        strong_field_events=strong_field_events, strong_field_weight=strong_field_weight,
        bridge_coeffs=bridge_coeffs,
        bridge_forward_calibration_rounds=bridge_forward_calibration_rounds,
        bridge_forward_calibration_reg=bridge_forward_calibration_reg,
    )
    if bool(coarse_init):
        gamma_grid = [0.0015, 0.0025, 0.004, 0.006, 0.01, 0.016, 0.025, 0.04, 0.07]
        v0_grid = [0.003, 0.0045, 0.0065, 0.0095, 0.014, 0.021, 0.032, 0.048]
        for cg in gamma_grid:
            for cv in v0_grid:
                c2_try, p_try = _theta_total_chi2(
                    cg, cv, cc=cc, sn=sn, bao=bao, sn_cov=sn_cov, alpha=alpha,
                    omegas=omegas,
                    strong_field_events=strong_field_events, strong_field_weight=strong_field_weight,
                    bridge_coeffs=bridge_coeffs,
                    bridge_forward_calibration_rounds=bridge_forward_calibration_rounds,
                    bridge_forward_calibration_reg=bridge_forward_calibration_reg,
                )
                if c2_try < best_c2:
                    g, v = float(cg), float(cv)
                    best_c2, best_parts = c2_try, p_try
    best = (best_c2, g, v, best_parts)

    step_g = 0.45
    step_v = 0.45
    accept = 0
    for i in range(max(200, nsamples)):
        cg = max(1e-3, min(10.0, g * math.exp(step_g * rng.gauss(0.0, 1.0))))
        cv = max(1e-4, min(10.0, v * math.exp(step_v * rng.gauss(0.0, 1.0))))
        c2, parts = _theta_total_chi2(
            cg, cv, cc=cc, sn=sn, bao=bao, sn_cov=sn_cov, alpha=alpha,
            omegas=omegas,
            strong_field_events=strong_field_events, strong_field_weight=strong_field_weight,
            bridge_coeffs=bridge_coeffs,
            bridge_forward_calibration_rounds=bridge_forward_calibration_rounds,
            bridge_forward_calibration_reg=bridge_forward_calibration_reg,
        )
        if (c2 <= best_c2) or (rng.random() < math.exp((best_c2 - c2) / 6.0)):
            g, v = cg, cv
            accept += 1
            if c2 < best_c2:
                best_c2, best_parts = c2, parts
                best = (c2, cg, cv, parts)
        if (i + 1) % 40 == 0:
            r = accept / 40.0
            if r > 0.45:
                step_g *= 1.12
                step_v *= 1.12
            elif r < 0.20:
                step_g /= 1.12
                step_v /= 1.12
            step_g = min(max(step_g, 0.04), 1.8)
            step_v = min(max(step_v, 0.04), 1.8)
            accept = 0
    c2, bg, bv, parts = best
    return c2, {"gamma": bg, "V0": bv, **parts}


def main(argv: Sequence[str] | None = None) -> int:
    args = parse_args(argv)
    cc = load_cc_dimless(args.cc)
    sn = load_sn_dimless(args.sn)
    sn_cov = None
    if args.sn_cov:
        import numpy as np

        sn_cov = np.loadtxt(args.sn_cov, dtype=float)
        if sn_cov.shape != (len(sn.z), len(sn.z)):
            raise ValueError(f"SN cov shape {sn_cov.shape} does not match SN length {len(sn.z)}")
    bao = load_bao_vector_from_desi(args.bao_mean, args.bao_cov)
    sf_events = None
    if bool(args.include_strong_field):
        sf_events = load_forward_event_catalog_json(args.forward_event_catalog_json)
    bridge_coeffs: dict | None = None
    bc_path = Path(args.bridge_coeffs_json)
    if bc_path.exists():
        bridge_coeffs = json.loads(bc_path.read_text(encoding="utf-8"))

    # alpha = c/(H0*rd) for *_over_rs BAO observables
    c_km_s = 299792.458
    alpha = c_km_s / (args.bao_h0_km_s_mpc * args.bao_rd_mpc)
    omegas = Omegas(
        H0=1.0,
        Omega_r=float(args.omega_r),
        Omega_b=float(args.omega_b),
        Omega_c=float(args.omega_c),
        Omega_L=float(args.omega_l),
    )

    chi2_theta, theta = _fit_theta(
        cc=cc,
        sn=sn,
        bao=bao,
        sn_cov=sn_cov,
        alpha=alpha,
        omegas=omegas,
        seed=args.seed,
        nsamples=args.theta_samples,
        strong_field_events=sf_events,
        strong_field_weight=float(args.strong_field_weight),
        bridge_coeffs=bridge_coeffs,
        bridge_forward_calibration_rounds=int(args.bridge_forward_calibration_rounds),
        bridge_forward_calibration_reg=float(args.bridge_forward_calibration_reg),
        coarse_init=bool(args.theta_coarse_init),
    )
    if str(args.comparison_mode) == "gr_only":
        chi2_lcdm, lcdm = _fit_gr(
            cc=cc,
            sn=sn,
            bao=bao,
            sn_cov=sn_cov,
            alpha=alpha,
            omegas=omegas,
            strong_field_events=sf_events,
            strong_field_weight=float(args.strong_field_weight),
        )
        baseline_label = "gr_only"
    else:
        chi2_lcdm, lcdm = _fit_lcdm(
            cc=cc,
            sn=sn,
            bao=bao,
            sn_cov=sn_cov,
            alpha=alpha,
            omegas=omegas,
            strong_field_events=sf_events,
            strong_field_weight=float(args.strong_field_weight),
            fit_omega_m=(str(args.comparison_mode) != "same_content"),
        )
        baseline_label = "lcdm_gr"
    bridge_diag: dict | None = None
    if bool(args.diagnose_bridge_dependency) and bool(args.include_strong_field) and int(args.bridge_forward_calibration_rounds) > 0:
        chi2_nb, parts_nb = _theta_total_chi2(
            float(theta["gamma"]),
            float(theta["V0"]),
            cc=cc,
            sn=sn,
            bao=bao,
            sn_cov=sn_cov,
            alpha=alpha,
            omegas=omegas,
            strong_field_events=sf_events,
            strong_field_weight=float(args.strong_field_weight),
            bridge_coeffs=bridge_coeffs,
            bridge_forward_calibration_rounds=0,
            bridge_forward_calibration_reg=float(args.bridge_forward_calibration_reg),
        )
        bridge_diag = {
            "chi2_total_no_bridge": float(chi2_nb),
            "chi2_sf_no_bridge": float(parts_nb["chi2_sf"]),
            "delta_chi2_total_no_bridge_minus_bridge": float(chi2_nb - chi2_theta),
            "delta_chi2_sf_no_bridge_minus_bridge": float(parts_nb["chi2_sf"] - theta["chi2_sf"]),
        }

    n_obs = len(cc.z) + len(sn.z) + len(bao.z) + (len(sf_events) if sf_events is not None else 0)
    # For gr_only we compare raw fit quality only, without model-selection penalty.
    # For the LCDM-fit / same_content branches we keep the existing parameter counts.
    if str(args.comparison_mode) == "gr_only":
        k_theta = 0
        k_lcdm = 0
    else:
        k_theta = 2
        # No latent bridge penalty in fixed-coeff mode:
        # compact-link coefficients are calibrated once outside this likelihood run.
        k_lcdm = 0 if (str(args.comparison_mode) == "same_content") else 1
    aic_theta = chi2_theta + 2 * k_theta
    aic_lcdm = chi2_lcdm + 2 * k_lcdm
    bic_theta = chi2_theta + k_theta * math.log(n_obs)
    bic_lcdm = chi2_lcdm + k_lcdm * math.log(n_obs)

    out_md = Path(args.output_md)
    out_md.parent.mkdir(parents=True, exist_ok=True)
    lines = []
    title = "Theta-Psi vs GR-only" if str(args.comparison_mode) == "gr_only" else "Theta-Psi vs flat LCDM/GR"
    lines.append(f"# Head-to-Head: {title}\n\n")
    lines.append(f"- N_obs: {n_obs}\n")
    lines.append(f"- SN covariance: {'full matrix' if sn_cov is not None else 'diagonal (sigma_mu)'}\n")
    lines.append(f"- BAO vector: {len(bao.z)} points (DESI mean+cov)\n\n")
    lines.append(f"- comparison_mode: {args.comparison_mode}\n")
    if str(args.comparison_mode) == "gr_only":
        lines.append("- model_selection_penalty: off (raw chi2 only)\n")
    lines.append(
        f"- shared_omegas: Or={float(args.omega_r):.6g}, "
        f"Ob={float(args.omega_b):.6g}, Oc={float(args.omega_c):.6g}, Ol={float(args.omega_l):.6g}\n\n"
    )
    lines.append(f"- include_strong_field: {'on' if bool(args.include_strong_field) else 'off'}\n")
    lines.append(f"- strong_field_weight: {float(args.strong_field_weight):.6g}\n")
    lines.append(f"- bridge_forward_calibration_rounds: {int(args.bridge_forward_calibration_rounds)}\n")
    lines.append(f"- bridge_forward_calibration_reg: {float(args.bridge_forward_calibration_reg):.6g}\n")
    lines.append(f"- bridge_coeffs_json: {args.bridge_coeffs_json}\n")
    lines.append("- bridge_inner_event_calibration: off (fixed global coefficients mode)\n")
    lines.append(f"- penalize_bridge_latent: {'on' if bool(args.penalize_bridge_latent) else 'off'}\n")
    lines.append("\n")
    lines.append("## Theta-Psi\n")
    lines.append(f"- gamma: {theta['gamma']:.8g}\n")
    lines.append(f"- V0: {theta['V0']:.8g}\n")
    lines.append(f"- mu0_SN(profile): {theta['mu0']:.8g}\n")
    lines.append(f"- chi2_total: {chi2_theta:.6f}\n")
    lines.append(f"- chi2_cc: {theta['chi2_cc']:.6f}\n")
    lines.append(f"- chi2_sn: {theta['chi2_sn']:.6f}\n")
    lines.append(f"- chi2_bao: {theta['chi2_bao']:.6f}\n")
    lines.append(f"- chi2_sf: {theta['chi2_sf']:.6f}\n")
    lines.append(f"- AIC: {aic_theta:.6f}\n")
    lines.append(f"- BIC: {bic_theta:.6f}\n\n")
    lines.append(f"## {baseline_label}\n")
    lines.append(f"- Omega_m: {lcdm['omega_m']:.8g}\n")
    lines.append(f"- mu0_SN(profile): {lcdm['mu0']:.8g}\n")
    lines.append(f"- chi2_total: {chi2_lcdm:.6f}\n")
    lines.append(f"- chi2_cc: {lcdm['chi2_cc']:.6f}\n")
    lines.append(f"- chi2_sn: {lcdm['chi2_sn']:.6f}\n")
    lines.append(f"- chi2_bao: {lcdm['chi2_bao']:.6f}\n")
    lines.append(f"- chi2_sf: {lcdm['chi2_sf']:.6f}\n")
    lines.append(f"- AIC: {aic_lcdm:.6f}\n")
    lines.append(f"- BIC: {bic_lcdm:.6f}\n\n")
    delta_label = "GR" if str(args.comparison_mode) == "gr_only" else "LCDM"
    lines.append(f"## Delta (Theta-Psi - {delta_label})\n")
    lines.append(f"- Delta chi2: {chi2_theta - chi2_lcdm:.6f}\n")
    lines.append(f"- Delta AIC: {aic_theta - aic_lcdm:.6f}\n")
    lines.append(f"- Delta BIC: {bic_theta - bic_lcdm:.6f}\n")
    if bridge_diag is not None:
        lines.append("\n## Bridge Dependency Diagnostic\n")
        lines.append(f"- chi2_total_no_bridge: {bridge_diag['chi2_total_no_bridge']:.6f}\n")
        lines.append(f"- chi2_sf_no_bridge: {bridge_diag['chi2_sf_no_bridge']:.6f}\n")
        lines.append(
            f"- delta_chi2_total_no_bridge_minus_bridge: "
            f"{bridge_diag['delta_chi2_total_no_bridge_minus_bridge']:.6f}\n"
        )
        lines.append(
            f"- delta_chi2_sf_no_bridge_minus_bridge: "
            f"{bridge_diag['delta_chi2_sf_no_bridge_minus_bridge']:.6f}\n"
        )
    out_md.write_text("".join(lines), encoding="utf-8")

    out_js = Path(args.output_json)
    out_js.parent.mkdir(parents=True, exist_ok=True)
    out_js.write_text(
        json.dumps(
            {
                "n_obs": n_obs,
                "sn_covariance": ("full" if sn_cov is not None else "diag"),
                "bao_points": len(bao.z),
                "comparison_mode": str(args.comparison_mode),
                "model_selection_penalty": (False if str(args.comparison_mode) == "gr_only" else True),
                "shared_omegas": {
                    "Omega_r": float(args.omega_r),
                    "Omega_b": float(args.omega_b),
                    "Omega_c": float(args.omega_c),
                    "Omega_L": float(args.omega_l),
                },
                "theta_psi": {
                    "gamma": theta["gamma"],
                    "V0": theta["V0"],
                    "mu0": theta["mu0"],
                    "chi2_total": chi2_theta,
                    "chi2_cc": theta["chi2_cc"],
                    "chi2_sn": theta["chi2_sn"],
                    "chi2_bao": theta["chi2_bao"],
                    "chi2_sf": theta["chi2_sf"],
                    "AIC": aic_theta,
                    "BIC": bic_theta,
                },
                baseline_label: {
                    "omega_m": lcdm["omega_m"],
                    "mu0": lcdm["mu0"],
                    "chi2_total": chi2_lcdm,
                    "chi2_cc": lcdm["chi2_cc"],
                    "chi2_sn": lcdm["chi2_sn"],
                    "chi2_bao": lcdm["chi2_bao"],
                    "chi2_sf": lcdm["chi2_sf"],
                    "AIC": aic_lcdm,
                    "BIC": bic_lcdm,
                },
                "delta_theta_minus_lcdm": {
                    "chi2": chi2_theta - chi2_lcdm,
                    "AIC": aic_theta - aic_lcdm,
                    "BIC": bic_theta - bic_lcdm,
                },
                "bridge_dependency_diagnostic": bridge_diag,
            },
            ensure_ascii=False,
            indent=2,
        )
        + "\n",
        encoding="utf-8",
    )
    print(f"Wrote {out_md}")
    print(f"Wrote {out_js}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
