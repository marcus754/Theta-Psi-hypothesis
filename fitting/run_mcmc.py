# -*- coding: utf-8 -*-
"""
Lightweight Metropolis sampler for (gamma, V0[, Q]) using dimensionless CC H(z) data.

CLI:
  python -m fitting.run_mcmc --data fitting/data/cc_demo_dimless.csv --steps 200
"""
from __future__ import annotations

import argparse
from dataclasses import dataclass
import json
from pathlib import Path
from typing import Sequence, Tuple
import math
import random
import statistics

from fitting.core_api import (
    BackgroundParams,
    CompactStarParams,
    Omegas,
    integrate_background_lnN,
    omegas_from_km_s_Mpc,
    solve_compact_star_profile,
)
from fitting.data_io import (
    load_cc_dimless,
    load_sn_dimless,
    load_bao_dimless,
    load_compact_targets_json,
    load_compact_star_targets_json,
    load_forward_targets_json,
    load_forward_event_catalog_json,
    CCData,
    SNData,
    BAOData,
)
from fitting.model import Hz_dimless, mu_dimless, DV_dimless
from fitting.likelihoods import loglike_gaussian
from fitting.likelihoods_compact import (
    CompactTargets,
    compact_metrics_model,
    loglike_compact,
    default_compact_targets,
)
from fitting.likelihoods_compact_star import (
    CompactStarTargets,
    compact_star_metrics_model,
    loglike_compact_star,
    default_compact_star_targets,
)
from fitting.strong_field_link import derive_compact_star_params_from_background
from fitting.forward_observables import (
    forward_observables_from_compact_metrics,
    calibrate_v1_coefficients,
)
from fitting.likelihoods_forward import ForwardTargets, loglike_forward
from fitting.forward_transfer import explicit_forward_observables
from fitting.likelihoods_forward_events import loglike_forward_events

_PRIOR_BOUNDS = {
    "gamma_min": 1e-3,
    "gamma_max": 10.0,
    "V0_min": 1e-4,
    "V0_max": 10.0,
    "Q_abs_max": 10.0,
}


@dataclass
class Params:
    gamma: float
    V0: float
    Q: float = 0.0

    def as_tuple(self) -> Tuple[float, float, float]:
        return (self.gamma, self.V0, self.Q)

def _stride_cc(data: CCData, stride: int) -> CCData:
    s = max(1, int(stride))
    if s == 1:
        return data
    return CCData(z=data.z[::s], H=data.H[::s], sigma=data.sigma[::s])


def _stride_sn(data: SNData | None, stride: int) -> SNData | None:
    if data is None:
        return None
    s = max(1, int(stride))
    if s == 1:
        return data
    return SNData(z=data.z[::s], mu=data.mu[::s], sigma=data.sigma[::s])


def _stride_bao(data: BAOData | None, stride: int) -> BAOData | None:
    if data is None:
        return None
    s = max(1, int(stride))
    if s == 1:
        return data
    return BAOData(z=data.z[::s], DV=data.DV[::s], sigma=data.sigma[::s])


def _resolve_compact_star_params(
    p: Params,
    *,
    link_model: bool,
    link_background: bool,
    theta_c: float,
    m2: float,
    lam: float,
    theta_c_from_gamma: float,
    m2_from_v0: float,
    lam_from_q2: float,
    bg_theta_quantile: float,
    bg_m2_scale: float,
    bg_lam_scale: float,
    forward_events_for_bridge: list[dict] | None = None,
    bridge_forward_calibration_rounds: int = 0,
    bridge_forward_calibration_reg: float = 0.0,
    om: Omegas | None = None,
) -> tuple[float, float, float]:
    if link_background:
        d = derive_compact_star_params_from_background(
            gamma=p.gamma,
            V0=p.V0,
            Q=p.Q,
            omegas=om,
            theta_quantile=bg_theta_quantile,
            m2_scale=bg_m2_scale,
            lam_scale=bg_lam_scale,
            forward_events=forward_events_for_bridge,
            calibrate_forward_rounds=bridge_forward_calibration_rounds,
            calibrate_forward_reg_strength=bridge_forward_calibration_reg,
            )
        return d["theta_c"], d["m2"], d["lam"]
    if not link_model:
        return theta_c, m2, lam
    # Simple phenomenological linking to make strong-field block responsive
    # to sampled model parameters.
    theta_c_eff = max(1e-6, theta_c_from_gamma * abs(p.gamma))
    m2_eff = max(1e-6, m2_from_v0 * p.V0)
    lam_eff = max(0.0, lam_from_q2 * (p.Q * p.Q) + lam)
    return theta_c_eff, m2_eff, lam_eff


def _weighted_sn_mu0(model_mu: list[float], data_mu: list[float], sigma_mu: list[float]) -> float:
    """Best-fit additive SN offset mu0 for fixed model shape."""
    wsum = 0.0
    num = 0.0
    for m, d, s in zip(model_mu, data_mu, sigma_mu):
        ss = max(float(s), 1e-18)
        w = 1.0 / (ss * ss)
        wsum += w
        num += (d - m) * w
    if wsum <= 0.0:
        return 0.0
    return num / wsum


def logprior(p: Params) -> float:
    # Broad, conservative priors in code units
    if not (_PRIOR_BOUNDS["gamma_min"] <= p.gamma <= _PRIOR_BOUNDS["gamma_max"]):
        return -math.inf
    if not (_PRIOR_BOUNDS["V0_min"] <= p.V0 <= _PRIOR_BOUNDS["V0_max"]):
        return -math.inf
    if not (abs(p.Q) <= _PRIOR_BOUNDS["Q_abs_max"]):
        return -math.inf
    # Log-uniform-ish preference via Jacobian terms (optional); keep flat for now
    return 0.0


def logposterior_cc(p: Params, cc: CCData, om: Omegas) -> float:
    lp = logprior(p)
    if not math.isfinite(lp):
        return -math.inf
    model = Hz_dimless(
        cc.z, gamma=p.gamma, V0=p.V0, Q=p.Q, omegas=om, a_min=1e-3, dN=0.02,
    )
    ll = loglike_gaussian(model, cc.H, cc.sigma)
    return lp + ll


def logposterior_combined(
    p: Params,
    *,
    om: Omegas,
    cc: CCData | None = None,
    sn: SNData | None = None,
    bao: BAOData | None = None,
    compact: CompactTargets | None = None,
    compact_theta_scale: float = 1.0,
    compact_weight: float = 1.0,
    compact_star: CompactStarTargets | None = None,
    compact_star_theta_c: float = 1.5,
    compact_star_m2: float = 1.0,
    compact_star_lam: float = 0.1,
    compact_star_r_max: float = 20.0,
    compact_star_dr: float = 0.02,
    compact_star_theta_scale: float = 1.0,
    compact_star_psi_grad_weight: float = 0.5,
    compact_star_theta_psi_coupling: float = 0.25,
    compact_star_refractive_strength: float = 0.5,
    compact_star_refractive_j_star: float = 1.0,
    compact_star_kappa_refr: float = 1.0,
    compact_star_weight: float = 1.0,
    compact_star_link_model: bool = False,
    compact_star_link_background: bool = False,
    compact_star_theta_c_from_gamma: float = 1.5,
    compact_star_m2_from_v0: float = 10.0,
    compact_star_lam_from_q2: float = 0.1,
    compact_star_bg_theta_quantile: float = 0.95,
    compact_star_bg_m2_scale: float = 10.0,
    compact_star_bg_lam_scale: float = 0.1,
    forward: ForwardTargets | None = None,
    forward_weight: float = 1.0,
    forward_mode: str = "explicit",
    forward_k_ring_m87: float = 3.0,
    forward_k_ring_sgra: float = 3.6,
    forward_k_echo_delay: float = 1.0,
    forward_ref_theta_c: float = 1.5,
    forward_ref_m2: float = 1.0,
    forward_ref_lam: float = 0.1,
    forward_ring_barrier_r_code: float = 3.0,
    forward_ring_band_index: int = 2,
    forward_ring_n_vac_tol: float = 0.01,
    forward_ring_r_min_code: float = 2.0,
    forward_echo_mass_msun: float = 60.0,
    forward_echo_barrier_r_code: float = 3.0,
    forward_echo_log_enhancement_coeff: float = 0.6366197723675814,
    forward_events: list[dict] | None = None,
    forward_events_weight: float = 1.0,
    sn_fit_mu0: bool = True,
    bridge_forward_calibration_rounds: int = 0,
    bridge_forward_calibration_reg: float = 0.0,
    **_legacy_kwargs,
) -> float:
    lp = logprior(p)
    if not math.isfinite(lp):
        return -math.inf
    total = lp
    if cc is not None:
        Hmodel = Hz_dimless(
            cc.z, gamma=p.gamma, V0=p.V0, Q=p.Q, omegas=om, a_min=1e-3, dN=0.02,
        )
        total += loglike_gaussian(Hmodel, cc.H, cc.sigma)
    if sn is not None:
        mu_model = mu_dimless(
            sn.z, gamma=p.gamma, V0=p.V0, Q=p.Q, omegas=om, a_min=1e-3, dN=0.02, mu0=0.0,
        )
        if sn_fit_mu0:
            mu0_fit = _weighted_sn_mu0(mu_model, sn.mu, sn.sigma)
            mu_model = [m + mu0_fit for m in mu_model]
        total += loglike_gaussian(mu_model, sn.mu, sn.sigma)
    if bao is not None:
        DV_model = DV_dimless(
            bao.z, gamma=p.gamma, V0=p.V0, Q=p.Q, omegas=om, a_min=1e-3, dN=0.02,
        )
        total += loglike_gaussian(DV_model, bao.DV, bao.sigma)
    if compact is not None:
        cmet = compact_metrics_model(
            gamma=p.gamma, V0=p.V0, Q=p.Q, omegas=om, theta_scale=compact_theta_scale
        )
        total += compact_weight * loglike_compact(cmet, compact)
    csmet = None
    csprof = None
    if compact_star is not None or forward is not None:
        ctheta, cm2, clam = _resolve_compact_star_params(
            p,
            link_model=compact_star_link_model,
            link_background=compact_star_link_background,
            theta_c=compact_star_theta_c,
            m2=compact_star_m2,
            lam=compact_star_lam,
            theta_c_from_gamma=compact_star_theta_c_from_gamma,
            m2_from_v0=compact_star_m2_from_v0,
            lam_from_q2=compact_star_lam_from_q2,
            bg_theta_quantile=compact_star_bg_theta_quantile,
            bg_m2_scale=compact_star_bg_m2_scale,
            bg_lam_scale=compact_star_bg_lam_scale,
            forward_events_for_bridge=forward_events,
            bridge_forward_calibration_rounds=bridge_forward_calibration_rounds,
            bridge_forward_calibration_reg=bridge_forward_calibration_reg,
            om=om,
        )
        cp = CompactStarParams(
            theta_c=ctheta,
            m2=cm2,
            lam=clam,
            r_max=compact_star_r_max,
            dr=compact_star_dr,
            theta_scale=compact_star_theta_scale,
            psi_grad_weight=compact_star_psi_grad_weight,
            theta_psi_coupling=compact_star_theta_psi_coupling,
            refractive_strength=compact_star_refractive_strength,
            refractive_j_star=compact_star_refractive_j_star,
            kappa_refr=compact_star_kappa_refr,
        )
        csprof = solve_compact_star_profile(cp)
        csmet = {
            "z_surface": float(csprof["z_surface"]),
            "delay_proxy": float(csprof["delay_proxy"]),
            "n_peak": float(csprof["n_peak"]),
        }
    if compact_star is not None and csmet is not None:
        total += compact_star_weight * loglike_compact_star(csmet, compact_star)
    if forward is not None and csmet is not None and csprof is not None:
        mode = str(forward_mode).strip().lower()
        if mode == "explicit":
            fobs = explicit_forward_observables(
                csprof,
                ring_barrier_r_code=forward_ring_barrier_r_code,
                ring_band_index=forward_ring_band_index,
                ring_n_vac_tol=forward_ring_n_vac_tol,
                ring_r_min_code=forward_ring_r_min_code,
                echo_mass_msun=forward_echo_mass_msun,
                echo_barrier_r_code=forward_echo_barrier_r_code,
                echo_log_enhancement_coeff=forward_echo_log_enhancement_coeff,
            )
        elif mode == "v1_calibrated":
            ref_prof = solve_compact_star_profile(
                CompactStarParams(
                    theta_c=forward_ref_theta_c,
                    m2=forward_ref_m2,
                lam=forward_ref_lam,
                r_max=compact_star_r_max,
                dr=compact_star_dr,
                theta_scale=compact_star_theta_scale,
                    psi_grad_weight=compact_star_psi_grad_weight,
                    theta_psi_coupling=compact_star_theta_psi_coupling,
                    refractive_strength=compact_star_refractive_strength,
                    refractive_j_star=compact_star_refractive_j_star,
                    kappa_refr=compact_star_kappa_refr,
                )
            )
            ref_met = {
                "z_surface": float(ref_prof["z_surface"]),
                "delay_proxy": float(ref_prof["delay_proxy"]),
                "n_peak": float(ref_prof["n_peak"]),
            }
            kk = calibrate_v1_coefficients(
                reference_metrics=ref_met,
                target_ring_m87_uas=forward.ring_m87_uas,
                target_ring_sgra_uas=forward.ring_sgra_uas,
                target_echo_delay_ms=forward.echo_delay_ms,
            )
            fobs = forward_observables_from_compact_metrics(csmet, **kk)
        else:
            fobs = forward_observables_from_compact_metrics(
                csmet,
                k_ring_m87=forward_k_ring_m87,
                k_ring_sgra=forward_k_ring_sgra,
                k_echo_delay=forward_k_echo_delay,
            )
        total += forward_weight * loglike_forward(fobs, forward)
    if forward_events is not None and csprof is not None:
        total += forward_events_weight * loglike_forward_events(
            csprof,
            forward_events,
            ring_barrier_r_code=forward_ring_barrier_r_code,
            ring_band_index=forward_ring_band_index,
            ring_n_vac_tol=forward_ring_n_vac_tol,
            ring_r_min_code=forward_ring_r_min_code,
            echo_barrier_r_code=forward_echo_barrier_r_code,
            echo_log_enhancement_coeff=forward_echo_log_enhancement_coeff,
        )
    return total


def safe_logposterior_combined(*, p: Params, **kwargs) -> float:
    """Return finite/inf-safe posterior; invalid trajectories map to -inf."""
    try:
        return logposterior_combined(p, **kwargs)
    except Exception:
        return -math.inf


def loglikes_components(
    p: Params,
    *,
    om: Omegas,
    cc: CCData | None = None,
    sn: SNData | None = None,
    bao: BAOData | None = None,
    compact: CompactTargets | None = None,
    compact_theta_scale: float = 1.0,
    compact_weight: float = 1.0,
    compact_star: CompactStarTargets | None = None,
    compact_star_theta_c: float = 1.5,
    compact_star_m2: float = 1.0,
    compact_star_lam: float = 0.1,
    compact_star_r_max: float = 20.0,
    compact_star_dr: float = 0.02,
    compact_star_theta_scale: float = 1.0,
    compact_star_psi_grad_weight: float = 0.5,
    compact_star_theta_psi_coupling: float = 0.25,
    compact_star_refractive_strength: float = 0.5,
    compact_star_refractive_j_star: float = 1.0,
    compact_star_kappa_refr: float = 1.0,
    compact_star_weight: float = 1.0,
    compact_star_link_model: bool = False,
    compact_star_link_background: bool = False,
    compact_star_theta_c_from_gamma: float = 1.5,
    compact_star_m2_from_v0: float = 10.0,
    compact_star_lam_from_q2: float = 0.1,
    compact_star_bg_theta_quantile: float = 0.95,
    compact_star_bg_m2_scale: float = 10.0,
    compact_star_bg_lam_scale: float = 0.1,
    forward: ForwardTargets | None = None,
    forward_weight: float = 1.0,
    forward_mode: str = "explicit",
    forward_k_ring_m87: float = 3.0,
    forward_k_ring_sgra: float = 3.6,
    forward_k_echo_delay: float = 1.0,
    forward_ref_theta_c: float = 1.5,
    forward_ref_m2: float = 1.0,
    forward_ref_lam: float = 0.1,
    forward_ring_barrier_r_code: float = 3.0,
    forward_ring_band_index: int = 2,
    forward_ring_n_vac_tol: float = 0.01,
    forward_ring_r_min_code: float = 2.0,
    forward_echo_mass_msun: float = 60.0,
    forward_echo_barrier_r_code: float = 3.0,
    forward_echo_log_enhancement_coeff: float = 0.6366197723675814,
    forward_events: list[dict] | None = None,
    forward_events_weight: float = 1.0,
    sn_fit_mu0: bool = True,
    bridge_forward_calibration_rounds: int = 0,
    bridge_forward_calibration_reg: float = 0.0,
    **_legacy_kwargs,
) -> dict:
    comps = {}
    if cc is not None:
        Hmodel = Hz_dimless(
            cc.z, gamma=p.gamma, V0=p.V0, Q=p.Q, omegas=om, a_min=1e-3, dN=0.02,
        )
        comps["cc"] = loglike_gaussian(Hmodel, cc.H, cc.sigma)
    if sn is not None:
        mu_model = mu_dimless(
            sn.z, gamma=p.gamma, V0=p.V0, Q=p.Q, omegas=om, a_min=1e-3, dN=0.02, mu0=0.0,
        )
        if sn_fit_mu0:
            mu0_fit = _weighted_sn_mu0(mu_model, sn.mu, sn.sigma)
            mu_model = [m + mu0_fit for m in mu_model]
        comps["sn"] = loglike_gaussian(mu_model, sn.mu, sn.sigma)
    if bao is not None:
        DV_model = DV_dimless(
            bao.z, gamma=p.gamma, V0=p.V0, Q=p.Q, omegas=om, a_min=1e-3, dN=0.02,
        )
        comps["bao"] = loglike_gaussian(DV_model, bao.DV, bao.sigma)
    if compact is not None:
        cmet = compact_metrics_model(
            gamma=p.gamma, V0=p.V0, Q=p.Q, omegas=om, theta_scale=compact_theta_scale
        )
        comps["compact"] = compact_weight * loglike_compact(cmet, compact)
    csmet = None
    csprof = None
    if compact_star is not None or forward is not None:
        ctheta, cm2, clam = _resolve_compact_star_params(
            p,
            link_model=compact_star_link_model,
            link_background=compact_star_link_background,
            theta_c=compact_star_theta_c,
            m2=compact_star_m2,
            lam=compact_star_lam,
            theta_c_from_gamma=compact_star_theta_c_from_gamma,
            m2_from_v0=compact_star_m2_from_v0,
            lam_from_q2=compact_star_lam_from_q2,
            bg_theta_quantile=compact_star_bg_theta_quantile,
            bg_m2_scale=compact_star_bg_m2_scale,
            bg_lam_scale=compact_star_bg_lam_scale,
            forward_events_for_bridge=forward_events,
            bridge_forward_calibration_rounds=bridge_forward_calibration_rounds,
            bridge_forward_calibration_reg=bridge_forward_calibration_reg,
            om=om,
        )
        cp = CompactStarParams(
            theta_c=ctheta,
            m2=cm2,
            lam=clam,
            r_max=compact_star_r_max,
            dr=compact_star_dr,
            theta_scale=compact_star_theta_scale,
            psi_grad_weight=compact_star_psi_grad_weight,
            theta_psi_coupling=compact_star_theta_psi_coupling,
            refractive_strength=compact_star_refractive_strength,
            refractive_j_star=compact_star_refractive_j_star,
            kappa_refr=compact_star_kappa_refr,
        )
        csprof = solve_compact_star_profile(cp)
        csmet = {
            "z_surface": float(csprof["z_surface"]),
            "delay_proxy": float(csprof["delay_proxy"]),
            "n_peak": float(csprof["n_peak"]),
        }
    if compact_star is not None and csmet is not None:
        comps["compact_star"] = compact_star_weight * loglike_compact_star(csmet, compact_star)
    if forward is not None and csmet is not None and csprof is not None:
        mode = str(forward_mode).strip().lower()
        if mode == "explicit":
            fobs = explicit_forward_observables(
                csprof,
                ring_barrier_r_code=forward_ring_barrier_r_code,
                ring_band_index=forward_ring_band_index,
                ring_n_vac_tol=forward_ring_n_vac_tol,
                ring_r_min_code=forward_ring_r_min_code,
                echo_mass_msun=forward_echo_mass_msun,
                echo_barrier_r_code=forward_echo_barrier_r_code,
                echo_log_enhancement_coeff=forward_echo_log_enhancement_coeff,
            )
        elif mode == "v1_calibrated":
            ref_prof = solve_compact_star_profile(
                CompactStarParams(
                    theta_c=forward_ref_theta_c,
                    m2=forward_ref_m2,
                lam=forward_ref_lam,
                r_max=compact_star_r_max,
                dr=compact_star_dr,
                theta_scale=compact_star_theta_scale,
                    psi_grad_weight=compact_star_psi_grad_weight,
                    theta_psi_coupling=compact_star_theta_psi_coupling,
                    refractive_strength=compact_star_refractive_strength,
                    refractive_j_star=compact_star_refractive_j_star,
                    kappa_refr=compact_star_kappa_refr,
                )
            )
            ref_met = {
                "z_surface": float(ref_prof["z_surface"]),
                "delay_proxy": float(ref_prof["delay_proxy"]),
                "n_peak": float(ref_prof["n_peak"]),
            }
            kk = calibrate_v1_coefficients(
                reference_metrics=ref_met,
                target_ring_m87_uas=forward.ring_m87_uas,
                target_ring_sgra_uas=forward.ring_sgra_uas,
                target_echo_delay_ms=forward.echo_delay_ms,
            )
            fobs = forward_observables_from_compact_metrics(csmet, **kk)
        else:
            fobs = forward_observables_from_compact_metrics(
                csmet,
                k_ring_m87=forward_k_ring_m87,
                k_ring_sgra=forward_k_ring_sgra,
                k_echo_delay=forward_k_echo_delay,
            )
        comps["forward"] = forward_weight * loglike_forward(fobs, forward)
    if forward_events is not None and csprof is not None:
        comps["forward_events"] = forward_events_weight * loglike_forward_events(
            csprof,
            forward_events,
            ring_barrier_r_code=forward_ring_barrier_r_code,
            ring_band_index=forward_ring_band_index,
            ring_n_vac_tol=forward_ring_n_vac_tol,
            ring_r_min_code=forward_ring_r_min_code,
            echo_barrier_r_code=forward_echo_barrier_r_code,
            echo_log_enhancement_coeff=forward_echo_log_enhancement_coeff,
        )
    return comps


def step_proposal(
    p: Params,
    rng: random.Random,
    scale: Tuple[float, float, float] = (0.1, 0.1, 0.0),
) -> Params:
    # Propose in log-space for positive parameters
    g = max(1e-6, p.gamma) * math.exp(scale[0] * rng.gauss(0.0, 1.0))
    v = max(1e-6, p.V0) * math.exp(scale[1] * rng.gauss(0.0, 1.0))
    q = p.Q + scale[2] * rng.gauss(0.0, 1.0)
    return Params(gamma=g, V0=v, Q=q)


def _looks_like_planck_H_values(cc: CCData) -> bool:
    if not cc.H:
        return False
    med_abs_h = statistics.median(abs(float(h)) for h in cc.H)
    return med_abs_h < 1e-20


def _validate_h0_and_cc_scale(cc: CCData, *, h0_km_s_mpc: float | None, h0_dimless: float) -> None:
    planck_like = _looks_like_planck_H_values(cc)
    if planck_like and h0_km_s_mpc is None and h0_dimless >= 1e-6:
        raise ValueError(
            "CC данные выглядят как планковские единицы (H ~ 1e-60), "
            "но задан demo-H0. Укажи --H0_km_s_Mpc (например 70)."
        )
    if (not planck_like) and h0_km_s_mpc is not None:
        raise ValueError(
            "CC данные не выглядят как планковские (H не ~1e-60), "
            "но указан --H0_km_s_Mpc. Убери --H0_km_s_Mpc для demo/dimless режима."
        )


def _initialize_valid_state(
    rng: random.Random,
    *,
    common_kwargs: dict,
    preferred: Params | None = None,
    max_tries: int = 512,
) -> tuple[Params, float]:
    seeds: list[Params] = []
    if preferred is not None:
        seeds.append(preferred)
    seeds.extend([
        Params(gamma=1.0, V0=0.1, Q=0.0),
        Params(gamma=0.3, V0=0.03, Q=0.0),
        Params(gamma=4.0, V0=1e-4, Q=0.0),
    ])
    best_p = seeds[0]
    best_lp = -math.inf
    for p in seeds:
        lp = safe_logposterior_combined(p=p, **common_kwargs)
        if math.isfinite(lp):
            return p, lp
        if lp > best_lp:
            best_lp = lp
            best_p = p
    for _ in range(max_tries):
        p = Params(
            gamma=10 ** rng.uniform(-3.0, 1.0),
            V0=10 ** rng.uniform(-4.0, 1.0),
            Q=0.0,
        )
        lp = safe_logposterior_combined(p=p, **common_kwargs)
        if math.isfinite(lp):
            return p, lp
        if lp > best_lp:
            best_lp = lp
            best_p = p
    return best_p, best_lp


def parse_args(argv: Sequence[str] | None = None) -> argparse.Namespace:
    p = argparse.ArgumentParser(description="Metropolis fit on CC/SN/BAO demo (dimless)")
    p.add_argument("--data", type=str, default=str(Path(__file__).resolve().parent / "data" / "cc_demo_dimless.csv"))
    p.add_argument("--H0_dimless", type=float, default=1.0, help="Dimensionless H0 for Omegas (use conversion helper for real H0)")
    p.add_argument("--H0_km_s_Mpc", type=float, default=None, help="If set, convert this H0 (km/s/Mpc) to code units and override --H0_dimless")
    p.add_argument("--sn", type=str, default=str(Path(__file__).resolve().parent / "data" / "sn_demo_dimless.csv"))
    p.add_argument("--bao", type=str, default=str(Path(__file__).resolve().parent / "data" / "bao_demo_dimless.csv"))
    p.add_argument("--cc_stride", type=int, default=1, help="Use every N-th CC point to speed up heavy runs")
    p.add_argument("--sn_stride", type=int, default=1, help="Use every N-th SN point to speed up heavy runs")
    p.add_argument("--bao_stride", type=int, default=1, help="Use every N-th BAO point to speed up heavy runs")
    p.add_argument("--sn_fit_mu0", action=argparse.BooleanOptionalAction, default=True, help="Fit additive SN magnitude offset mu0 analytically at each likelihood call")
    p.add_argument("--prop_scale_gamma", type=float, default=0.1)
    p.add_argument("--prop_scale_v0", type=float, default=0.1)
    p.add_argument("--prop_scale_q", type=float, default=0.0)
    p.add_argument("--prior_gamma_min", type=float, default=1e-3)
    p.add_argument("--prior_gamma_max", type=float, default=10.0)
    p.add_argument("--prior_V0_min", type=float, default=1e-4)
    p.add_argument("--prior_V0_max", type=float, default=10.0)
    p.add_argument("--prior_Q_abs_max", type=float, default=10.0)
    p.add_argument("--priors_json", type=str, default=None, help="Optional JSON with prior bounds overriding scalar prior flags")
    p.add_argument("--target_accept", type=float, default=0.30, help="Target acceptance for adaptive proposal during burn-in")
    p.add_argument("--adapt_interval", type=int, default=50, help="Adapt proposal scales every N steps during burn-in")
    p.add_argument("--adapt_until_frac", type=float, default=0.5, help="Fraction of chain used for adaptation (0..1)")
    p.add_argument("--steps", type=int, default=200)
    p.add_argument("--seed", type=int, default=7)
    p.add_argument("--output", type=str, default=None, help="Optional path to write a one-line summary")
    p.add_argument("--save_csv", type=str, default=None, help="Optional path to write CSV chain: step,gamma,V0,Q,logpost,accept")
    p.add_argument("--report_md", type=str, default=None, help="Optional path to write a Markdown report with best-fit and component loglikes")
    p.add_argument("--use_compact", action="store_true", help="Include compact-object component in loglike")
    p.add_argument("--compact_theta_scale", type=float, default=1.0, help="theta scale used in compact metrics")
    p.add_argument("--compact_target_z", type=float, default=None, help="Target z_surface")
    p.add_argument("--compact_sigma_z", type=float, default=None, help="Sigma for z_surface")
    p.add_argument("--compact_target_delay", type=float, default=None, help="Target delay")
    p.add_argument("--compact_sigma_delay", type=float, default=None, help="Sigma for delay")
    p.add_argument("--compact_target_n_peak", type=float, default=None, help="Target n_peak (demo if omitted)")
    p.add_argument("--compact_sigma_n_peak", type=float, default=None, help="Sigma for n_peak")
    p.add_argument("--compact_targets_json", type=str, default=None, help="Path to compact targets JSON (overrides scalar compact target args)")
    p.add_argument("--compact_weight", type=float, default=1.0, help="Weight for compact likelihood")
    p.add_argument("--use_compact_star", action="store_true", help="Include stationary compact-star component in loglike")
    p.add_argument("--compact_star_targets_json", type=str, default=None, help="Path to compact-star targets JSON")
    p.add_argument("--compact_star_theta_c", type=float, default=1.5)
    p.add_argument("--compact_star_m2", type=float, default=1.0)
    p.add_argument("--compact_star_lam", type=float, default=0.1)
    p.add_argument("--compact_star_link_model", action="store_true", help="Link compact-star profile params to sampled (gamma,V0,Q)")
    p.add_argument("--compact_star_link_background", action="store_true", help="Derive compact-star params from background trajectory stats")
    p.add_argument("--compact_star_theta_c_from_gamma", type=float, default=1.5)
    p.add_argument("--compact_star_m2_from_v0", type=float, default=10.0)
    p.add_argument("--compact_star_lam_from_q2", type=float, default=0.1)
    p.add_argument("--compact_star_bg_theta_quantile", type=float, default=0.95)
    p.add_argument("--compact_star_bg_m2_scale", type=float, default=10.0)
    p.add_argument("--compact_star_bg_lam_scale", type=float, default=0.1)
    p.add_argument("--compact_star_r_max", type=float, default=20.0)
    p.add_argument("--compact_star_dr", type=float, default=0.02)
    p.add_argument("--compact_star_theta_scale", type=float, default=1.0)
    p.add_argument("--compact_star_psi_grad_weight", type=float, default=0.5)
    p.add_argument("--compact_star_theta_psi_coupling", type=float, default=0.25)
    p.add_argument("--compact_star_refractive_strength", type=float, default=0.5)
    p.add_argument("--compact_star_refractive_j_star", type=float, default=1.0)
    p.add_argument("--compact_star_kappa_refr", type=float, default=1.0)
    p.add_argument("--compact_star_target_z", type=float, default=None)
    p.add_argument("--compact_star_sigma_z", type=float, default=None)
    p.add_argument("--compact_star_target_delay", type=float, default=None)
    p.add_argument("--compact_star_sigma_delay", type=float, default=None)
    p.add_argument("--compact_star_weight", type=float, default=1.0, help="Weight for stationary compact-star likelihood")
    p.add_argument("--use_forward", action="store_true", help="Include forward EHT/GW component in loglike")
    p.add_argument("--forward_targets_json", type=str, default="fitting/data/forward_observational_proxy.json")
    p.add_argument("--forward_weight", type=float, default=1.0)
    p.add_argument("--forward_mode", type=str, default="explicit", choices=["explicit", "v1_calibrated", "v1_fixed"])
    p.add_argument("--forward_k_ring_m87", type=float, default=3.0)
    p.add_argument("--forward_k_ring_sgra", type=float, default=3.6)
    p.add_argument("--forward_k_echo_delay", type=float, default=1.0)
    p.add_argument("--forward_ref_theta_c", type=float, default=1.5)
    p.add_argument("--forward_ref_m2", type=float, default=1.0)
    p.add_argument("--forward_ref_lam", type=float, default=0.1)
    p.add_argument("--forward_ring_barrier_r_code", type=float, default=3.0)
    p.add_argument("--forward_ring_band_index", type=int, default=2)
    p.add_argument("--forward_ring_n_vac_tol", type=float, default=0.01)
    p.add_argument("--forward_ring_r_min_code", type=float, default=2.0)
    p.add_argument("--forward_echo_mass_msun", type=float, default=60.0)
    p.add_argument("--forward_echo_barrier_r_code", type=float, default=3.0)
    p.add_argument("--forward_echo_log_enhancement_coeff", type=float, default=0.6366197723675814)
    p.add_argument("--use_forward_events", action="store_true", help="Include event-level forward component")
    p.add_argument("--forward_event_catalog_json", type=str, default="fitting/data/forward_event_catalog.json")
    p.add_argument("--forward_events_weight", type=float, default=1.0)
    p.add_argument(
        "--bridge_forward_calibration_rounds",
        type=int,
        default=0,
        help="Local rounds of compact-bridge calibration against forward events during background linking",
    )
    p.add_argument(
        "--bridge_forward_calibration_reg",
        type=float,
        default=0.0,
        help="Regularization strength to keep bridge calibration near background-prior params",
    )
    return p.parse_args(argv)


def main(argv: Sequence[str] | None = None) -> int:
    args = parse_args(argv)
    if args.priors_json:
        pdata = json.loads(Path(args.priors_json).read_text(encoding="utf-8"))
        for key, arg_name in (
            ("gamma_min", "prior_gamma_min"),
            ("gamma_max", "prior_gamma_max"),
            ("V0_min", "prior_V0_min"),
            ("V0_max", "prior_V0_max"),
            ("Q_abs_max", "prior_Q_abs_max"),
        ):
            if key in pdata:
                setattr(args, arg_name, float(pdata[key]))
    _PRIOR_BOUNDS["gamma_min"] = float(args.prior_gamma_min)
    _PRIOR_BOUNDS["gamma_max"] = float(args.prior_gamma_max)
    _PRIOR_BOUNDS["V0_min"] = float(args.prior_V0_min)
    _PRIOR_BOUNDS["V0_max"] = float(args.prior_V0_max)
    _PRIOR_BOUNDS["Q_abs_max"] = float(args.prior_Q_abs_max)
    rng = random.Random(args.seed)
    cc = load_cc_dimless(args.data)
    # Optional SN/BAO demo datasets
    sn = load_sn_dimless(args.sn) if args.sn else None
    bao = load_bao_dimless(args.bao) if args.bao else None
    cc = _stride_cc(cc, args.cc_stride)
    sn = _stride_sn(sn, args.sn_stride)
    bao = _stride_bao(bao, args.bao_stride)
    _validate_h0_and_cc_scale(cc, h0_km_s_mpc=args.H0_km_s_Mpc, h0_dimless=args.H0_dimless)
    if args.H0_km_s_Mpc is not None:
        om = omegas_from_km_s_Mpc(args.H0_km_s_Mpc)
    else:
        om = Omegas(H0=args.H0_dimless)
    compact = None
    if args.use_compact:
        if args.compact_targets_json:
            cj = load_compact_targets_json(args.compact_targets_json)
            compact = CompactTargets(
                mode=str(cj.get("mode", "gaussian")),
                z_surface=cj.get("z_surface"),
                sigma_z_surface=cj.get("sigma_z_surface"),
                delay_proxy=cj.get("delay_proxy"),
                sigma_delay_proxy=cj.get("sigma_delay_proxy"),
                n_peak=cj.get("n_peak"),
                sigma_n_peak=cj.get("sigma_n_peak"),
                z_surface_min=cj.get("z_surface_min"),
                sigma_z_surface_min=cj.get("sigma_z_surface_min"),
                delay_proxy_min=cj.get("delay_proxy_min"),
                sigma_delay_proxy_min=cj.get("sigma_delay_proxy_min"),
                n_peak_max=cj.get("n_peak_max"),
                sigma_n_peak_max=cj.get("sigma_n_peak_max"),
            )
        elif (
            args.compact_target_z is None
            or args.compact_sigma_z is None
            or args.compact_target_delay is None
            or args.compact_sigma_delay is None
            or args.compact_target_n_peak is None
            or args.compact_sigma_n_peak is None
        ):
            compact = default_compact_targets(omegas=om, theta_scale=args.compact_theta_scale)
        else:
            compact = CompactTargets(
                z_surface=float(args.compact_target_z),
                sigma_z_surface=float(args.compact_sigma_z),
                delay_proxy=float(args.compact_target_delay),
                sigma_delay_proxy=float(args.compact_sigma_delay),
                n_peak=float(args.compact_target_n_peak),
                sigma_n_peak=float(args.compact_sigma_n_peak),
            )
    compact_star = None
    if args.use_compact_star:
        if args.compact_star_targets_json:
            csj = load_compact_star_targets_json(args.compact_star_targets_json)
            compact_star = CompactStarTargets(
                z_surface=csj["z_surface"],
                sigma_z_surface=csj["sigma_z_surface"],
                delay_proxy=csj["delay_proxy"],
                sigma_delay_proxy=csj["sigma_delay_proxy"],
            )
        elif (
            args.compact_star_target_z is None
            or args.compact_star_sigma_z is None
            or args.compact_star_target_delay is None
            or args.compact_star_sigma_delay is None
        ):
            compact_star = default_compact_star_targets(
                theta_c=args.compact_star_theta_c,
                m2=args.compact_star_m2,
                lam=args.compact_star_lam,
                r_max=args.compact_star_r_max,
                dr=args.compact_star_dr,
                theta_scale=args.compact_star_theta_scale,
                psi_grad_weight=args.compact_star_psi_grad_weight,
                theta_psi_coupling=args.compact_star_theta_psi_coupling,
                refractive_strength=args.compact_star_refractive_strength,
                refractive_j_star=args.compact_star_refractive_j_star,
                kappa_refr=args.compact_star_kappa_refr,
            )
        else:
            compact_star = CompactStarTargets(
                z_surface=float(args.compact_star_target_z),
                sigma_z_surface=float(args.compact_star_sigma_z),
                delay_proxy=float(args.compact_star_target_delay),
                sigma_delay_proxy=float(args.compact_star_sigma_delay),
            )
    forward = None
    if args.use_forward:
        fj = load_forward_targets_json(args.forward_targets_json)
        forward = ForwardTargets(
            ring_m87_uas=fj["ring_m87_uas"],
            sigma_ring_m87_uas=fj["sigma_ring_m87_uas"],
            ring_sgra_uas=fj["ring_sgra_uas"],
            sigma_ring_sgra_uas=fj["sigma_ring_sgra_uas"],
            echo_delay_ms=fj.get("echo_delay_ms"),
            sigma_echo_delay_ms=fj.get("sigma_echo_delay_ms"),
        )
    forward_events = None
    if args.use_forward_events:
        forward_events = load_forward_event_catalog_json(args.forward_event_catalog_json)
    # Initialize near unity values.
    p0 = Params(gamma=1.0, V0=0.1, Q=0.0)
    common_kwargs = dict(
        om=om, cc=cc, sn=sn, bao=bao, compact=compact,
        compact_theta_scale=args.compact_theta_scale,
        compact_weight=args.compact_weight,
        compact_star=compact_star,
        compact_star_theta_c=args.compact_star_theta_c,
        compact_star_m2=args.compact_star_m2,
        compact_star_lam=args.compact_star_lam,
        compact_star_r_max=args.compact_star_r_max,
        compact_star_dr=args.compact_star_dr,
        compact_star_theta_scale=args.compact_star_theta_scale,
        compact_star_psi_grad_weight=args.compact_star_psi_grad_weight,
        compact_star_theta_psi_coupling=args.compact_star_theta_psi_coupling,
        compact_star_refractive_strength=args.compact_star_refractive_strength,
        compact_star_refractive_j_star=args.compact_star_refractive_j_star,
        compact_star_kappa_refr=args.compact_star_kappa_refr,
        compact_star_weight=args.compact_star_weight,
        compact_star_link_model=args.compact_star_link_model,
        compact_star_link_background=args.compact_star_link_background,
        compact_star_theta_c_from_gamma=args.compact_star_theta_c_from_gamma,
        compact_star_m2_from_v0=args.compact_star_m2_from_v0,
        compact_star_lam_from_q2=args.compact_star_lam_from_q2,
        compact_star_bg_theta_quantile=args.compact_star_bg_theta_quantile,
        compact_star_bg_m2_scale=args.compact_star_bg_m2_scale,
        compact_star_bg_lam_scale=args.compact_star_bg_lam_scale,
        forward=forward,
        forward_weight=args.forward_weight,
        forward_mode=args.forward_mode,
        forward_k_ring_m87=args.forward_k_ring_m87,
        forward_k_ring_sgra=args.forward_k_ring_sgra,
        forward_k_echo_delay=args.forward_k_echo_delay,
        forward_ref_theta_c=args.forward_ref_theta_c,
        forward_ref_m2=args.forward_ref_m2,
        forward_ref_lam=args.forward_ref_lam,
        forward_ring_barrier_r_code=args.forward_ring_barrier_r_code,
        forward_ring_band_index=args.forward_ring_band_index,
        forward_ring_n_vac_tol=args.forward_ring_n_vac_tol,
        forward_ring_r_min_code=args.forward_ring_r_min_code,
        forward_echo_mass_msun=args.forward_echo_mass_msun,
        forward_echo_barrier_r_code=args.forward_echo_barrier_r_code,
        forward_echo_log_enhancement_coeff=args.forward_echo_log_enhancement_coeff,
        forward_events=forward_events,
        forward_events_weight=args.forward_events_weight,
        bridge_forward_calibration_rounds=int(args.bridge_forward_calibration_rounds),
        bridge_forward_calibration_reg=float(args.bridge_forward_calibration_reg),
        sn_fit_mu0=bool(args.sn_fit_mu0),
    )
    p, logp = _initialize_valid_state(rng, common_kwargs=common_kwargs, preferred=p0)
    if not math.isfinite(logp):
        summary = (
            "posterior_loglike=-inf params={gamma:1, V0:0.1, Q:0} accept=0.00 "
            "best={logp:-inf, gamma:1, V0:0.1, Q:0} error={init_failed:true}"
        )
        print(summary)
        if args.output:
            outp = Path(args.output)
            outp.parent.mkdir(parents=True, exist_ok=True)
            outp.write_text(summary + "\n", encoding="utf-8")
        return 1
    accept = 0
    best = (logp, p)
    csv_path = Path(args.save_csv) if args.save_csv else None
    if csv_path is not None:
        csv_path.parent.mkdir(parents=True, exist_ok=True)
        with csv_path.open("w", newline="", encoding="utf-8") as fh:
            fh.write("step,gamma,V0,Q,logpost,accept\n")
    s_gamma = float(args.prop_scale_gamma)
    s_v0 = float(args.prop_scale_v0)
    s_q = float(args.prop_scale_q)
    w_accept = 0
    adapt_until = int(max(0, min(args.steps, round(float(args.adapt_until_frac) * args.steps))))
    for i in range(args.steps):
        cand = step_proposal(
            p,
            rng,
            scale=(s_gamma, s_v0, s_q),
        )
        logc = safe_logposterior_combined(p=cand, **common_kwargs)
        acc = 0
        if math.isfinite(logc) and (logc - logp >= math.log(rng.random())):
            p, logp = cand, logc
            accept += 1
            w_accept += 1
            acc = 1
            if logp > best[0]:
                best = (logp, p)
        # Lightweight Robbins-Monro style adaptation in burn-in window.
        if (
            args.adapt_interval > 0
            and (i + 1) <= adapt_until
            and ((i + 1) % int(args.adapt_interval) == 0)
        ):
            w_rate = w_accept / float(args.adapt_interval)
            if w_rate > float(args.target_accept) + 0.05:
                s_gamma *= 1.15
                s_v0 *= 1.15
                if s_q > 0.0:
                    s_q *= 1.10
            elif w_rate < float(args.target_accept) - 0.05:
                s_gamma /= 1.15
                s_v0 /= 1.15
                if s_q > 0.0:
                    s_q /= 1.10
            s_gamma = min(max(s_gamma, 1e-4), 2.0)
            s_v0 = min(max(s_v0, 1e-4), 2.0)
            s_q = min(max(s_q, 0.0), 2.0)
            w_accept = 0
        if csv_path is not None:
            with csv_path.open("a", newline="", encoding="utf-8") as fh:
                fh.write(f"{i},{p.gamma},{p.V0},{p.Q},{logp},{acc}\n")
    acc_rate = accept / max(args.steps, 1)
    summary = (
        f"posterior_loglike={logp:.3f} "
        f"params={{gamma:{p.gamma:.4g}, V0:{p.V0:.4g}, Q:{p.Q:.4g}}} "
        f"accept={acc_rate:.2f} "
        f"scales={{sg:{s_gamma:.4g}, sv:{s_v0:.4g}, sq:{s_q:.4g}}} "
        f"best={{logp:{best[0]:.3f}, gamma:{best[1].gamma:.4g}, V0:{best[1].V0:.4g}, Q:{best[1].Q:.4g}}}"
    )
    print(summary)
    if args.output:
        outp = Path(args.output)
        outp.parent.mkdir(parents=True, exist_ok=True)
        outp.write_text(summary + "\n", encoding="utf-8")
    if args.report_md:
        bp = best[1]
        try:
            comps = loglikes_components(
                bp, om=om, cc=cc, sn=sn, bao=bao, compact=compact,
                compact_theta_scale=args.compact_theta_scale,
                compact_weight=args.compact_weight,
                compact_star=compact_star,
                compact_star_theta_c=args.compact_star_theta_c,
                compact_star_m2=args.compact_star_m2,
                compact_star_lam=args.compact_star_lam,
                compact_star_r_max=args.compact_star_r_max,
                compact_star_dr=args.compact_star_dr,
                compact_star_theta_scale=args.compact_star_theta_scale,
                compact_star_psi_grad_weight=args.compact_star_psi_grad_weight,
                compact_star_theta_psi_coupling=args.compact_star_theta_psi_coupling,
                compact_star_refractive_strength=args.compact_star_refractive_strength,
                compact_star_refractive_j_star=args.compact_star_refractive_j_star,
                compact_star_kappa_refr=args.compact_star_kappa_refr,
                compact_star_weight=args.compact_star_weight,
                compact_star_link_model=args.compact_star_link_model,
                compact_star_link_background=args.compact_star_link_background,
                compact_star_theta_c_from_gamma=args.compact_star_theta_c_from_gamma,
                compact_star_m2_from_v0=args.compact_star_m2_from_v0,
                compact_star_lam_from_q2=args.compact_star_lam_from_q2,
                compact_star_bg_theta_quantile=args.compact_star_bg_theta_quantile,
                compact_star_bg_m2_scale=args.compact_star_bg_m2_scale,
                compact_star_bg_lam_scale=args.compact_star_bg_lam_scale,
                forward=forward,
                forward_weight=args.forward_weight,
                forward_mode=args.forward_mode,
                forward_k_ring_m87=args.forward_k_ring_m87,
                forward_k_ring_sgra=args.forward_k_ring_sgra,
                forward_k_echo_delay=args.forward_k_echo_delay,
                forward_ref_theta_c=args.forward_ref_theta_c,
                forward_ref_m2=args.forward_ref_m2,
                forward_ref_lam=args.forward_ref_lam,
                forward_ring_barrier_r_code=args.forward_ring_barrier_r_code,
                forward_ring_band_index=args.forward_ring_band_index,
                forward_ring_n_vac_tol=args.forward_ring_n_vac_tol,
                forward_ring_r_min_code=args.forward_ring_r_min_code,
                forward_echo_mass_msun=args.forward_echo_mass_msun,
                forward_echo_barrier_r_code=args.forward_echo_barrier_r_code,
                forward_echo_log_enhancement_coeff=args.forward_echo_log_enhancement_coeff,
                forward_events=forward_events,
                forward_events_weight=args.forward_events_weight,
                bridge_forward_calibration_rounds=int(args.bridge_forward_calibration_rounds),
                bridge_forward_calibration_reg=float(args.bridge_forward_calibration_reg),
                sn_fit_mu0=bool(args.sn_fit_mu0),
            )
        except Exception:
            comps = {}
        lines = []
        lines.append(f"# Fit Report\n")
        lines.append(f"\n")
        lines.append(f"- Steps: {args.steps}\n")
        lines.append(f"- Seed: {args.seed}\n")
        if args.H0_km_s_Mpc is not None:
            lines.append(f"- H0: {args.H0_km_s_Mpc} km/s/Mpc → dimless {om.H0:.3e}\n")
        else:
            lines.append(f"- H0 (dimless): {om.H0:.3e}\n")
        lines.append(f"- Accept rate: {acc_rate:.2f}\n")
        lines.append(f"- bridge_forward_calibration_rounds: {int(args.bridge_forward_calibration_rounds)}\n")
        lines.append(f"- bridge_forward_calibration_reg: {float(args.bridge_forward_calibration_reg):.6g}\n")
        lines.append(f"\n## Best parameters\n")
        lines.append(f"- logpost: {best[0]:.3f}\n")
        lines.append(f"- gamma: {bp.gamma:.6g}\n")
        lines.append(f"- V0: {bp.V0:.6g}\n")
        lines.append(f"- Q: {bp.Q:.6g}\n")
        lines.append(f"\n## Component loglikes\n")
        if cc is not None:
            lines.append(f"- CC (H): {comps.get('cc', 0.0):.3f}  (N={len(cc.z)})\n")
        if sn is not None:
            lines.append(f"- SN (mu): {comps.get('sn', 0.0):.3f}  (N={len(sn.z)})\n")
        if bao is not None:
            lines.append(f"- BAO (DV): {comps.get('bao', 0.0):.3f}  (N={len(bao.z)})\n")
        if compact is not None:
            lines.append(f"- Compact: {comps.get('compact', 0.0):.3f}\n")
        if compact_star is not None:
            lines.append(f"- CompactStar (stationary): {comps.get('compact_star', 0.0):.3f}\n")
        if forward is not None:
            lines.append(f"- Forward (EHT/GW): {comps.get('forward', 0.0):.3f}\n")
        if forward_events is not None:
            lines.append(f"- ForwardEvents (event-level): {comps.get('forward_events', 0.0):.3f}\n")
        lines.append(f"\n## Component chi2\n")
        for k, v in comps.items():
            lines.append(f"- {k}: {-2.0*float(v):.3f}\n")
        if comps:
            # Most negative loglike -> largest chi2 pull
            worst = min(comps.items(), key=lambda kv: float(kv[1]))
            lines.append(f"\n## Dominant Tension\n")
            lines.append(f"- component: {worst[0]}\n")
            lines.append(f"- loglike: {float(worst[1]):.3f}\n")
            lines.append(f"- chi2: {-2.0*float(worst[1]):.3f}\n")
        rep = Path(args.report_md)
        rep.parent.mkdir(parents=True, exist_ok=True)
        rep.write_text("".join(lines), encoding="utf-8")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
