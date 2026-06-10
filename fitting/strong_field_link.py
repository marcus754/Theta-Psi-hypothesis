# -*- coding: utf-8 -*-
"""
Background-driven linking from (gamma, V0, Q) to compact-star parameters.

This is a provisional mapping layer for inference and scans.
"""
from __future__ import annotations

from typing import Dict, List
import math
import random

from fitting.core_api import (
    BackgroundParams,
    CompactStarParams,
    Omegas,
    integrate_background_lnN,
    solve_compact_star_profile,
)
from fitting.likelihoods_forward_events import loglike_forward_events


def _quantile(xs: List[float], q: float) -> float:
    if not xs:
        return 0.0
    ys = sorted(xs)
    qq = min(1.0, max(0.0, float(q)))
    pos = qq * (len(ys) - 1)
    i0 = int(math.floor(pos))
    i1 = int(math.ceil(pos))
    if i0 == i1:
        return ys[i0]
    w = pos - i0
    return (1.0 - w) * ys[i0] + w * ys[i1]


def resolve_bridge_scales(
    *,
    gamma: float,
    V0: float,
    bridge_coeffs: dict | None,
) -> dict:
    """
    Resolve (theta_c_scale, m2_scale, lam_scale) from bridge config.

    Supported payloads:
    - constant:
      {"theta_c_scale":..., "m2_scale":..., "lam_scale":...}
    - powerlaw2d:
      {
        "bridge_model":"powerlaw2d",
        "gamma_ref":..., "V0_ref":...,
        "theta_c": {"a0":..., "ag":..., "av":...},
        "m2": {"a0":..., "ag":..., "av":...},
        "lam": {"a0":..., "ag":..., "av":...}
      }
    """
    cfg = bridge_coeffs or {}
    model = str(cfg.get("bridge_model", "constant")).strip().lower()
    if model != "powerlaw2d":
        return {
            "theta_c_scale": float(cfg.get("theta_c_scale", 1.0)),
            "m2_scale": float(cfg.get("m2_scale", 10.0)),
            "lam_scale": float(cfg.get("lam_scale", 0.1)),
        }

    gref = max(1e-12, float(cfg.get("gamma_ref", 1.0)))
    vref = max(1e-12, float(cfg.get("V0_ref", 1.0)))
    lg = math.log(max(float(gamma), 1e-12) / gref)
    lv = math.log(max(float(V0), 1e-12) / vref)

    def _scale(block: dict, default: float) -> float:
        a0 = math.log(max(float(block.get("a0", default)), 1e-12))
        ag = float(block.get("ag", 0.0))
        av = float(block.get("av", 0.0))
        x = a0 + ag * lg + av * lv
        x = min(30.0, max(-30.0, x))
        return float(math.exp(x))

    return {
        "theta_c_scale": _scale(dict(cfg.get("theta_c", {})), 1.0),
        "m2_scale": _scale(dict(cfg.get("m2", {})), 10.0),
        "lam_scale": _scale(dict(cfg.get("lam", {})), 0.1),
    }


def derive_compact_star_params_from_background(
    *,
    gamma: float,
    V0: float,
    Q: float = 0.0,
    omegas: Omegas | None = None,
    theta_quantile: float = 0.95,
    theta_c_scale: float = 1.0,
    m2_scale: float = 10.0,
    lam_scale: float = 0.1,
    a0: float = 1e-3,
    x0: float = 0.1,
    xdot0: float = 0.0,
    psi0: float = 0.05,
    psidot0: float = 0.0,
    dN: float = 0.02,
    forward_events: list[dict] | None = None,
    calibrate_forward_rounds: int = 0,
    calibrate_theta_c: bool = False,
    calibrate_forward_reg_strength: float = 0.0,
    **_legacy_kwargs,
) -> Dict[str, float]:
    """
    Map background trajectory statistics to compact-star profile parameters.

    Strong-field sensitivity is derived from the direct two-field refractive
    invariant track (`J_refr`, `I_grad`) and the reduced time field `theta`.
    """
    if omegas is None:
        omegas = Omegas()
    p = BackgroundParams(gamma=gamma, V0=V0, Q=Q, omegas=omegas)
    out = integrate_background_lnN(
        a0, x0, xdot0, psi0, psidot0, N1=0.0, p=p, hN=dN,
    )
    theta = [abs(float(t)) for t in out["theta"]]
    hvals = [max(abs(float(h)), 1e-18) for h in out["H"]]
    dtheta_dN = [float(td) / h for td, h in zip(out["theta_dot"], hvals)]
    j_refr = [max(0.0, float(v)) for v in out.get("J_refr", [])]
    i_grad = [max(0.0, float(v)) for v in out.get("I_grad", [])]

    # Physical scale derivation instead of quantile-based matching.
    # From docs/17: rc^2 ~ gamma^2 / V0 (approx scale of the bridge).
    # We use theta_c_scale as a fixed theoretical multiplier, not a fit knob.
    rc2 = max(1e-12, (gamma * gamma) / max(abs(V0), 1e-12))
    rc = math.sqrt(rc2)
    
    # theta_c (value at center) is pinned to the background theta_0 (1.0 in normalized units)
    # plus a correction from the potential depth.
    theta_c = max(1e-6, float(theta_c_scale) * (1.0 + 0.1 * rc))
    
    # m2 (mass scale) is derived from the field stiffness
    m2 = max(1e-6, float(m2_scale) * (1.0 / rc2 + abs(V0)))
    
    # lam (non-linearity) is derived from the phase charge and potential
    lam = max(0.0, float(lam_scale) * (abs(V0) + 0.01 * Q * Q))

    out_map = {
        "theta_c": theta_c,
        "m2": m2,
        "lam": lam,
        "rc_derived": rc,
    }

    # Optional local calibration of compact-profile parameters directly against
    # forward event-level residuals with regularization to the background prior.
    if forward_events and int(calibrate_forward_rounds) > 0:
        tc0 = float(theta_c)
        m20 = float(m2)
        lam0 = float(lam)
        reg = max(0.0, float(calibrate_forward_reg_strength))

        def _chi2_forward(tc_try: float, m2_try: float, lam_try: float) -> float:
            prof = solve_compact_star_profile(
                CompactStarParams(
                    theta_c=max(1e-6, tc_try),
                    m2=max(1e-6, m2_try),
                    lam=max(0.0, lam_try),
                    r_max=20.0,
                    dr=0.02,
                    theta_scale=1.0,
                )
            )
            chi2_f = -2.0 * float(loglike_forward_events(prof, forward_events))
            if reg <= 0.0:
                return chi2_f
            d_tc = math.log(max(tc_try, 1e-12) / max(tc0, 1e-12))
            d_m2 = math.log(max(m2_try, 1e-12) / max(m20, 1e-12))
            d_lm = math.log(max(lam_try + 1e-12, 1e-12) / max(lam0 + 1e-12, 1e-12))
            return chi2_f + reg * (d_tc * d_tc + d_m2 * d_m2 + d_lm * d_lm)

        # Optimize in unconstrained coordinates:
        #   a = ln theta_c, b = ln m2, c = ln(lam + eps)
        a = math.log(max(tc0, 1e-12))
        b = math.log(max(m20, 1e-12))
        c = math.log(max(lam0 + 1e-12, 1e-12))

        def _decode(a1: float, b1: float, c1: float) -> tuple[float, float, float]:
            # Numerical guard: local optimizer may propose extreme values.
            # Clamp log-parameters to keep exp() in a safe range.
            a1 = min(50.0, max(-50.0, float(a1)))
            b1 = min(50.0, max(-50.0, float(b1)))
            c1 = min(50.0, max(-50.0, float(c1)))
            return math.exp(a1), math.exp(b1), max(0.0, math.exp(c1) - 1e-12)

        def _obj(a1: float, b1: float, c1: float) -> float:
            tc_try, m2_try, lam_try = _decode(a1, b1, c1)
            if not calibrate_theta_c:
                tc_try = tc0
            return _chi2_forward(tc_try, m2_try, lam_try)

        def _grad(a1: float, b1: float, c1: float) -> tuple[float, float, float]:
            h = 5e-3
            f0 = _obj(a1, b1, c1)
            ga = (_obj(a1 + h, b1, c1) - _obj(a1 - h, b1, c1)) / (2.0 * h)
            gb = (_obj(a1, b1 + h, c1) - _obj(a1, b1 - h, c1)) / (2.0 * h)
            gc = (_obj(a1, b1, c1 + h) - _obj(a1, b1, c1 - h)) / (2.0 * h)
            return ga, gb, gc

        lr = 0.15
        rounds = max(1, int(calibrate_forward_rounds))
        best_val = _obj(a, b, c)
        best_abc = (a, b, c)
        for _ in range(rounds * 12):
            ga, gb, gc = _grad(a, b, c)
            # Backtracking line-search
            improved = False
            for _ls in range(8):
                an = a - lr * ga
                bn = b - lr * gb
                cn = c - lr * gc
                fn = _obj(an, bn, cn)
                if fn < best_val:
                    a, b, c = an, bn, cn
                    best_val = fn
                    best_abc = (a, b, c)
                    improved = True
                    lr *= 1.05
                    break
                lr *= 0.5
            if not improved:
                # Small exploratory coordinate move
                for da, db, dc in ((0.01, 0, 0), (-0.01, 0, 0), (0, 0.01, 0), (0, -0.01, 0), (0, 0, 0.01), (0, 0, -0.01)):
                    fn = _obj(a + da, b + db, c + dc)
                    if fn < best_val:
                        a, b, c = a + da, b + db, c + dc
                        best_val = fn
                        best_abc = (a, b, c)
                        improved = True
                        break
            if lr < 1e-5:
                break

        tc, mm, ll = _decode(*best_abc)
        if not calibrate_theta_c:
            tc = tc0
        out_map["theta_c"] = float(tc)
        out_map["m2"] = float(mm)
        out_map["lam"] = float(ll)
        out_map["bridge_forward_chi2"] = float(best_val)

    return out_map


def calibrate_global_bridge_coeffs(
    *,
    gamma: float,
    V0: float,
    events: list[dict],
    Q: float = 0.0,
    rounds: int = 32,
    reg_strength: float = 0.2,
    theta_c_scale0: float = 1.0,
    m2_scale0: float = 10.0,
    lam_scale0: float = 0.1,
    **_legacy_kwargs,
) -> Dict[str, float]:
    """
    One-shot global bridge calibration.
    Returns fixed coefficients to be reused in likelihood calls.
    """
    tc0 = max(1e-6, float(theta_c_scale0))
    m20 = max(1e-6, float(m2_scale0))
    lm0 = max(1e-9, float(lam_scale0))
    reg = max(0.0, float(reg_strength))

    def _decode(a: float, b: float, c: float) -> tuple[float, float, float]:
        a = min(30.0, max(-30.0, float(a)))
        b = min(30.0, max(-30.0, float(b)))
        c = min(30.0, max(-30.0, float(c)))
        return math.exp(a), math.exp(b), math.exp(c)

    def _obj(a: float, b: float, c: float) -> float:
        tcs, m2s, lms = _decode(a, b, c)
        d = derive_compact_star_params_from_background(
            gamma=gamma,
            V0=V0,
            Q=Q,
            theta_c_scale=tcs,
            m2_scale=m2s,
            lam_scale=lms,
            # Explicitly disable event-level inner fit; this is global-only mode.
            forward_events=None,
            calibrate_forward_rounds=0,
        )
        prof = solve_compact_star_profile(
            CompactStarParams(
                theta_c=max(1e-6, float(d["theta_c"])),
                m2=max(1e-6, float(d["m2"])),
                lam=max(0.0, float(d["lam"])),
                r_max=20.0,
                dr=0.02,
                theta_scale=1.0,
            )
        )
        chi2_f = -2.0 * float(loglike_forward_events(prof, events))
        if reg <= 0.0:
            return chi2_f
        da = math.log(max(tcs, 1e-12) / max(tc0, 1e-12))
        db = math.log(max(m2s, 1e-12) / max(m20, 1e-12))
        dc = math.log(max(lms, 1e-12) / max(lm0, 1e-12))
        return chi2_f + reg * (da * da + db * db + dc * dc)

    a = math.log(tc0)
    b = math.log(m20)
    c = math.log(lm0)
    best = _obj(a, b, c)
    best_abc = (a, b, c)
    # Coarse pre-scan to avoid poor local starts.
    tc_grid = [0.1, 0.3, 1.0, 3.0, 10.0]
    m2_grid = [0.1, 1.0, 10.0, 100.0, 1000.0]
    lm_grid = [0.001, 0.01, 0.1, 1.0, 10.0]
    for tcs in tc_grid:
        for m2s in m2_grid:
            for lms in lm_grid:
                val = _obj(math.log(tcs), math.log(m2s), math.log(lms))
                if val < best:
                    best = val
                    best_abc = (math.log(tcs), math.log(m2s), math.log(lms))
    a, b, c = best_abc
    lr = 0.2
    nsteps = max(8, int(rounds) * 14)

    for _ in range(nsteps):
        h = 5e-3
        ga = (_obj(a + h, b, c) - _obj(a - h, b, c)) / (2.0 * h)
        gb = (_obj(a, b + h, c) - _obj(a, b - h, c)) / (2.0 * h)
        gc = (_obj(a, b, c + h) - _obj(a, b, c - h)) / (2.0 * h)
        improved = False
        for _ls in range(8):
            an = a - lr * ga
            bn = b - lr * gb
            cn = c - lr * gc
            fn = _obj(an, bn, cn)
            if fn < best:
                a, b, c = an, bn, cn
                best = fn
                best_abc = (a, b, c)
                lr *= 1.05
                improved = True
                break
            lr *= 0.5
        if not improved and lr < 1e-5:
            break

    tcs, m2s, lms = _decode(*best_abc)
    out = {
        "theta_c_scale": float(tcs),
        "m2_scale": float(m2s),
        "lam_scale": float(lms),
        "calibrated_chi2_proxy": float(best),
        "reg_strength": float(reg),
        "calibration_gamma": float(gamma),
        "calibration_V0": float(V0),
    }
    return out


def calibrate_global_bridge_coeffs_multi(
    *,
    points: list[dict],
    events: list[dict],
    rounds: int = 36,
    reg_strength: float = 0.2,
    theta_c_scale0: float = 1.0,
    m2_scale0: float = 10.0,
    lam_scale0: float = 0.1,
    **_legacy_kwargs,
) -> Dict[str, float]:
    """
    Multi-point global bridge calibration.

    `points` item format:
      {
        "gamma": float,
        "V0": float,
        "Q": float (optional, default 0),
        "weight": float (optional, default 1),
      }
    """
    if not points:
        raise ValueError("points is empty")

    tc0 = max(1e-6, float(theta_c_scale0))
    m20 = max(1e-6, float(m2_scale0))
    lm0 = max(1e-9, float(lam_scale0))
    reg = max(0.0, float(reg_strength))

    def _decode(a: float, b: float, c: float) -> tuple[float, float, float]:
        a = min(30.0, max(-30.0, float(a)))
        b = min(30.0, max(-30.0, float(b)))
        c = min(30.0, max(-30.0, float(c)))
        return math.exp(a), math.exp(b), math.exp(c)

    def _obj(a: float, b: float, c: float) -> float:
        tcs, m2s, lms = _decode(a, b, c)
        chi2_sum = 0.0
        w_sum = 0.0
        for p in points:
            gamma = float(p["gamma"])
            v0 = float(p["V0"])
            q = float(p.get("Q", 0.0))
            w = max(0.0, float(p.get("weight", 1.0)))
            if w <= 0.0:
                continue
            d = derive_compact_star_params_from_background(
                gamma=gamma,
                V0=v0,
                Q=q,
                theta_c_scale=tcs,
                m2_scale=m2s,
                lam_scale=lms,
                forward_events=None,
                calibrate_forward_rounds=0,
            )
            prof = solve_compact_star_profile(
                CompactStarParams(
                    theta_c=max(1e-6, float(d["theta_c"])),
                    m2=max(1e-6, float(d["m2"])),
                    lam=max(0.0, float(d["lam"])),
                    r_max=20.0,
                    dr=0.02,
                    theta_scale=1.0,
                )
            )
            chi2_f = -2.0 * float(loglike_forward_events(prof, events))
            chi2_sum += w * chi2_f
            w_sum += w
        if w_sum <= 0.0:
            return 1e18
        val = chi2_sum / w_sum
        if reg > 0.0:
            da = math.log(max(tcs, 1e-12) / max(tc0, 1e-12))
            db = math.log(max(m2s, 1e-12) / max(m20, 1e-12))
            dc = math.log(max(lms, 1e-12) / max(lm0, 1e-12))
            val += reg * (da * da + db * db + dc * dc)
        return val

    a = math.log(tc0)
    b = math.log(m20)
    c = math.log(lm0)
    best = _obj(a, b, c)
    best_abc = (a, b, c)
    tc_grid = [0.1, 0.3, 1.0, 3.0, 10.0]
    m2_grid = [0.1, 1.0, 10.0, 100.0, 1000.0]
    lm_grid = [0.001, 0.01, 0.1, 1.0, 10.0]
    for tcs in tc_grid:
        for m2s in m2_grid:
            for lms in lm_grid:
                val = _obj(math.log(tcs), math.log(m2s), math.log(lms))
                if val < best:
                    best = val
                    best_abc = (math.log(tcs), math.log(m2s), math.log(lms))
    a, b, c = best_abc
    lr = 0.2
    nsteps = max(8, int(rounds) * 14)
    for _ in range(nsteps):
        h = 5e-3
        ga = (_obj(a + h, b, c) - _obj(a - h, b, c)) / (2.0 * h)
        gb = (_obj(a, b + h, c) - _obj(a, b - h, c)) / (2.0 * h)
        gc = (_obj(a, b, c + h) - _obj(a, b, c - h)) / (2.0 * h)
        improved = False
        for _ls in range(8):
            an = a - lr * ga
            bn = b - lr * gb
            cn = c - lr * gc
            fn = _obj(an, bn, cn)
            if fn < best:
                a, b, c = an, bn, cn
                best = fn
                best_abc = (a, b, c)
                lr *= 1.05
                improved = True
                break
            lr *= 0.5
        if not improved and lr < 1e-5:
            break

    tcs, m2s, lms = _decode(*best_abc)
    return {
        "theta_c_scale": float(tcs),
        "m2_scale": float(m2s),
        "lam_scale": float(lms),
        "calibrated_chi2_proxy": float(best),
        "reg_strength": float(reg),
        "n_points": int(len(points)),
    }


def calibrate_bridge_powerlaw2d_multi(
    *,
    points: list[dict],
    events: list[dict],
    rounds: int = 40,
    reg_strength: float = 0.05,
    gamma_ref: float = 0.0020339622,
    V0_ref: float = 0.012364227,
    theta0: float = 1.0,
    m20: float = 10.0,
    lam0: float = 0.1,
) -> Dict[str, float]:
    """
    Calibrate functional bridge scales:
      s(g,V) = s0 * (g/g_ref)^ag * (V/V_ref)^av
    for theta_c, m2 and lam.
    """
    if not points:
        raise ValueError("points is empty")

    reg = max(0.0, float(reg_strength))
    gref = max(1e-12, float(gamma_ref))
    vref = max(1e-12, float(V0_ref))
    max_slope = 4.0
    # params: [lt0, atg, atv, lm0, amg, amv, ll0, alg, alv]
    p = [
        math.log(max(theta0, 1e-12)),
        0.0,
        0.0,
        math.log(max(m20, 1e-12)),
        0.0,
        0.0,
        math.log(max(lam0, 1e-12)),
        0.0,
        0.0,
    ]

    def _decode(pp: list[float], gamma: float, V0: float) -> tuple[float, float, float, tuple[float, float, float, float, float, float]]:
        lg = math.log(max(gamma, 1e-12) / gref)
        lv = math.log(max(V0, 1e-12) / vref)
        lt0 = min(8.0, max(-6.0, pp[0]))
        lm0 = min(8.0, max(-6.0, pp[3]))
        ll0 = min(8.0, max(-6.0, pp[6]))
        atg = max_slope * math.tanh(pp[1])
        atv = max_slope * math.tanh(pp[2])
        amg = max_slope * math.tanh(pp[4])
        amv = max_slope * math.tanh(pp[5])
        alg = max_slope * math.tanh(pp[7])
        alv = max_slope * math.tanh(pp[8])
        lt = lt0 + atg * lg + atv * lv
        lm = lm0 + amg * lg + amv * lv
        ll = ll0 + alg * lg + alv * lv
        lt = min(30.0, max(-30.0, lt))
        lm = min(30.0, max(-30.0, lm))
        ll = min(30.0, max(-30.0, ll))
        return math.exp(lt), math.exp(lm), math.exp(ll), (lt0, atg, atv, lm0, amg, amv, ll0, alg, alv)

    def _obj(pp: list[float]) -> float:
        chi2_sum = 0.0
        w_sum = 0.0
        for pt in points:
            g = float(pt["gamma"])
            v = float(pt["V0"])
            q = float(pt.get("Q", 0.0))
            w = max(0.0, float(pt.get("weight", 1.0)))
            if w <= 0.0:
                continue
            tcs, m2s, lms, _tr = _decode(pp, g, v)
            d = derive_compact_star_params_from_background(
                gamma=g,
                V0=v,
                Q=q,
                theta_c_scale=tcs,
                m2_scale=m2s,
                lam_scale=lms,
                forward_events=None,
                calibrate_forward_rounds=0,
            )
            prof = solve_compact_star_profile(
                CompactStarParams(
                    theta_c=max(1e-6, float(d["theta_c"])),
                    m2=max(1e-6, float(d["m2"])),
                    lam=max(0.0, float(d["lam"])),
                    r_max=20.0,
                    dr=0.02,
                    theta_scale=1.0,
                )
            )
            chi2_f = -2.0 * float(loglike_forward_events(prof, events))
            chi2_sum += w * chi2_f
            w_sum += w
        if w_sum <= 0.0:
            return 1e18
        val = chi2_sum / w_sum
        if reg > 0.0:
            # weak L2 prior around ag=av=0 and base scales (theta0,m20,lam0)
            _tcs, _m2s, _lms, tr = _decode(pp, gref, vref)
            lt0, atg, atv, lm0, amg, amv, ll0, alg, alv = tr
            val += reg * (
                (lt0 - math.log(max(theta0, 1e-12))) ** 2
                + (lm0 - math.log(max(m20, 1e-12))) ** 2
                + (ll0 - math.log(max(lam0, 1e-12))) ** 2
                + atg ** 2
                + atv ** 2
                + amg ** 2
                + amv ** 2
                + alg ** 2
                + alv ** 2
            )
        return val

    # Fast SPSA optimizer: 2 objective calls per step instead of 2*N.
    best_p = p[:]
    best = _obj(best_p)
    cur = best_p[:]
    cur_val = best
    rng = random.Random(42)
    nsteps = max(16, int(rounds) * 6)
    for k in range(1, nsteps + 1):
        ak = 0.12 / (k ** 0.602)
        ck = 0.08 / (k ** 0.101)
        delta = [1.0 if rng.random() < 0.5 else -1.0 for _ in range(len(cur))]
        up = [c + ck * d for c, d in zip(cur, delta)]
        dn = [c - ck * d for c, d in zip(cur, delta)]
        fu = _obj(up)
        fd = _obj(dn)
        ghat = [((fu - fd) / (2.0 * ck * d)) for d in delta]
        trial = [c - ak * g for c, g in zip(cur, ghat)]
        ftrial = _obj(trial)
        if ftrial < cur_val:
            cur, cur_val = trial, ftrial
        if cur_val < best:
            best_p, best = cur[:], cur_val

    _tcs, _m2s, _lms, trb = _decode(best_p, gref, vref)
    lt0, atg, atv, lm0, amg, amv, ll0, alg, alv = trb
    def _exp_safe(x: float) -> float:
        x = min(8.0, max(-6.0, float(x)))
        return float(math.exp(x))

    return {
        "bridge_model": "powerlaw2d",
        "gamma_ref": float(gref),
        "V0_ref": float(vref),
        "theta_c": {"a0": _exp_safe(lt0), "ag": float(atg), "av": float(atv)},
        "m2": {"a0": _exp_safe(lm0), "ag": float(amg), "av": float(amv)},
        "lam": {"a0": _exp_safe(ll0), "ag": float(alg), "av": float(alv)},
        "calibrated_chi2_proxy": float(best),
        "reg_strength": float(reg),
        "n_points": int(len(points)),
    }
