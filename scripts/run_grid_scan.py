# -*- coding: utf-8 -*-
"""
Grid scan of (gamma, V0) with health checks.

Использует фоновой интегратор `frw_background` и набор проверок из `checks/`
для фильтрации параметров, которые нарушают устойчивость или простые EFT-ограничения.

Usage (CLI):
    python scripts/run_grid_scan.py --gammas 0.5 1.0 1.5 --v0 0.05 0.1
"""
from __future__ import annotations

import argparse
import csv
from dataclasses import dataclass, field
from typing import Iterable, List, Sequence
import math

from fitting.core_api import BackgroundParams, Omegas, integrate_background, integrate_background_lnN, integrate_background_rk45, omegas_from_km_s_Mpc
from checks import cutoff_ok, has_superluminal, has_superluminal_adm, finite_k_stable_along_traj
from src.linear_modes import stability_conditions
from checks.health import evaluate_health as _eval_health
from checks.compact_report import compute_compact_metrics, write_compact_markdown

# Optional ln a integrator
try:
    from fitting.core_api import integrate_background_lnN as _integrate_lnN
    _HAS_LNN = True
except Exception:
    _HAS_LNN = False


@dataclass
class InitialConditions:
    a0: float = 1e-3
    x0: float = 0.1
    xdot0: float = 0.0
    psi0: float = 0.05
    psidot0: float = 0.0

    def as_tuple(self):
        return (self.a0, self.x0, self.xdot0, self.psi0, self.psidot0)


@dataclass
class ScanConfig:
    h: float = 0.01
    nsteps: int = 200
    ctheta2: float = 1.0
    cpsi2: float = 1.0
    cutoff: float = 1e-1
    # Background cosmology
    omegas: Omegas = field(default_factory=Omegas)
    rk45: bool = False
    lnN: bool = False
    use_health: bool = False
    check_no_horizon: bool = False
    theta_scale: float = 1.0
    # Phase charge of complex Ψ; Q=0 reduces to real-amplitude case
    Q: float = 0.0
    ic: InitialConditions = field(default_factory=InitialConditions)


def scan_point(gamma: float, V0: float, cfg: ScanConfig) -> dict:
    params = BackgroundParams(gamma=gamma, V0=V0, Q=cfg.Q, omegas=cfg.omegas)
    try:
        if cfg.lnN and _HAS_LNN:
            # Integrate over N = ln a with step h treated as ΔN
            import math as _m
            N0 = _m.log(cfg.ic.a0)
            N1 = N0 + cfg.h * cfg.nsteps
            out = _integrate_lnN(cfg.ic.a0, cfg.ic.x0, cfg.ic.xdot0, cfg.ic.psi0, cfg.ic.psidot0,
                                  N1=N1, p=params, hN=cfg.h)
        elif cfg.rk45:
            t1 = cfg.h * cfg.nsteps
            out = integrate_background_rk45(0.0, cfg.ic.as_tuple(), t1=t1, p=params, h0=cfg.h)
        else:
            out = integrate_background(
                0.0, cfg.ic.as_tuple(), h=cfg.h, nsteps=cfg.nsteps, p=params
            )
    except Exception as exc:
        return {
            "gamma": gamma,
            "V0": V0,
            "ok": False,
            "error": str(exc),
        }

    theta = out["theta"]
    theta_dot = out["theta_dot"]
    H_list = out["H"]

    use_health = cfg.use_health
    adm_opts = getattr(cfg, "_adm", {})
    if use_health:
        hres = _eval_health(
            out,
            gamma=gamma,
            V0=V0,
            Q=cfg.Q,
            ctheta2=cfg.ctheta2,
            cpsi2=cfg.cpsi2,
            use_adm=adm_opts.get("use", False),
            cR2=adm_opts.get("cR2", 1.0),
            cphi2=adm_opts.get("cphi2", 1.0),
            epsKtr=adm_opts.get("epsKtr", 0.0),
            epsKtphi=adm_opts.get("epsKtphi", 0.0),
            epsKRphi=adm_opts.get("epsKRphi", 0.0),
            epsGtr=adm_opts.get("epsGtr", 0.0),
            epsGtphi=adm_opts.get("epsGtphi", 0.0),
            epsGRphi=adm_opts.get("epsGRphi", 0.0),
            alpha_k2_theta=adm_opts.get("alpha_k2_theta", 0.0),
            alpha_k2_R=adm_opts.get("alpha_k2_R", 0.0),
            alpha_k2_phi=adm_opts.get("alpha_k2_phi", 0.0),
            finite_k=adm_opts.get("finite_k", False),
            kmin=adm_opts.get("kmin", 1e-2),
            kmax=adm_opts.get("kmax", 1.0),
            nk=adm_opts.get("nk", 5),
            cutoff=cfg.cutoff,
            strict_first_principles=True,
            gravM=adm_opts.get("gravM", False),
            check_no_horizon=cfg.check_no_horizon,
            theta_scale=cfg.theta_scale,
        )
        ok_lin = hres["ok_linear"]
        ok_sub = hres["ok_sub"]
        ok_cutoff = hres["ok_cutoff"]
        ok_no_horizon = hres.get("ok_no_horizon", True)
        ok_n_bounds = ok_no_horizon
        ok = hres["ok"]
    else:
        # Legacy path
        ok_lin = True
        for th in theta:
            cond = stability_conditions(th, gamma, cfg.ctheta2, cfg.cpsi2)
            if not (cond["no_ghost"] and cond["no_gradient_instability"]):
                ok_lin = False
                break
        ok_sub = True
        use_adm = adm_opts.get("use", False)
        if use_adm:
            cR2 = adm_opts.get("cR2", 1.0)
            cphi2 = adm_opts.get("cphi2", 1.0)
            epsKtr = adm_opts.get("epsKtr", 0.0)
            epsKtphi = adm_opts.get("epsKtphi", 0.0)
            epsKRphi = adm_opts.get("epsKRphi", 0.0)
            epsGtr = adm_opts.get("epsGtr", 0.0)
            epsGtphi = adm_opts.get("epsGtphi", 0.0)
            epsGRphi = adm_opts.get("epsGRphi", 0.0)
            alpha_k2_theta = adm_opts.get("alpha_k2_theta", 0.0)
            alpha_k2_R = adm_opts.get("alpha_k2_R", 0.0)
            alpha_k2_phi = adm_opts.get("alpha_k2_phi", 0.0)
            for th, R, a in zip(theta, out["psi"], out["a"]):
                if has_superluminal_adm(
                    th, R, a, 1.0, gamma, cfg.Q, V0,
                    ctheta2=cfg.ctheta2, cR2=cR2, cphi2=cphi2,
                    epsK_theta_R=epsKtr, epsK_theta_phi=epsKtphi, epsK_R_phi=epsKRphi,
                    epsG_theta_R=epsGtr, epsG_theta_phi=epsGtphi, epsG_R_phi=epsGRphi,
                    alpha_k2_theta=alpha_k2_theta, alpha_k2_R=alpha_k2_R, alpha_k2_phi=alpha_k2_phi,
                    adm_strict=adm_opts.get("strict", False),
                ):
                    ok_sub = False
                    break
        else:
            for th, H in zip(theta, H_list):
                if has_superluminal(th, H, gamma, V0, cfg.ctheta2, cfg.cpsi2):
                    ok_sub = False
                    break
        # Optional finite-k stability (ADM)
        if ok_sub and use_adm and adm_opts.get("finite_k", False):
            cR2 = adm_opts.get("cR2", 1.0)
            cphi2 = adm_opts.get("cphi2", 1.0)
            kmin = adm_opts.get("kmin", 1e-2)
            kmax = adm_opts.get("kmax", 1.0)
            nk = adm_opts.get("nk", 5)
            fk_ok = finite_k_stable_along_traj(
                theta,
                out["psi"],
                out["a"],
                gamma,
                cfg.Q,
                V0,
                ctheta2=cfg.ctheta2,
                cR2=cR2,
                cphi2=cphi2,
                kmin=kmin,
                kmax=kmax,
                nk=nk,
                epsK_theta_R=adm_opts.get("epsKtr", 0.0),
                epsK_theta_phi=adm_opts.get("epsKtphi", 0.0),
                epsK_R_phi=adm_opts.get("epsKRphi", 0.0),
                epsG_theta_R=adm_opts.get("epsGtr", 0.0),
                epsG_theta_phi=adm_opts.get("epsGtphi", 0.0),
                epsG_R_phi=adm_opts.get("epsGRphi", 0.0),
                alpha_k2_theta=adm_opts.get("alpha_k2_theta", 0.0),
                alpha_k2_R=adm_opts.get("alpha_k2_R", 0.0),
                alpha_k2_phi=adm_opts.get("alpha_k2_phi", 0.0),
                gravM=adm_opts.get("gravM", False),
                adm_strict=adm_opts.get("strict", False),
                H_list=out["H"],
                theta_dot_list=out["theta_dot"],
                R_dot_list=out["psi_dot"],
                Hdot_list=out.get("Hdot", []),
            )
            if not fk_ok:
                ok_sub = False
        ok_cutoff = cutoff_ok(theta, theta_dot, gamma, cfg.cutoff)
        ok_no_horizon = True
        ok_n_bounds = ok_no_horizon
        ok = ok_lin and ok_sub and ok_cutoff

    theta_star = math.sqrt(3.0) / abs(gamma) if gamma != 0.0 else float("inf")
    cmetrics = compute_compact_metrics(theta, theta_scale=cfg.theta_scale, n_obs=1.0)

    return {
        "gamma": gamma,
        "V0": V0,
        "ok": ok,
        "ok_linear": ok_lin,
        "ok_sub": ok_sub,
        "ok_cutoff": ok_cutoff,
        "ok_n_bounds": ok_n_bounds,
        "ok_no_horizon": ok_no_horizon,
        "H_end": H_list[-1],
        "theta_end": theta[-1],
        "theta_ratio_max": max(abs(th) for th in theta) / theta_star,
        "n_surface": cmetrics["n_surface"],
        "n_peak": cmetrics["n_peak"],
        "z_surface": cmetrics["z_surface"],
        "z_surface_proxy": cmetrics["z_surface_proxy"],
        "delay": cmetrics["delay"],
        "delay_proxy": cmetrics["delay_proxy"],
    }


def run_scan(gammas: Sequence[float], V0_values: Sequence[float], cfg: ScanConfig) -> List[dict]:
    rows: List[dict] = []
    for g in gammas:
        for v in V0_values:
            rows.append(scan_point(g, v, cfg))
    return rows


def write_csv(rows: Iterable[dict], path: str) -> None:
    rows = list(rows)
    if not rows:
        return
    fieldnames = sorted(rows[0].keys())
    with open(path, "w", newline="", encoding="utf-8") as fh:
        writer = csv.DictWriter(fh, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def parse_args(argv: Sequence[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Run parameter grid scan")
    parser.add_argument("--gammas", type=float, nargs="+", required=True)
    parser.add_argument("--v0", type=float, nargs="+", required=True)
    parser.add_argument("--output", type=str, default="scan_results.csv")
    parser.add_argument("--compact_report", type=str, default=None, help="Optional Markdown report for compact no-horizon metrics")
    parser.add_argument("--nsteps", type=int, default=200)
    parser.add_argument("--step", type=float, default=0.01, help="Integration step h")
    parser.add_argument("--cutoff", type=float, default=1e-1)
    parser.add_argument("--H0_km_s_Mpc", type=float, default=None, help="If set, convert this H0 (km/s/Mpc) to code units and use in background Omegas")
    parser.add_argument("--Q", type=float, default=0.0, help="Conserved phase charge for complex Ψ; Q=0 disables phase dynamics")
    parser.add_argument("--ctheta2", type=float, default=1.0, help="Sound speed squared for θ-sector in the linear-mode scaffold")
    parser.add_argument("--cpsi2", type=float, default=1.0, help="Sound speed squared for Ψ-sector in the linear-mode scaffold")
    parser.add_argument("--rk45", action="store_true", help="Use adaptive RK45 for background evolution")
    parser.add_argument("--lnN", action="store_true", help="Integrate background over N=ln a (step is ΔN)")
    parser.add_argument("--health", action="store_true", help="Use unified health checks (aggregated vetoes)")
    parser.add_argument("--check_no_horizon", action="store_true", help="Enable explicit finite-branch no-horizon check")
    parser.add_argument("--theta_scale", type=float, default=1.0, help="Scale for finite n(theta) profile")
    parser.add_argument("--adm", action="store_true", help="Use ADM 3×3 superluminality check (θ,R,φ)")
    parser.add_argument("--adm-strict", action="store_true", help="Use strict σ-model elimination (luminal G=K, k^0 mass corrections)")
    parser.add_argument("--cR2", type=float, default=1.0, help="Sound speed squared for R-sector (ADM 3×3)")
    parser.add_argument("--cphi2", type=float, default=1.0, help="Sound speed squared for φ-sector (ADM 3×3)")
    parser.add_argument("--epsKtr", type=float, default=0.0, help="K mixing θ–R")
    parser.add_argument("--epsKtphi", type=float, default=0.0, help="K mixing θ–φ")
    parser.add_argument("--epsKRphi", type=float, default=0.0, help="K mixing R–φ")
    parser.add_argument("--epsGtr", type=float, default=0.0, help="G mixing θ–R")
    parser.add_argument("--epsGtphi", type=float, default=0.0, help="G mixing θ–φ")
    parser.add_argument("--epsGRphi", type=float, default=0.0, help="G mixing R–φ")
    parser.add_argument("--finite_k", action="store_true", help="Also require finite-k stability (ω^2>0) along trajectory")
    parser.add_argument("--kmin", type=float, default=1e-2)
    # duplicate line cleanup handled; keep only one set
    parser.add_argument("--kmax", type=float, default=1.0)
    parser.add_argument("--nk", type=int, default=5)
    parser.add_argument("--gravM", action="store_true", help="Use gravitational k^0 mass corrections in finite-k stability")
    # removed: --elim-defaults (heuristic elimination couplings are deprecated)
    return parser.parse_args(argv)


def main(argv: Sequence[str] | None = None) -> int:
    args = parse_args(argv)
    # Build background Omegas (with optional H0 conversion)
    if args.H0_km_s_Mpc is not None:
        om = omegas_from_km_s_Mpc(args.H0_km_s_Mpc)
    else:
        om = Omegas()
    cfg = ScanConfig(
        h=args.step,
        nsteps=args.nsteps,
        cutoff=args.cutoff,
        omegas=om,
        Q=args.Q,
    )
    cfg.ctheta2 = args.ctheta2
    cfg.cpsi2 = args.cpsi2
    cfg.rk45 = args.rk45
    cfg.lnN = args.lnN
    cfg.use_health = args.health
    cfg.check_no_horizon = args.check_no_horizon
    cfg.theta_scale = args.theta_scale
    # Attach ADM choices to cfg to avoid globals
    cfg._adm = {
        "use": args.adm,
        "strict": args.adm_strict,
        "cR2": args.cR2,
        "cphi2": args.cphi2,
        "epsKtr": args.epsKtr,
        "epsKtphi": args.epsKtphi,
        "epsKRphi": args.epsKRphi,
        "epsGtr": args.epsGtr,
        "epsGtphi": args.epsGtphi,
        "epsGRphi": args.epsGRphi,
        "finite_k": args.finite_k,
        "kmin": args.kmin,
        "kmax": args.kmax,
        "nk": args.nk,
        "gravM": args.gravM,
        # no elim_defaults
    }
    rows = run_scan(args.gammas, args.v0, cfg)
    write_csv(rows, args.output)
    if args.compact_report:
        write_compact_markdown(rows, args.compact_report)
    print(f"Wrote {len(rows)} rows to {args.output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
