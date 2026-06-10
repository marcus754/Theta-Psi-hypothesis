# -*- coding: utf-8 -*-
"""
Head-to-head likelihood helpers for Theta-Psi vs flat LCDM/GR.

Implements:
- SN profile chi2 with optional full covariance and analytic mu0.
- BAO vector chi2 for DESI-like mean/cov tables in *_over_rs units.
- Shared prediction helpers for Theta-Psi (via model.py) and flat LCDM.
"""
from __future__ import annotations

from dataclasses import dataclass
from bisect import bisect_left
from pathlib import Path
import math
from typing import Iterable

import numpy as np

from fitting.core_api import Omegas
from fitting.model import Hz_dimless, mu_dimless


@dataclass
class BAOVectorData:
    z: list[float]
    obs: list[str]  # one of {"DV_over_rs", "DM_over_rs", "DH_over_rs"}
    y: np.ndarray
    cov: np.ndarray
    cov_inv: np.ndarray


def load_bao_vector_from_desi(mean_path: str | Path, cov_path: str | Path) -> BAOVectorData:
    z: list[float] = []
    obs: list[str] = []
    y: list[float] = []
    with Path(mean_path).open("r", encoding="utf-8") as f:
        for ln in f:
            t = ln.strip()
            if (not t) or t.startswith("#"):
                continue
            p = t.split()
            if len(p) < 3:
                continue
            z.append(float(p[0]))
            y.append(float(p[1]))
            obs.append(str(p[2]))
    cov = np.loadtxt(str(cov_path), dtype=float)
    if cov.ndim == 1:
        n = int(round(math.sqrt(float(cov.size))))
        cov = cov.reshape((n, n))
    n = len(y)
    if cov.shape != (n, n):
        raise ValueError(f"BAO cov shape {cov.shape} does not match mean length {n}")
    cov_inv = np.linalg.inv(cov)
    return BAOVectorData(z=z, obs=obs, y=np.asarray(y, dtype=float), cov=cov, cov_inv=cov_inv)


def _interp(x: list[float], y: list[float], xq: float) -> float:
    if xq <= x[0]:
        return float(y[0])
    if xq >= x[-1]:
        return float(y[-1])
    j = bisect_left(x, xq)
    x0, x1 = x[j - 1], x[j]
    y0, y1 = y[j - 1], y[j]
    t = (xq - x0) / max(x1 - x0, 1e-18)
    return float(y0 + t * (y1 - y0))


def _cumtrapz(y: list[float], x: list[float]) -> list[float]:
    out = [0.0] * len(x)
    s = 0.0
    for i in range(1, len(x)):
        dx = x[i] - x[i - 1]
        s += 0.5 * (y[i] + y[i - 1]) * dx
        out[i] = s
    return out


def _build_dc_from_E(E_fn, zmax: float, nint: int = 3000) -> tuple[list[float], list[float]]:
    n = max(64, int(nint))
    zs = [zmax * i / n for i in range(n + 1)]
    Es = [float(E_fn(z)) for z in zs]
    invE = [1.0 / max(e, 1e-18) for e in Es]
    dc = _cumtrapz(invE, zs)  # H0*Dc/c
    return zs, dc


def profile_mu0_diag(mu_model: Iterable[float], mu_data: Iterable[float], sigma: Iterable[float]) -> float:
    wsum = 0.0
    num = 0.0
    for m, d, s in zip(mu_model, mu_data, sigma):
        ss = max(float(s), 1e-18)
        w = 1.0 / (ss * ss)
        wsum += w
        num += (d - m) * w
    return 0.0 if wsum <= 0.0 else (num / wsum)


def profile_mu0_cov(mu_model: np.ndarray, mu_data: np.ndarray, cov_inv: np.ndarray) -> float:
    ones = np.ones_like(mu_model)
    r = mu_data - mu_model
    den = float(ones @ cov_inv @ ones)
    if den <= 0.0:
        return 0.0
    num = float(ones @ cov_inv @ r)
    return num / den


def sn_chi2(
    mu_model: list[float],
    mu_data: list[float],
    sigma_diag: list[float],
    *,
    cov: np.ndarray | None = None,
    fit_mu0: bool = True,
) -> tuple[float, float]:
    m = np.asarray(mu_model, dtype=float)
    d = np.asarray(mu_data, dtype=float)
    if cov is None:
        s = np.asarray(sigma_diag, dtype=float)
        mu0 = profile_mu0_diag(m, d, s) if fit_mu0 else 0.0
        r = (m + mu0 - d) / np.maximum(s, 1e-18)
        return float(np.sum(r * r)), float(mu0)
    cov_use = np.array(cov, dtype=float, copy=True)
    cov_inv = np.linalg.inv(cov_use)
    mu0 = profile_mu0_cov(m, d, cov_inv) if fit_mu0 else 0.0
    r = m + mu0 - d
    chi2 = float(r @ cov_inv @ r)
    return chi2, float(mu0)


def theta_model_E_fn(
    gamma: float,
    V0: float,
    Q: float = 0.0,
    **_legacy_kwargs,
):
    def _E(z: float) -> float:
        return float(
            Hz_dimless(
                [z], gamma=gamma, V0=V0, Q=Q, omegas=Omegas(H0=1.0), a_min=1e-3, dN=0.02,
            )[0]
        )

    return _E


def lcdm_E_fn(omega_m: float):
    om = float(omega_m)
    ol = 1.0 - om

    def _E(z: float) -> float:
        return math.sqrt(max(om * (1.0 + z) ** 3 + ol, 1e-18))

    return _E


def gr_E_fn(omega_m: float, omega_r: float = 8.4e-5):
    om = float(omega_m)
    orad = float(omega_r)

    def _E(z: float) -> float:
        return math.sqrt(max(orad * (1.0 + z) ** 4 + om * (1.0 + z) ** 3, 1e-18))

    return _E


def bao_vector_model(
    data: BAOVectorData,
    *,
    E_fn,
    alpha_c_over_h0_rd: float,
) -> np.ndarray:
    zmax = max(data.z)
    zs, dc = _build_dc_from_E(E_fn, zmax=zmax, nint=3500)
    out = []
    for z, o in zip(data.z, data.obs):
        Ez = float(E_fn(z))
        Dm = _interp(zs, dc, z)  # H0*Dm/c
        if o == "DM_over_rs":
            out.append(alpha_c_over_h0_rd * Dm)
        elif o == "DH_over_rs":
            out.append(alpha_c_over_h0_rd / max(Ez, 1e-18))
        elif o == "DV_over_rs":
            dv_dimless = ((z / max(Ez, 1e-18)) * (Dm * Dm)) ** (1.0 / 3.0)
            out.append(alpha_c_over_h0_rd * dv_dimless)
        else:
            raise ValueError(f"Unsupported BAO observable: {o}")
    return np.asarray(out, dtype=float)


def bao_vector_chi2(y_model: np.ndarray, data: BAOVectorData) -> float:
    r = y_model - data.y
    return float(r @ data.cov_inv @ r)


def theta_predictions_for_sn(
    sn_z: list[float],
    *,
    gamma: float,
    V0: float,
    Q: float = 0.0,
    **_legacy_kwargs,
) -> list[float]:
    return mu_dimless(
        sn_z, gamma=gamma, V0=V0, Q=Q, omegas=Omegas(H0=1.0), a_min=1e-3, dN=0.02, mu0=0.0,
    )


def lcdm_predictions_for_sn(sn_z: list[float], *, omega_m: float) -> list[float]:
    E = lcdm_E_fn(omega_m)
    zmax = max(sn_z) if sn_z else 0.0
    zs, dc = _build_dc_from_E(E, zmax=zmax, nint=3500)
    out = []
    for z in sn_z:
        Dl = (1.0 + z) * _interp(zs, dc, z)
        out.append(5.0 * math.log10(max(Dl, 1e-18)))
    return out


def gr_predictions_for_sn(sn_z: list[float], *, omega_m: float, omega_r: float = 8.4e-5) -> list[float]:
    E = gr_E_fn(omega_m, omega_r=omega_r)
    zmax = max(sn_z) if sn_z else 0.0
    zs, dc = _build_dc_from_E(E, zmax=zmax, nint=3500)
    out = []
    for z in sn_z:
        Dl = (1.0 + z) * _interp(zs, dc, z)
        out.append(5.0 * math.log10(max(Dl, 1e-18)))
    return out


def theta_predictions_for_cc(
    cc_z: list[float],
    *,
    gamma: float,
    V0: float,
    Q: float = 0.0,
    **_legacy_kwargs,
) -> list[float]:
    return Hz_dimless(
        cc_z, gamma=gamma, V0=V0, Q=Q, omegas=Omegas(H0=1.0), a_min=1e-3, dN=0.02,
    )


def lcdm_predictions_for_cc(cc_z: list[float], *, omega_m: float) -> list[float]:
    E = lcdm_E_fn(omega_m)
    return [E(z) for z in cc_z]


def gr_predictions_for_cc(cc_z: list[float], *, omega_m: float, omega_r: float = 8.4e-5) -> list[float]:
    E = gr_E_fn(omega_m, omega_r=omega_r)
    return [E(z) for z in cc_z]
