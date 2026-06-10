# -*- coding: utf-8 -*-
"""
Background integrator with explicit scale-factor dynamics.

State vector: (a, x, xdot, psi, psidot), где
    a      — масштабный фактор,
    x      — регуляризованная переменная, x = artanh(γ θ / √3),
    psi    — reduced amplitude representative R of the complex energy sector
             Ψ = R e^{iφ}.

Интегратор использует простой RK4-шаг и суммирует вклад материи, излучения и Λ,
заданные в `background.components.Omegas`. Основное предназначение — сканы по
параметрам и проверка того, что reduced amplitude sector и phase charge не
портят стандартную космологию.
"""
from __future__ import annotations

from dataclasses import dataclass, field
from typing import Dict, List, Tuple
import math

from background.components import Omegas, rho_m, rho_r, rho_L
from theory.optical_metric import diagnostic_refractive_index_from_theta
from background.frw_symbolic import (
    energy_density,
    theta_from_x,
    theta_dot_from_x,
    x_eom_rhs,
)


@dataclass
class Params:
    gamma: float = 1.0
    V0: float = 0.1
    # Conserved phase charge for complex Ψ = R e^{iφ}:
    # Q = a^3 R^2 φ̇. When Q ≠ 0, the amplitude R=psi obeys an extra
    # repulsive term and the energy gets a contribution ½ Q^2/(a^6 R^2).
    Q: float = 0.0
    # Optical index controls the local weak-field response.
    theta_scale: float = 1.0
    n_kappa_wf: float = 2.0
    omegas: Omegas = field(default_factory=Omegas)


def rho_fields(x: float, xdot: float, psi: float, psidot: float, p: Params) -> float:
    theta = theta_from_x(x, p.gamma)
    thetadot = theta_dot_from_x(x, xdot, p.gamma)
    return energy_density(theta, thetadot, psi, psidot, p.gamma, p.V0)


def H_total(
    a: float,
    x: float,
    xdot: float,
    psi: float,
    psidot: float,
    p: Params,
) -> float:
    # Base field energy from minisuperspace L_eff (θ, Ψ amplitude only)
    rhof = rho_fields(x, xdot, psi, psidot, p)
    # Add phase energy if complex Ψ carries nonzero conserved charge Q
    if p.Q != 0.0:
        if abs(psi) < 1e-12:
            raise ValueError("R amplitude too small for nonzero Q (divergent phase energy)")
        rhof += 0.5 * (p.Q * p.Q) / (a ** 6 * (psi ** 2))
    rhot = rhof + rho_m(a, p.omegas) + rho_r(a, p.omegas) + rho_L(a, p.omegas)
    if rhot <= 0.0:
        raise ValueError("Total energy density non-positive")
    return math.sqrt(rhot / 3.0)


def Hdot_total(
    a: float,
    x: float,
    xdot: float,
    psi: float,
    psidot: float,
    H: float,
    p: Params,
) -> float:
    """
    Compute Ḣ using the ln a expression but via t-variables (xdot, psidot):
    yN = xdot/H and uN = psidot/H.
    """
    if H == 0.0:
        return 0.0
    yN = xdot / H
    uN = psidot / H
    return _Hdot_from_N_state(a, x, yN, psi, uN, H, p)


# --- Log-scale (N = ln a) integrator helpers --------------------------------
def _sech2(x: float) -> float:
    c = math.cosh(x)
    return 1.0 / (c * c)


def _ln_n_of_theta(theta: float, p: Params) -> float:
    nval = diagnostic_refractive_index_from_theta(
        theta,
        theta_scale=max(float(p.theta_scale), 1e-12),
        kappa_wf=float(p.n_kappa_wf),
        profile="asinh",
    )
    return math.log(max(nval, 1e-300))


def _ln_n_of_theta_eff(theta_eff: float, p: Params) -> float:
    """Alias for ln n(θ_eff)."""
    return _ln_n_of_theta(theta_eff, p)


def _lnn_theta_derivatives(theta: float, p: Params) -> tuple[float, float]:
    """Numerical derivatives d(ln n)/dtheta and d2(ln n)/dtheta2."""
    h = max(1e-8, 1e-5 * max(abs(theta), 1.0))
    l0 = _ln_n_of_theta(theta, p)
    lp = _ln_n_of_theta(theta + h, p)
    lm = _ln_n_of_theta(theta - h, p)
    d1 = (lp - lm) / (2.0 * h)
    d2 = (lp - 2.0 * l0 + lm) / (h * h)
    return d1, d2


def _lnn_theta_eff_derivatives(theta_eff: float, p: Params) -> tuple[float, float]:
    """Alias for derivatives of ln n with respect to θ_eff."""
    return _lnn_theta_derivatives(theta_eff, p)


def _rho_fields_coeff_H2(x: float, yN: float, uN: float, *, gamma: float) -> float:
    """
    Coefficient A such that ρ_fields = A H^2 + 1/2 V0 R^2 (phase term handled elsewhere),
    when using variables yN = dx/dN and uN = dpsi/dN.
    """
    # θ̇/θ = (sech^2 x / tanh x) * ẋ,   ẋ = yN * H
    tanh = math.tanh(x)
    if abs(tanh) < 1e-15:
        # Avoid division by zero near x~0 -> use series tanh x ≈ x
        tanh = x if x != 0.0 else 1e-15
    s2 = _sech2(x)
    a1 = 1.5 * (s2 * s2) / (tanh * tanh) * (yN * yN)           # from 3/2 (θ̇/θ)^2
    a2 = 0.5 * (uN * uN)                                        # from 1/2 Ṙ^2
    a3 = 1.5 * (s2 * s2) * (yN * yN)                            # from 1/2 γ^2 θ̇^2 with γ-cancelled
    return a1 + a2 + a3


def _H2_from_N_state(
    a: float,
    x: float,
    yN: float,
    psi: float,
    uN: float,
    p: Params,
) -> float:
    """
    Solve H^2 from Friedmann with variables defined in N=ln a.

    Direct implementation without artificial tanh-saturation. 
    If A >= 3, the system is in a non-physical or ghost regime.
    """
    A = _rho_fields_coeff_H2(x, yN, uN, gamma=p.gamma)
    rhobg = rho_m(a, p.omegas) + rho_r(a, p.omegas) + rho_L(a, p.omegas)
    Vterm = 0.5 * p.V0 * (psi * psi)
    phase = 0.0
    if p.Q != 0.0:
        if abs(psi) < 1e-12:
            raise ValueError("R amplitude too small for nonzero Q (divergent phase energy)")
        phase = 0.5 * (p.Q * p.Q) / (a ** 6 * (psi ** 2))
    
    denom = 3.0 - A
    if denom <= 0.0:
        raise ValueError(f"Friedmann singularity: A={A} >= 3.0. Kinetic energy of fields dominates.")
        
    return (rhobg + Vterm + phase) / denom


def _Hdot_from_N_state(
    a: float,
    x: float,
    yN: float,
    psi: float,
    uN: float,
    H: float,
    p: Params,
) -> float:
    """
    Compute Ḣ using `Hdot = -1/2 (ρ + p)_tot`.

    For fields: (ρ+p)_fields = (γ^2 θ̇^2 + Ṙ^2 + 3 (θ̇/θ)^2) = B H^2, with
    B = uN^2 + 3 s2^2/tanh^2 yN^2 + 3 s2^2 yN^2, where s2 = sech^2 x.
    For phase: (ρ+p)_phase = Q^2/(a^6 psi^2).
    For matter/radiation/Λ: ρ_m, 4/3 ρ_r, 0 respectively.
    """
    s2 = _sech2(x)
    tanh = math.tanh(x)
    if abs(tanh) < 1e-15:
        tanh = x if x != 0.0 else 1e-15
    B = (uN * uN) + 3.0 * (s2 * s2) * (yN * yN) * (1.0 + 1.0 / (tanh * tanh))
    fields = 0.5 * B * (H * H)
    matter = 0.5 * rho_m(a, p.omegas)
    rad = 0.5 * (4.0 / 3.0) * rho_r(a, p.omegas)
    phase = 0.0
    if p.Q != 0.0:
        if abs(psi) < 1e-12:
            raise ValueError("R amplitude too small for nonzero Q (divergent phase energy)")
        phase = 0.5 * (p.Q * p.Q) / (a ** 6 * (psi ** 2))
    return -(fields + matter + rad + phase)


def integrate_background_lnN(
    a0: float,
    x0: float,
    xdot0: float,
    psi0: float,
    psidot0: float,
    N1: float,
    p: Params,
    hN: float = 1e-2,
) -> Dict[str, List[float]]:
    """
    Integrate background using N = ln a as the independent variable (RK4).

    State in N: (x, y, psi, u, t) where y = dx/dN, u = dpsi/dN.
    Outputs arrays for N,a,t,H,x,theta,theta_dot,theta_eff,theta_eff_dot,psi.
    """
    N = math.log(a0)
    # Compute initial H from the t-based initial derivatives
    H0 = H_total(a0, x0, xdot0, psi0, psidot0, p)
    y = xdot0 / H0
    u = psidot0 / H0
    x = x0
    psi = psi0
    t = 0.0

    out: Dict[str, List[float]] = {
        "N": [],
        "a": [],
        "t": [],
        "H": [],
        "Hdot": [],
        "x": [],
        "theta": [],
        "theta_dot": [],
        "theta_eff": [],
        "theta_eff_dot": [],
        "psi": [],
        "psi_dot": [],
        "J_refr": [],
        "I_grad": [],
    }

    def rhs_N(N: float, state: Tuple[float, float, float, float, float]):
        x, y, psi, u, t = state
        a = math.exp(N)
        H2 = _H2_from_N_state(a, x, y, psi, u, p)
        H = math.sqrt(H2)
        Hdot = _Hdot_from_N_state(a, x, y, psi, u, H, p)
        xdot = y * H
        psidot = u * H
        xdd = x_eom_rhs(x, xdot, H, p.gamma)
        psidd = -3.0 * H * psidot - p.V0 * psi
        if p.Q != 0.0:
            if abs(psi) < 1e-12:
                raise ValueError("R amplitude too small for nonzero Q (divergent R-equation)")
            psidd += (p.Q * p.Q) / (a ** 6 * (psi ** 3))
        dy = (xdd - y * Hdot) / (H * H)
        # Convert to N-derivatives
        du = (psidd - u * Hdot) / (H * H)
        dtdN = 1.0 / H
        return dy, du, dtdN, H

    nsteps = int(abs((N1 - N) / hN)) + 1
    h = math.copysign(abs(hN), N1 - N)
    for _ in range(nsteps):
        a = math.exp(N)
        # Output sample at the beginning of the step
        H2 = _H2_from_N_state(a, x, y, psi, u, p)
        H = math.sqrt(H2)
        out["N"].append(N)
        out["a"].append(a)
        out["t"].append(t)
        out["H"].append(H)
        out["x"].append(x)
        theta = theta_from_x(x, p.gamma)
        out["theta"].append(theta)
        theta_dot = theta_dot_from_x(x, y * H, p.gamma)
        out["theta_dot"].append(theta_dot)
        out["theta_eff"].append(theta)
        out["theta_eff_dot"].append(theta_dot)
        out["psi"].append(psi)
        out["psi_dot"].append(u * H)
        out["I_grad"].append(0.0)
        out["J_refr"].append(0.0)
        out["Hdot"].append(_Hdot_from_N_state(a, x, y, psi, u, H, p))

        # RK4 step in N for (x,y,psi,u,t)
        # We pack derivatives explicitly to keep code clear
        def pack(x, y, psi, u, t):
            return (x, y, psi, u, t)

        # k1
        dy1, du1, dt1, H1 = rhs_N(N, pack(x, y, psi, u, t))
        k1 = (y, dy1, u, du1, dt1)

        # k2
        x2 = x + 0.5 * h * k1[0]
        y2 = y + 0.5 * h * k1[1]
        psi2 = psi + 0.5 * h * k1[2]
        u2 = u + 0.5 * h * k1[3]
        t2 = t + 0.5 * h * k1[4]
        dy2, du2, dt2, H2_mid = rhs_N(N + 0.5 * h, pack(x2, y2, psi2, u2, t2))
        k2 = (y2, dy2, u2, du2, dt2)

        # k3
        x3 = x + 0.5 * h * k2[0]
        y3 = y + 0.5 * h * k2[1]
        psi3 = psi + 0.5 * h * k2[2]
        u3 = u + 0.5 * h * k2[3]
        t3 = t + 0.5 * h * k2[4]
        dy3, du3, dt3, H3_mid = rhs_N(N + 0.5 * h, pack(x3, y3, psi3, u3, t3))
        k3 = (y3, dy3, u3, du3, dt3)

        # k4
        x4 = x + h * k3[0]
        y4 = y + h * k3[1]
        psi4 = psi + h * k3[2]
        u4 = u + h * k3[3]
        t4 = t + h * k3[4]
        dy4, du4, dt4, H4 = rhs_N(N + h, pack(x4, y4, psi4, u4, t4))
        k4 = (y4, dy4, u4, du4, dt4)

        # Update
        x += (h / 6.0) * (k1[0] + 2 * k2[0] + 2 * k3[0] + k4[0])
        y += (h / 6.0) * (k1[1] + 2 * k2[1] + 2 * k3[1] + k4[1])
        psi += (h / 6.0) * (k1[2] + 2 * k2[2] + 2 * k3[2] + k4[2])
        u += (h / 6.0) * (k1[3] + 2 * k2[3] + 2 * k3[3] + k4[3])
        t += (h / 6.0) * (k1[4] + 2 * k2[4] + 2 * k3[4] + k4[4])
        N += h

    return out


def rhs(
    t: float, y: Tuple[float, float, float, float, float], p: Params
) -> Tuple[float, float, float, float, float]:
    a, x, xdot, psi, psidot = y
    H = H_total(a, x, xdot, psi, psidot, p)
    xdd = x_eom_rhs(x, xdot, H, p.gamma)
    # Ψ amplitude (R): R̈ + 3H Ṙ - Q^2/(a^6 R^3) + V0 R = 0
    psidd = -3.0 * H * psidot - p.V0 * psi
    if p.Q != 0.0:
        if abs(psi) < 1e-12:
            raise ValueError("R amplitude too small for nonzero Q (divergent R-equation)")
        psidd += (p.Q * p.Q) / (a ** 6 * (psi ** 3))
    adot = a * H
    return adot, xdot, xdd, psidot, psidd


def rk4_step(t: float, y: Tuple[float, ...], h: float, p: Params) -> Tuple[float, Tuple[float, ...]]:
    def add(vec, deriv, coeff=1.0):
        return tuple(vi + coeff * di for vi, di in zip(vec, deriv))

    k1 = rhs(t, y, p)
    k2 = rhs(t + 0.5 * h, add(y, k1, 0.5 * h), p)
    k3 = rhs(t + 0.5 * h, add(y, k2, 0.5 * h), p)
    k4 = rhs(t + h, add(y, k3, h), p)
    y_next = tuple(
        yi + (h / 6.0) * (k1i + 2 * k2i + 2 * k3i + k4i)
        for yi, k1i, k2i, k3i, k4i in zip(y, k1, k2, k3, k4)
    )
    return t + h, y_next


def rk45_step(t: float, y: Tuple[float, ...], h: float, p: Params, atol: float = 1e-9, rtol: float = 1e-6):
    """Single adaptive RK45 (Cash–Karp) step for the background system."""
    # Cash–Karp coefficients
    a2 = 1 / 5; a3 = 3 / 10; a4 = 3 / 5; a5 = 1.0; a6 = 7 / 8
    b21 = 1 / 5
    b31 = 3 / 40; b32 = 9 / 40
    b41 = 3 / 10; b42 = -9 / 10; b43 = 6 / 5
    b51 = -11 / 54; b52 = 5 / 2; b53 = -70 / 27; b54 = 35 / 27
    b61 = 1631 / 55296; b62 = 175 / 512; b63 = 575 / 13824; b64 = 44275 / 110592; b65 = 253 / 4096

    c1 = 37 / 378; c3 = 250 / 621; c4 = 125 / 594; c6 = 512 / 1771
    c1s = 2825 / 27648; c3s = 18575 / 48384; c4s = 13525 / 55296; c5s = 277 / 14336; c6s = 1 / 4

    def add(vec, deriv, coeff=1.0):
        return tuple(vi + coeff * di for vi, di in zip(vec, deriv))

    k1 = rhs(t, y, p)
    y2 = add(y, k1, h * b21)
    k2 = rhs(t + a2 * h, y2, p)

    y3 = tuple(yi + h * (b31 * k1i + b32 * k2i) for yi, k1i, k2i in zip(y, k1, k2))
    k3 = rhs(t + a3 * h, y3, p)

    y4 = tuple(yi + h * (b41 * k1i + b42 * k2i + b43 * k3i) for yi, k1i, k2i, k3i in zip(y, k1, k2, k3))
    k4 = rhs(t + a4 * h, y4, p)

    y5 = tuple(yi + h * (b51 * k1i + b52 * k2i + b53 * k3i + b54 * k4i) for yi, k1i, k2i, k3i, k4i in zip(y, k1, k2, k3, k4))
    k5 = rhs(t + a5 * h, y5, p)

    y6 = tuple(yi + h * (b61 * k1i + b62 * k2i + b63 * k3i + b64 * k4i + b65 * k5i) for yi, k1i, k2i, k3i, k4i, k5i in zip(y, k1, k2, k3, k4, k5))
    k6 = rhs(t + a6 * h, y6, p)

    y4th = tuple(yi + h * (c1 * k1i + c3 * k3i + c4 * k4i + c6 * k6i) for yi, k1i, k3i, k4i, k6i in zip(y, k1, k3, k4, k6))
    y5th = tuple(yi + h * (c1s * k1i + c3s * k3i + c4s * k4i + c5s * k5i + c6s * k6i) for yi, k1i, k3i, k4i, k5i, k6i in zip(y, k1, k3, k4, k5, k6))

    err = max(abs(a - b) for a, b in zip(y5th, y4th))
    tol = atol + rtol * max(max(abs(val) for val in y), max(abs(val) for val in y5th))
    accept = (err <= tol) or (h <= 1e-16)
    safety = 0.9
    if err == 0.0:
        h_new = h * 2.0
    else:
        h_new = h * safety * (tol / err) ** 0.2
        h_new = max(min(h_new, 2.5 * h), 0.1 * h)
    return y5th, err, accept, h_new


def integrate_background_rk45(
    t0: float,
    y0: Tuple[float, float, float, float, float],
    t1: float,
    p: Params,
    h0: float = 1e-2,
    atol: float = 1e-9,
    rtol: float = 1e-6,
) -> Dict[str, List[float]]:
    """Adaptive background integration using RK45 (Cash–Karp)."""
    t = t0
    h = math.copysign(abs(h0), (t1 - t0))
    y = y0
    out: Dict[str, List[float]] = {
        "t": [],
        "a": [],
        "H": [],
        "Hdot": [],
        "x": [],
        "theta": [],
        "theta_dot": [],
        "theta_eff": [],
        "theta_eff_dot": [],
        "psi": [],
        "psi_dot": [],
    }
    steps = 0
    while (t < t1 and h > 0) or (t > t1 and h < 0):
        a, x, xdot, psi, psidot = y
        H = H_total(a, x, xdot, psi, psidot, p)
        out["t"].append(t)
        out["a"].append(a)
        out["H"].append(H)
        out["x"].append(x)
        theta = theta_from_x(x, p.gamma)
        out["theta"].append(theta)
        theta_dot = theta_dot_from_x(x, xdot, p.gamma)
        out["theta_dot"].append(theta_dot)
        out["theta_eff"].append(theta)
        out["theta_eff_dot"].append(theta_dot)
        out["psi"].append(psi)
        out["psi_dot"].append(psidot)
        out["Hdot"].append(Hdot_total(a, x, xdot, psi, psidot, H, p))

        y_new, err, accept, h_new = rk45_step(t, y, h, p, atol=atol, rtol=rtol)
        if accept:
            t = t + h
            y = y_new
        h = h_new
        steps += 1
        if steps > 100000:
            raise RuntimeError("Exceeded max steps in adaptive background integrator")
    return out


def integrate_background(
    t0: float,
    y0: Tuple[float, float, float, float, float],
    h: float,
    nsteps: int,
    p: Params,
) -> Dict[str, List[float]]:
    t = t0
    y = y0
    out: Dict[str, List[float]] = {
        "t": [],
        "a": [],
        "H": [],
        "Hdot": [],
        "x": [],
        "theta": [],
        "theta_dot": [],
        "theta_eff": [],
        "theta_eff_dot": [],
        "psi": [],
        "psi_dot": [],
    }
    for _ in range(nsteps):
        a, x, xdot, psi, psidot = y
        H = H_total(a, x, xdot, psi, psidot, p)
        out["t"].append(t)
        out["a"].append(a)
        out["H"].append(H)
        out["x"].append(x)
        theta = theta_from_x(x, p.gamma)
        out["theta"].append(theta)
        theta_dot = theta_dot_from_x(x, xdot, p.gamma)
        out["theta_dot"].append(theta_dot)
        out["theta_eff"].append(theta)
        out["theta_eff_dot"].append(theta_dot)
        out["psi"].append(psi)
        out["psi_dot"].append(psidot)
        out["Hdot"].append(Hdot_total(a, x, xdot, psi, psidot, H, p))
        t, y = rk4_step(t, y, h, p)
    return out


def demo() -> Dict[str, float]:
    """
    Быстрый прогон: стартуем в ранней радиационной эпохе и смотрим,
    как поле ведёт себя на фоне стандартной ΛCDM смеси.
    """
    params = Params()
    y0 = (1e-3, 0.1, 0.0, 0.01, 0.0)  # (a, x, xdot, psi, psidot)
    out = integrate_background(0.0, y0, h=0.01, nsteps=200, p=params)
    return {
        "a_end": out["a"][-1],
        "H_end": out["H"][-1],
        "theta_end": out["theta"][-1],
    }


if __name__ == "__main__":
    print(demo())
