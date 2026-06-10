
"""
FRW integrator for the Θ–Ψ model — robust & stable version
- Canonical sign for R:  +½ Ṙ² - ½ V0 R²
- Variable transform to eliminate the pole at 1 - γ² θ² / 3 = 0:
    θ = (√3 / γ) * tanh x      (regularization variable only)
    θ̇ = (√3 / γ) * sech²(x) * ẋ
- Dynamics in x:
    ẍ + 3H ẋ = ẋ² ( 2 tanh x + coth x )
    (derived from θ̈ + 3H θ̇ = θ̇² / [ θ (1 - γ² θ² / 3) ])
- Background H(t) may be supplied by user; default H ≡ 0 for testing.

Integrator: adaptive Runge–Kutta Cash–Karp (RK45 embedded)
with step-size control. No artificial caps are imposed on x; the tanh-map is a coordinate regularization, not a physical constraint.
"""

from __future__ import annotations
from typing import Callable, Dict, Tuple, List, Optional
import math
import csv

from background.frw_symbolic import (
    theta_from_x as _theta_from_x,
    theta_dot_from_x as _theta_dot_from_x,
    x_eom_rhs,
)

SQRT3 = math.sqrt(3.0)

# Model parameters container
class Params:
    def __init__(self, gamma: float, V0: float, H_of_t: Optional[Callable[[float], float]] = None):
        assert gamma != 0.0, "gamma must be non-zero"
        assert V0 >= 0.0, "V0 must be ≥ 0"
        self.gamma = float(gamma)
        self.V0 = float(V0)
        self.H_of_t = H_of_t or (lambda t: 0.0)

    @property
    def theta_star(self) -> float:
        return SQRT3/abs(self.gamma)

# State helpers --------------------------------------------------------------
def pack_state(x, xdot, R, Rdot) -> List[float]:
    return [x, xdot, R, Rdot]

def unpack_state(y: List[float], p: Params):
    x, xdot, R, Rdot = y
    return x, xdot, R, Rdot

def theta_from_x(x: float, p: Params) -> float:
    return _theta_from_x(x, p.gamma)

def theta_dot_from_x(x: float, xdot: float, p: Params) -> float:
    return _theta_dot_from_x(x, xdot, p.gamma)

# RHS of ODEs ---------------------------------------------------------------
def rhs(t: float, y: List[float], p: Params) -> List[float]:
    x, xdot, R, Rdot = unpack_state(y, p)
    H = p.H_of_t(t)

    xddot = x_eom_rhs(x, xdot, H, p.gamma)

    # R: R̈ + 3H Ṙ = - V0 R   (canonical sign)
    Rddot = -3.0*H*Rdot - p.V0*R

    return [xdot, xddot, Rdot, Rddot]

# Optional first-principles RHS in θ-variables (no coordinate transform) -----
def rhs_theta(t: float, y: list[float], p: Params) -> list[float]:
    """
    Direct evolution in (theta, thetadot, R, Rdot) without x-regularization.

    Equation (first principles):
        θ̈ + 3H θ̇ = θ̇² / [ θ (1 - γ² θ² / 3) ].
    R remains canonical: R̈ + 3H Ṙ + V0 R = 0.
    """
    theta, thetadot, R, Rdot = y
    H = p.H_of_t(t)
    denom = (1.0 - (p.gamma * p.gamma) * (theta * theta) / 3.0)
    thetaddot = -3.0 * H * thetadot + (thetadot * thetadot) / (theta * denom)
    Rddot = -3.0 * H * Rdot - p.V0 * R
    return [thetadot, thetaddot, Rdot, Rddot]


# Adaptive RK45 (Cash–Karp) -------------------------------------------------
def rk45_step(fun, t, y, h, p, atol=1e-9, rtol=1e-6):
    # Cash–Karp coefficients
    a2=1/5; a3=3/10; a4=3/5; a5=1.0; a6=7/8
    b21=1/5
    b31=3/40; b32=9/40
    b41=3/10; b42=-9/10; b43=6/5
    b51=-11/54; b52=5/2; b53=-70/27; b54=35/27
    b61=1631/55296; b62=175/512; b63=575/13824; b64=44275/110592; b65=253/4096

    c1=37/378; c3=250/621; c4=125/594; c6=512/1771
    c1s=2825/27648; c3s=18575/48384; c4s=13525/55296; c5s=277/14336; c6s=1/4

    k1 = fun(t, y, p)
    y2 = [yi + h*b21*k1i for yi,k1i in zip(y,k1)]
    k2 = fun(t + a2*h, y2, p)

    y3 = [yi + h*(b31*k1i + b32*k2i) for yi,k1i,k2i in zip(y,k1,k2)]
    k3 = fun(t + a3*h, y3, p)

    y4 = [yi + h*(b41*k1i + b42*k2i + b43*k3i) for yi,k1i,k2i,k3i in zip(y,k1,k2,k3)]
    k4 = fun(t + a4*h, y4, p)

    y5 = [yi + h*(b51*k1i + b52*k2i + b53*k3i + b54*k4i) for yi,k1i,k2i,k3i,k4i in zip(y,k1,k2,k3,k4)]
    k5 = fun(t + a5*h, y5, p)

    y6 = [yi + h*(b61*k1i + b62*k2i + b63*k3i + b64*k4i + b65*k5i) for yi,k1i,k2i,k3i,k4i,k5i in zip(y,k1,k2,k3,k4,k5)]
    k6 = fun(t + a6*h, y6, p)

    y4th = [yi + h*(c1*k1i + c3*k3i + c4*k4i + c6*k6i) for yi,k1i,k3i,k4i,k6i in zip(y,k1,k3,k4,k6)]
    y5th = [yi + h*(c1s*k1i + c3s*k3i + c4s*k4i + c5s*k5i + c6s*k6i) for yi,k1i,k3i,k4i,k5i,k6i in zip(y,k1,k3,k4,k5,k6)]

    # Error estimate
    err = max(abs(a-b) for a,b in zip(y5th, y4th))
    # Tolerance
    tol = atol + rtol * max(max(abs(val) for val in y), max(abs(val) for val in y5th))
    # Step acceptance
    accept = err <= tol or h <= 1e-16
    # New step size
    safety = 0.9
    if err == 0.0:
        h_new = h*2.0
    else:
        h_new = h * safety * (tol/err) ** 0.2  # 1/(order+1) with order=4
        h_new = max(min(h_new, 2.5*h), 0.1*h)
    return y5th, err, accept, h_new

def integrate(t0: float, t1: float, y0: List[float], params: Params,
              h0: float=1e-2, atol: float=1e-9, rtol: float=1e-6, max_steps: int=100000):
    t = t0
    h = math.copysign(abs(h0), t1-t0)
    y = list(y0)
    traj = [(t, y.copy())]
    steps = 0

    while (t < t1 and h > 0) or (t > t1 and h < 0):
        y_new, err, accept, h_new = rk45_step(rhs, t, y, h, params, atol, rtol)
        if accept:
            t = t + h
            y = y_new
            traj.append((t, y.copy()))
        h = h_new
        steps += 1
        if steps > max_steps:
            raise RuntimeError("Exceeded max_steps — likely stiffness or divergence.")
    return traj

# Optional integrator directly in θ ------------------------------------------
def integrate_theta(
    t0: float,
    t1: float,
    y0: list[float],
    params: Params,
    *,
    h0: float = 1e-2,
    atol: float = 1e-9,
    rtol: float = 1e-6,
    max_steps: int = 100000,
):
    """
    Adaptive RK45 integration in (theta, thetadot, R, Rdot) directly.

    This variant imposes no coordinate regularization. Any divergence near
    1 - γ² θ² / 3 → 0 arises dynamically from the equation itself.
    """
    t = t0
    h = (h0 if t1>=t0 else -abs(h0))
    y = list(y0)
    traj = [(t, y.copy())]
    steps = 0
    while (t < t1 and h > 0) or (t > t1 and h < 0):
        y_new, err, accept, h_new = rk45_step(rhs_theta, t, y, h, params, atol, rtol)
        if accept:
            t = t + h
            y = y_new
            traj.append((t, y.copy()))
        h = h_new
        steps += 1
        if steps > max_steps:
            raise RuntimeError("Exceeded max_steps — likely stiffness or divergence.")
    return traj

# Convenience routine --------------------------------------------------------
def run_demo():
    # Constant H as a toy background
    H0 = 0.05
    params = Params(gamma=1.0, V0=0.1, H_of_t=lambda t: H0)
    # Initial conditions in x, not theta
    x0 = 0.1
    xdot0 = 0.01
    R0 = 0.05
    Rdot0 = 0.0
    y0 = pack_state(x0, xdot0, R0, Rdot0)
    traj = integrate(0.0, 10.0, y0, params, h0=1e-2)

    # Export CSV
    try:
        with open("frw_demo_trajectory.csv", "w", newline="") as f:
            w = csv.writer(f)
            w.writerow(["t", "theta", "theta_dot", "R", "R_dot", "H"])
            for t, state in traj:
                x, xdot, R, Rdot = unpack_state(state, params)
                th = theta_from_x(x, params)
                thdot = theta_dot_from_x(x, xdot, params)
                w.writerow([t, th, thdot, R, Rdot, params.H_of_t(t)])
        print("Saved demo trajectory to frw_demo_trajectory.csv")
    except Exception as e:
        print("CSV save failed:", e)

if __name__ == "__main__":
    run_demo()
