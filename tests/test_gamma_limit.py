import math

import background.frw_background as bg


def _x_from_theta(theta: float, gamma: float) -> float:
    # Clamp argument for atanh to the open interval (-1,1)
    arg = gamma * theta / math.sqrt(3.0)
    arg = max(min(arg, 1.0 - 1e-15), -1.0 + 1e-15)
    return 0.5 * math.log((1 + arg) / (1 - arg))


def test_background_insensitive_when_theta_fixed_and_static():
    # If we set identical (theta, theta_dot, R, Rdot) conditions via x-variables,
    # and keep R=0 and theta_dot=0, short-run H(a) should match regardless of gamma.
    a0 = 1e-3
    theta0 = 1e-6
    theta_dot0 = 0.0
    R0 = 0.0
    Rdot0 = 0.0

    # Two gammas: O(1) and very small
    g1 = 1.0
    g2 = 1e-6

    # Map desired (theta0, theta_dot0) to x-variables for each gamma
    x1 = _x_from_theta(theta0, g1)
    x2 = _x_from_theta(theta0, g2)
    # theta_dot = (sqrt(3)/gamma) * sech^2(x) * xdot ⇒ xdot = 0 for theta_dot=0
    xdot1 = 0.0
    xdot2 = 0.0

    p1 = bg.Params(gamma=g1, V0=0.1)
    p2 = bg.Params(gamma=g2, V0=0.1)

    N0 = math.log(a0)
    N1 = N0 + 0.05
    out1 = bg.integrate_background_lnN(a0, x1, xdot1, R0, Rdot0, N1=N1, p=p1, hN=0.01)
    out2 = bg.integrate_background_lnN(a0, x2, xdot2, R0, Rdot0, N1=N1, p=p2, hN=0.01)

    # Compare H trajectories (relative difference)
    H1, H2 = out1["H"], out2["H"]
    assert len(H1) == len(H2)
    rel_diffs = [abs(h1 - h2) / max(1e-18, abs(h1)) for h1, h2 in zip(H1, H2)]
    assert max(rel_diffs) < 1e-9
