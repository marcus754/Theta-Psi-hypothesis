import math

import pytest

from background.frw_symbolic import (
    chi_dot_from_theta,
    chi_from_theta,
    lagrangian,
    lagrangian_chi,
    theta_dot_from_chi,
    theta_from_chi,
)


def test_chi_theta_round_trip():
    theta = 0.42
    thetadot = 0.013
    chi = chi_from_theta(theta)
    chidot = chi_dot_from_theta(theta, thetadot)

    assert math.isclose(theta_from_chi(chi), theta, rel_tol=1e-12)
    assert math.isclose(theta_dot_from_chi(chi, chidot), thetadot, rel_tol=1e-12)


def test_lagrangian_chi_matches_theta():
    theta = 0.33
    thetadot = 0.025
    chi = chi_from_theta(theta)
    chidot = chi_dot_from_theta(theta, thetadot)
    psi = 0.2
    psidot = -0.07
    gamma = 0.7
    V0 = 0.2

    lag_theta = lagrangian(theta, thetadot, psi, psidot, gamma, V0)
    lag_chi = lagrangian_chi(chi, chidot, psi, psidot, gamma, V0)
    assert math.isclose(lag_theta, lag_chi, rel_tol=1e-12)


def test_chi_from_theta_requires_positive():
    with pytest.raises(ValueError):
        chi_from_theta(0.0)
    with pytest.raises(ValueError):
        chi_dot_from_theta(-0.1, 0.5)
