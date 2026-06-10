# -*- coding: utf-8 -*-
"""Simple EFT-style cutoff checks."""
from __future__ import annotations
from typing import Iterable


def cutoff_ok(
    theta_list: Iterable[float],
    thetadot_list: Iterable[float],
    gamma: float,
    cutoff: float = 1e-1,
) -> bool:
    """
    Enforce |θ̇/θ| < Λ and |γ θ̇| < Λ along a trajectory.

    Parameters
    ----------
    theta_list, thetadot_list
        Iterables with the background values sampled along the trajectory.
    gamma
        Coupling parameter of the model.
    cutoff
        Crude EFT cutoff Λ; default 0.1 in Planck units.
    """
    for theta, thetadot in zip(theta_list, thetadot_list):
        if abs(theta) < 1e-12:
            # If theta crosses zero the minisuperspace reduction already broke down.
            return False
        if abs(thetadot / theta) >= cutoff:
            return False
        if abs(gamma * thetadot) >= cutoff:
            return False
    return True

