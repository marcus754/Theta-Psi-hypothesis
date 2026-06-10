# -*- coding: utf-8 -*-
"""
Thin access layer from observational code to the physical kernel.

The goal of this module is architectural, not physical: `fitting/` and
`scripts/` should depend on one stable API surface instead of importing
`background/` and `checks/` internals all over the tree.
"""
from __future__ import annotations

from background.components import H0_from_km_s_Mpc, Omegas, omegas_from_km_s_Mpc
from background.frw_background import Params as BackgroundParams
from background.frw_background import integrate_background, integrate_background_lnN, integrate_background_rk45
from checks.compact_star import CompactStarParams, scan_compact_star_grid, solve_compact_star_profile

__all__ = [
    "BackgroundParams",
    "CompactStarParams",
    "H0_from_km_s_Mpc",
    "Omegas",
    "integrate_background",
    "integrate_background_lnN",
    "integrate_background_rk45",
    "omegas_from_km_s_Mpc",
    "scan_compact_star_grid",
    "solve_compact_star_profile",
]
