# -*- coding: utf-8 -*-
"""
Standard FRW background components.

The helper functions below keep track of the usual radiation + matter + Λ
content and expose utilities that operate purely on the scale factor `a`.

Units: all quantities are assumed to be in a consistent system where the
Friedmann equation reads `H^2 = ρ_tot / 3` (e.g. reduced Planck units with
8πG = 1 and c = 1). The parameter `H0` must be supplied in these units by the
caller. The default value is a dimensionless placeholder for toy runs; if you
use km/s/Mpc in applications, convert to the chosen unit system before passing
`H0` here.
"""
from __future__ import annotations
from dataclasses import dataclass
import math


@dataclass
class Omegas:
    """
    Background density parameters at `a = 1`.

    Parameters are given as fractional energy densities.  Separate baryon and
    cold-matter pieces are allowed for future bookkeeping, but for the simple
    functions below only their sum matters.
    """

    H0: float = 1.0                     # dimensionless H0 in the chosen unit system (user must convert if needed)
    Omega_r: float = 8.4e-5             # effective radiation today (photons + neutrinos)
    Omega_b: float = 0.049              # baryons
    Omega_c: float = 0.251              # cold dark matter
    Omega_L: float = 0.7                # cosmological constant / dark energy

    @property
    def Omega_m(self) -> float:
        return self.Omega_b + self.Omega_c


def E2(a: float, om: Omegas) -> float:
    """Return E(a)^2 = H(a)^2 / H0^2 for the given scale factor."""
    am3 = a ** -3
    am4 = a ** -4
    return om.Omega_r * am4 + om.Omega_m * am3 + om.Omega_L


def H_of_a(a: float, om: Omegas) -> float:
    """Hubble parameter H(a) using the supplied background."""
    return om.H0 * math.sqrt(E2(a, om))


def rho_r(a: float, om: Omegas) -> float:
    """Radiation energy density ρ_r(a) in reduced Planck units."""
    return 3.0 * (om.H0 ** 2) * om.Omega_r * (a ** -4)


def rho_m(a: float, om: Omegas) -> float:
    """Matter energy density ρ_m(a) (baryons + cold matter)."""
    return 3.0 * (om.H0 ** 2) * om.Omega_m * (a ** -3)


def rho_L(a: float, om: Omegas) -> float:
    """Vacuum energy density ρ_Λ(a) (independent of a)."""
    return 3.0 * (om.H0 ** 2) * om.Omega_L


# ---------------------------------------------------------------------------
# Unit conversion helpers (optional, for convenience)

# 1 Mpc in meters (CODATA 2018 precision is sufficient for our needs)
_MPC_IN_M: float = 3.085677581491367e22
# Reduced Planck time t̄_P = √(8π) * t_P ≈ 2.701177... × 10^-43 s
# We embed a constant to avoid carrying fundamental constants around.
_T_REDUCED_PLANCK_S: float = 2.70117733e-43


def H0_from_km_s_Mpc(H0_km_s_Mpc: float) -> float:
    """
    Convert H0 given in km/s/Mpc to the dimensionless H0 used by this module
    (units with H^2 = ρ/3, i.e. 8πG = 1, c = 1).

    Steps: km/s/Mpc → s^-1, then multiply by reduced Planck time t̄_P.
    Example: 70 km/s/Mpc → ~6.13e-61 in these units.
    """
    H0_per_s = (H0_km_s_Mpc * 1000.0) / _MPC_IN_M
    return H0_per_s * _T_REDUCED_PLANCK_S


def omegas_from_km_s_Mpc(H0_km_s_Mpc: float = 70.0, *, Omega_r: float = 8.4e-5, Omega_b: float = 0.049, Omega_c: float = 0.251, Omega_L: float = 0.7) -> Omegas:
    """
    Factory returning `Omegas` with H0 converted from km/s/Mpc to the
    dimensionless units used in this module (8πG=1, c=1).

    Example:
        om = omegas_from_km_s_Mpc(70.0)  # close to ΛCDM defaults
    """
    H0_dimless = H0_from_km_s_Mpc(H0_km_s_Mpc)
    return Omegas(H0=H0_dimless, Omega_r=Omega_r, Omega_b=Omega_b, Omega_c=Omega_c, Omega_L=Omega_L)
