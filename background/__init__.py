"""FRW background module for the reduced Θ–Ψ cosmology.

This package belongs to the physical kernel. It should not know about
head-to-head comparisons, falsifiability policy or external baselines.
"""

from .components import Omegas, rho_r, rho_m, rho_L, H_of_a, E2

__all__ = ["Omegas", "rho_r", "rho_m", "rho_L", "H_of_a", "E2"]
