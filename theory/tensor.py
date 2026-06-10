# -*- coding: utf-8 -*-
"""
Tensor-sector bookkeeping.

Baseline statement (optical metric, docs/14, 15):
- A single universal index n(θ) induces an effective (optical) metric \tilde g.
- The low-energy universal response includes an Einstein–Hilbert term for \tilde g,
  which yields two transverse-traceless tensor polarizations with luminal speed:
  c_T = 1. On FRW backgrounds the TT equation reads (in conformal time):
    h''_{ij} + 2 ℋ h'_{ij} + k^2 h_{ij} = 0.

This module keeps a minimal API (c_T^2, optional damping) so downstream code
has a single place to query tensor properties if/when the sector is extended.
"""
from dataclasses import dataclass


@dataclass
class TensorProps:
    c_T2: float = 1.0  # gravitational-wave speed squared
    nu: float = 0.0    # running/attenuation parameter (unused placeholder)


def tensor_properties() -> TensorProps:
    """Return default tensor-sector properties."""
    return TensorProps()
