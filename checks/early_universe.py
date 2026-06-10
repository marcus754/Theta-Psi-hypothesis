# -*- coding: utf-8 -*-
"""Early-Universe sanity checks."""
from __future__ import annotations
from typing import Callable


def bbn_safe(H_of_a: Callable[[float], float]) -> bool:
    """
    Crude requirement for Big-Bang-nucleosynthesis safety.

    We demand that the expansion rate H(a) monotonically decreases between
    a = 1e-10 and a = 1e-6, which is a minimal expectation for a radiation
    dominated epoch.
    """
    a_vals = [10 ** (-10 + i * (4.0 / 20.0)) for i in range(21)]  # log grid
    H_vals = [H_of_a(a) for a in a_vals]
    return all(H_vals[i] >= H_vals[i + 1] for i in range(len(H_vals) - 1))

