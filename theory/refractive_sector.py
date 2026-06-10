# -*- coding: utf-8 -*-
"""
Minimal helpers for the refractive sector L_refr = M_*^4 F(J_refr).

The canonical normalized representative is

    F_v1(J) = 1 - exp(-J),

which satisfies
    F(0)=0, F'(J)>0, F''(J)<=0, F(+inf)=1.

Any explicit scale `J_*` is a diagnostic normalization choice that can be
absorbed into the definition of `J_refr`. The current action derives the
argument `J_refr` and admissibility constraints on `F`, but not the remaining
shape beyond those constraints.
"""
from __future__ import annotations

import math


def f_v1(j_refr: float) -> float:
    """Canonical normalized saturating representative F_v1(J)."""
    j = max(0.0, float(j_refr))
    return 1.0 - math.exp(-j)


def f_v1_with_scale(j_refr: float, *, j_star: float = 1.0) -> float:
    """Diagnostic saturating representative with an explicit scale."""
    js = max(float(j_star), 1e-12)
    j = max(0.0, float(j_refr))
    return 1.0 - math.exp(-j / js)


def f_v1_prime(j_refr: float) -> float:
    """First derivative of the canonical normalized representative."""
    j = max(0.0, float(j_refr))
    return math.exp(-j)


def f_v1_prime_with_scale(j_refr: float, *, j_star: float = 1.0) -> float:
    """Diagnostic first derivative for the scaled representative."""
    js = max(float(j_star), 1e-12)
    j = max(0.0, float(j_refr))
    return math.exp(-j / js) / js


def f_v1_double_prime(j_refr: float) -> float:
    """Second derivative of the canonical normalized representative."""
    j = max(0.0, float(j_refr))
    return -math.exp(-j)


def f_v1_double_prime_with_scale(j_refr: float, *, j_star: float = 1.0) -> float:
    """Diagnostic second derivative for the scaled representative."""
    js = max(float(j_star), 1e-12)
    j = max(0.0, float(j_refr))
    return -math.exp(-j / js) / (js * js)
