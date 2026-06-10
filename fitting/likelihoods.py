# -*- coding: utf-8 -*-
"""
Gaussian likelihoods for simple background datasets.
"""
from __future__ import annotations

from typing import List
import math


def chi2_gaussian(model: List[float], data: List[float], sigma: List[float]) -> float:
    assert len(model) == len(data) == len(sigma)
    return sum(((m - d) / s) ** 2 for m, d, s in zip(model, data, sigma))


def loglike_gaussian(model: List[float], data: List[float], sigma: List[float]) -> float:
    return -0.5 * chi2_gaussian(model, data, sigma)

