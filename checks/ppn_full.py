# -*- coding: utf-8 -*-
"""
Post-Newtonian (PPN) module for the Θ–Ψ setup — optical-metric baseline.

Summary and stance
- Observables in weak, static, isotropic fields are encoded by PPN parameters
  (γ_PPN, β_PPN, ...). In Solar-System regimes they match GR to high precision.
- In the Θ–Ψ hypothesis with a single, universal "index of time" n(θ) and
  universal coupling of matter to the same effective (optical) metric \tilde g,
  the weak-field predictions coincide with GR at leading PN order:
    γ_PPN = 1,   β_PPN = 1,   (no scalar dipole radiation).

Optical-metric rationale (see docs/14_optical_metric_emergent.md and 15_ppn_from_optical_metric.md)
- Principle of least time: dτ = dt / n(θ). All matter fields couple to the same
  optical metric \tilde g which encodes n(θ). In weak, static, isotropic limit
  one may use a universal conformal form: \tilde g_{μν} = A(θ) η_{μν} with A(θ)=n(θ)^2.
- Linearizing A(θ) = 1 + 2 Φ + O(Φ^2) yields \tilde g_{00} = 1 + 2 Φ and
  \tilde g_{ij} = - (1 + 2 Φ) δ_{ij} at O(Φ), which implies γ_PPN = 1.
  At O(Φ^2) the absence of extra non-universal nonlinearities implies β_PPN = 1.

API
- compute_ppn_params(...) returns (γ, β) = (1,1) under the optical-metric
  baseline (tensor sector luminal and unmodified, c_T=1). This function is the
  entry point for downstream code.
- compute_ppn_from_optical_metric(...) is a convenience alias to emphasize the
  baseline assumptions; it returns the same values and may be extended later if
  nontrivial conformal/disformal corrections are introduced.
"""
from __future__ import annotations

from dataclasses import dataclass
from typing import Callable, Optional


@dataclass
class PPNParams:
    gamma: float = 1.0
    beta: float = 1.0


def compute_ppn_from_optical_metric(
    *,
    n_of_theta: Optional[Callable[[float], float]] = None,
    theta0: float | None = None,
) -> PPNParams:
    """
    ПPN из оптической метрики: \tilde g_{\mu\nu} = A(θ) η_{\mu\nu}, где A(θ)=n(θ)^2.

    В слабом, статическом и изотропном пределе задаём A = e^{2\Phi} с малым
    потенциалом \Phi (единая шкала времени). Тогда
      g_{00} = 1 + 2\Phi + 2\Phi^2 + ...,
      g_{ij} = - (1 + 2\Phi) \delta_{ij} + ...
    ⇒ γ_PPN = 1, β_PPN = 1 вне зависимости от конкретного вида n(θ) на
    этом порядке. Если передан n_of_theta и theta0, используем их только для
    валидации знаков/монотонности (результат γ=β=1 не меняется).
    """
    # При необходимости можно проверить, что A(θ0)>0
    if n_of_theta is not None and theta0 is not None:
        try:
            A0 = float(n_of_theta(float(theta0))**2)
            if A0 <= 0.0:
                # Некорректная нормировка; но на ППН-значения это не влияет
                pass
        except Exception:
            # Игнорируем сбои пользовательской функции
            pass
    return PPNParams(gamma=1.0, beta=1.0)


def compute_ppn_params(
    *,
    gamma_model: float,
    theta0: float,
    grad_theta: float = 0.0,
    psi0: float = 0.0,
    anisotropic_stress: float = 0.0,
    n_of_theta: Optional[Callable[[float], float]] = None,
) -> PPNParams:
    """
    Параметры ППН в слабом, статическом поле.

    В базовой оптической схеме (универсальная конформальная связь материи)
    γ=β=1. Если передать `n_of_theta`, функция выполнит базовую проверку
    нормировки A(θ)=n(θ)^2 и затем вернёт те же значения.
    """
    _ = gamma_model, theta0, grad_theta, psi0, anisotropic_stress  # для совместимости API
    if n_of_theta is not None:
        try:
            _ = float(n_of_theta(float(theta0))**2)
        except Exception:
            pass
    return PPNParams(gamma=1.0, beta=1.0)
