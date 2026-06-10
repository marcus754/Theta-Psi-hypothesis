# -*- coding: utf-8 -*-
from src.linear_modes_adm import cs2_eigs_adm, dispersion_eigs_adm, stability_finite_k_adm


def test_sigma_mode_extends_spectrum():
    base = cs2_eigs_adm(0.2, 0.1, 1.0, 1.0, 1e-3, 0.0, 0.1)
    ext = cs2_eigs_adm(0.2, 0.1, 1.0, 1.0, 1e-3, 0.0, 0.1, include_sigma=True, csigma2=0.8)
    assert len(base) == 3
    assert len(ext) == 4
    assert min(ext) <= 0.8 + 1e-12


def test_sigma_mode_finite_k_stability():
    w = dispersion_eigs_adm(
        0.2, 0.1, 1.0, 1.0, 1e-3, 0.0, 0.1,
        include_sigma=True, csigma2=0.9, msigma2=0.02,
    )
    assert len(w) == 4
    assert all(val > -1e-12 for val in w)

    res = stability_finite_k_adm(
        0.2, 0.1, 1.0, 1.0, 1e-3, 0.0, 0.1,
        include_sigma=True, csigma2=0.9, msigma2=0.02,
    )
    assert res["no_ghost"]
    assert res["no_tachyon"]
