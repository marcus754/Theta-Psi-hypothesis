from pathlib import Path

from fitting.run_mcmc import Params, logposterior_combined, loglikes_components
from fitting.data_io import load_cc_dimless, load_sn_dimless, load_bao_dimless
from background.components import Omegas
from fitting.likelihoods_compact import CompactTargets, loglike_compact
from fitting.likelihoods_compact_star import CompactStarTargets
from fitting.likelihoods_forward import ForwardTargets
from fitting.strong_field_link import derive_compact_star_params_from_background


def test_loglikes_components_sum_matches_total():
    cc = load_cc_dimless(Path("fitting/data/cc_demo_dimless.csv"))
    sn = load_sn_dimless(Path("fitting/data/sn_demo_dimless.csv"))
    bao = load_bao_dimless(Path("fitting/data/bao_demo_dimless.csv"))
    om = Omegas(H0=1.0)
    p = Params(gamma=1.0, V0=0.1, Q=0.0)
    total = logposterior_combined(p, om=om, cc=cc, sn=sn, bao=bao)
    comps = loglikes_components(p, om=om, cc=cc, sn=sn, bao=bao)
    # Prior=0 in allowed region, so total equals sum of components
    s = sum(comps.values())
    assert abs(total - s) < 1e-9


def test_loglikes_components_sum_matches_total_with_compact():
    cc = load_cc_dimless(Path("fitting/data/cc_demo_dimless.csv"))
    sn = load_sn_dimless(Path("fitting/data/sn_demo_dimless.csv"))
    bao = load_bao_dimless(Path("fitting/data/bao_demo_dimless.csv"))
    om = Omegas(H0=1.0)
    p = Params(gamma=1.0, V0=0.1, Q=0.0)
    compact = CompactTargets(
        z_surface=20.0,
        sigma_z_surface=10.0,
        delay_proxy=200.0,
        sigma_delay_proxy=100.0,
        n_peak=21.0,
        sigma_n_peak=10.0,
    )
    total = logposterior_combined(
        p, om=om, cc=cc, sn=sn, bao=bao, compact=compact,
        compact_theta_scale=1.0
    )
    comps = loglikes_components(
        p, om=om, cc=cc, sn=sn, bao=bao, compact=compact,
        compact_theta_scale=1.0
    )
    s = sum(comps.values())
    assert abs(total - s) < 1e-9


def test_loglikes_components_sum_matches_total_with_compact_star():
    cc = load_cc_dimless(Path("fitting/data/cc_demo_dimless.csv"))
    sn = load_sn_dimless(Path("fitting/data/sn_demo_dimless.csv"))
    bao = load_bao_dimless(Path("fitting/data/bao_demo_dimless.csv"))
    om = Omegas(H0=1.0)
    p = Params(gamma=1.0, V0=0.1, Q=0.0)
    cstar = CompactStarTargets(
        z_surface=350.0,
        sigma_z_surface=200.0,
        delay_proxy=700.0,
        sigma_delay_proxy=400.0,
    )
    total = logposterior_combined(
        p, om=om, cc=cc, sn=sn, bao=bao, compact_star=cstar
    )
    comps = loglikes_components(
        p, om=om, cc=cc, sn=sn, bao=bao, compact_star=cstar
    )
    s = sum(comps.values())
    assert abs(total - s) < 1e-9


def test_compact_star_component_changes_with_params_when_linked():
    om = Omegas(H0=1.0)
    cstar = CompactStarTargets(
        z_surface=350.0,
        sigma_z_surface=200.0,
        delay_proxy=700.0,
        sigma_delay_proxy=400.0,
    )
    p1 = Params(gamma=0.5, V0=0.05, Q=0.0)
    p2 = Params(gamma=2.0, V0=0.5, Q=0.0)
    c1 = loglikes_components(
        p1, om=om, compact_star=cstar, compact_star_link_model=True
    )["compact_star"]
    c2 = loglikes_components(
        p2, om=om, compact_star=cstar, compact_star_link_model=True
    )["compact_star"]
    assert c1 != c2


def test_compact_constraints_loglike_is_one_sided():
    targets = CompactTargets(
        mode="constraints",
        z_surface_min=20.0,
        sigma_z_surface_min=10.0,
        delay_proxy_min=200.0,
        sigma_delay_proxy_min=100.0,
    )
    ok = {"z_surface": 30.0, "delay_proxy": 250.0, "n_peak": 31.0}
    bad = {"z_surface": 10.0, "delay_proxy": 150.0, "n_peak": 11.0}
    assert loglike_compact(ok, targets) == 0.0
    assert loglike_compact(bad, targets) < 0.0


def test_background_link_returns_positive_params():
    d = derive_compact_star_params_from_background(gamma=1.0, V0=0.1, Q=0.0)
    assert d["theta_c"] > 0.0
    assert d["m2"] > 0.0
    assert d["lam"] >= 0.0
    assert d["mean_j_refr"] >= 0.0
    assert d["mean_i_grad"] >= 0.0
    assert d["j_refr_q95"] >= d["j_refr_q50"]


def test_loglikes_components_sum_matches_total_with_forward():
    cc = load_cc_dimless(Path("fitting/data/cc_demo_dimless.csv"))
    sn = load_sn_dimless(Path("fitting/data/sn_demo_dimless.csv"))
    bao = load_bao_dimless(Path("fitting/data/bao_demo_dimless.csv"))
    om = Omegas(H0=1.0)
    p = Params(gamma=1.0, V0=0.1, Q=0.0)
    cstar = CompactStarTargets(
        z_surface=350.0,
        sigma_z_surface=200.0,
        delay_proxy=700.0,
        sigma_delay_proxy=400.0,
    )
    fwd = ForwardTargets(
        ring_m87_uas=42.0,
        sigma_ring_m87_uas=3.0,
        ring_sgra_uas=51.8,
        sigma_ring_sgra_uas=2.3,
    )
    total = logposterior_combined(
        p, om=om, cc=cc, sn=sn, bao=bao, compact_star=cstar, forward=fwd,
        compact_star_link_background=True
    )
    comps = loglikes_components(
        p, om=om, cc=cc, sn=sn, bao=bao, compact_star=cstar, forward=fwd,
        compact_star_link_background=True
    )
    s = sum(comps.values())
    assert abs(total - s) < 1e-9
