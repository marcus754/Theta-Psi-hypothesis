from pathlib import Path

from fitting.data_io import load_sn_dimless, load_bao_dimless
from fitting.model import mu_dimless, DV_dimless
from background.components import Omegas
from fitting.likelihoods import chi2_gaussian


def test_sn_pipeline_dimless_runs():
    data_path = Path("fitting/data/sn_demo_dimless.csv")
    sn = load_sn_dimless(data_path)
    om = Omegas(H0=1.0)
    mu_model = mu_dimless(sn.z, gamma=1.0, V0=0.1, Q=0.0, omegas=om, a_min=1e-3, dN=0.05, mu0=0.0)
    chi2 = chi2_gaussian(mu_model, sn.mu, sn.sigma)
    assert chi2 >= 0.0


def test_bao_pipeline_dimless_runs():
    data_path = Path("fitting/data/bao_demo_dimless.csv")
    bao = load_bao_dimless(data_path)
    om = Omegas(H0=1.0)
    DV_model = DV_dimless(bao.z, gamma=1.0, V0=0.1, Q=0.0, omegas=om, a_min=1e-3, dN=0.05)
    chi2 = chi2_gaussian(DV_model, bao.DV, bao.sigma)
    assert chi2 >= 0.0

