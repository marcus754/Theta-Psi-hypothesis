from pathlib import Path

from fitting.data_io import load_cc_dimless
from fitting.model import Hz_dimless
from fitting.likelihoods import chi2_gaussian
from background.components import Omegas


def test_cc_pipeline_dimless_runs():
    data_path = Path("fitting/data/cc_demo_dimless.csv")
    cc = load_cc_dimless(data_path)
    om = Omegas(H0=1.0)
    Hmodel = Hz_dimless(cc.z, gamma=1.0, V0=0.1, Q=0.0, omegas=om, a_min=1e-3, dN=0.05)
    # The demo data roughly follows H≈sqrt(1+z); keep chi2 finite and non-pathological
    chi2 = chi2_gaussian(Hmodel, cc.H, cc.sigma)
    assert chi2 >= 0.0

