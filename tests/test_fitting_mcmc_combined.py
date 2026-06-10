from pathlib import Path

from fitting.run_mcmc import Params, logposterior_combined
from fitting.data_io import load_cc_dimless, load_sn_dimless, load_bao_dimless
from background.components import Omegas


def test_logposterior_combined_runs_on_demos():
    cc = load_cc_dimless(Path("fitting/data/cc_demo_dimless.csv"))
    sn = load_sn_dimless(Path("fitting/data/sn_demo_dimless.csv"))
    bao = load_bao_dimless(Path("fitting/data/bao_demo_dimless.csv"))
    om = Omegas(H0=1.0)
    p = Params(gamma=1.0, V0=0.1, Q=0.0)
    lp = logposterior_combined(p, om=om, cc=cc, sn=sn, bao=bao)
    assert lp <= 0.0  # Gaussian loglike is non-positive
    assert lp > -1e12

