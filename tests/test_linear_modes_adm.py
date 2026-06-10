import math
import random

from src.linear_modes_adm import stability_adm, cs2_eigs_adm


def test_adm_modes_stability_no_mixing():
    random.seed(7)
    gamma = 1.0
    V0 = 0.2
    Q = 0.0
    a = 1.0
    k = 0.1
    for _ in range(20):
        x = random.uniform(-1.0, 1.0)
        theta = (math.sqrt(3.0)/gamma) * math.tanh(x)
        R = random.uniform(0.05, 0.5)
        res = stability_adm(theta, R, a, k, gamma, Q, V0, ctheta2=0.8, cR2=0.9, cphi2=0.95)
        assert res["no_ghost"]
        assert res["no_gradient_instability"]
        for s in res["cs2_eigs"]:
            assert s > 0.0 and s <= 1.0

