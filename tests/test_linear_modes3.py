import math
import random

from src.linear_modes3 import stability_3, cs2_eigs_3


def test_linear_modes3_basic_stability_no_mixing():
    random.seed(3)
    gamma = 1.0
    for _ in range(20):
        x = random.uniform(-1.0, 1.0)
        theta = (math.sqrt(3.0)/gamma) * math.tanh(x)
        R = random.uniform(0.05, 0.5)
        res = stability_3(theta, gamma, R, ctheta2=0.8, cR2=0.9, cphi2=0.95)
        assert res["no_ghost"]
        assert res["no_gradient_instability"]
        s = res["cs2_eigs"]
        assert all(si > 0 for si in s)

