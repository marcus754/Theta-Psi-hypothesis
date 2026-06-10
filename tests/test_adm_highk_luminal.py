import math
import random

from src.linear_modes_adm import cs2_eigs_adm


def test_adm_highk_luminal_default_params():
    random.seed(42)
    gamma = 1.0
    V0 = 0.1
    Q = 0.0
    a = 1.0
    k = 0.3
    for _ in range(20):
        x = random.uniform(-1.0, 1.0)
        theta = (math.sqrt(3.0)/gamma) * math.tanh(x)
        R = random.uniform(0.05, 0.5)
        c = cs2_eigs_adm(theta, R, a, k, gamma, Q, V0)
        # All three modes should be luminal at high-k in the sigma-model
        for ci in c:
            assert abs(ci - 1.0) < 1e-12

