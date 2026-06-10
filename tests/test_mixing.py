
import math, random
from src.linear_modes import stability, cs2_eigs

def test_mixing_keeps_stability_for_small_eps():
    random.seed(42)
    gamma = 1.0
    for _ in range(50):
        x = random.uniform(-1.0, 1.0)
        theta = (math.sqrt(3.0)/gamma)*math.tanh(x)
        ctheta2 = random.uniform(0.3, 1.5)
        cpsi2 = random.uniform(0.7, 1.5)
        epsK = random.uniform(-0.1, 0.1)
        epsG = random.uniform(-0.1, 0.1)
        cond = stability(theta, gamma, ctheta2, cpsi2, epsK, epsG)
        assert cond["no_ghost"]
        assert cond["no_gradient_instability"]
        s1, s2 = cond["cs2_eigs"]
        assert s1 > 0 and s2 > 0
