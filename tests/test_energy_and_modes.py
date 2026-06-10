
import math
import random

from background.frw_symbolic import energy_density
from src.linear_modes import (
    kinetic_matrix,
    gradient_matrix,
    cs2_eigenvalues,
    stability_conditions,
)

def test_stability_and_energy_positive():
    random.seed(1)
    gamma = 1.0
    V0 = 0.3
    for _ in range(50):
        x = random.uniform(-1.2,1.2)
        theta = (math.sqrt(3.0)/gamma)*math.tanh(x)
        theta = max(min(theta, 0.99*math.sqrt(3.0)/abs(gamma)), -0.99*math.sqrt(3.0)/abs(gamma))
        theta_dot = random.uniform(-0.2, 0.2)
        psi = random.uniform(-0.5, 0.5)
        psi_dot = random.uniform(-0.5, 0.5)

        cond = stability_conditions(theta, gamma, ctheta2=0.7, cpsi2=1.0)
        assert cond["no_ghost"] and cond["no_gradient_instability"]

        rho = energy_density(theta, theta_dot, psi, psi_dot, gamma, V0)
        # We allow small negative values near theta~0 due to truncations? Use soft bound:
        assert rho > -1e-12
