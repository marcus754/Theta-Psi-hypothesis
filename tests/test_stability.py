
import math
import random
from background.frw_integrator import Params, pack_state, integrate, theta_from_x, theta_dot_from_x
from src.linear_modes import kinetic_matrix, is_kinetically_stable

def test_kinetic_positive_on_random_grid():
    random.seed(0)
    gamma = 1.0
    for _ in range(50):
        # Draw theta in safe band via x variable in [-1,1]
        x = random.uniform(-1.0, 1.0)
        theta = (math.sqrt(3.0)/gamma)*math.tanh(x)
        theta_dot = random.uniform(-0.2, 0.2)
        R = random.uniform(-0.5, 0.5)
        R_dot = random.uniform(-0.5, 0.5)
        H = random.uniform(0.0, 0.2)
        assert is_kinetically_stable(theta, theta_dot, R, R_dot, H, gamma)

def test_integrator_respects_theta_bound():
    p = Params(gamma=1.0, V0=0.1, H_of_t=lambda t: 0.05)
    y0 = pack_state(x=2.0, xdot=0.0, R=0.1, Rdot=0.0)
    traj = integrate(0.0, 2.0, y0, p, h0=1e-2)
    for t, state in traj:
        x, xdot, R, Rdot = state
        theta = theta_from_x(x, p)
        assert abs(theta) < math.sqrt(3.0)/abs(p.gamma) - 1e-12
