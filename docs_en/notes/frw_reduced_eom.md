FRW REDUCED EOM (V6)

Short working version of the background equations in the form actually used in background/frw_background.py.

VARIABLES

- N = ln a.
- x — regularization of the time field: theta = (sqrt(3)/gamma) tanh(x).
- y = dx/dN.
- R — reduced amplitude Psi = R e^{iφ}, u = dR/dN.

LIMITED KINETIC LOAD

A_raw(x, y, u) = 3/2 * (sech^4 x / tanh^2 x) * y^2 + 1/2 * u^2 + 3/2 * sech^4 x * y^2
A = 3 * tanh^2( sqrt(max(A_raw, 0)/3) ), 0 <= A < 3

FRIEDMANN

In N-form:
H^2 = [ rho_m + rho_r + rho_L + 1/2 V0 R^2 + 1/2 Q^2/(a^6 R^2) ] / [ 3 - A ] .

On a homogeneous FRW, the spatial source of the refractive regime is absent: I_grad=0, therefore J_refr=0 and the strong-field/refractive sector does not enter the background equations.

EQUATION FOR Hdot

Hdot = -1/2 * [ gamma^2 theta_dot^2 + R_dot^2 + 3 (theta_dot/theta)^2 + Q^2/(a^6 R^2) + rho_m + 4/3 rho_r ] .

STRONG-FIELD SECTOR ON FRW

On a homogeneous background, the spatial source of the regime is absent:
I_grad = 0 (on a homogeneous background)

Result:
J_refr = Phi_eff(theta_eff) * I_grad = 0 .

There is no longer a separate sigma regime field in the canon. Old sigma* parameters in the code are allowed only as legacy-compatibility and should not be used for deriving physical formulas.
