11. Derivation Line of Theta–Psi: Tensor Time as Canon

Goal: to fix without ambiguity that in the canonical version, the fundamental object of time is the Theta_munu tensor, and the scalar theta is only its reduced projection for current FRW calculations.

1) Canon Postulates

1. The Theta_munu tensor (symmetric tensor of the 2nd rank) is accepted as the fundamental object of time.
2. The energy sector is modeled by a complex Psi field.
3. It is assumed that the observed dynamics of matter are determined by an effective vacuum response, which is determined by invariants of Theta_munu.
4. In the weak field limit, the working version of the theory must reproduce verified observational tests.
5. Since time in the hypothesis is understood as a physical object, permissible dynamics must not translate the temporal regime into zero or a negative branch.

2) Basic Covariant Action (structure)

S = ∫ d^4x √(-g) [ (M_pl^2/2)R + L_m + L_Theta(Theta_munu, g) + L_Psi(Psi, g) + L_int(Theta,Psi,g) ].

2.1) Unified reduced action for FRW and strong field

We fix a single reduced action that reproduces both working limits:
S_red = ∫ d^4x √(-g) [
  (M_pl²/2)R
  + L_Theta
  + L_R
  - U_red(theta_eff,R)
  + M_*^4 F(J_refr)
].

In this canon:
L_Theta = 3/2 · ∇_mu ln theta_eff ∇^mu ln theta_eff + 1/2 · gamma² ∇_mu theta_eff ∇^mu theta_eff
L_R = 1/2 · ∇_mu R ∇^mu R

With a minimal potential:
U_red(theta_eff,R) = 1/2 m_theta² theta_eff² + λ_theta/4 theta_eff⁴ + 1/2 V0 R² + λ_ThetaPsi/2 theta_eff² R².

Meaning of this record:
1. On the FRW background ∇_i theta_eff = ∇_i Psi = 0, so J_refr = 0 and a homogeneous Lagrangian remains.
2. In a stationary inhomogeneous regime, the gradients of time/energy turn on F(J_refr), giving the answer in a strong field from the same action.
3. This is a single reduced covariant action, not two independent Lagrangians.

2.2) Current register of covariant action

In the repository, this transition is now supported not only by text, but also by a separate symbolic register module:
- theory/covariant_action.py

Its role:
1. to fix theta_eff² = Theta_munu u^mu u^nu;
2. to store the canonical Phi_eff^can = ln(theta_eff/theta_0);
3. to hold a single symbolic source for S_red, stationary J_refr, and reduced energy-momentum blocks;
4. to allow verification that FRW and stationary strong field are indeed read as two limits of a single reduced action, rather than two independent substitutions.

This is the canonical action register: it exists, is tested, and sets the working canon S_red.

2.3) Equations of motion and reduced T_munu

From S_red for a homogeneous FRW slice, we fix a single minisuperspace register:
E_theta = (3/theta² + gamma²) thetä - 3 thetȧ² / theta³ = 0,
E_Psi = Psï + V0 Psi = 0.

Corresponding energy density and pressure:
ρ_red = 3/2 · (thetȧ/theta)² + 1/2 · Ṙ² + 1/2 · V0 R² + 1/2 · gamma² thetȧ²,
p_red = 3/2 · (thetȧ/theta)² + 1/2 · Ṙ² - 1/2 · V0 R² + 1/2 · gamma² thetȧ².

Consequently:
T^(red)_00 = ρ_red,
T^(red)_ii = a² p_red (without summation over i).

In the reduced canon, this is the minimal closed set: action, variation by fields, and summary of T_munu in one source.

2.4) Full covariant closure

The full covariant object of the canon is fixed as

S_cov = S_red + S_m[g_hat, matter],
g_hat = g + (n²(theta_eff,Theta,Psi) - 1) u ⊗ u.

Here g_hat is not a new ontological geometry, but a compact record of the refractive and clock response of the vacuum. S_m[g_hat, matter] means that matter and radiation receive the same propagation law through this response. This is the full covariant closure for the current canon: a reduced core plus a single matter coupling rule.

In this record, g and g_hat have the status of an accounting language, not primary objects of the hypothesis. Theta_munu and Psi remain primary. The metric form is only needed to write a local law of process propagation without different rules for light, clocks, and massive matter.

Variation by metric:
G_munu = 8πG · (T^(m)_munu + T^(Theta)_munu + T^(Psi)_munu + T^(int)_munu).

Critically important: observable tensor perturbations should be collective perturbations of the effective propagation law, rather than a consequence of a single scalar reduction theta_eff.

3) Tensor Lagrangian L_Theta: minimal form

3.1) Requirements for L_Theta

1. Diffeomorphism invariance: scalar relative to coordinate transformations.
2. Minimal number of derivatives: no higher than 2nd order for hyperbolicity.
3. Positive energy: ghost-free.
4. Correct inert limit: at Theta_munu → const, the effective vacuum response ceases to introduce additional inhomogeneity.

3.2) Minimal Lagrangian on the space of tensors

The most general ghost-free form for a symmetric tensor of the 2nd rank:
L_Theta = 1/2 · G^munuρσ(Theta) · ∇_mu Theta_nuρ · ∇_σ Theta_αbeta · g^αbeta - U(Theta_munuTheta^munu),

where G^munuρσ is a metric on the Theta configuration space.

For an isotropic FRW background with a highlighted time u^mu = (1,0,0,0) in the comoving frame:
Theta_munu = diag(theta²(t), 0, 0, 0) (in synchronous calibration),

then the only non-trivial invariant:
I_1 = Theta_munu u^mu u^nu = theta²(t).

3.3) Reduction to the scalar theta_eff

Define:
theta_eff ≡ √(Theta_munu u^mu u^nu), u^mu u_mu = -1.

In the FRW metric ds² = -dt² + a²(t) dx²:
u^mu = (1, 0, 0, 0), u_mu = (-1, 0, 0, 0),
Theta_00 = theta²(t), Theta_ij = 0,
theta_eff = √(Theta_00 u^0 u^0) = √(theta² · 1 · 1) = theta(t).

3.4) Kinetic term derivation

From the tensor Lagrangian:
L_Theta ⊃ 1/2 · G^0000(Theta) · ∂_0 Theta_00 · ∂_0 Theta_00.

For Theta_00 = theta²:
∂_0 Theta_00 = 2theta thetȧ,
L_Theta ⊃ 1/2 · G^0000 · (2theta thetȧ)² = 2 G^0000 · theta² thetȧ².

Requirement for canonical normalization (see 08_foundations.md):
L_eff ⊃ 3/2 · (thetȧ/theta)² = 3/2 · theta^-2 · thetȧ².

Comparing:
2 G^0000 · theta² = 3/2 · theta^-2 ⇒ G^0000 = 3/4 · theta^-4.

Conclusion: configuration space metric on Theta:
G^munuρσ(Theta) = 3/4 · (Theta_αbeta u^α u^beta)^-2 · P^munuρσ,

where P^munuρσ is the projector onto the temporal component.

3.5) Mixing term gamma² thetȧ²

This term arises from the interaction of Theta and Psi:
L_int ⊃ gamma²/2 · (∇_mu Theta_nuρ) (∇^mu Psi) (∇^nu Psi) (∇^ρ Psi) / (Psī Psi).

In the FRW limit with Psi = Psi(t):
L_int ⊃ gamma²/2 · thetȧ² · Psi̇² / |Psi|².

For a homogeneous background Psi̇²/|Psi|² ≈ 1 (field normalization), what remains is:
L_int ⊃ gamma²/2 · thetȧ².

We assume that on a homogeneous background, the Psi field is normalized such that Psi̇²/|Psi|² = 1. Under this condition, gamma²/2 · thetȧ² follows from L_int.

Physical meaning: gamma is a dimensionless coupling between the time gradient and the energy gradient.

4) Final effective Lagrangian

After the reduction Theta_munu → theta_eff:
L_eff = 3/2 · (thetȧ/theta)² + 1/2 · Ṙ² - 1/2 · V0 · R² + 1/2 · gamma² · thetȧ².

4.1) Limits check

| Limit | L_eff | Physics |
|--------|-------|--------|
| gamma → 0 | 3/2 (thetȧ/theta)² + 1/2 Ṙ² - 1/2 V0 R² | Decoupled Theta and Psi |
| theta → const | 1/2 Ṙ² - 1/2 V0 R² | inert limit of Psi-sector |
| Psi → 0 | 3/2 (thetȧ/theta)² + 1/2 gamma² thetȧ² | Pure Theta-field (stiff matter) |
| thetȧ → 0 | 1/2 Ṙ² - 1/2 V0 R² | static temporal regime |

4.2) Energy and pressure

ρ = (∂L/∂thetȧ)thetȧ + (∂L/∂Psi̇)Psi̇ - L
  = 3/2 (thetȧ/theta)² + 1/2 Ṙ² + 1/2 V0 R² + 1/2 gamma² thetȧ²,

p = L = 3/2 (thetȧ/theta)² + 1/2 Ṙ² - 1/2 V0 R² + 1/2 gamma² thetȧ².

Equation of state:
w = p/ρ = [3/2 (thetȧ/theta)² + 1/2 Ṙ² - 1/2 V0 R² + 1/2 gamma² thetȧ²] / [3/2 (thetȧ/theta)² + 1/2 Ṙ² + 1/2 V0 R² + 1/2 gamma² thetȧ²].

Limits:
- Kinetic dominance: w → +1 (stiff matter)
- Potential dominance: w → -1 (dark energy)
- thetȧ = 0, Psi̇ = 0: w = -1 (cosmological constant)

5) Optical metric from Theta_munu tensor

5.1) Effective metric for matter

Matter moves along the geodesics of the effective metric:
g_hat_munu = g_munu + (n²(theta_eff) - 1) · u_mu u_nu,

where u_mu = ∂_mu theta_eff / √(∂theta·∂theta) — highlighted temporal direction.

In the comoving frame u^mu = (1,0,0,0):
g_hat_00 = n²(theta_eff), g_hat_ij = -delta_ij.

5.1.1) Minimal working class of L_refr(Theta,Psi)

The link

Theta, Psi -> J_refr(Theta,Psi) -> n(theta_eff,Theta,Psi) -> g_hat_munu

fixes the minimal working class of the refractive sector. Its task is to define how the universal refractive response of the vacuum arises from the temporal object and the energy substrate.

Requirements for L_refr:
1. It depends on Theta and Psi only through invariants compatible with the canon.
2. It is universal for matter and radiation.
3. On the FRW background, its contribution is turned off automatically at I_grad = 0.
4. In the weak field limit, it gives a linear correspondence ln n κ Phi_eff.
5. In the strong field regime, it saturates and does not drive n into pathology.

Minimal record:
L_refr = M_*^4 · F(J_refr),
J_refr = Phi_eff(theta_eff) · I_grad(Theta,Psi),
Phi_eff(theta_eff) = ln(theta_eff / theta_0).

Meaning:
- Phi_eff — canonical weak field indicator of the temporal regime;
- I_grad measures the strength of the spatial inhomogeneity of the Theta and Psi interaction;
- F defines a universal vacuum response, the same for light, clocks, and massive bodies through a common effective metric g_hat_munu.

Then the effective index is not introduced separately "from above," but is defined as a universal vacuum response to J_refr:
n - 1 = 𝒩(J_refr),

where 𝒩(0)=0, 𝒩'(0)=κ, 𝒩(J)>=0, 𝒩(J) is monotonic and has finite saturation.

Equivalent working form:
n(theta_eff,Theta,Psi)=1+𝒩(J_refr).

This form fixes three properties at once:
1. at I_grad=0 we get J_refr=0 and n=1, that is, FRW does not turn on the refractive sector;
2. at small J_refr we have a weak field law n ≈ 1 + κ J_refr;
3. at large J_refr, the response saturates, and numerical iteration does not go into an unstable region.

5.1.1a) Status of Phi_eff: canon vs numerical regularization

To remove ambiguity:
1. Canonical definition:
Phi_eff^can(theta_eff) = ln(theta_eff / theta_0), theta_eff > 0.

2. Numerical regularization in solvers:
Phi_eff^num(theta_eff) = smooth positive approximation(theta_eff/theta_0),

which coincides with the canonical form in weak field alignment and is used only where a smooth positive function without logarithmic rigidity at machine theta_eff 0 is needed.

Consequence: the hypothesis canon does not depend on a specific smooth approximation given three conditions:
- Phi_eff^num >= 0,
- Phi_eff^num -> 0 as theta_eff -> theta_0,
- in the working weak field Phi_eff^num = Phi_eff^can + O((theta_eff-theta_0)^2).

5.1.2) How L_refr enters the matter action

L_refr restructures the effective local clock sector of the vacuum. Working record:
S_m -> S_m[g_hat_munu, matter],
g_hat_munu = g_munu + (n²(theta_eff,Theta,Psi)-1) u_mu u_nu.

5.1.2a) Energy-momentum tensor of the reduced sector

For reduced action, we define:
T_munu^(red) = -2 / √(-g) · delta(√(-g)L_red) / deltag^munu
               = T_munu^(theta) + T_munu^(Psi) + T_munu^(refr) - g_munuU_red.

Minimal pieces:
T_munu^(theta) =
  3 · ∇_mu ln theta_eff ∇_nu ln theta_eff
  + gamma² ∇_mu theta_eff ∇_nu theta_eff
  - g_munu[ 3/2 (∇lntheta_eff)^2 + 1/2 gamma² (∇theta_eff)^2 ],

T_munu^(Psi) =
  ∇_mu Psi ∇_nu Psi
  - g_munu[ 1/2 (∇Psi)^2 ],

T_munu^(refr) =
  -g_munu M_*^4 F(J_refr)
  + 2 M_*^4 F'(J_refr) · deltaJ_refr/deltag^munu.

On the FRW background:
I_grad = 0 => J_refr = 0 => T_munu^(refr) = 0,

and the energy-momentum tensor reduces to a homogeneous Theta–Psi sector. In stationary regime, T_munu^(refr) plays the role of an additional stress-sector, which creates a finite optical response without introducing a separate switch field.

5.1.2b) Linear stability of the reduced sector

Minimal criteria before full ADM-closure:
1. FRW no-ghost:
K_FRW = diag( 3/theta_eff² + gamma², 1 ) ≻ 0.

2. FRW no-gradient-instability:
G_FRW = diag( c_theta² K_theta, c_Psi² ) ≻ 0.

3. Stationary ellipticity for radial two-field solver:
M_11 > 0, M_22 > 0, det(M) > 0,

where M is a local quasilinear matrix from item 5.1.3.

4. No extra modes:
no_extra_modes = no_ghost ∧ no_gradient_instability ∧ elliptic.

This gives a simple and verifiable criterion: the working parameter range must be simultaneously no-ghost on FRW, elliptic in stationary strong field sector and without signs of a hidden additional mode.

Note: the code implements a numerically stable estimate of eigenvalues via Cholesky reduction (see src/linear_modes_adm.py). Symbolic framework for future full elimination of Lapse/Shift is in theory/adm_symbolic.py.

Lapse/Shift elimination (c = (d_N, B)):
L^(2) ⊃ − dqᵀ B c − 1/2 · cᵀ C c ⇒ c* = − C⁻¹ Bᵀ dq,
A_eff = A − B C⁻¹ Bᵀ, A = G (k²/a²) + M.

From here, G_eff and M_eff are recovered (terms at k²/a² and constants). In the numerical API, passing elim=... includes this step (see src/linear_modes_adm.py).

5.1.3) Minimal stationary sector and Euler–Lagrange equations

To ensure the strong field solver is Lagrange-closed, we fix a minimal stationary radial functional for the reduced two-entity sector:
L_stat = 1/2 · I_grad - U(theta,Psi) + α · F(J_refr),
I_grad = theta'^2 + w_Psi Psi'^2,
J_refr = Phi_eff(theta) · I_grad,
Phi_eff(theta) = ln(theta/theta_0) or its smooth positive approximation.

Here:
- α — strength of the refractive feedback;
- w_Psi > 0 — relative weight of the Psi gradient contribution;
- prime denotes derivative with respect to the stationary radial coordinate r.

Minimal working stationary potential, consistent with the current computational solver:
U(theta,Psi) =
  1/2 m² theta² + λ/4 theta⁴
  - 1/2 m² Psi² + c_ThetaPsi λ/2 · theta² Psi².

Then:
U_,theta = m² theta + λ theta³ + c_ThetaPsi λ theta Psi²,
U_,Psi = -m² Psi + c_ThetaPsi λ theta² Psi.

For a spherically symmetric stationary profile, the Euler–Lagrange equations are:
1/r² · d/dr [ r² A theta' ] + U_,theta - α F'(J_refr) Phi'_(eff) I_grad = 0,
1/r² · d/dr [ r² A w_Psi Psi' ] + U_,Psi = 0,

where

A(theta,theta',Psi') = 1 + 2 α Phi_eff(theta) F'(J_refr).

Since A depends on gradients through J_refr, the system is quasilinear. It is convenient to rewrite it as a local matrix system for second derivatives:
M_11 theta'' + M_12 Psi'' = R_theta,
M_21 theta'' + M_22 Psi'' = R_Psi,

with coefficients

A_,theta = 2 α Phi'_(eff) [F' + Phi_eff I_grad F''],
A_,theta' = 4 α Phi_eff² F'' theta',
A_,Psi' = 4 α Phi_eff² F'' w_Psi Psi',

M_11 = A + A_,theta' theta',
M_12 = A_,Psi' theta',
M_21 = A_,theta' w_Psi Psi',
M_22 = w_Psi (A + A_,Psi' Psi').

Right side:
R_theta = -A_,theta theta'^2 - 2Atheta'/r - U_,theta + αF' Phi'_(eff) I_grad,
R_Psi = -A_,theta w_Psi theta'Psi' - 2Aw_PsiPsi'/r - U_,Psi.

This is the minimal form in which the current two-field solver can already be considered Lagrange-motivated, rather than just an ad hoc closure.

5.1.4) Weak field and FRW limits of the stationary sector

If I_grad -> 0, then:
J_refr -> 0,
F'(J_refr) -> F'(0),
A -> 1,

and the radial system reduces to two ordinary weakly coupled scalar profile equations. On the FRW background, there are no spatial gradients, so the stationary refractive block is automatically turned off.

That is, L_refr defines a universal scalar response n, and already through g_hat this response is read by all matter sectors equally.

This is the substantial reason why:
- light is lensed,
- clocks slow down,
- massive bodies experience acceleration,

but at the same time, three different mechanisms are not introduced for three different phenomena.

5.1.3) FRW-limit and strong field activation

On a homogeneous FRW background:
I_grad = 0 => J_refr = 0.

Therefore, the FRW background itself does not turn on the refractive sector.

In a strong field or quasistationary regime:
I_grad > 0 => J_refr = Phi_eff I_grad > 0

and the strong-field regime turns on directly, without an intermediate regime field.

This means that the strong field regime is turned on like this:
1. Theta, Psi gradients create I_grad;
2. I_grad together with Phi_eff(theta_eff) creates J_refr;
3. J_refr turns on the finite refractive response n(theta_eff,Theta,Psi).

5.1.4) How saturation of n is obtained from this

For small J_refr:
F(J) = J + O(J²),

therefore

n - 1 = 𝒩(J_refr) ≈ κ J_refr = κ Phi_eff I_grad.

That is, in the weak strong field regime, the response grows analytically and controllably.

For large J_refr, the function 𝒩(J) must saturate:
𝒩(J) -> 𝒩_∞ < ∞ as J -> +∞.

Consequently, even with the growth of I_grad, the temporal response does not go into a pathological branch. This is the mathematical mechanism of how a finite strong-field regime is born from the direct interaction of Theta and Psi without a manual switch.

5.1.5) Canonical class of F(J) functions

To ensure the regime sector is not arbitrary, we fix the minimal allowed properties of the functions even before fitting and strong-field tests.

For F(J), we canonically require:
F(0) = 0,
F'(J) > 0 for J >= 0,
F''(J) <= 0 for J >= 0,
F(J) -> F_∞ < ∞ as J -> +∞.

Meaning:
- F(0)=0 guarantees the absence of artificial refractive shift on the background;
- F'(J)>0 guarantees that strengthening the regime does not reduce the vacuum response;
- F''(J)<=0 makes the response saturable rather than explosive;
- finite F_∞ excludes the escape into an infinite optical branch.

Canonical normalized family-v1:
F_v1(J) = 1 - exp(-J).

Any explicit scale J_* can be removed by normalizing J_refr to a dimensionless variable; after that, it is not a physical parameter.

To make the exponent itself a derivation, a separate principle is needed for the response remainder:
S(J) = 1 - F(J),
S(J1 + J2) = S(J1) S(J2).

That is, independent increments of J act without memory and multiplicatively on the remaining unused part of the response. With continuity, such a law gives

S(J) = exp(-J/J_*),
F(J) = 1 - exp(-J/J_*).

It is this step that makes the exponent a consequence of the semigroup principle for the response remainder, rather than a separate choice of form.

5.1.6) Why exactly this class is considered canonical

For family-v1 we get:
F_v1'(J) = exp(-J) > 0,
F_v1''(J) = -exp(-J) < 0,
F_v1(∞) = 1.

That is, all the requirements above are fulfilled automatically.

In addition, for small J:
F_v1(J) = J + O(J²),

and the weak field regime is obtained analytical and controllable. For large J, the response saturates, which makes family-v1 a minimal smooth candidate without kinks and without additional thresholds.

5.1.7) Connection of family-v1 with the working form for n

If we introduce a numerically regularized family

𝒩(J) = saturating_response(J),

then the family for n takes the form

n(theta_eff,Theta,Psi) = 1 + 𝒩(J_refr),
J_refr = Phi_eff(theta_eff) I_grad(Theta,Psi).

This makes the structure uniform:
- F(J) controls how the refractive sector enters the effective action;
- 𝒩(J) controls how this sector is read as an observable index in the solver implementation;
- both functions belong to the same class: monotonic, concave, saturating response without threshold and without branch change.

5.1.8) What is allowed as an alternative and what is not

Allowed as a blind-test family if fixed in advance:
F_alt,1(J) = J / (1 + J),
F_alt,2(J) = tanh(J),

with the same requirements of monotonicity and saturation.

Not allowed in the canon:
1. functions with F'(J) < 0 on part of the domain;
2. functions with unlimited growth of F(J);
3. piecewise thresholds of the form "if J>J_c, then turn on the regime";
4. redefinition of the family after viewing the data.

5.1.9) Canon Lemmas

Below are fixed minimal statements that must be true for any allowed F(J_refr). F_v1 in normalized record 1-exp(-J) is a canonical representative-form; explicit J_* is only normalization.

Lemma A. Disconnection of the refractive sector on FRW

Prerequisites:
I_grad = 0,
F(0)=0.

Then

J_refr = 0, n = 1.

Meaning: FRW does not turn on the refractive sector.

Lemma B. Monotonicity and finiteness of the refractive response

Prerequisites:
J >= 0,
𝒩(0)=0,
𝒩'(J) > 0,
𝒩''(J) <= 0,
lim_J->∞ 𝒩(J) < ∞.

Then

1 <= n(J) = 1 + 𝒩(J),
dn/dJ > 0.

Consequently, strengthening the regime monotonically increases the index but does not allow it to escape to the infinite optical branch on the physical branch of the solution; numerical regularization serves only for solver stability.

Lemma C. Weak field analyticity

If F(J) and 𝒩(J) are analytical at zero and

F(J) = J + O(J²),
𝒩(J) = κ J + O(J²),

then for small I_grad and weak Phi_eff:
n - 1 = κ Phi_eff I_grad + O((Phi_eff I_grad)²).

Consequently, the weak field response does not contain a hidden threshold and turns on smoothly.

Lemma D. Inaccessibility of the true horizon in the reduced circuit

Prerequisites:
ρ(t) <= ρ_max < ∞,
L_eff ⊃ 3/2 (thetȧ/theta)²,
theta(t0) > 0.

Then

|d ln theta / dt| <= sqrt(2ρ_max/3),

and after integration:
theta(t) >= theta(t0) exp(-K |t-t0|) > 0

for any finite t.

Consequently, the temporal regime does not reach physical vanishing in a finite course, and in the final optical branch this means finite redshift and delay.

Lemma E. Equivalence as a consequence of universal g_hat_munu

If matter enters only through

S_m[g_hat_munu, matter],
g_hat_munu = g_munu + (n²-1) u_mu u_nu,

and the function n is the same for all matter fields, then in the non-relativistic limit we get

L ≈ (m/2n) v² - m/n,
a = -∇ ln n,

and the mass of the test body cancels out of the equation of motion.

Consequently, universality of free fall here is not a separate postulate, but a consequence of the universal vacuum response.

5.1.10) What should be checked numerically

Minimal verification criteria follow from the lemmas above:
1. FRW runs with I_grad=0 must give J_refr=0 and not turn on the refractive sector.
2. The chosen family 𝒩(J) must remain monotonic and finite on the entire working parameter range.
3. In strong field scans, a finite J_refr exists without branch change and without blow-up of n.
4. For finite energy, numerically accessible theta=0 must not appear.

5.2) Connection of n(theta_eff,Theta,Psi) with the tensor

The refractive index is a function of the Theta invariant. Important here: the regularized form is not a primary ontological postulate, but only a working record of the finite physical branch and numerical solver stability:
n(theta_eff,Theta,Psi) = 1 + 𝒩(J_refr).

where Phi_eff(theta_eff) — weak field proxy:
Phi_eff(theta_eff) = ln(theta_eff / theta_0) (logarithmic link).

As theta_eff → theta_0 and a fixed weak regime:
n ≈ 1 + κ · J_refr ≈ 1 + κ · Phi_eff I_grad.

The canonical choice κ = 2 sets the weak field normalization on the branch where J_refr is already reduced to an effective response coefficient:
n ≈ 1 + 2 Phi_N, where Phi_N is the Newtonian potential.

5.2.1) Where gravity comes from in this language

If Theta_munu is inhomogeneous in space, then theta_eff is also inhomogeneous, and along with it, the vacuum response n(theta_eff) is also inhomogeneous. This means that the vacuum sets a non-uniform local passage of processes in different points.

For a test body:
S_probe = -m ∫ dtau, dtau = dt / n(theta_eff).

In the non-relativistic limit:
L ≈ (m/2n) v² - m/n,

and the equation of motion gives

a = -∇ ln n(theta_eff).

That is, gravity in this frame arises not as a separate force entity, but as a response of bodies to a spatial gradient of causal-temporal vacuum conductivity.

5.2.2) What "vacuum refraction" means

This expression should be understood broadly. It is not a single photon and not the "speed of light itself" that is refracted, but the regime of causality transfer in the vacuum determined by Theta.

- For photons, this gives ray bending and signal delay.
- For massive bodies, this gives trajectory restructuring through the same vacuum response.
- For local clocks, this gives a slowing or acceleration of the passage rate.

Therefore, "vacuum refraction" in Theta–Psi means a change in how the vacuum conducts causal processes, and light is just one of the observed probes of this regime.

5.3) Tensor waves from g_hat_munu

Variation of the action by g_hat_munu gives an induced Einstein term (Sakharov):
Γ_eff[g_hat] ⊃ ∫ d^4x √(-g_hat) · (M_ind²/2) · R(g_hat).

Linearization g_hat_munu = eta_munu + h_munu gives the tensor wave equation:
□ h_ij^TT = 0 ⇒ c_T = 1.

Conclusion: tensor waves are TT-modes of the effective metric g_hat_munu(Theta), not of Theta_munu itself.

Important: there are no fundamental gravitons in Theta–Psi. Tensor waves are an emergent phenomenon, oscillations of the induced metric g_hat_munu. They are not quantized within classical Theta–Psi (the quantum limit requires a separate derivation).

5.4) Consequence for the horizon

If the temporal regime is required to remain positive, then the reduced optical record cannot lead to a physically reachable vanishing of time in a finite course. In the final form, it looks like the finiteness of n(theta_eff) and, consequently, the finiteness of redshift/delay.

Semantic priority is not on the upper limit on n, but on the positive status of time as an object; finite n here is a working computational form of this requirement.

6) What is already fixed in v2 and what is not

Fixed:
1. Time canon in documentation: Theta_munu.
2. Effective FRW circuit on theta_eff as an approximation.
3. Minimal working class L_refr = M_*^4 F(J_refr) with J_refr = Phi_eff(theta_eff) I_grad(Theta,Psi) fixed as a canonical bridge between Theta, Psi and n(theta_eff,Theta,Psi).
4. Response class n(theta_eff,Theta,Psi) is set in advance of fit.
5. Tensor derivation L_Theta → L_eff (this document).
6. Diagnostic combined stability_register check.

Not fixed (requires adjustment):
1. Full high-k diagnostics for tensor modes in the extended sector.
2. Unique Psi signature outside the background fitting.
3. Microphysical origin of the F(J_refr) response function beyond minimal class.
4. Quantum limit (dq_k → HUP → unitarity).

6.1) Four closure steps before fit

This list sets the working sequence and does not allow the Theta–Psi to disperse into a set of independent fittings.

Step 1. Strict covariant derivation

What is already there:
- theory/covariant_action.py as a symbolic register;
- reduced FRW/strong field action form in documents;
- canonical bridge theta_eff, Phi_eff, J_refr.

What must be finalized:
1. one covariant action S[g, Theta, Psi];
2. explicit variation by g_munu, Theta_munu, Psi;
3. clear separation of kinetics, potential, coupling, and matter;
4. matter coupling rule through a single propagation law encoded by g_hat_munu;
5. explicit indication of which members are fixed by the canon and which are working definitions;
6. list of symmetries and constraints that fix the number of degrees of freedom.

Acceptance criteria:
- no hidden additional fields masquerading as parameters;
- FRW and stationary branches are indeed limits of one action;
- stress-energy is obtained by one variation, rather than manually assembled blocks.

Step 2. Stability and absence of ghost/extra modes

What is already there:
- FRW no-ghost/no-gradient diagnostics in theory/twofield_stability.py;
- stationary ellipticity diagnostics for radial strong field solver;
- finite n(theta_eff) as a restriction on the physical branch.

What must be finalized:
1. explicit quadratic form for perturbations around FRW;
2. explicit quadratic form for perturbations around stationary background;
3. checking signs of kinetic matrix and gradient matrix;
4. checking tachyon/gradient/ghost in the working region;
5. checking that constraints do not produce a hidden third mode;
6. high-k behavior and causal structure for tensor perturbations.

Acceptance criteria:
- kinetic matrix is positive definite in the working area;
- gradient matrix is positive definite;
- stationary operator remains elliptical on the used range;
- no runaway or blow-up of a hidden mode occurs during numerical evolution.

Step 3. Reduce to 1-2 distinctive observable predictions

What is already there:
- strong-field packages of observational checks;
- background CC+SN+BAO;
- event-level direct observables.

What must be finalized:
1. select one main strong-field discriminator;
2. select at most one additional discriminator;
3. leave all other quantities as auxiliary diagnostics;
4. declare in advance which observables are primary and which are secondary;
5. fix the mapping from theory output to data-facing measurement without manual search.

For the current version, a reasonable choice is:
- primary 1: saturation of compact/strong field curve z_surface and delay together with ok_no_horizon as a single theoretical strong field-test;
- primary 2: ring_uas and echo_delay_ms as downstream event-level comparator for direct data-facing check.

Important restriction:
- CC+SN+BAO should not be the only "proof" of the hypothesis;
- they are needed as a background that is preserved when the strong-field regime is turned on.

Acceptance criteria:
- the number of primary predictions does not grow from run to run;
- each primary prediction has a defined sign/direction of deviation;
- each primary prediction has a direct data-facing comparator.

Step 4. Only then fit and comparison with observational baselines

What is already there:
- head-to-head and validation pipeline;
- k-fold/falsifiability protocols;
- reports on AIC/BIC and χ² components.

What must be finalized:
1. fit only after closing steps 1-3;
2. no hidden tuning on test;
3. one train/test protocol for all models;
4. comparison not only by χ², but also by AIC/BIC and out-of-samp <= stability;
5. ban on reassembling observab <= list after viewing the result.

Acceptance criteria:
- fit does not change the canon, but only evaluates parameters of an already closed model;
- test set does not participate in the choice of bridge/observab <= family;
- if the improvement disappears on holdout, the model is not considered the winner.

6.2) Strict proof for the current canon

Below is fixed not a promotional text, but a working proof-skeleton for the already written reduced canon. It is strict in the sense that each conclusion rests on explicitly written prerequisites. It does not invent new microphysics; it closes what is already contained in S_red, L_refr, stability_register and theory/adm_symbolic.py.

Theorem 1. Covariant closure of reduced canon

Let the following conditions be fulfilled:
1. theta_eff^2 = Theta_munu u^mu u^nu and Phi_eff = ln(theta_eff/theta_0) are scalars;
2. S_red = ∫ √(-g)[ EH + L_Theta + L_Psi - U_red + M_*^4 F(J_refr)] d^4x;
3. J_refr = Phi_eff · I_grad, I_grad = l_theta^2|∇_⊥theta|^2 + l_Psi^2|∇_⊥Psi|^2;
4. matter couples only through the single propagation law encoded by g_hat_munu = g_munu + (n^2-1)u_mu u_nu;
5. F depends only on J_refr and not on gauge-dependent objects.

Then:
1. each member of S_red is diffeomorphism-invariant;
2. variation by g_munu gives one consistent T^(red)_munu;
3. FRW background and stationary strong field sector are two limit realizations of the same action, rather than different theories;
4. the reduced canon is closed at the level of field invariants theta_eff, Phi_eff, J_refr, g_hat_munu.

Proof:
1. EH, √(-g), ∇_mu and g_munu are standard covariant objects.
2. theta_eff and Phi_eff are built from a full convolution by u^mu, therefore they are scalars; I_grad and J_refr are composed of scalar norms of spatial gradients and therefore are also scalars.
3. U_red depends only on theta_eff and Psi, and L_refr depends only on J_refr; consequently, no member requires an external threshold, manual switch, or hidden field.
4. On the FRW background ∇_i theta_eff = ∇_i Psi = 0, so I_grad = 0, J_refr = 0, L_refr = 0. On a stationary inhomogeneous background I_grad > 0, and the same term is turned on without a change of law.
5. Consequently, the covariant reduced canon is closed.

Theorem 2. Absence of ghost and extra modes in reduced canon

Let the following conditions be fulfilled:
1. theta_eff != 0 in the considered region;
2. gamma >= 0, V0 >= 0;
3. ctheta2 > 0, cpsi2 > 0;
4. F(0)=0, F'(J)>0, F''(J)<=0 for J>=0;
5. J_refr >= 0;
6. Lapse/Shift do not have their own kinetic terms and remain algebraic constraints;
7. for stationary operator M_11 > 0, M_22 > 0, det(M) > 0 is fulfilled.

Then:
1. FRW kinetic matrix is positive definite;
2. FRW gradient matrix is positive definite;
3. stationary radial operator is elliptical;
4. additional propagating modes are absent;
5. ghost/gradient instabilities are excluded in the working region of the reduced canon.

Proof:
1. FRW kinetic matrix in the current canon is equal to K_FRW = diag(3/theta_eff² + gamma², 1). At theta_eff != 0 and gamma >= 0 diagonal elements are strictly positive, meaning K_FRW ≻ 0.
2. FRW gradient matrix is given as G_FRW = diag(c_theta² K_theta, c_Psi²). At c_theta² > 0, c_Psi² > 0 and K_theta > 0 it is also positive definite. Consequently, gradient instability is absent.
3. F'(J)>0 and F''(J)<=0 mean that refractive response is monotonic and saturating, rather than accelerating. This excludes sign flip and runaway in the refractive sector.
4. M_11 > 0, M_22 > 0, det(M) > 0 mean positive definiteness of the local quasilinear operator. Thus stationary PDE is elliptical, and hidden propagating branch does not arise.
5. Lapse/shift in ADM-skeleton do not receive kinetic terms, therefore they do not add new physical modes; they remain constraints rather than propagating degrees of freedom.
6. Combining points 1-5, we get absence of ghost, gradient and extra modes in the reduced canon on the entire working area where listed prerequisites are fulfilled.

Consequence

If the reduced canon passes stability_register(...).healthy == True, then it passed:
1. covariant closure at the level of current S_red;
2. checking for ghost/gradient/tachyon/extra-mode in canonical working area;
3. baseline-condition for transition to 1-2 primary observables and only then to fit.

Theorem 3. Full perturbation register closes ADM-canon

Let the conditions of Theorem 1 and Theorem 2 be fulfilled, as well as:
1. for quadratic sector matrices K, G, M are given in the form L^(2) = 1/2 dq̇ᵀ K dq̇ - 1/2 dqᵀ[ G(k²/a²) + M ]dq;
2. Lapse/Shift are eliminated via Schur complement;
3. direct and Schur channels give the same result for K, G, M;
4. stability_register gives healthy = True on the working area.

Then:
1. ADM-reduction is consistent with covariant action register;
2. quadratic perturbation sector does not contain a hidden split between direct and Schur completion;
3. working perturbation proof is closed within the limits of the current canon.

For this purpose, a single register exists in the repository:
- theory/perturbation_proof.py

Its outputs:
1. S_cov as a full covariant completion;
2. K_direct/G_direct/M_direct and K_schur/G_schur/M_schur;
3. routes_match;
4. health;
5. proof_closed.

6.3) Minimal derived S_micro

S_micro -> S_red -> S_cov is fixed as a working hierarchy. In order not to start the canon with S_cov as a postulate, below is set a minimal derived microphysics layer, derived from the Theta_munu, Psi basis. It is needed only as a local layer from which current S_red is obtained by reduction, and S_cov by adding matter completion. It does not introduce new geometry.

6.3.1) Fields and structure

Fundamental fields of derived layer:
1. Theta_munu - clock sector;
2. Psi - energy sector;
3. χ_a - material fields;
4. local deviation scalar R_micro, fixed as the minimal local quadratic deviation from the canonical weak field point;
5. derived preferred-slicing volume form measure_dt_d3x = dt d^3x;
6. local response coefficients Z_t(R_micro) and Z_s(R_micro).

Minimal local action in preferred slicing:
S_micro = ∫ dt d^3x [ L_Theta^micro + L_Psi^micro + L_matter^micro + L_int^micro ].

Composition:
R_micro = (theta_eff/theta_0 - 1)^2 + Psi^2/M_pl^2 + chi^2/M_pl^2 + |nabla theta_eff|^2/Lambda_theta^4 + |nabla Psi|^2/Lambda_Psi^4 + |nabla chi|^2/Lambda_chi^4,
A_Theta = 1 + R_micro,
B_Theta = gamma^2 + R_micro,
A_Psi = 1 + R_micro,
B_Psi = 1 + R_micro,
Z_t = 1 + R_micro,
Z_s = 1 + 2 * R_micro,
V_Theta_micro = 1/2 * m_theta^2 * theta_eff^2 + 1/4 * lambda_theta * theta_eff^4,
V_Psi_micro = 1/2 * V0 * R^2 + 1/4 * lambda_Psi * R^4,
V_matter_micro = 1/2 * m_m^2 * chi^2 + 1/4 * lambda_m * chi^4,
L_Theta_micro = 3/2 * M_pl^2 * A_Theta * |nabla ln theta_eff|^2 + 1/2 * B_Theta * |nabla theta_eff|^2 - V_Theta_micro,
L_Psi_micro = 1/2 * A_Psi * |nabla R|^2 - V_Psi_micro,
L_matter_micro = 1/2 * Z_t * chi_dot^2 - 1/2 * Z_s * |nabla chi|^2 - V_matter_micro,
L_int_micro = (R_micro / M_pl^2) * (theta_eff^2 * R^2 + theta_eff^2 * chi^2).

6.3.2) Weak field alignment

In weak field limit:
R_micro → 0,
A_Theta → 1, B_Theta → gamma², A_Psi → 1, B_Psi → 1,
Z_t → 1, Z_s → 1.

Alignment-rule for refractive sector:
n² = Z_s / Z_t.

Consequently:
1. n → 1 in homogeneous weak field limit;
2. the difference from weak field lives in local response coefficients, not in a new ontology of geometry;
3. S_red is a reduced/covariant completion of this local action after integrating heavy modes;
4. S_cov is obtained as S_red + S_m[g_hat, matter], i.e., through adding matter completion on top of reduced layer.
5. In code the bridge is fixed as coincidence of S_red_from_micro = S_red and S_cov_from_micro = S_cov at the level of canonical registers.
6. The transition from S_micro to S_red is now algebraic: it goes through R_micro -> 0, rather than through free functions of coefficients or a separate knob for response force.

6.3.2.1) Closing items 1-3: propagation law

This sub-item fixes why g/g_hat, universality of coupling and weak field alignment are not separate postulates on top of Theta–Psi, but are read as one chain of derivation.

Item 1: g/g_hat status.

In S_micro primary remain only:
Theta_munu, Psi, χ_a, R_micro, Z_t(R_micro), Z_s(R_micro).

For matter-field χ local quadratic block has the form:
L_matter^micro =
  1/2 Z_t(R_micro) χ̇²
  - 1/2 Z_s(R_micro) |∇χ|²
  - V_matter(χ).

The ratio of coefficients sets the local propagation cone:
c_matter² = Z_s/Z_t.

Therefore the record

n² = Z_s/Z_t

is not a new geometry, but a compressed way to say: local clocks, light signals and matter-processes use the same response ratio. g_hat is introduced only after this as bookkeeping-record:
S_matter[χ; Z_t,Z_s] -> S_matter[χ; g_hat],
g_hat = g + (n² - 1) u ⊗ u.

Consequently, g and g_hat are not a third ontological entity. They are allowed only if they do not carry independent degrees of freedom beyond Theta_munu, Psi and response coefficients.

Item 2: universality.

Universality of coupling means:
∀χ_a: Z_t^a(R_micro) = Z_t(R_micro),
      Z_s^a(R_micro) = Z_s(R_micro).

Then for any matter-sectors:
n_a² = Z_s^a/Z_t^a = Z_s/Z_t = n².

If species-dependent coefficients appear

Z_t^a != Z_t^b or Z_s^a != Z_s^b,

then it is already not the current Theta–Psi canon: different substances see different local propagation laws, and the principle of unified clock-sector response breaks. Therefore, universality of Z_t, Z_s is an mandatory acceptance criterion for derived micro layer.

Item 3: weak field alignment.

Weak field point is defined as:
R_micro -> 0,
Z_t -> 1,
Z_s -> 1,
n² -> 1,
n -> 1.

The first physically significant correction is recorded through weak field scalar Phi_eff:
n(theta_eff) = 1 + 2 Phi_eff + O(Phi_eff²).

This is not a new ontology, but a normalization of observab <= layer. It fixes one coefficient so that one and the same n gives consistent weak field effects:
Δnu/nu ≈ Δ ln n,
t_signal = ∫ n dl,
α_lens ≈ ∫ ∇_⊥ ln n dl,
a_probe = -∇ ln n.

In other words, weak field alignment is closed in the canon if and only if:
1. R_micro -> 0 returns inert local regime;
2. all matter-sectors have the same n;
3. linear normalization n = 1 + 2 Phi_eff + O(Phi_eff²) is fixed before fit;
4. higher-order corrections do not introduce composition-dependent propagation.

Total chain:

local Theta–Psi basis
  -> R_micro
  -> universal Z_t,Z_s
  -> n² = Z_s/Z_t
  -> g_hat as bookkeeping only
  -> weak field observables

This is the closure of items 1-3 in the current canon. In code it is secured by propagation_law_closure_register() function in theory/micro_action.py and test_propagation_law_closure_pins_points_1_2_3 test:
1. g_hat_independent_dof = 0;
2. Z_t_chi_a = Z_t_chi_b, Z_s_chi_a = Z_s_chi_b;
3. R_micro -> 0 gives Z_t -> 1, Z_s -> 1, n -> 1;
4. linear observable-normalization n = 1 + 2 Phi_eff + O(Phi_eff²) is fixed before fit.

Consequently, items 1-3 are considered closed at the level of current canon, register and tests. Further work may refine higher-order corrections, but must not change this chain without explicit canon change.

6.3.2.1a) Anti-fit guardrails

To prevent derived micro layer from turning into a set of knobs to get required numbers, explicit separation is introduced:
canonical = fixed before fit or derived from the action,
diagnostic = useful for scans, but not a physical claim.

In current canon are fixed:
primary ontology = Theta_munu, Psi,
theta_eff² = Theta_munu u^mu u^nu,
Phi_eff = ln(theta_eff/theta_0),
weak field slope n ≈ 1 + 2Phi_eff,
weak field profile n(Phi_eff) = 1 + 2 asinh(Phi_eff),
n² = Z_s/Z_t,
EPR singlet rule E(a,b) = -a·b,
particle = stable localized Theta-projection of a Psi-sector.

Into diagnostic-only layer are moved:
weak field profile alternatives: linear/exp/tanh2,
legacy explicit profile selector,
F_v1(J) = 1 - exp(-J),
ADM eps*/alpha* mixings,
legacy/exploratory sigma extension,
strong field bridge calibration coefficients,
inverse-BAO required multipliers.

They cannot be used as physical derivation. Any such element can be promoted to canon only if simultaneously fulfilled:
1. derivation from action or from a smaller set of postulates;
2. fixation before viewing data-facing fit;
3. absence of choice by test set;
4. absence of new light modes without observational price;
5. independent observational prediction.

Forbidden moves:
choose_profile_after_data,
retune_bridge_on_test,
promote_inverse_diagnostic_to_correction,
use_adm_eps_without_derivation,
add_sigma_as_physical_field_without_new_canon.

Code fixation: anti_fit_register() and canonical_promotion_check() in theory/anti_fit.py, tests in tests/test_anti_fit.py.

6.3.2.1b) What is already derived from action

From current minimal Theta–Psi action for q = (d_theta, d_R, deltaφ) follows diagonal high-k quadratic block:
K = diag(3/theta² + gamma², 1, R²),
G = diag(3/theta² + gamma², 1, R²),
M = diag(0, V0 + 3Q²/(a⁶R⁴), 0).

Consequence:
epsK* = 0,
epsG* = 0,
alpha_k2* = 0,
c_s² = (1, 1, 1).

This is not an additional fit, but a negative derivation from action: if in action are no non-diagonal kinetic/gradient-operators and higher-spatial derivative operators, then corresponding mixings in quadratic action do not appear. Non-zero epsK*, epsG*, alpha_k2* can enter canon only together with explicitly recorded new operators in S_cov.

For refractive sector current action derives the argument

J_refr = Phi_eff · I_grad,
L_refr = M_*⁴ F(J_refr),

and requirements for the function:
F(0)=0,
F'(J) ≥ 0,
finite response for the no-horizon branch.

But choice of exponential form itself from current action does not follow. It follows only that J_* is not a separate physical parameter and can be absorbed into J_refr normalization.

Strong field bridge (theta_c_scale, m2_scale, lam_scale, powerlaw2d) also do not follow from current action. For their derivation one needs to solve sourced static Euler-Lagrange equations with boundary conditions, rather than calibrate coefficients by EHT/GW target-set.

Code fixation: theory/action_derivation.py, tests/test_action_derivation.py.

Principle which makes exponent derivation rather than choice of form:
S(J) = 1 - F(J),
S(J1 + J2) = S(J1) S(J2).

That is independent increments of J act without memory and multiplicatively on the remaining unused part of response. With continuity such law gives

S(J) = exp(-J/J_*),
F(J) = 1 - exp(-J/J_*).

In this record J_* is not physical knob anymore, but only normalization of response remainder.

6.3.2.2) Vacuum offset source-map

Vacuum catastrophe check in current canon is formulated not as statement "zero energy is absent". Zero energies of matter-modes can exist and manifest in boundary-differences, for example in Casimir effect. Something else is checked: absolute additive constant in effective matter action is not a direct source for Theta–Psi.

Minimal record:
Γ_matter^eff = Γ_response[Z_t,Z_s,n,boundary data] + ∫ dt d³x Λ_vac_offset.

Source map for Theta_munu, Psi basis is taken only from variations by local response data, i.e. by those quantities which change propagation regime:
source_ThetaPsi deltaΓ_eff/d_theta_eff, deltaΓ_eff/d_R, deltaΓ_eff/d_R_micro.

For absolute offset:
deltaΛ_vac_offset/d_theta_eff = 0,
deltaΛ_vac_offset/d_R = 0,
deltaΛ_vac_offset/d_R_micro = 0.

Meaning such offset does not enter Theta–Psi source as huge cosmological constant. At same time boundary/response-contributions remain physical:
Γ_boundary ΔΓ_boundary · n(R_micro),
Γ_response Δresponse · n(R_micro),
deltaΓ_boundary/d_R_micro != 0,
deltaΓ_response/d_R_micro != 0.

In other words:
1. absolute sum 1/2 sum ħω is not an observable Theta–Psi response itself;
2. Casimir is preserved because it is difference of spectra at different boundary data;
3. local QFT preserves zero modes as part of effective action, but gravitating/source-contribution goes through response map, rather than through arbitrary additive constant.

Code fixation: vacuum_offset_source_register() and test_vacuum_offset_source_register_* tests in tests/test_micro_action.py.

6.3.3) Status

This is derived layer from the Theta_munu, Psi basis, not separate ontology. It is needed to:
1. not start canon with S_cov as postulate;
2. keep explicit n² = Z_s/Z_t rule;
3. show where exactly matter response appears.

6.3.3.1) Closing stationary strong field branch

Stationary strong field branch is closed by separate register:
strongfield_branch_closure_register()

This register closes stationary profile. Its scope:
closure_scope = stationary_strongfield_branch_only,
dynamic_collapse_closed = handled_by_dynamic_collapse_closure_register.

Criteria of stationary branch closure:
1. theta_positive = 1: theta_eff remains on positive branch;
2. n_finite = 1: n profile is finite;
3. finite_redshift = 1: z_surface is finite;
4. finite_delay = 1: optical delay is finite;
5. outgoing_characteristic_exists = 1: outgoing characteristic does not disappear;
6. finite_energy = 1: radial integral of energy is finite;
7. radial_operator_elliptic = 1: stationary radial operator does not lose ellipticity;
8. regular_center = 1: center of profile is regular.

Boundary conditions of stationary branch:
theta_eff(0) > 0,
theta_eff'(0) = 0,
R'(0) = 0,
theta_eff(r -> infinity) = theta_0,
R(r -> infinity) = R_0.

Direct map is fixed as mandatory chain:
n_profile
  -> z_surface
  -> delay
  -> ray_map_from_n_profile
  -> transfer_function_from_finite_delay.

Branch is falsified if at least one is fulfilled:
theta_eff reaches zero,
n diverges,
energy diverges,
radial operator is not elliptic,
no outgoing characteristic,
direct observables incompatible.

In code this is secured by test_strongfield_branch_closure_pins_horizonless_criteria test.

6.3.3.2) Closing dynamic collapse as evolution contract

Dynamic collapse is closed not as numerical PDE-solver, but as canonical evolution contract:
dynamic_collapse_closure_register()

Scope:
closure_scope = dynamic_collapse_register_contract,
dynamic_collapse_closed = 1,
requires_numeric_pde_implementation = 1.

Meaning: solver can be different, but it is obliged to implement the same contract. Evolution variables:
theta_eff(t,r),
R(t,r),
n(t,r),
E_total(t),
min theta_eff(t),
max n(t).

Criteria of closure:
1. initial_finite_energy = 1;
2. theta_eff_positive_all_t = 1;
3. n_finite_all_t = 1;
4. energy_finite_all_t = 1;
5. no_blow_up_before_branch_decision = 1;
6. outgoing_characteristic_exists_all_t = 1;
7. hyperbolic_evolution_operator = 1;
8. stationary_or_dispersive_outcome = 1.

Allowed outcomes:
settle_to_stationary_branch,
disperse_to_weak_field.

Forbidden outcomes:
theta_eff_zero,
n_infinite,
energy_blow_up,
loss_of_hyperbolicity,
no_outgoing_characteristic.

Total chain:
finite_energy_initial_data
  -> evolve_Theta_Psi
  -> track_theta_min_and_n_max
  -> preserve_outgoing_characteristic
  -> avoid_blow_up
  -> settle_or_disperse
  -> reject_invalid_outcomes.

In code this is secured by test_dynamic_collapse_closure_contract_pins_evolution_criteria test.

6.3.4) Reduction check

To ensure the bridge is not declarative but computable, we fix a weak-field reduction test:
1. A_Theta -> 1, B_Theta -> 1, A_Psi -> 1, B_Psi -> 1;
2. Z_t -> 1, Z_s -> 1;
3. n² -> 1 and consequently n -> 1;
4. L_Theta^micro and L_Psi^micro return to canonical quadratic blocks;
5. L_matter^micro returns to canonical local matter kinetic block;
6. heavy modes are integrated and do not remain separate propagating sector.

This test does not replace physical derivation from Theta_munu, Psi basis. It fixes that derived microphysics layer does not break weak field alignment and indeed leads to already accepted reduced canon.

6.3.5) Micro perturbation check

Separately fix minimal perturbation sanity for derived micro layer:
1. A_Theta > 0, B_Theta > 0, A_Psi > 0, B_Psi > 0;
2. Z_t > 0, Z_s > 0;
3. c_Theta² = B_Theta / A_Theta > 0, c_Psi² = B_Psi / A_Psi > 0;
4. n² = Z_s / Z_t > 0;
5. no_new_modes == True on working branch.

This is not a full UV-proof at the level of fundamental microphysics. This is an honest micro-perturbation barrier: derived layer must be healthy already at the local quadratic level, otherwise it cannot be raised higher.

6.3.6) Derivation ladder

To prevent derived layer from dispersing into separate crutches, explicit derivation ladder is fixed in code and canon:
1. S_micro;
2. weak_field_projection;
3. integrate_heavy_modes;
4. assemble_S_red;
5. add_S_matter_completion;
6. S_cov.

This ladder does not claim UV-uniqueness. It fixes that current derived microphysics layer is not just "compatible" with Theta_munu, Psi basis canon, but leads to it in the same reduction order.

6.3.7) Strong field consequence

From derived micro layer we preserve only one strong-field output:
1. finite redshift on compact branch;
2. finite delay on compact branch;
3. ok_no_horizon == True.

This is not a second channel with separate physics and not a set of proxies. It is one strong field consequence read through micro_strongfield_primary_register() and micro_to_strongfield_primary_bridge(). It fixes that derived strong field consequence does not disperse by alternative strong-field indicators, but reduces to one observable class: compact_no_horizon.

6.3.8) Symbol audit

For public analysis each micro-layer entity gets a status:
1. derived - derived from projection, alignment or reduction from the Theta–Psi basis;
2. assumed - minimal object without which reduced completion does not assemble;
3. diagnostic - final check, not basic physics.

In current derived micro layer this is distributed like this:
derived: theta_eff, measure_dt_d3x, R_micro, A_Theta, B_Theta, A_Psi, B_Psi, Z_t, Z_s,
         n, n^2, S_micro, S_red, S_cov, bridge_to_covariant
assumed: R, chi
diagnostic: finite redshift, finite delay, ok_no_horizon

This is honest public form of answer to question "why exactly these coefficients": most symbols are not chosen arbitrarily, but are marked as derived from reduction rules; measure_dt_d3x = dt d^3x - derived preferred-slicing volume form, and R_micro - minimal local quadratic deviation scalar derived from locality, parity, positivity and weak field normalization.

7) Canonical n(theta_eff) class

Family-v1 (current working class):
n(theta_eff; p, thetas) = 1 + 𝒩(p · thetas / theta_eff),
p>0, thetas>0.

Properties of 𝒩:
1. 𝒩(0)=0;
2. 𝒩 is monotonic on the working branch;
3. on the physical branch n remains finite without introducing a separate ceiling as a physical entity.

Rules:
1. The family is fixed before fitting.
2. Changing the family after viewing the result is prohibited.
3. Alternatives are allowed only as a pre-announced blind-test.

8) Hard rejection criteria

The model/branch is rejected if at least one item is stably violated:
1. Weak field tests (Mercury/Shapiro/redshift/light deflection).
2. Stability (ghost/gradient) in the working parameter range.
3. Joint CC+SN+BAO+strong field circuit.
4. Agreement is achieved only by hidden tuning on the test set.
5. There is no stable solution for a fixed n(theta_eff) class.

9) Psi Status (mandatory nested-test)

The Psi sector is considered canonical only upon passing:
1. M0: Theta-only (without dynamic Psi);
2. M1: Theta+Psi;
3. single train/test protocol;
4. acceptance criterion: ΔBIC = BIC(M1)-BIC(M0) < -6 stably on k-fold.

If ΔBIC >= 0 in most folds, Psi is excluded from the canon.

10) Link with code

- Full tensor version: canonical physical layer of documentation.
- Current scripts background/*, fitting/*: reduced theta_eff layer.
- In case of divergence, priority is given to the canon: first the mathematical formulation is corrected, then the code.
- Symbolic formula verification: theory/tensor_symbolic.py — exploratory, not a canon source.

10.1) tensor_symbolic.py API

This API is used as a symbolic sanity-check for reduced formulas. It does not set the canon and does not replace background/frw_symbolic.py.

from theory.tensor_symbolic import (
    kinetic_term_from_tensor, # 3/2 (thetȧ/theta)² from Theta_munu
    mixing_term_gamma, # gamma²/2 thetȧ² from L_int
    full_effective_lagrangian, # full L_eff
    energy_density, # ρ(theta,Psi)
    pressure, # p(theta,Psi)
    equation_of_state, # w = p/ρ
    weak_field_limit, # n ≈ 1 + κ Phi_eff
    tensor_wave_equation, # h'' + 2ℋh' + k²h = 0
    check_dimensions, # dimensions check
)

Tests: tests/test_tensor_symbolic.py (14 tests, all PASS).

11) Appendices

A. Symbol Table

| Symbol | Definition | Dimension (8πG=1, c=1) |
|--------|-------------|--------------------------|
| Theta_munu | Time tensor | [mass]² |
| theta_eff | √(Theta_munu u^mu u^nu) | [mass] |
| Psi | Complex field | [mass] |
| gamma | Theta-Psi coupling | dimensionless |
| V0 | Psi mass | [mass]² |
| Q | Phase charge | [mass]⁰ |
| n(theta_eff) | Refractive index | dimensionless |
| g_hat_munu | Effective metric | [mass]⁰ |

B. L_eff dimension check

In units 8πG=1, c=1:
[L_eff] = [mass]⁴ (energy density).
[M_pl^2 * (theta_dot/theta)^2] = [mass]^2 * [mass]^2 = [mass]⁴.
[R_dot^2] = [mass]⁴.
[V0 R^2] = [mass]^2 * [mass]^2 = [mass]⁴.
[gamma^2 theta_dot^2] = [mass]^0 * ([mass]^2)^2 = [mass]⁴.
All terms have correct dimension [mass]⁴. ✓
