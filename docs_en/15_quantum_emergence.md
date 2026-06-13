15. EMERGENCE OF QUANTUM AND ELECTRODYNAMIC EFFECTS

Section Status: This is a consistent reconstruction rather than a strict first-principles derivation from the action. It is shown that known structures (quantization, electrodynamics, the Schrödinger equation) allow for a natural interpretation within the Theta–Psi framework but are not derived from it in a unique way at this stage.

The EPR/CHSH check in the current version is considered closed:
- the quantum Tsirelson bound is reproduced;
- the absence of superluminal signal transmission is preserved;
- the correlation sector does not use free fitting parameters.

However, a full derivation of quantum mechanics from the base action is not yet complete.

TIME CANON

- The fundamental object of time in the hypothesis is the Theta_munu tensor (see docs/11_derivation_path.md).
- The formulas in this section use the reduced FRW variable theta ≡ theta_eff, where theta_eff^2 = Theta_munu u^mu u^nu, u^mu u_mu = -1.
- All notations of the form n(theta) in this file should be read as n(theta_eff, Theta, Psi).
- This section describes the computational reduced circuit rather than the full tensor dynamics of Theta_munu.

CANONICAL FORMULATION (MANDATORY FOR ALL SECTIONS)

In the project, Theta_munu and Psi are taken as the base; the reduced circuit does not require an external threshold and does not use manual branch cutting.
S = ∫ d^4x √(-g) [ (M_pl^2/2)R + L_m + L_ThetaPsi + M_*^4 F(J_refr) ],
J_refr(Theta,Psi) = Phi_eff(theta_eff) * I_grad(Theta,Psi).

On a homogeneous FRW background: I_grad = 0, J_refr = 0.

The document fixes a possible mechanism for the emergence of quantum and electrodynamic effects from the Theta–Psi basis. It is not claimed here that all of quantum theory has already been strictly derived from first principles. It is a matter of a consistent reconstruction in which known structures appear as natural limits of the Theta–Psi dynamics.

1. EMERGENCE OF PLANCK'S CONSTANT hbar

Step 1.1. Phase Charge
In the hypothesis, energy is represented by a complex field:
Psi = R · exp(i φ)
Global symmetry with respect to phase shift
φ -> φ + const
generates a conserved Noether current and a corresponding charge.

In the coordinates of preferred time slicing, the charge has the form:
Q = ∫ d³x √(-g) · R² · φ̇
Phase charge has the dimension of action, i.e., energy multiplied by time.

Step 1.2. Topological Sectors
Since the phase φ is defined on the circle S¹, the field allows for topologically non-trivial configurations.
For a closed loop around a defect where R = 0, the condition is met:
∮ dφ = 2πn, where n is an integer.
This defines the discrete winding sectors of the Psi field.

For stationary configurations with finite energy, the phase charge is linked to the topological number of the sector:
Q_n = n · q0
where q0 is the minimum stable portion of phase action.

Thus, quantization arises as a consequence of the topological structure of the Psi phase sector.

Step 1.3. Identification with hbar
In Theta–Psi, the minimum quantum of phase action is identified with the observed Planck's constant:
hbar = q0
In this interpretation, hbar is not an external fundamental constant but represents the minimum stable quantum of Psi phase dynamics.

2. PHYSICAL NATURE OF DE BROGLIE WAVES

Consider a localized soliton-like Psi configuration moving with velocity v:
Psi(x, t) = R(x − vt) · exp(i[φ₀(x − vt) + kx − ωt])
The phase part of the action defines the canonical momentum:
p = ħ · ∇φ
For the plane component:
grad(phi) = k = 2 pi / lambda
Consequently:
p = ħ · 2π/λ = h/λ
Thus, the de Broglie wavelength arises as the spatial period of the Psi phase structure.

3. EMERGENCE OF ELECTRODYNAMICS

Step 3.1. Local Phase Covariance
The Theta tensor sets the local structure of causal propagation. If the Psi phase dynamics are to remain locally consistent, global U(1) symmetry naturally extends to local symmetry:
Psi(x) -> exp(iα(x)) · Psi(x)
To maintain local covariance, a compensating coupling field is required.

Step 3.2. Gauge Field
The minimal locally invariant derivative has the form:
Dmu = ∂mu − i(e/ħ) Amu
The field A_mu arises as a gauge-connectivity ensuring the consistency of the Psi local phase structure.

Step 3.3. Dynamics and Propagation
The dynamics of the gauge sector are set by the standard invariant:
Fmunu = ∂mu Anu − ∂nu Amu
The low-energy action of the gauge sector has the form:
L_gauge = −1/4 · √(−g_hat) · Fmunu F^munu
All fields propagate in the effective metric:
g_hat_mu_nu = g_mu_nu + (n^2 - 1) · u_mu u_nu
where n(theta) is the universal refractive index of the Theta sector.

Variation of the action with respect to Amu gives:
∇mu F^munu = J^nu
Thus, electrodynamics arises as a low-energy gauge sector of Psi phase dynamics in a medium of variable causal structure.

4. NON-RELATIVISTIC LIMIT AND THE SCHRÖDINGER EQUATION

Step 4.1. Relativistic Psi Sector
The relativistic dynamics of Psi are set by the Klein–Gordon equation:
(Dmu D^mu + m²) Psi = 0

Step 4.2. Non-relativistic Limit
Let's isolate the fast rest phase:
Psi = psi · exp(-i m t / hbar)
In the limit of low velocities v << c, the slow envelope psi satisfies:
i hbar · dpsi/dt = [ -hbar^2 / (2m) · grad^2 + V_ext ] psi
Thus, the Schrödinger equation arises as the non-relativistic limit of the Psi sector.

Step 4.3. Gravitational Potential
In Theta–Psi, an external potential is a consequence of the inhomogeneity of the event-frequency field:
V_ext = m · Phi_N
where Phi_N is the Newtonian limit of the Theta sector.

5. UNIFIED LAYER

The full low-energy action is written as:
S = ∫ d4x sqrt(-g) · [ (M_Pl^2 / 2) R + L_Theta + L_gauge + L_Psi + L_refr - U_tot ]
where:
L_Theta = 3/2 · (thetȧ/theta)² + 1/2 · gamma² · thetȧ²
(L_Theta here is the kinetic part of the Theta sector; the potential term −1/2 V0 R² is included in L_Psi of the canonical L_eff)
L_gauge = −1/4 · F²
L_Psi = |DmuPsi|² − 1/2 V0 R²
L_refr = M_*⁴ · F(J_refr)

The refractive index is defined as:
n(theta) = 1 + 𝒩(J_refr)

CONCLUSION

In the framework of Theta–Psi, quantization is linked to the topological structure of the Psi phase sector. Electrodynamics arises as locally-covariant phase dynamics. The Schrödinger equation appears as the non-relativistic limit of the Psi sector. Gravity manifests via the inhomogeneity of the causal/event-frequency structure of Theta.

The model does not claim to be a complete derivation of all quantum theory from first principles, but it sets an internally consistent emergent framework in which the basic structures of modern physics arise as consistent limits of a single Theta–Psi dynamics.
