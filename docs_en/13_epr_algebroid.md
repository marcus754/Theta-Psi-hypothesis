13. EPR/CHSH AND ALGEBROID LAYER

GOAL

This section fixes the minimal verifiable formulation of the image:
Theta slices/layers Psi as a local causal history.

The verification should not turn into a free model of correlations. Therefore, the first test is chosen to be as strict as possible: a singlet pair and CHSH.

MINIMAL STRUCTURE

A -> M
E -> M
ρ_Theta : A -> TM
Psi_AB ∈ Γ(E)

Meaning:
- M — the observable spacetime of local events;
- E — the state-bundle in which Psi is recorded;
- A — a layer of allowed internal/state transformations;
- ρ_Theta — the anchor set by Theta;
- ker ρ_Theta contains correlation generators: they physically set the connectivity of Psi but do not project into a controllable signal in TM.

Canonical phrase:
Psi carries the indecomposability of the state.
Theta sets the causal anchor through which this indecomposability manifests as a local history without a superluminal signal.

CHSH SANITY-CHECK

For the first test, free parameters are prohibited:
free_fit_parameters = 0
free_correlation_function = 0

The state is fixed as a singlet:
Psi_AB = singlet.

The only correlation rule:
E(a, b) = - a · b.

For optimal directions:
a0 = (1, 0, 0),
a1 = (0, 1, 0),
b0 = (1, 1, 0) / sqrt(2),
b1 = (1, -1, 0) / sqrt(2).

CHSH:
S = E(a0, b0) + E(a0, b1) + E(a1, b0) - E(a1, b1),
|S| = 2 sqrt(2).

Consequently:
|S| > 2 violates the local hidden-variable boundary,
|S| <= 2√2 does not exceed the Tsirelson bound.

NO-SIGNALLING

Joint probabilities are taken as standard for a singlet:
P(A=α, B=beta | a, b) = 1/4 [1 + αbeta E(a, b)], α, beta ∈ {-1, +1}.

Then the marginals:
P(A=α | a, b) = 1/2,
P(B=beta | a, b) = 1/2.

They do not depend on the remote setting. This means entanglement remains a physical Psi-connectivity but does not become a channel for information transmission.

ANCHOR CONDITION

The algebroid formulation separates correlation and signal:
ρ_Theta(C_AB) = 0,
ρ_Theta(S_signal) ⊂ causal cone(TM).

C_AB — correlation generator. It is responsible for the indecomposability of Psi_AB, but its anchor is zero: it does not provide a controllable displacement/signal in the observed spacetime.

S_signal — signal generator. If there is a signal, its projection via ρ_Theta must remain within the causal cone.

ACCEPTANCE CRITERIA

The model passes the first EPR sanity-check only if simultaneously:
1. |S_CHSH| > 2;
2. |S_CHSH| <= 2√2;
3. no-signalling marginals are equal to 1/2;
4. ρ_Theta(correlation generator) = 0;
5. signal anchor is causal;
6. there is no free correlation function;
7. there are no fit parameters for CHSH.

WHAT IS EXPLICITLY REJECTED

No-signalling by itself is insufficient. A PR-box has:
E00 = 1,
E01 = 1,
E10 = 1,
E11 = -1,
S_CHSH = 4.

It does not provide a controllable signal but exceeds the quantum boundary:
4 > 2√2.

Therefore, the acceptance rule:
no-signalling = required but not sufficient,
|S_CHSH| > 2√2 -> reject.

This protects the layer from turning into a model that allows any no-signalling correlations.

FOLIATION INDEPENDENCE

If Theta sets the slicing/anchor, the observed probabilities should not depend on which description event A stands before B or B before A:
P_AB(α, beta | a, b; A-before-B) = P_AB(α, beta | a, b; B-before-A).

For a singlet, this is satisfied because the joint probability is set by a single Psi_AB, not by a dynamic signal from the first measurement to the second.

GHZ/MERMIN CONUALITY

CHSH shows the violation of the local hidden-variable boundary statistically. The next sanity-check is GHZ/Mermin, where the contradiction with global pre-set values becomes algebraic.

For GHZ constraints:
XYY = +1,
YXY = +1,
YYX = +1,
XXX = -1.

The product of the quantum right-hand sides is -1. But any nonconual hidden valuation with pre-set X_i, Y_i = ±1 gives a product of +1, because each local factor enters twice.

Result:
global_hidden_valuation_exists = 0,
conuality_required = 1.

This fixes that Psi is not a container of pre-made local values. It sets contextual state-connectivity.

BRIDGE TO THE QUANTUM LAYER

docs/04_quantum.md and this section do not introduce two different Psi. The bridge:
Psi_state_layer
  -> two_subsystem_entanglement
  -> singlet_projection
  -> fixed_CHSH_correlations
  -> Theta_anchor_no_signalling.

That is, the EPR layer is a special two-subsystem check of the same state layer described in the quantum section via states, covariances, and modes.

NO-FIT EXPERIMENTAL PROFILE

For any set of coplanar angles a0, a1, b0, b1, the prediction is set immediately:
Eij = -cos(ai - bj),
S = E00 + E01 + E10 - E11.

There are no fit parameters. A real experimental table can be substituted as a set of angles and measured correlations, but the theoretical profile is already fixed.

Code fixation:
theory.epr_algebroid.epr_algebroid_registry()
tests/test_epr_algebroid.py

WHAT THIS GIVES THE HYPOTHESIS

This layer does not prove the entire Theta–Psi hypothesis. It fixes that the image of Theta as a causal slicing/anchor for Psi can reproduce the known EPR picture without arbitrary knobs:
Bell violation: yes
Tsirelson bound: yes
no-signalling: yes
PR-box rejected: yes
GHZ global valuation rejected: yes
foliation-order dependence: no
free parameters: no

That is, the algebroid here is not a decoration but a way to distinguish the physical connectivity of Psi from the signal causality of Theta.
