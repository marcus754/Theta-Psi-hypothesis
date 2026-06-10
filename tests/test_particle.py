import pytest

from theory.particle import (
    LocalizationProfile,
    QuantumLabels,
    classical_limit_ledger,
    conserved_label_condition,
    detector_event_map,
    localization_condition,
    multiparticle_state_ledger,
    particle_ledger,
)


def test_particle_ledger_is_derived_not_primary_ontology():
    d = particle_ledger()
    assert d["new_ontology"] == 0
    assert d["primary_ontology"][0].name == "Theta_mu_nu"
    assert d["primary_ontology"][1].name == "Psi"
    assert d["definition"].name == "stable_localized_Theta_projection_of_Psi_sector"


def test_localization_condition_requires_finite_profile():
    ok = localization_condition(LocalizationProfile(norm=1.0, width=2.0, energy=3.0))
    assert ok["localized"] is True

    bad = localization_condition(LocalizationProfile(norm=1.0, width=0.0, energy=3.0))
    assert bad["localized"] is False
    assert bad["finite_width"] is False


def test_conserved_labels_define_repeatable_particle_sector():
    electron = QuantumLabels(mass=1.0, charge=-1.0, spin=0.5)
    same = QuantumLabels(mass=1.0, charge=-1.0, spin=0.5)
    changed = QuantumLabels(mass=1.0, charge=0.0, spin=0.5)

    assert conserved_label_condition(electron, same)["labels_conserved"] is True
    assert conserved_label_condition(electron, changed)["labels_conserved"] is False
    assert conserved_label_condition(electron, changed)["charge_conserved"] is False


def test_detector_event_map_thresholds_local_theta_event():
    click = detector_event_map(coupling=2.0, threshold=1.0, local_density=0.75)
    assert click["click"] is True
    assert click["event_type"].name == "Theta_local_detector_event"

    no_click = detector_event_map(coupling=0.5, threshold=1.0, local_density=1.0)
    assert no_click["click"] is False

    with pytest.raises(ValueError):
        detector_event_map(threshold=0.0)


def test_multiparticle_state_distinguishes_product_and_entangled_sectors():
    entangled = multiparticle_state_ledger(entangled=True)
    product = multiparticle_state_ledger(entangled=False)

    assert entangled["independent_complete_states"] == 0
    assert product["independent_complete_states"] == 1
    assert entangled["local_theta_events"] == 2
    assert product["local_theta_events"] == 2
    assert "!=" in entangled["factorization"].name
    assert "=" in product["factorization"].name


def test_classical_limit_requires_stationary_phase_and_decoherence():
    ok = classical_limit_ledger(action_scale=1000.0, hbar_scale=1.0, decoherence_strength=1.0)
    assert ok["classical_trajectory_limit"] is True

    too_quantum = classical_limit_ledger(action_scale=1.0, hbar_scale=1.0, decoherence_strength=1.0)
    assert too_quantum["classical_trajectory_limit"] is False

    coherent = classical_limit_ledger(action_scale=1000.0, hbar_scale=1.0, decoherence_strength=0.0)
    assert coherent["classical_trajectory_limit"] is False

    with pytest.raises(ValueError):
        classical_limit_ledger(hbar_scale=0.0)


def test_particle_acceptance_pins_all_required_checks():
    d = particle_ledger()
    acc = d["acceptance"]
    assert acc["localized"] is True
    assert acc["labels_conserved"] is True
    assert acc["detector_click_possible"] is True
    assert acc["classical_limit_exists"] is True
    assert d["entangled_pair_state"]["independent_complete_states"] == 0
