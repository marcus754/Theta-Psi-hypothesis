from checks.strongfield_reference import strongfield_public_reference_check


def test_strongfield_public_reference_is_within_one_sigma():
    ref = strongfield_public_reference_check()
    assert ref.n_events >= 3
    assert ref.within_tolerance
    assert ref.max_abs_pull_sigma < 1.0
