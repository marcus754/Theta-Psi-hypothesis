from theory.tensor import tensor_properties


def test_tensor_speed_is_luminal():
    tp = tensor_properties()
    assert tp.c_T2 == 1.0

