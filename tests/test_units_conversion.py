from background.components import H0_from_km_s_Mpc


def test_H0_conversion_reasonable_scale():
    # 70 km/s/Mpc should result in a very small dimensionless H0 ~ 1e-61
    H0 = H0_from_km_s_Mpc(70.0)
    assert H0 > 0.0
    assert 1e-62 < H0 < 1e-60

