"""Tests for chemistrykit.analytical.systems.chromatography against closed-form/known results."""

import numpy as np
import pytest

from chemistrykit.analytical.systems.chromatography import (
    kovats_retention_index,
    minimum_plate_height,
    optimum_flow_velocity,
    plate_height,
    purnell_resolution,
    resolution,
    retention_factor,
    selectivity_factor,
    theoretical_plates,
    van_deemter_H,
)


def test_theoretical_plates_base_width_formula():
    assert theoretical_plates(10.0, 0.5, width_type="base") == pytest.approx(6400.0)


def test_theoretical_plates_base_and_half_height_agree():
    sigma = 0.2
    N_base = theoretical_plates(10.0, 4.0 * sigma, width_type="base")
    N_half = theoretical_plates(10.0, 2.0 * np.sqrt(2.0 * np.log(2.0)) * sigma, width_type="half_height")
    # 5.545 is a rounded textbook constant (exactly 8*ln2 = 5.5452...), so
    # the two width conventions agree only to that rounding's precision.
    assert N_base == pytest.approx(N_half, rel=1e-3)


def test_theoretical_plates_rejects_bad_width_type():
    with pytest.raises(ValueError):
        theoretical_plates(10.0, 0.5, width_type="bogus")


def test_plate_height_matches_definition():
    assert plate_height(25.0, 6400.0) == pytest.approx(25.0 / 6400.0)


def test_van_deemter_H_matches_formula():
    assert van_deemter_H(2.0, A=1.0, B=2.0, C=0.05) == pytest.approx(1.0 + 1.0 + 0.1)


@pytest.mark.parametrize(("A", "B", "C"), [(1.0, 2.0, 0.05), (0.5, 5.0, 0.02), (2.0, 1.0, 0.1)])
def test_van_deemter_H_minimum_matches_dH_du_zero(A, B, C):
    u_opt = optimum_flow_velocity(B, C)
    u = np.linspace(0.01, 20.0, 2_000_000)
    H = van_deemter_H(u, A, B, C)
    u_numeric_min = u[np.argmin(H)]
    assert u_numeric_min == pytest.approx(u_opt, rel=1e-3)
    assert van_deemter_H(u_opt, A, B, C) == pytest.approx(minimum_plate_height(A, B, C))


def test_van_deemter_H_at_optimum_is_global_minimum():
    A, B, C = 1.0, 2.0, 0.05
    u_opt = optimum_flow_velocity(B, C)
    H_min = van_deemter_H(u_opt, A, B, C)
    for u_test in [u_opt * 0.5, u_opt * 0.9, u_opt * 1.1, u_opt * 2.0]:
        assert van_deemter_H(u_test, A, B, C) > H_min


def test_retention_factor_matches_definition():
    assert retention_factor(12.0, 2.0) == pytest.approx(5.0)


def test_selectivity_factor_matches_definition():
    assert selectivity_factor(2.0, 5.0) == pytest.approx(2.5)


def test_resolution_matches_definition():
    assert resolution(9.0, 10.0, 0.5, 0.5) == pytest.approx(2.0)


def test_resolution_increases_with_peak_separation():
    r_close = resolution(9.5, 10.0, 0.5, 0.5)
    r_far = resolution(8.0, 10.0, 0.5, 0.5)
    assert r_far > r_close


def test_kovats_index_of_alkanes_is_100_times_carbon_number():
    t0, t8, t9 = 0.8, 4.0, 7.5
    assert kovats_retention_index(t8, t8, t9, n=8, dead_time=t0) == pytest.approx(800.0)
    assert kovats_retention_index(t9, t8, t9, n=8, dead_time=t0) == pytest.approx(900.0)


def test_kovats_index_follows_log_linear_alkane_series():
    t0 = 1.0
    t_adj = lambda c: 0.05 * 1.8**c  # noqa: E731 -- log t' linear in carbon number
    tx = t0 + t_adj(9.37)
    assert kovats_retention_index(tx, t0 + t_adj(9), t0 + t_adj(10), n=9, dead_time=t0) == pytest.approx(937.0)
    assert kovats_retention_index(tx, t0 + t_adj(8), t0 + t_adj(12), n=8, N=12, dead_time=t0) == pytest.approx(937.0)


def test_purnell_matches_direct_resolution_for_equal_widths():
    for N, k1, alpha in [(2500.0, 2.0, 1.05), (10000.0, 5.0, 1.2), (40000.0, 0.5, 1.02)]:
        k2 = alpha * k1
        t0 = 1.3
        tR1, tR2 = t0 * (1 + k1), t0 * (1 + k2)
        w = 4 * tR2 / np.sqrt(N)
        assert purnell_resolution(N, alpha, k2) == pytest.approx(resolution(tR1, tR2, w, w))
