"""Tests for chemistrykit.analytical.systems.chromatography against closed-form/known results."""

import numpy as np
import pytest

from chemistrykit.analytical.systems.chromatography import (
    minimum_plate_height,
    optimum_flow_velocity,
    plate_height,
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
