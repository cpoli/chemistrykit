"""Tests for chemistrykit.crystal.systems.crystal_chemistry: radius-ratio limits and the tolerance factor."""

import numpy as np
import pytest

from chemistrykit.crystal.systems.crystal_chemistry import (
    RADIUS_RATIO_LIMITS,
    goldschmidt_tolerance_factor,
    radius_ratio_coordination,
)


def test_radius_ratio_limits_are_closed_form_touching_ratios():
    limits = [lim[0] for lim in RADIUS_RATIO_LIMITS]
    expected = [0.0, 2.0 / np.sqrt(3.0) - 1.0, np.sqrt(1.5) - 1.0, np.sqrt(2.0) - 1.0, np.sqrt(3.0) - 1.0]
    np.testing.assert_allclose(limits, expected)


@pytest.mark.parametrize(
    ("ratio", "cn"),
    [(0.10, 2), (0.20, 3), (0.30, 4), (0.50, 6), (0.90, 8), (1.2, 8)],
)
def test_radius_ratio_coordination_bands(ratio, cn):
    assert radius_ratio_coordination(ratio * 100.0, 100.0).coordination_number == cn


def test_radius_ratio_exact_boundary_goes_to_higher_coordination():
    assert radius_ratio_coordination(np.sqrt(2.0) - 1.0, 1.0).coordination_number == 6


def test_radius_ratio_rejects_nonpositive_radii():
    with pytest.raises(ValueError):
        radius_ratio_coordination(0.0, 1.0)


def test_tolerance_factor_is_one_for_ideal_geometry():
    r_b, r_x = 60.0, 140.0
    r_a = np.sqrt(2.0) * (r_b + r_x) - r_x
    assert goldschmidt_tolerance_factor(r_a, r_b, r_x) == pytest.approx(1.0)


def test_tolerance_factor_closed_form_and_vectorized():
    r_a = np.array([134.0, 144.0, 161.0])
    t = goldschmidt_tolerance_factor(r_a, 60.5, 140.0)
    np.testing.assert_allclose(t, (r_a + 140.0) / (np.sqrt(2.0) * 200.5))
    assert np.all(np.diff(t) > 0)
