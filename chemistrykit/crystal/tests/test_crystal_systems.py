"""Tests for chemistrykit.crystal.systems.crystal_systems against closed-form/known results."""

import numpy as np
import pytest

from chemistrykit.crystal.systems.crystal_systems import classify_crystal_system, unit_cell_volume


@pytest.mark.parametrize(
    ("a", "b", "c", "alpha", "beta", "gamma", "expected"),
    [
        (5.0, 5.0, 5.0, 90.0, 90.0, 90.0, "cubic"),
        (5.0, 5.0, 7.0, 90.0, 90.0, 90.0, "tetragonal"),
        (5.0, 6.0, 7.0, 90.0, 90.0, 90.0, "orthorhombic"),
        (5.0, 5.0, 7.0, 90.0, 90.0, 120.0, "hexagonal"),
        (5.0, 5.0, 5.0, 80.0, 80.0, 80.0, "trigonal"),
        (5.0, 6.0, 7.0, 90.0, 95.0, 90.0, "monoclinic"),
        (5.0, 6.0, 7.0, 80.0, 85.0, 95.0, "triclinic"),
    ],
)
def test_classify_crystal_system_all_seven(a, b, c, alpha, beta, gamma, expected):
    assert classify_crystal_system(a, b, c, alpha, beta, gamma) == expected


def test_unit_cell_volume_right_angle_cell_is_abc():
    assert unit_cell_volume(2.0, 3.0, 4.0, 90.0, 90.0, 90.0) == pytest.approx(24.0)


def test_unit_cell_volume_hexagonal_matches_shortcut_formula():
    a, c = 2.46, 6.71
    general = unit_cell_volume(a, a, c, 90.0, 90.0, 120.0)
    shortcut = (np.sqrt(3.0) / 2.0) * a**2 * c
    assert general == pytest.approx(shortcut)


def test_unit_cell_volume_vectorized():
    a = np.array([1.0, 2.0, 3.0])
    result = unit_cell_volume(a, a, a, 90.0, 90.0, 90.0)
    np.testing.assert_allclose(result, a**3)
