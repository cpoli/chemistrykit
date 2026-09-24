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


def test_bravais_lattices_total_fourteen_over_seven_systems():
    from chemistrykit.crystal.systems.crystal_systems import BRAVAIS_LATTICES

    assert len(BRAVAIS_LATTICES) == 7
    assert sum(len(c) for c in BRAVAIS_LATTICES.values()) == 14


@pytest.mark.parametrize(("centering", "n_points", "per_cell"), [("P", 8, 1), ("I", 9, 2), ("F", 14, 4)])
def test_cubic_lattice_points_counts(centering, n_points, per_cell):
    from chemistrykit.crystal.systems.crystal_systems import cubic_lattice_points

    assert len(cubic_lattice_points(centering)) == n_points
    # A large block contains (n+1)^3-ish boundary points; the per-cell density
    # is exact for the half-open block [0, n)^3.
    n = 4
    pts = cubic_lattice_points(centering, a=2.0, n_cells=n)
    interior = np.all(pts < 2.0 * n - 1e-9, axis=1)
    assert interior.sum() == per_cell * n**3


def test_cubic_lattice_points_nearest_neighbor_distances():
    from chemistrykit.crystal.systems.crystal_systems import cubic_lattice_points

    for centering, d_nn in (("P", 1.0), ("I", np.sqrt(3.0) / 2.0), ("F", 1.0 / np.sqrt(2.0))):
        pts = cubic_lattice_points(centering, a=1.0, n_cells=2)
        diffs = np.linalg.norm(pts[:, None, :] - pts[None, :, :], axis=-1)
        assert diffs[diffs > 1e-9].min() == pytest.approx(d_nn)


@pytest.mark.parametrize(
    ("intercepts", "expected"),
    [((1, 2, 3), (6, 3, 2)), ((1, 1, 1), (1, 1, 1)), ((0.5, 1, float("inf")), (2, 1, 0)), ((-1, 1, float("inf")), (-1, 1, 0)), ((2, 3, 6), (3, 2, 1))],
)
def test_miller_indices_from_intercepts(intercepts, expected):
    from chemistrykit.crystal.systems.crystal_systems import miller_indices_from_intercepts

    assert miller_indices_from_intercepts(*intercepts) == expected


def test_miller_indices_rejects_degenerate_faces():
    from chemistrykit.crystal.systems.crystal_systems import miller_indices_from_intercepts

    with pytest.raises(ValueError):
        miller_indices_from_intercepts(0, 1, 1)
    with pytest.raises(ValueError):
        miller_indices_from_intercepts(float("inf"), float("inf"), float("inf"))


def test_interplanar_angle_cubic_known_values():
    from chemistrykit.crystal.systems.crystal_systems import interplanar_angle_cubic

    assert interplanar_angle_cubic((1, 0, 0), (0, 0, 1)) == pytest.approx(90.0)
    assert interplanar_angle_cubic((1, 1, 1), (1, 1, -1)) == pytest.approx(np.degrees(np.arccos(1.0 / 3.0)))
    assert interplanar_angle_cubic((1, 0, 0), (1, 1, 0)) == pytest.approx(45.0)
    # Scale invariance: (222) is parallel to (111).
    assert interplanar_angle_cubic((1, 1, 1), (2, 2, 2)) == pytest.approx(0.0, abs=1e-6)
