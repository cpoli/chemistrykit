"""Tests for chemistrykit.structure.core.base_system.Molecule against closed-form geometry."""

import numpy as np
import pytest

from chemistrykit.structure.core.base_system import Molecule, angle_between, unit_vector


def _water():
    angle = np.radians(104.5)
    r = 0.958
    return Molecule(
        symbols=["O", "H", "H"],
        coordinates=[
            [0.0, 0.0, 0.0],
            [r * np.sin(angle / 2), 0.0, r * np.cos(angle / 2)],
            [-r * np.sin(angle / 2), 0.0, r * np.cos(angle / 2)],
        ],
        bonds=[(0, 1), (0, 2)],
    )


def test_bond_length_matches_construction():
    water = _water()
    assert water.bond_length(0, 1) == pytest.approx(0.958)
    assert water.bond_length(0, 2) == pytest.approx(0.958)


def test_bond_angle_matches_construction():
    water = _water()
    assert water.bond_angle(1, 0, 2) == pytest.approx(104.5)


def test_n_atoms():
    assert _water().n_atoms == 3


def test_rejects_mismatched_coordinate_shape():
    with pytest.raises(ValueError):
        Molecule(symbols=["H", "H"], coordinates=[[0, 0, 0]])


def test_rejects_bond_index_out_of_range():
    with pytest.raises(ValueError):
        Molecule(symbols=["H", "H"], coordinates=[[0, 0, 0], [1, 0, 0]], bonds=[(0, 5)])


def test_adjacency_is_symmetric():
    water = _water()
    adjacency = water.adjacency()
    assert adjacency[0] == [1, 2]
    assert adjacency[1] == [0]
    assert adjacency[2] == [0]


def test_center_of_mass_is_between_atoms_weighted_toward_heavier():
    water = _water()
    com = water.center_of_mass()
    centroid = water.centroid()
    # Oxygen (at the origin) is much heavier than hydrogen, so the mass-weighted
    # center of mass sits closer to the origin (z=0) than the unweighted centroid.
    assert abs(com[2]) < abs(centroid[2])


def test_is_linear_true_for_collinear_atoms():
    co2 = Molecule(symbols=["O", "C", "O"], coordinates=[[0, 0, -1.16], [0, 0, 0], [0, 0, 1.16]])
    assert co2.is_linear() is True


def test_is_linear_false_for_bent_molecule():
    assert _water().is_linear() is False


def test_unit_vector_normalizes():
    v = unit_vector([3.0, 4.0, 0.0])
    assert np.linalg.norm(v) == pytest.approx(1.0)


def test_unit_vector_rejects_zero_vector():
    with pytest.raises(ValueError):
        unit_vector([0.0, 0.0, 0.0])


def test_angle_between_orthogonal_vectors():
    assert angle_between([1, 0, 0], [0, 1, 0]) == pytest.approx(90.0)


def test_angle_between_parallel_vectors_is_zero():
    assert angle_between([1, 0, 0], [2, 0, 0]) == pytest.approx(0.0)
