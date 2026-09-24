"""Tests for the Kekule, ring-strain and dipole-moment modules of chemistrykit.structure."""

import numpy as np
import pytest

from chemistrykit.structure.core.base_system import angle_between
from chemistrykit.structure.systems.dipole import E_ANGSTROM_IN_DEBYE, bond_dipole_sum, dipole_moment
from chemistrykit.structure.systems.kekule import count_kekule_structures, kekule_structures
from chemistrykit.structure.systems.ring_strain import TETRAHEDRAL_ANGLE, baeyer_angle_strain, chair_cyclohexane_coordinates, planar_ring_angle


def _polyacene(n_rings):
    """Linear acene skeleton: 4n+2 carbons (benzene, naphthalene, anthracene, ...)."""
    top = list(range(0, 2 * n_rings + 1))
    bottom = list(range(2 * n_rings + 1, 4 * n_rings + 2))
    bonds = [(top[i], top[i + 1]) for i in range(len(top) - 1)]
    bonds += [(bottom[i], bottom[i + 1]) for i in range(len(bottom) - 1)]
    bonds += [(top[2 * k], bottom[2 * k]) for k in range(n_rings + 1)]
    return 4 * n_rings + 2, bonds


@pytest.mark.parametrize("n_rings", [1, 2, 3, 4, 5])
def test_acene_kekule_count_is_n_plus_one(n_rings):
    n, bonds = _polyacene(n_rings)
    assert count_kekule_structures(n, bonds) == n_rings + 1


def test_phenanthrene_has_five_kekule_structures():
    bonds = [(0, 1), (1, 2), (2, 3), (3, 4), (4, 5), (5, 0), (4, 6), (6, 7), (7, 8), (8, 9), (9, 3), (8, 10), (10, 11), (11, 12), (12, 13), (13, 9)]
    assert count_kekule_structures(14, bonds) == 5


def test_every_kekule_structure_is_a_perfect_matching():
    n, bonds = _polyacene(3)
    bondset = {tuple(sorted(b)) for b in bonds}
    for structure in kekule_structures(n, bonds):
        atoms = [a for pair in structure for a in pair]
        assert sorted(atoms) == list(range(n))
        assert set(structure) <= bondset


def test_odd_rings_and_polyene_chains():
    assert count_kekule_structures(5, [(i, (i + 1) % 5) for i in range(5)]) == 0
    assert count_kekule_structures(8, [(i, i + 1) for i in range(7)]) == 1


def test_baeyer_strain_closed_form_and_signs():
    for n in range(3, 12):
        assert baeyer_angle_strain(n) == pytest.approx(0.5 * (np.degrees(np.arccos(-1 / 3)) - 180.0 * (n - 2) / n))
    assert baeyer_angle_strain(3) > baeyer_angle_strain(4) > baeyer_angle_strain(5) > 0 > baeyer_angle_strain(6)
    assert planar_ring_angle(4) == 90.0
    with pytest.raises(ValueError):
        planar_ring_angle(2)


def test_chair_cyclohexane_is_strain_free():
    x = chair_cyclohexane_coordinates(1.53)
    for k in range(6):
        a, b, c = x[k - 1], x[k], x[(k + 1) % 6]
        assert np.linalg.norm(b - a) == pytest.approx(1.53)
        assert angle_between(a - b, c - b) == pytest.approx(TETRAHEDRAL_ANGLE)


def test_dipole_unit_conversion_and_origin_independence():
    assert E_ANGSTROM_IN_DEBYE == pytest.approx(4.80320, rel=1e-5)
    rng = np.random.default_rng(0)
    x = rng.normal(size=(5, 3))
    q = rng.normal(size=5)
    q -= q.mean()
    assert np.allclose(dipole_moment(x, q), dipole_moment(x + [3.0, -1.0, 2.0], q))


def test_symmetric_bond_dipoles_cancel_and_water_matches_closed_form():
    tetra = np.array([[1, 1, 1], [1, -1, -1], [-1, 1, -1], [-1, -1, 1]], dtype=float)
    assert np.allclose(bond_dipole_sum(tetra, [1.3] * 4), 0.0)
    half = np.radians(52.25)
    u = [[np.sin(half), 0, np.cos(half)], [-np.sin(half), 0, np.cos(half)]]
    assert np.linalg.norm(bond_dipole_sum(u, [1.5, 1.5])) == pytest.approx(3.0 * np.cos(half))
    with pytest.raises(ValueError):
        bond_dipole_sum([[0, 0, 0]], [1.0])
