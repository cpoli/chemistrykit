"""Tests for chemistrykit.structure.systems.bonding against closed-form/known results."""

import numpy as np
import pytest

from chemistrykit.quantum.systems.huckel import HuckelSystem
from chemistrykit.structure.systems.bonding import bond_length_from_order, bond_order_from_length, coulson_pi_bond_order


def test_bond_order_from_length_is_one_at_reference_length():
    assert bond_order_from_length(1.54, 1.54) == pytest.approx(1.0)


def test_bond_order_from_length_and_bond_length_from_order_are_inverses():
    D1, order_in = 1.54, 1.8
    length = bond_length_from_order(D1, order_in)
    order_out = bond_order_from_length(D1, length)
    assert order_out == pytest.approx(order_in)


def test_higher_order_gives_shorter_bond():
    lengths = [bond_length_from_order(1.54, n) for n in (1, 2, 3)]
    assert lengths[0] > lengths[1] > lengths[2]


def test_bond_length_from_order_rejects_nonpositive_order():
    with pytest.raises(ValueError):
        bond_length_from_order(1.54, 0.0)


def test_benzene_bond_order_between_single_and_double():
    n = bond_order_from_length(single_bond_length=1.54, observed_length=1.397)
    assert 1.0 < n < 2.0


def test_coulson_bond_order_ethene_is_one():
    ethene = HuckelSystem(n_atoms=2, bonds=[(0, 1)])
    result = ethene.solve()
    order = np.argsort(result.energies)
    occupations = np.zeros(2)
    occupations[order[0]] = 2.0
    assert coulson_pi_bond_order(result.coefficients, occupations, 0, 1) == pytest.approx(1.0)


def test_coulson_bond_order_benzene_is_two_thirds():
    benzene = HuckelSystem.cyclic_polyene(n_atoms=6)
    result = benzene.solve()
    order = np.argsort(result.energies)
    occupations = np.zeros(6)
    occupations[order[:3]] = 2.0  # 6 pi electrons fill the 3 lowest (bonding) MOs
    for i in range(6):
        j = (i + 1) % 6
        assert coulson_pi_bond_order(result.coefficients, occupations, i, j) == pytest.approx(2.0 / 3.0, abs=1e-6)


def test_coulson_bond_order_zero_for_unbonded_atoms_in_polyene():
    """Butadiene's 1,4 (non-adjacent) carbons have a much weaker pi bond order than the 1,2 (bonded) pair."""
    butadiene = HuckelSystem.linear_polyene(n_atoms=4)
    result = butadiene.solve()
    order = np.argsort(result.energies)
    occupations = np.zeros(4)
    occupations[order[:2]] = 2.0
    p12 = coulson_pi_bond_order(result.coefficients, occupations, 0, 1)
    p14 = coulson_pi_bond_order(result.coefficients, occupations, 0, 3)
    assert abs(p12) > abs(p14)
