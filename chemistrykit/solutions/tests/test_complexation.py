"""Tests for chemistrykit.solutions.systems.complexation."""

import numpy as np
import pytest

from chemistrykit.solutions.systems.complexation import average_ligand_number, complex_fractions, cumulative_formation_constants, solve_complexation

BETA_AG_NH3 = [10**3.31, 10**7.23]


def test_fractions_sum_to_one_and_follow_beta():
    L = np.logspace(-6, 0, 13)
    alpha = complex_fractions(L, BETA_AG_NH3)
    np.testing.assert_allclose(alpha.sum(axis=0), 1.0)
    np.testing.assert_allclose(alpha[2] / alpha[0], BETA_AG_NH3[1] * L**2)


def test_single_step_formation_function():
    L = np.logspace(-5, -1, 9)
    np.testing.assert_allclose(average_ligand_number(L, [1e3]), 1e3 * L / (1 + 1e3 * L))


def test_cumulative_constants():
    np.testing.assert_allclose(cumulative_formation_constants([10**3.31, 10**3.92]), [10**3.31, 10**7.23])


def test_mass_balances_hold():
    eq = solve_complexation(M_total=0.01, L_total=0.015, beta=BETA_AG_NH3)
    assert eq.species.sum() == pytest.approx(0.01)
    bound = np.dot(np.arange(3), eq.species)
    assert eq.free_ligand + bound == pytest.approx(0.015, rel=1e-9)
    assert eq.species[1] / (eq.free_metal * eq.free_ligand) == pytest.approx(BETA_AG_NH3[0], rel=1e-9)


def test_invalid_totals_raise():
    with pytest.raises(ValueError):
        solve_complexation(1e-3, 0.0, BETA_AG_NH3)
