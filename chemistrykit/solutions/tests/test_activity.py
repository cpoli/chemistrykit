"""Tests for chemistrykit.solutions.systems.activity against closed-form results."""

import numpy as np
import pytest

from chemistrykit.solutions.systems.activity import (
    activity_coefficient_debye_huckel_extended,
    activity_coefficient_debye_huckel_limiting,
    ionic_strength,
)


def test_ionic_strength_of_1_1_electrolyte_equals_concentration():
    assert ionic_strength([0.05, 0.05], [1, -1]) == pytest.approx(0.05)


def test_ionic_strength_of_2_1_electrolyte():
    """CaCl2: I = 0.5*(c*4 + 2c*1) = 3c."""
    c = 0.02
    assert ionic_strength([c, 2 * c], [2, -1]) == pytest.approx(3.0 * c)


def test_debye_huckel_limiting_law_gamma_is_one_at_zero_ionic_strength():
    assert activity_coefficient_debye_huckel_limiting(z=1, I=0.0) == pytest.approx(1.0)


def test_debye_huckel_limiting_law_gamma_decreases_with_ionic_strength():
    gammas = [activity_coefficient_debye_huckel_limiting(z=1, I=I) for I in (0.0, 0.001, 0.01, 0.1)]
    assert all(g1 > g2 for g1, g2 in zip(gammas[:-1], gammas[1:], strict=True))


def test_debye_huckel_limiting_law_higher_charge_deviates_more():
    gamma_1 = activity_coefficient_debye_huckel_limiting(z=1, I=0.01)
    gamma_2 = activity_coefficient_debye_huckel_limiting(z=2, I=0.01)
    gamma_3 = activity_coefficient_debye_huckel_limiting(z=3, I=0.01)
    assert gamma_1 > gamma_2 > gamma_3


def test_debye_huckel_limiting_law_matches_manual_formula():
    A, z, I = 0.509, 2, 0.005
    expected = 10.0 ** (-A * z**2 * np.sqrt(I))
    assert activity_coefficient_debye_huckel_limiting(z, I, A=A) == pytest.approx(expected)


def test_debye_huckel_extended_reduces_to_limiting_at_low_ionic_strength():
    I = 1e-8
    extended = activity_coefficient_debye_huckel_extended(z=1, I=I)
    limiting = activity_coefficient_debye_huckel_limiting(z=1, I=I)
    assert extended == pytest.approx(limiting, rel=1e-3)


def test_debye_huckel_extended_predicts_higher_gamma_than_limiting_at_moderate_ionic_strength():
    extended = activity_coefficient_debye_huckel_extended(z=1, I=0.1)
    limiting = activity_coefficient_debye_huckel_limiting(z=1, I=0.1)
    assert extended > limiting


def test_debye_huckel_extended_matches_manual_formula():
    A, z, I, Ba = 0.509, 1, 0.05, 1.5
    expected = 10.0 ** (-A * z**2 * np.sqrt(I) / (1.0 + Ba * np.sqrt(I)))
    assert activity_coefficient_debye_huckel_extended(z, I, A=A, Ba=Ba) == pytest.approx(expected)


def test_davies_matches_manual_formula():
    from chemistrykit.solutions.systems.activity import activity_coefficient_davies

    A, z, I, b = 0.509, 2, 0.1, 0.3
    expected = 10.0 ** (-A * z**2 * (np.sqrt(I) / (1 + np.sqrt(I)) - b * I))
    assert activity_coefficient_davies(z, I, A=A, b=b) == pytest.approx(expected)


def test_davies_reduces_to_limiting_law_at_low_ionic_strength():
    from chemistrykit.solutions.systems.activity import activity_coefficient_davies

    I = 1e-8
    assert activity_coefficient_davies(z=2, I=I) == pytest.approx(activity_coefficient_debye_huckel_limiting(z=2, I=I), rel=1e-6)


def test_davies_log_gamma_has_minimum_where_derivative_vanishes():
    """d/dI [sqrt(I)/(1+sqrt(I)) - bI] = 0  <=>  1/(2 sqrt(I)(1+sqrt(I))^2) = b."""
    from scipy.optimize import brentq

    from chemistrykit.solutions.systems.activity import activity_coefficient_davies

    b = 0.3
    I_min = brentq(lambda I: 1.0 / (2 * np.sqrt(I) * (1 + np.sqrt(I)) ** 2) - b, 1e-3, 5.0)
    I_grid = np.linspace(0.01, 2.0, 20001)
    gammas = np.array([activity_coefficient_davies(1, I, b=b) for I in I_grid])
    assert I_grid[np.argmin(gammas)] == pytest.approx(I_min, abs=1e-3)
