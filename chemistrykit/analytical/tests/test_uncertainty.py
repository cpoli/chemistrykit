"""Tests for chemistrykit.analytical.systems.uncertainty against closed-form propagation results."""

import numpy as np
import pytest

from chemistrykit.analytical.systems.uncertainty import (
    propagate_power,
    propagate_product,
    propagate_sum,
    propagate_uncertainty,
)


def test_propagate_sum_two_equal_uncertainties():
    assert propagate_sum([0.02, 0.02]) == pytest.approx(0.02 * np.sqrt(2.0))


def test_propagate_sum_reduces_to_single_uncertainty_for_one_term():
    assert propagate_sum([0.05]) == pytest.approx(0.05)


def test_propagate_product_matches_relative_quadrature_formula():
    y = 2.0 * 3.0
    sigma = propagate_product([2.0, 3.0], [0.1, 0.2])
    expected_relative = np.sqrt((0.1 / 2.0) ** 2 + (0.2 / 3.0) ** 2)
    assert sigma == pytest.approx(y * expected_relative)


def test_propagate_power_doubles_relative_uncertainty_for_square():
    x, sigma_x = 3.0, 0.1
    sigma_y = propagate_power(x, sigma_x, 2.0)
    assert sigma_y / x**2 == pytest.approx(2.0 * sigma_x / x)


def test_propagate_power_with_exponent_one_is_identity():
    assert propagate_power(5.0, 0.3, 1.0) == pytest.approx(0.3)


def test_propagate_uncertainty_matches_closed_form_product():
    sigma_general = propagate_uncertainty(lambda a, b: a * b, [2.0, 3.0], [0.1, 0.2])
    sigma_closed_form = propagate_product([2.0, 3.0], [0.1, 0.2])
    assert sigma_general == pytest.approx(sigma_closed_form, rel=1e-4)


def test_propagate_uncertainty_matches_closed_form_sum():
    sigma_general = propagate_uncertainty(lambda a, b: a + b, [2.0, 3.0], [0.1, 0.2])
    sigma_closed_form = propagate_sum([0.1, 0.2])
    assert sigma_general == pytest.approx(sigma_closed_form, rel=1e-4)


def test_propagate_uncertainty_matches_closed_form_power():
    sigma_general = propagate_uncertainty(lambda x: x**3, [4.0], [0.2])
    sigma_closed_form = propagate_power(4.0, 0.2, 3.0)
    assert sigma_general == pytest.approx(sigma_closed_form, rel=1e-4)
