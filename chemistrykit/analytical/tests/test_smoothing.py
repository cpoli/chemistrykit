"""Tests for chemistrykit.analytical.systems.smoothing: Savitzky-Golay filter."""

import numpy as np
import pytest
from scipy.signal import savgol_filter

from chemistrykit.analytical.systems.smoothing import savitzky_golay, savitzky_golay_coefficients


def test_classic_table_coefficients():
    assert savitzky_golay_coefficients(5, 2) * 35 == pytest.approx([-3, 12, 17, 12, -3])
    assert savitzky_golay_coefficients(7, 2) * 21 == pytest.approx([-2, 3, 6, 7, 6, 3, -2])
    # first-derivative, 5-point quadratic: (-2, -1, 0, 1, 2)/10
    assert savitzky_golay_coefficients(5, 2, deriv=1) * 10 == pytest.approx([-2, -1, 0, 1, 2])


def test_smoothing_weights_sum_to_one():
    assert savitzky_golay_coefficients(11, 4).sum() == pytest.approx(1.0)


@pytest.mark.parametrize(("window", "order", "deriv"), [(5, 2, 0), (11, 3, 0), (9, 4, 1), (15, 2, 2)])
def test_matches_scipy_savgol_filter_interp_mode(window, order, deriv):
    rng = np.random.default_rng(3)
    y = np.sin(np.linspace(0, 6, 120)) + rng.normal(scale=0.1, size=120)
    ours = savitzky_golay(y, window, order, deriv=deriv, delta=0.05)
    ref = savgol_filter(y, window, order, deriv=deriv, delta=0.05, mode="interp")
    assert ours == pytest.approx(ref, abs=1e-9)


def test_derivative_of_polynomial_is_exact():
    x = np.linspace(0.0, 2.0, 41)
    y = x**3 - 2 * x
    dy = savitzky_golay(y, 7, 3, deriv=1, delta=x[1] - x[0])
    assert dy == pytest.approx(3 * x**2 - 2, abs=1e-9)


def test_smoothing_reduces_noise():
    rng = np.random.default_rng(4)
    x = np.linspace(0, 1, 400)
    clean = np.exp(-((x - 0.5) ** 2) / 0.005)
    noisy = clean + rng.normal(scale=0.05, size=x.size)
    smoothed = savitzky_golay(noisy, 21, 3)
    assert np.std(smoothed - clean) < 0.5 * np.std(noisy - clean)
