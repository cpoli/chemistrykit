"""Tests for chemistrykit.photochem.systems.fluorescence."""

import numpy as np
import pytest

from chemistrykit.constants import K_B
from chemistrykit.photochem.systems.fluorescence import perrin_anisotropy, rotational_correlation_time, stokes_shift


def test_stokes_shift_closed_form_and_sign():
    assert stokes_shift(400.0, 500.0) == pytest.approx(5000.0)
    assert stokes_shift(500.0, 500.0) == pytest.approx(0.0)
    assert np.all(stokes_shift(np.array([350.0, 450.0]), np.array([420.0, 520.0])) > 0)


def test_perrin_limits_and_linear_perrin_plot():
    assert perrin_anisotropy(0.4, tau=1e-12, theta=1.0) == pytest.approx(0.4)
    assert perrin_anisotropy(0.4, tau=4.0, theta=4.0) == pytest.approx(0.2)
    # 1/r is linear in T/eta (Perrin plot), with intercept 1/r0
    T = np.linspace(280.0, 340.0, 5)
    eta, V, tau = 1e-3, 1e-27, 4e-9
    theta = rotational_correlation_time(eta, V, T)
    inv_r = 1.0 / perrin_anisotropy(0.4, tau, theta)
    slope, intercept = np.polyfit(T / eta, inv_r, 1)
    assert intercept == pytest.approx(2.5)
    assert slope == pytest.approx(2.5 * tau * K_B / V)
