"""Tests for chemistrykit.surface.systems.dubinin against closed-form results."""

import numpy as np
import pytest

from chemistrykit.constants import R
from chemistrykit.surface.systems.dubinin import (
    DubininRadushkevichIsotherm,
    dubinin_radushkevich_loading,
    fit_dubinin_radushkevich,
    polanyi_potential,
)


def test_polanyi_potential_closed_form_and_zero_at_saturation():
    assert polanyi_potential(0.1, P0=1.0, T=77.0) == pytest.approx(R * 77.0 * np.log(10.0))
    assert polanyi_potential(2.0, P0=1.0, T=77.0) == 0.0


def test_dr_loading_is_W0_over_e_when_A_equals_E():
    T, E = 250.0, 8.0e3
    P = np.exp(-E / (R * T))
    assert dubinin_radushkevich_loading(P, W0=3.0, E=E, P0=1.0, T=T) == pytest.approx(3.0 / np.e)


def test_characteristic_curve_is_temperature_invariant():
    A = np.linspace(100.0, 20e3, 30)
    curves = []
    for T, P0 in [(250.0, 0.3), (300.0, 1.0), (350.0, 3.0)]:
        P = P0 * np.exp(-A / (R * T))
        curves.append(dubinin_radushkevich_loading(P, W0=0.5, E=10e3, P0=P0, T=T))
    np.testing.assert_allclose(curves[0], curves[1], rtol=1e-12)
    np.testing.assert_allclose(curves[1], curves[2], rtol=1e-12)


def test_fit_recovers_parameters_from_noisy_data():
    rng = np.random.default_rng(3)
    P = np.logspace(-5, -0.5, 12)
    W = dubinin_radushkevich_loading(P, W0=0.4, E=9e3, P0=1.0, T=77.0) * (1 + rng.normal(scale=0.002, size=P.shape))
    fit = fit_dubinin_radushkevich(P, W, P0=1.0, T=77.0)
    assert fit.W0 == pytest.approx(0.4, rel=1e-2)
    assert fit.E == pytest.approx(9e3, rel=1e-2)


def test_isotherm_saturates_at_W0():
    iso = DubininRadushkevichIsotherm(W0=0.5, E=12e3, P0=1.0, T=300.0)
    assert iso.loading(1.0) == pytest.approx(0.5)
    with pytest.raises(ValueError):
        DubininRadushkevichIsotherm(W0=-1.0, E=1.0, P0=1.0, T=1.0)
