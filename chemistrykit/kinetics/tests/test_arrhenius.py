"""Tests for chemistrykit.kinetics.systems.arrhenius."""

import numpy as np
import pytest

from chemistrykit.constants import R
from chemistrykit.kinetics.systems.arrhenius import arrhenius_rate_constant, fit_arrhenius


def test_arrhenius_rate_constant_reduces_to_A_when_Ea_is_zero():
    """With no activation energy, k = A regardless of temperature."""
    A = 1e6
    k = arrhenius_rate_constant(A=A, Ea=0.0, T=np.array([250.0, 300.0, 500.0]))
    np.testing.assert_allclose(k, A)


def test_arrhenius_rate_constant_increases_with_temperature():
    T = np.array([280.0, 300.0, 320.0])
    k = arrhenius_rate_constant(A=1e10, Ea=40e3, T=T)
    assert np.all(np.diff(k) > 0)


def test_fit_arrhenius_recovers_exact_parameters():
    """Fitting noiseless data generated from known (A, Ea) should recover them exactly."""
    A_true, Ea_true = 2.5e11, 55_000.0
    T = np.linspace(280.0, 360.0, 8)
    k = arrhenius_rate_constant(A=A_true, Ea=Ea_true, T=T)

    fit = fit_arrhenius(T, k)

    assert fit.Ea == pytest.approx(Ea_true, rel=1e-6)
    assert fit.A == pytest.approx(A_true, rel=1e-6)
    assert fit.r_squared == pytest.approx(1.0, abs=1e-9)


def test_fit_arrhenius_predict_matches_original_data():
    A_true, Ea_true = 8e12, 62_000.0
    T = np.linspace(290.0, 350.0, 6)
    k = arrhenius_rate_constant(A=A_true, Ea=Ea_true, T=T)
    fit = fit_arrhenius(T, k)
    np.testing.assert_allclose(fit.predict(T), k, rtol=1e-6)


def test_fit_arrhenius_default_gas_constant_matches_module_constant():
    T = np.linspace(280.0, 360.0, 5)
    k = arrhenius_rate_constant(A=1e12, Ea=50e3, T=T)
    fit = fit_arrhenius(T, k)
    assert fit.R_gas == R
