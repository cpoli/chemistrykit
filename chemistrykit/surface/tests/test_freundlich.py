"""Tests for chemistrykit.surface.systems.freundlich against closed-form/known results."""

import numpy as np
import pytest

from chemistrykit.surface.systems.freundlich import FreundlichIsotherm, fit_freundlich, freundlich_loading


def test_freundlich_loading_equals_Kf_at_unit_pressure():
    assert freundlich_loading(Kf=3.7, n=2.2, P=1.0) == pytest.approx(3.7)


def test_freundlich_loading_is_linear_when_n_equals_one():
    P = np.linspace(0.1, 10.0, 20)
    q = freundlich_loading(Kf=2.0, n=1.0, P=P)
    np.testing.assert_allclose(q, 2.0 * P)


def test_freundlich_isotherm_matches_function():
    iso = FreundlichIsotherm(Kf=1.5, n=3.0)
    P = np.array([0.5, 1.0, 2.0])
    np.testing.assert_allclose(iso.loading(P), freundlich_loading(1.5, 3.0, P))


def test_freundlich_isotherm_rejects_nonpositive_parameters():
    with pytest.raises(ValueError):
        FreundlichIsotherm(Kf=0.0, n=2.0)
    with pytest.raises(ValueError):
        FreundlichIsotherm(Kf=1.0, n=-1.0)


def test_freundlich_isotherm_has_no_saturation_coverage():
    iso = FreundlichIsotherm(Kf=1.0, n=2.0)
    with pytest.raises(NotImplementedError):
        iso.fractional_coverage(1.0)


def test_fit_freundlich_recovers_known_parameters_from_exact_data():
    Kf_true, n_true = 4.0, 2.5
    P = np.array([0.1, 0.5, 1.0, 2.0, 5.0, 10.0])
    q = freundlich_loading(Kf_true, n_true, P)
    fit = fit_freundlich(P, q)
    assert fit.Kf == pytest.approx(Kf_true, rel=1e-6)
    assert fit.n == pytest.approx(n_true, rel=1e-6)
    assert fit.r_squared == pytest.approx(1.0, abs=1e-9)
