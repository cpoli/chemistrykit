"""Tests for chemistrykit.surface.systems.langmuir against closed-form/known results."""

import numpy as np
import pytest

from chemistrykit.surface.systems.langmuir import LangmuirIsotherm, fit_langmuir, langmuir_coverage


def test_langmuir_coverage_is_half_at_half_saturation_pressure():
    K = 3.5
    assert langmuir_coverage(K, P=1.0 / K) == pytest.approx(0.5)


def test_langmuir_coverage_zero_at_zero_pressure():
    assert langmuir_coverage(K=2.0, P=0.0) == pytest.approx(0.0)


def test_langmuir_coverage_saturates_at_high_pressure():
    assert langmuir_coverage(K=2.0, P=1.0e8) == pytest.approx(1.0, abs=1e-6)


def test_langmuir_coverage_is_monotonic_increasing():
    P = np.linspace(0.0, 100.0, 50)
    theta = langmuir_coverage(K=1.5, P=P)
    assert np.all(np.diff(theta) > 0)


def test_langmuir_isotherm_loading_equals_qmax_times_coverage():
    iso = LangmuirIsotherm(K=2.0, qmax=10.0)
    P = 3.0
    assert iso.loading(P) == pytest.approx(iso.qmax * langmuir_coverage(iso.K, P))


def test_langmuir_isotherm_half_saturation_pressure_is_reciprocal_of_K():
    iso = LangmuirIsotherm(K=4.0, qmax=1.0)
    assert iso.half_saturation_pressure() == pytest.approx(0.25)
    assert iso.fractional_coverage(iso.half_saturation_pressure()) == pytest.approx(0.5)


def test_langmuir_isotherm_rejects_nonpositive_parameters():
    with pytest.raises(ValueError):
        LangmuirIsotherm(K=-1.0, qmax=1.0)
    with pytest.raises(ValueError):
        LangmuirIsotherm(K=1.0, qmax=0.0)


def test_fit_langmuir_recovers_known_parameters_from_exact_data():
    K_true, qmax_true = 3.0, 8.0
    P = np.array([0.1, 0.5, 1.0, 2.0, 5.0, 10.0, 20.0])
    q = LangmuirIsotherm(K=K_true, qmax=qmax_true).loading(P)
    fit = fit_langmuir(P, q)
    assert fit.K == pytest.approx(K_true, rel=1e-6)
    assert fit.qmax == pytest.approx(qmax_true, rel=1e-6)
    assert fit.r_squared == pytest.approx(1.0, abs=1e-9)


def test_langmuir_fit_predict_matches_to_isotherm():
    P = np.array([0.1, 0.5, 1.0, 2.0, 5.0])
    q = LangmuirIsotherm(K=2.0, qmax=6.0).loading(P)
    fit = fit_langmuir(P, q)
    np.testing.assert_allclose(fit.predict(P), fit.to_isotherm().loading(P))
