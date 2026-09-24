"""Tests for chemistrykit.surface.systems.temkin against closed-form results."""

import numpy as np
import pytest
from scipy.integrate import quad

from chemistrykit.constants import R
from chemistrykit.surface.systems.langmuir import langmuir_coverage
from chemistrykit.surface.systems.temkin import TemkinIsotherm, fit_temkin, temkin_loading, uniform_energy_coverage


def test_uniform_energy_coverage_matches_numerical_average_of_langmuir():
    K_max, f = 50.0, 6.0
    for P in [0.001, 0.1, 1.0, 30.0]:
        numeric = quad(lambda u: langmuir_coverage(K_max * np.exp(-u), P), 0.0, f)[0] / f
        assert uniform_energy_coverage(K_max, f, P) == pytest.approx(numeric, rel=1e-8)


def test_uniform_energy_coverage_reduces_to_langmuir_for_zero_spread():
    assert uniform_energy_coverage(2.0, 1e-9, 0.7) == pytest.approx(langmuir_coverage(2.0, 0.7), rel=1e-6)


def test_uniform_energy_coverage_mid_range_is_temkin_logarithm():
    K_max, f = 1e8, 25.0
    P = np.logspace(-6, -3, 10)
    np.testing.assert_allclose(uniform_energy_coverage(K_max, f, P), np.log(K_max * P) / f, atol=1e-3)


def test_temkin_loading_closed_form():
    assert temkin_loading(2.0, 1e3, 5.0, 300.0) == pytest.approx(R * 300.0 / 1e3 * np.log(10.0))
    assert TemkinIsotherm(2.0, 1e3, 300.0).loading(0.5) == pytest.approx(0.0, abs=1e-15)
    with pytest.raises(NotImplementedError):
        TemkinIsotherm(2.0, 1e3, 300.0).fractional_coverage(1.0)


def test_fit_temkin_recovers_parameters_from_noisy_data():
    rng = np.random.default_rng(5)
    P = np.logspace(-0.5, 1.5, 10)
    q = temkin_loading(3.0, 800.0, P, 300.0) + rng.normal(scale=1e-3, size=P.shape)
    fit = fit_temkin(P, q, T=300.0)
    assert fit.A_T == pytest.approx(3.0, rel=1e-2)
    assert fit.b_T == pytest.approx(800.0, rel=1e-2)
