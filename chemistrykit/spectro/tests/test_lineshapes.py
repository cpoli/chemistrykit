"""Tests for chemistrykit.spectro.utils.lineshapes against closed-form/known results."""

import numpy as np
import pytest

from chemistrykit.spectro.utils.lineshapes import broaden_stick_spectrum, gaussian, lorentzian, voigt


def test_gaussian_normalized():
    x = np.linspace(-30, 30, 300001)
    assert np.trapezoid(gaussian(x, 0.0, fwhm=2.0), x) == pytest.approx(1.0, abs=1e-3)


def test_gaussian_symmetric_about_center():
    assert gaussian(3.0, 0.0, fwhm=1.0) == pytest.approx(gaussian(-3.0, 0.0, fwhm=1.0))


def test_gaussian_half_max_at_half_fwhm():
    peak = gaussian(0.0, 0.0, fwhm=4.0)
    half = gaussian(2.0, 0.0, fwhm=4.0)
    assert half == pytest.approx(peak / 2.0, rel=1e-9)


def test_lorentzian_normalized():
    x = np.linspace(-3000, 3000, 3_000_001)
    assert np.trapezoid(lorentzian(x, 0.0, fwhm=2.0), x) == pytest.approx(1.0, abs=1e-3)


def test_lorentzian_half_max_at_half_fwhm():
    peak = lorentzian(0.0, 0.0, fwhm=3.0)
    half = lorentzian(1.5, 0.0, fwhm=3.0)
    assert half == pytest.approx(peak / 2.0, rel=1e-9)


def test_voigt_normalized():
    x = np.linspace(-1000, 1000, 400001)
    integral = np.trapezoid(voigt(x, 0.0, fwhm_gaussian=1.0, fwhm_lorentzian=1.0), x)
    assert integral == pytest.approx(1.0, abs=5e-3)


def test_voigt_reduces_to_gaussian_when_lorentzian_width_vanishes():
    x = np.linspace(-20, 20, 4001)
    voigt_values = voigt(x, 0.0, fwhm_gaussian=2.0, fwhm_lorentzian=1e-8)
    gaussian_values = gaussian(x, 0.0, fwhm=2.0)
    assert np.allclose(voigt_values, gaussian_values, atol=1e-3)


def test_voigt_reduces_to_lorentzian_when_gaussian_width_vanishes():
    x = np.linspace(-50, 50, 200001)
    voigt_values = voigt(x, 0.0, fwhm_gaussian=1e-8, fwhm_lorentzian=2.0)
    lorentzian_values = lorentzian(x, 0.0, fwhm=2.0)
    assert np.allclose(voigt_values, lorentzian_values, atol=1e-3)


def test_broaden_stick_spectrum_single_stick_matches_scaled_lineshape():
    x = np.linspace(-10, 10, 501)
    spectrum = broaden_stick_spectrum([1.0], [3.0], x, shape="gaussian", fwhm=1.0)
    assert np.allclose(spectrum, 3.0 * gaussian(x, 1.0, 1.0))


def test_broaden_stick_spectrum_is_additive_over_sticks():
    x = np.linspace(-10, 10, 501)
    combined = broaden_stick_spectrum([-2.0, 2.0], [1.0, 1.0], x, shape="lorentzian", fwhm=1.0)
    separate = lorentzian(x, -2.0, 1.0) + lorentzian(x, 2.0, 1.0)
    assert np.allclose(combined, separate)


def test_broaden_stick_spectrum_rejects_unknown_shape():
    with pytest.raises(ValueError):
        broaden_stick_spectrum([0.0], [1.0], np.array([0.0]), shape="triangular")


def test_broaden_stick_spectrum_voigt_requires_lorentzian_fwhm():
    with pytest.raises(ValueError):
        broaden_stick_spectrum([0.0], [1.0], np.array([0.0]), shape="voigt", fwhm=1.0)
