"""Smoke tests for chemistrykit.spectro.visualizers: right return type, no crash."""

import matplotlib

matplotlib.use("Agg")

import matplotlib.axes
import numpy as np

from chemistrykit.spectro.core.base_system import Spectrum
from chemistrykit.spectro.visualizers.spectro_plots import plot_beer_lambert_deviation, plot_broadened_spectrum, plot_stick_spectrum


def _spectrum():
    return Spectrum(positions=np.array([100.0, 200.0, 300.0]), intensities=np.array([1.0, 0.6, 0.2]))


def test_plot_stick_spectrum_returns_axes():
    ax = plot_stick_spectrum(_spectrum())
    assert isinstance(ax, matplotlib.axes.Axes)


def test_plot_broadened_spectrum_gaussian_returns_axes():
    x = np.linspace(0, 400, 500)
    ax = plot_broadened_spectrum(_spectrum(), x, shape="gaussian", fwhm=5.0)
    assert isinstance(ax, matplotlib.axes.Axes)
    assert len(ax.lines) == 1


def test_plot_broadened_spectrum_voigt_returns_axes():
    x = np.linspace(0, 400, 500)
    ax = plot_broadened_spectrum(_spectrum(), x, shape="voigt", fwhm=5.0, fwhm_lorentzian=3.0)
    assert isinstance(ax, matplotlib.axes.Axes)


def test_plot_beer_lambert_deviation_returns_axes():
    concentrations = np.linspace(1e-6, 1e-3, 50)
    ax = plot_beer_lambert_deviation(molar_absorptivity=5000.0, path_length=1.0, concentrations=concentrations, stray_light_fraction=0.001)
    assert isinstance(ax, matplotlib.axes.Axes)
    assert len(ax.lines) == 2
