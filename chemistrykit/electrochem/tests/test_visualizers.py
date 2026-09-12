"""Smoke tests for chemistrykit.electrochem.visualizers: right return type, no crash."""

import matplotlib

matplotlib.use("Agg")

import matplotlib.axes

from chemistrykit.electrochem.systems.battery import ConstantCurrentBattery
from chemistrykit.electrochem.systems.butler_volmer import butler_volmer_current_density, fit_tafel_plot
from chemistrykit.electrochem.visualizers.electrochem_plots import (
    plot_discharge_curve,
    plot_nernst_concentration_dependence,
    plot_tafel,
)


def test_plot_tafel_returns_axes():
    ax = plot_tafel(i0=1e-6, alpha=0.5, n=1)
    assert isinstance(ax, matplotlib.axes.Axes)


def test_plot_tafel_with_fit_overlay_returns_axes():
    import numpy as np

    i0, alpha, n = 1e-6, 0.5, 1
    eta = np.linspace(0.2, 0.4, 10)
    i = butler_volmer_current_density(i0, eta, alpha=alpha, n=n)
    fit = fit_tafel_plot(eta, i)
    ax = plot_tafel(i0=i0, alpha=alpha, n=n, fit=fit)
    assert isinstance(ax, matplotlib.axes.Axes)


def test_plot_nernst_concentration_dependence_returns_axes():
    ax = plot_nernst_concentration_dependence(E_standard=0.34, n=2)
    assert isinstance(ax, matplotlib.axes.Axes)


def test_plot_discharge_curve_returns_axes():
    battery = ConstantCurrentBattery(capacity_peukert=5.0, current=1.0, v_nominal=3.7)
    ax = plot_discharge_curve(battery)
    assert isinstance(ax, matplotlib.axes.Axes)
