"""Smoke tests for chemistrykit.surface.visualizers: right return type, no crash."""

import matplotlib

matplotlib.use("Agg")

import matplotlib.axes
import numpy as np

from chemistrykit.surface.systems.langmuir import LangmuirIsotherm
from chemistrykit.surface.systems.langmuir_hinshelwood import lh_rate_single_site
from chemistrykit.surface.visualizers.surface_plots import plot_isotherm, plot_lh_rate_vs_pressure, plot_linearization


def test_plot_isotherm_returns_axes():
    iso = LangmuirIsotherm(K=2.0, qmax=5.0)
    ax = plot_isotherm(iso, P_max=10.0)
    assert isinstance(ax, matplotlib.axes.Axes)


def test_plot_isotherm_with_data_overlay_returns_axes():
    iso = LangmuirIsotherm(K=2.0, qmax=5.0)
    P_data = np.array([0.5, 1.0, 2.0])
    q_data = iso.loading(P_data)
    ax = plot_isotherm(iso, P_max=10.0, P_data=P_data, q_data=q_data)
    assert isinstance(ax, matplotlib.axes.Axes)


def test_plot_linearization_returns_axes():
    x = np.array([1.0, 2.0, 3.0])
    y = np.array([2.0, 4.0, 6.0])
    ax = plot_linearization(x, y, fit=(2.0, 0.0))
    assert isinstance(ax, matplotlib.axes.Axes)


def test_plot_lh_rate_vs_pressure_returns_axes():
    P_A = np.linspace(0.01, 10.0, 50)
    rate = lh_rate_single_site(k=4.0, K_A=2.0, P_A=P_A)
    ax = plot_lh_rate_vs_pressure(P_A, rate)
    assert isinstance(ax, matplotlib.axes.Axes)
