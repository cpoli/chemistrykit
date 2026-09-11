"""Smoke tests for chemistrykit.solutions.visualizers: right return type, no crash."""

import matplotlib

matplotlib.use("Agg")

import matplotlib.axes
import numpy as np

from chemistrykit.solutions.systems.titration import StrongAcidStrongBaseTitration, WeakAcidStrongBaseTitration
from chemistrykit.solutions.visualizers.solutions_plots import plot_activity_coefficients, plot_titration_curve


def test_plot_titration_curve_returns_axes():
    titration = StrongAcidStrongBaseTitration(Ca=0.1, Va=0.05, Cb=0.1)
    Vb = np.linspace(1e-6, 0.09, 500)
    ax = plot_titration_curve(titration, Vb)
    assert isinstance(ax, matplotlib.axes.Axes)


def test_plot_titration_curve_weak_acid_returns_axes():
    titration = WeakAcidStrongBaseTitration(Ca=0.1, Va=0.05, Ka=1.8e-5, Cb=0.1)
    Vb = np.linspace(1e-6, 0.09, 200)
    ax = plot_titration_curve(titration, Vb, mark_equivalence=False)
    assert isinstance(ax, matplotlib.axes.Axes)


def test_plot_activity_coefficients_returns_axes():
    I_range = np.linspace(1e-4, 0.1, 50)
    ax = plot_activity_coefficients(I_range, z_values=[1, 2, 3])
    assert isinstance(ax, matplotlib.axes.Axes)
    assert len(ax.lines) == 3


def test_plot_activity_coefficients_extended_returns_axes():
    I_range = np.linspace(1e-4, 0.1, 50)
    ax = plot_activity_coefficients(I_range, z_values=[1], extended=True)
    assert isinstance(ax, matplotlib.axes.Axes)
