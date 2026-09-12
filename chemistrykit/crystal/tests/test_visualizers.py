"""Smoke tests for chemistrykit.crystal.visualizers: right return type, no crash."""

import matplotlib

matplotlib.use("Agg")

import matplotlib.axes
import numpy as np

from chemistrykit.crystal.systems.defects import schottky_defect_concentration
from chemistrykit.crystal.systems.madelung import madelung_constant_nacl
from chemistrykit.crystal.systems.xrd import powder_xrd_peaks
from chemistrykit.crystal.visualizers.crystal_plots import (
    plot_defect_concentration_vs_temperature,
    plot_madelung_convergence,
    plot_packing_fractions,
    plot_xrd_pattern,
)


def test_plot_packing_fractions_returns_axes():
    ax = plot_packing_fractions(["SC", "BCC", "FCC"], [0.52, 0.68, 0.74])
    assert isinstance(ax, matplotlib.axes.Axes)


def test_plot_madelung_convergence_returns_axes():
    n_shells = [2, 5, 10]
    values = [madelung_constant_nacl(n) for n in n_shells]
    ax = plot_madelung_convergence(n_shells, values, literature_value=1.747565)
    assert isinstance(ax, matplotlib.axes.Axes)


def test_plot_xrd_pattern_returns_axes():
    peaks = powder_xrd_peaks("FCC", a=408.6, wavelength=154.18, hkl_max=2)
    ax = plot_xrd_pattern(peaks)
    assert isinstance(ax, matplotlib.axes.Axes)


def test_plot_defect_concentration_vs_temperature_returns_axes():
    T = np.linspace(300.0, 1000.0, 20)
    n = schottky_defect_concentration(N=1e22, delta_h=2.0e-19, T=T)
    ax = plot_defect_concentration_vs_temperature(T, n, label="Schottky")
    assert isinstance(ax, matplotlib.axes.Axes)
