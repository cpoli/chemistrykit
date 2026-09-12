"""Smoke tests for chemistrykit.analytical.visualizers: right return type, no crash."""

import matplotlib

matplotlib.use("Agg")

import matplotlib.axes
import numpy as np

from chemistrykit.analytical.systems.calibration import fit_calibration
from chemistrykit.analytical.systems.chromatography import (
    minimum_plate_height,
    optimum_flow_velocity,
    simulate_chromatogram,
    van_deemter_H,
)
from chemistrykit.analytical.systems.titration import RedoxTitration
from chemistrykit.analytical.visualizers.analytical_plots import (
    plot_calibration_curve,
    plot_chromatogram,
    plot_titration_curve,
    plot_van_deemter,
)


def test_plot_titration_curve_returns_axes():
    titration = RedoxTitration(E1_standard=0.771, n1=1, E2_standard=1.72, n2=1, C_analyte=0.10, V_analyte=0.050, C_titrant=0.10)
    V = np.linspace(1e-6, 0.09, 200)
    curve = titration.curve(V)
    ax = plot_titration_curve(curve.V, curve.response, V_equiv=titration.equivalence_volume())
    assert isinstance(ax, matplotlib.axes.Axes)


def test_plot_van_deemter_returns_axes():
    A, B, C = 1.0, 2.0, 0.05
    u = np.linspace(0.1, 10.0, 200)
    H = van_deemter_H(u, A, B, C)
    u_opt = optimum_flow_velocity(B, C)
    H_min = minimum_plate_height(A, B, C)
    ax = plot_van_deemter(u, H, u_opt=u_opt, H_min=H_min)
    assert isinstance(ax, matplotlib.axes.Axes)


def test_plot_chromatogram_returns_axes():
    t = np.linspace(8.0, 12.0, 500)
    chrom = simulate_chromatogram(t, centers=[9.5, 10.5], N=5000.0)
    ax = plot_chromatogram(t, chrom)
    assert isinstance(ax, matplotlib.axes.Axes)


def test_plot_calibration_curve_returns_axes():
    conc = np.array([0.0, 1.0, 2.0, 3.0, 4.0])
    signal = 2.5 * conc + 0.3
    cal = fit_calibration(conc, signal)
    ax = plot_calibration_curve(conc, signal, calibration=cal, lod=cal.lod())
    assert isinstance(ax, matplotlib.axes.Axes)
