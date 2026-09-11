"""Smoke tests for chemistrykit.kinetics.visualizers: right return type, no crash."""

import matplotlib

matplotlib.use("Agg")

import matplotlib.axes
import numpy as np

from chemistrykit.kinetics.systems.arrhenius import arrhenius_rate_constant, fit_arrhenius
from chemistrykit.kinetics.systems.enzyme import fit_lineweaver_burk, michaelis_menten_rate
from chemistrykit.kinetics.systems.networks import StoichiometricNetwork
from chemistrykit.kinetics.systems.oscillators import Brusselator
from chemistrykit.kinetics.visualizers.kinetics_plots import (
    plot_arrhenius,
    plot_concentration_vs_time,
    plot_lineweaver_burk,
    plot_phase_portrait,
)


def test_plot_concentration_vs_time_returns_axes():
    net = StoichiometricNetwork.consecutive(k1=1.0, k2=0.3)
    result = net.integrate((0.0, 5.0), dt=1e-2, method="rk4")
    ax = plot_concentration_vs_time(result)
    assert isinstance(ax, matplotlib.axes.Axes)


def test_plot_concentration_vs_time_subset_of_species():
    net = StoichiometricNetwork.consecutive(k1=1.0, k2=0.3)
    result = net.integrate((0.0, 5.0), dt=1e-2, method="rk4")
    ax = plot_concentration_vs_time(result, species=["A", "C"])
    assert isinstance(ax, matplotlib.axes.Axes)
    assert len(ax.lines) == 2


def test_plot_phase_portrait_returns_axes():
    system = Brusselator(A=1.0, B=3.0)
    result = system.integrate((0.0, 50.0), dt=1e-2, method="rk4")
    ax = plot_phase_portrait(result, "X", "Y")
    assert isinstance(ax, matplotlib.axes.Axes)


def test_plot_arrhenius_returns_axes():
    T = np.linspace(280.0, 360.0, 8)
    k = arrhenius_rate_constant(A=1e12, Ea=50e3, T=T)
    fit = fit_arrhenius(T, k)
    ax = plot_arrhenius(T, k, fit=fit)
    assert isinstance(ax, matplotlib.axes.Axes)


def test_plot_lineweaver_burk_returns_axes():
    S = np.array([0.5, 1.0, 2.0, 4.0, 8.0])
    v = michaelis_menten_rate(S, Vmax=10.0, Km=2.0)
    fit = fit_lineweaver_burk(S, v)
    ax = plot_lineweaver_burk(S, v, fit=fit)
    assert isinstance(ax, matplotlib.axes.Axes)
