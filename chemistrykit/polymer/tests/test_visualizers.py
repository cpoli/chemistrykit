"""Smoke tests for chemistrykit.polymer.visualizers: right return type, no crash."""

import matplotlib

matplotlib.use("Agg")

import matplotlib.axes
import numpy as np

from chemistrykit.polymer.systems.chain_growth import free_radical_network
from chemistrykit.polymer.systems.chain_statistics import IdealChain, RealChain
from chemistrykit.polymer.systems.molecular_weight_distribution import flory_schulz_number_fraction
from chemistrykit.polymer.systems.step_growth import degree_of_polymerization
from chemistrykit.polymer.visualizers.polymer_plots import (
    plot_carothers_curve,
    plot_chain_scaling,
    plot_free_radical_kinetics,
    plot_molecular_weight_distribution,
)


def test_plot_chain_scaling_returns_axes():
    models = {"ideal": IdealChain(), "good": RealChain.good_solvent(), "poor": RealChain.poor_solvent()}
    n_values = np.logspace(1, 4, 20)
    ax = plot_chain_scaling(models, n_values, b=0.5)
    assert isinstance(ax, matplotlib.axes.Axes)


def test_plot_molecular_weight_distribution_returns_axes():
    x = np.arange(1, 200)
    fraction = flory_schulz_number_fraction(x, p=0.95)
    ax = plot_molecular_weight_distribution(x, fraction, label="number fraction")
    assert isinstance(ax, matplotlib.axes.Axes)


def test_plot_carothers_curve_returns_axes():
    p = np.linspace(0.0, 0.99, 100)
    Xn = degree_of_polymerization(p)
    ax = plot_carothers_curve(p, Xn)
    assert isinstance(ax, matplotlib.axes.Axes)


def test_plot_free_radical_kinetics_returns_axes():
    net = free_radical_network(kd=1.0e-5, f=0.5, kp=1.0e3, kt=1.0e7, I0=0.01, M0=5.0)
    result = net.integrate((0.0, 10.0), dt=1e-2, method="rk4")
    ax = plot_free_radical_kinetics(result)
    assert isinstance(ax, matplotlib.axes.Axes)
