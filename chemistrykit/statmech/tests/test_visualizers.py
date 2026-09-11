"""Smoke tests for chemistrykit.statmech.visualizers: right return type, no crash."""

import matplotlib

matplotlib.use("Agg")

import matplotlib.axes
import numpy as np

from chemistrykit.statmech.systems.lattice_gas import LatticeGasAdsorption
from chemistrykit.statmech.systems.maxwell_boltzmann import MaxwellBoltzmannSpeedDistribution
from chemistrykit.statmech.systems.partition_functions import VibrationalPartitionFunctionHarmonic
from chemistrykit.statmech.visualizers.statmech_plots import (
    plot_adsorption_isotherm,
    plot_heat_capacity_vs_temperature,
    plot_maxwell_boltzmann,
)


def test_plot_heat_capacity_vs_temperature_returns_axes():
    q = VibrationalPartitionFunctionHarmonic(frequency=8.7e13)
    T = np.linspace(50.0, 5000.0, 50)
    ax = plot_heat_capacity_vs_temperature(q, T)
    assert isinstance(ax, matplotlib.axes.Axes)


def test_plot_maxwell_boltzmann_returns_axes():
    dist = MaxwellBoltzmannSpeedDistribution(mass=6.63e-26, temperature=298.15)
    ax = plot_maxwell_boltzmann(dist)
    assert isinstance(ax, matplotlib.axes.Axes)


def test_plot_adsorption_isotherm_returns_axes():
    model = LatticeGasAdsorption(adsorption_energy=3.0e-20, mass=4.65e-26, T=300.0)
    P = np.logspace(0, 8, 100)
    ax = plot_adsorption_isotherm(model, P)
    assert isinstance(ax, matplotlib.axes.Axes)
    assert len(ax.lines) == 2  # the isotherm curve plus the theta=0.5 reference line
