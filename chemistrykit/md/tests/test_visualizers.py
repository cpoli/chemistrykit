"""Smoke tests for chemistrykit.md.visualizers: right return type, no crash."""

import matplotlib

matplotlib.use("Agg")

import matplotlib.axes
import numpy as np

from chemistrykit.md.systems.lj_fluid import LJFluid
from chemistrykit.md.visualizers.md_plots import (
    plot_energy_conservation,
    plot_radial_distribution_function,
    plot_speed_distribution,
    plot_trajectory_2d,
)


def test_plot_energy_conservation_returns_axes():
    fluid = LJFluid.from_lattice(n_per_side=3, cutoff=1.5, density=0.6, temperature=1.0, rng=0)
    result = fluid.run(dt=0.001, n_steps=50, sample_every=10)
    ax = plot_energy_conservation(result)
    assert isinstance(ax, matplotlib.axes.Axes)
    assert len(ax.lines) == 3


def test_plot_radial_distribution_function_returns_axes():
    fluid = LJFluid.from_lattice(n_per_side=3, cutoff=1.5, density=0.6, temperature=1.0, rng=1)
    r, g = fluid.radial_distribution_function()
    ax = plot_radial_distribution_function(r, g)
    assert isinstance(ax, matplotlib.axes.Axes)


def test_plot_trajectory_2d_returns_axes():
    fluid = LJFluid.from_lattice(n_per_side=3, cutoff=1.5, density=0.6, temperature=1.0, rng=2)
    result = fluid.run(dt=0.001, n_steps=50, sample_every=10)
    ax = plot_trajectory_2d(result, particle_indices=[0, 1])
    assert isinstance(ax, matplotlib.axes.Axes)
    assert len(ax.lines) == 2


def test_plot_speed_distribution_returns_axes():
    rng = np.random.default_rng(3)
    speeds = rng.normal(2.0, 0.3, size=500)
    speeds = np.abs(speeds)
    ax = plot_speed_distribution(speeds, mass=1.0, temperature=1.0)
    assert isinstance(ax, matplotlib.axes.Axes)
