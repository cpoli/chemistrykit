"""MD plotting: energy conservation, trajectories, g(r), and speed-distribution checks.

matplotlib-based, following the established chemistrykit visualizer
convention (:mod:`chemistrykit.kinetics.visualizers.kinetics_plots`):
thin functions that accept an optional ``ax`` (creating a new figure if
omitted) and return the ``Axes`` they drew on.
"""

from __future__ import annotations

import matplotlib.pyplot as plt
import numpy as np

__all__ = [
    "plot_energy_conservation",
    "plot_radial_distribution_function",
    "plot_trajectory_2d",
    "plot_speed_distribution",
]


def plot_energy_conservation(result, ax=None, **kwargs):
    """Plot kinetic, potential, and total energy against time.

    Parameters
    ----------
    result : chemistrykit.md.core.base_system.MDResult
    ax : matplotlib.axes.Axes, optional
        Axes to draw on; a new figure is created if omitted.
    **kwargs
        Forwarded to each ``ax.plot`` call.

    Returns
    -------
    matplotlib.axes.Axes
    """
    if ax is None:
        _, ax = plt.subplots()
    ax.plot(result.t, result.kinetic_energy, label="kinetic", **kwargs)
    ax.plot(result.t, result.potential_energy, label="potential", **kwargs)
    ax.plot(result.t, result.total_energy, label="total", **kwargs)
    ax.set_xlabel("t")
    ax.set_ylabel("energy")
    ax.set_title("MD energy conservation")
    ax.legend()
    return ax


def plot_radial_distribution_function(r, g, ax=None, **kwargs):
    """Plot a radial distribution function g(r).

    Parameters
    ----------
    r, g : array-like of float
        Output of :meth:`chemistrykit.md.systems.lj_fluid.LJFluid.radial_distribution_function`.
    ax : matplotlib.axes.Axes, optional
        Axes to draw on; a new figure is created if omitted.
    **kwargs
        Forwarded to ``ax.plot``.

    Returns
    -------
    matplotlib.axes.Axes
    """
    if ax is None:
        _, ax = plt.subplots()
    ax.plot(r, g, **kwargs)
    ax.axhline(1.0, color="gray", linestyle=":", linewidth=0.8)
    ax.set_xlabel("r")
    ax.set_ylabel("g(r)")
    ax.set_title("Radial distribution function")
    return ax


def plot_trajectory_2d(result, particle_indices=None, ax=None, **kwargs):
    """Plot an x-y projection of one or more particles' trajectories.

    Parameters
    ----------
    result : chemistrykit.md.core.base_system.MDResult
    particle_indices : sequence of int, optional
        Which particles to plot; defaults to all of them.
    ax : matplotlib.axes.Axes, optional
        Axes to draw on; a new figure is created if omitted.
    **kwargs
        Forwarded to each ``ax.plot`` call.

    Returns
    -------
    matplotlib.axes.Axes
    """
    if ax is None:
        _, ax = plt.subplots()
    n = result.positions.shape[1]
    indices = range(n) if particle_indices is None else particle_indices
    for i in indices:
        ax.plot(result.positions[:, i, 0], result.positions[:, i, 1], **kwargs)
    ax.set_xlabel("x")
    ax.set_ylabel("y")
    ax.set_title("Particle trajectories (x-y projection)")
    ax.set_aspect("equal", adjustable="datalim")
    return ax


def plot_speed_distribution(speeds, mass, temperature, ax=None, n_bins: int = 40, k_b: float = 1.0, **kwargs):
    """Plot a histogram of particle speeds against the predicted Maxwell-Boltzmann pdf.

    A direct cross-check between an MD trajectory's velocity ensemble and
    :mod:`chemistrykit.statmech`'s
    :class:`~chemistrykit.statmech.systems.maxwell_boltzmann.MaxwellBoltzmannSpeedDistribution`
    (imported lazily here to keep :mod:`chemistrykit.md` independently
    importable of :mod:`chemistrykit.statmech`).

    Parameters
    ----------
    speeds : array-like of float
        Sampled particle speeds (e.g. from :meth:`chemistrykit.md.core.base_system.MDResult.speeds`).
    mass : float
    temperature : float
        The distribution's temperature (e.g. the MD run's target/measured
        temperature), in the same unit system as `mass` and `k_b`.
    ax : matplotlib.axes.Axes, optional
        Axes to draw on; a new figure is created if omitted.
    n_bins : int, default 40
    k_b : float, default 1.0
    **kwargs
        Forwarded to ``ax.hist``.

    Returns
    -------
    matplotlib.axes.Axes
    """
    from chemistrykit.statmech.systems.maxwell_boltzmann import MaxwellBoltzmannSpeedDistribution

    if ax is None:
        _, ax = plt.subplots()
    speeds = np.asarray(speeds, dtype=np.float64).ravel()
    ax.hist(speeds, bins=n_bins, density=True, alpha=0.6, label="MD data", **kwargs)
    distribution = MaxwellBoltzmannSpeedDistribution(mass=mass, temperature=temperature, k_b=k_b)
    v = np.linspace(0.0, speeds.max() * 1.2, 300)
    ax.plot(v, distribution.pdf(v), color="crimson", label="Maxwell-Boltzmann")
    ax.set_xlabel("speed")
    ax.set_ylabel("probability density")
    ax.set_title("Speed distribution")
    ax.legend()
    return ax
