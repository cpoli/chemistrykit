"""Statmech plotting: partition-function thermodynamics, Maxwell-Boltzmann, adsorption isotherms.

matplotlib-based, following the established chemistrykit visualizer
convention (:mod:`chemistrykit.kinetics.visualizers.kinetics_plots`):
thin functions that accept an optional ``ax`` (creating a new figure if
omitted) and return the ``Axes`` they drew on.
"""

from __future__ import annotations

import matplotlib.pyplot as plt
import numpy as np

__all__ = [
    "plot_heat_capacity_vs_temperature",
    "plot_maxwell_boltzmann",
    "plot_adsorption_isotherm",
]


def plot_heat_capacity_vs_temperature(partition_function, T_range, ax=None, **kwargs):
    """Plot a partition function's heat capacity Cv against temperature.

    Parameters
    ----------
    partition_function : chemistrykit.statmech.core.base_system.PartitionFunction
    T_range : array-like of float
        Temperatures at which to evaluate Cv, in K.
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
    T_range = np.asarray(T_range, dtype=np.float64)
    Cv = np.array([partition_function.heat_capacity_v(T, N=1.0) for T in T_range])
    ax.plot(T_range, Cv, **kwargs)
    ax.set_xlabel("T (K)")
    ax.set_ylabel("Cv (J/K per molecule)")
    ax.set_title("Heat capacity vs. temperature")
    return ax


def plot_maxwell_boltzmann(distribution, v_max=None, ax=None, n_points: int = 300, **kwargs):
    """Plot a Maxwell-Boltzmann speed pdf, with markers for the characteristic speeds.

    Parameters
    ----------
    distribution : chemistrykit.statmech.systems.maxwell_boltzmann.MaxwellBoltzmannSpeedDistribution
    v_max : float, optional
        Upper speed limit to plot; defaults to 3x the rms speed.
    ax : matplotlib.axes.Axes, optional
        Axes to draw on; a new figure is created if omitted.
    n_points : int, default 300
    **kwargs
        Forwarded to ``ax.plot``.

    Returns
    -------
    matplotlib.axes.Axes
    """
    if ax is None:
        _, ax = plt.subplots()
    v_max = 3.0 * distribution.rms_speed() if v_max is None else v_max
    v = np.linspace(0.0, v_max, n_points)
    ax.plot(v, distribution.pdf(v), **kwargs)
    for speed, label in (
        (distribution.most_probable_speed(), "v_p"),
        (distribution.mean_speed(), "v_mean"),
        (distribution.rms_speed(), "v_rms"),
    ):
        ax.axvline(speed, linestyle=":", linewidth=0.8, label=label)
    ax.set_xlabel("speed")
    ax.set_ylabel("probability density")
    ax.set_title("Maxwell-Boltzmann speed distribution")
    ax.legend()
    return ax


def plot_adsorption_isotherm(model, P_range, ax=None, **kwargs):
    """Plot a lattice-gas adsorption isotherm, coverage vs. pressure.

    Parameters
    ----------
    model : chemistrykit.statmech.systems.lattice_gas.LatticeGasAdsorption
    P_range : array-like of float
        Pressures at which to evaluate coverage, in Pa.
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
    P_range = np.asarray(P_range, dtype=np.float64)
    ax.plot(P_range, model.coverage(P_range), **kwargs)
    ax.axhline(0.5, color="gray", linestyle=":", linewidth=0.8)
    ax.set_xlabel("P (Pa)")
    ax.set_ylabel("coverage (theta)")
    ax.set_title("Lattice-gas (Langmuir) adsorption isotherm")
    return ax
