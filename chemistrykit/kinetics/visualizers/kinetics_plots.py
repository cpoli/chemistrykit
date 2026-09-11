"""Kinetics plotting: concentration-vs-time, Arrhenius, and Lineweaver-Burk plots.

matplotlib-based, following physicskit's visualizer convention: thin
functions that accept an optional ``ax`` (creating a new figure if
omitted) and return the ``Axes`` they drew on.
"""

from __future__ import annotations

import matplotlib.pyplot as plt
import numpy as np

__all__ = [
    "plot_concentration_vs_time",
    "plot_phase_portrait",
    "plot_arrhenius",
    "plot_lineweaver_burk",
]


def plot_concentration_vs_time(result, species=None, ax=None, **kwargs):
    """Plot species concentration trajectories against time.

    Parameters
    ----------
    result : chemistrykit.kinetics.core.base_system.KineticsResult
        The result of a :meth:`~chemistrykit.kinetics.core.base_system.ReactionNetwork.integrate` call.
    species : sequence of str, optional
        Subset of ``result.species`` to plot; defaults to all of them.
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
    names = result.species if species is None else species
    for name in names:
        ax.plot(result.t, result.concentration(name), label=name, **kwargs)
    ax.set_xlabel("t")
    ax.set_ylabel("concentration")
    ax.set_title("Concentration vs. time")
    ax.legend()
    return ax


def plot_phase_portrait(result, x_species: str, y_species: str, ax=None, **kwargs):
    """Plot one species' concentration against another's (a 2-species phase portrait).

    Useful for e.g. visualizing a
    :class:`~chemistrykit.kinetics.systems.oscillators.Brusselator` limit
    cycle in its ``(X, Y)`` plane.

    Parameters
    ----------
    result : chemistrykit.kinetics.core.base_system.KineticsResult
    x_species, y_species : str
        Species names from ``result.species`` for the horizontal/vertical axes.
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
    ax.plot(result.concentration(x_species), result.concentration(y_species), **kwargs)
    ax.set_xlabel(x_species)
    ax.set_ylabel(y_species)
    ax.set_title(f"{y_species} vs. {x_species}")
    return ax


def plot_arrhenius(T, k, fit=None, ax=None, **kwargs):
    """Plot an Arrhenius plot: ``ln k`` against ``1/T``, with an optional fitted line.

    Parameters
    ----------
    T, k : array-like of float
        Temperatures and measured rate constants.
    fit : chemistrykit.kinetics.systems.arrhenius.ArrheniusFit, optional
        A fit from :func:`chemistrykit.kinetics.systems.arrhenius.fit_arrhenius`;
        if given, its predicted line is overlaid.
    ax : matplotlib.axes.Axes, optional
        Axes to draw on; a new figure is created if omitted.
    **kwargs
        Forwarded to the data ``ax.scatter`` call.

    Returns
    -------
    matplotlib.axes.Axes
    """
    if ax is None:
        _, ax = plt.subplots()
    T = np.asarray(T, dtype=np.float64)
    k = np.asarray(k, dtype=np.float64)
    ax.scatter(1.0 / T, np.log(k), label="data", **kwargs)
    if fit is not None:
        T_line = np.linspace(np.min(T), np.max(T), 100)
        ax.plot(1.0 / T_line, np.log(fit.predict(T_line)), color="crimson", label="fit")
        ax.legend()
    ax.set_xlabel("1/T")
    ax.set_ylabel("ln k")
    ax.set_title("Arrhenius plot")
    return ax


def plot_lineweaver_burk(S, v, fit=None, ax=None, **kwargs):
    """Plot a Lineweaver-Burk (double-reciprocal) plot: ``1/v`` against ``1/[S]``.

    Parameters
    ----------
    S, v : array-like of float
        Substrate concentrations and measured initial rates.
    fit : chemistrykit.kinetics.systems.enzyme.MichaelisMentenFit, optional
        A fit from :func:`chemistrykit.kinetics.systems.enzyme.fit_lineweaver_burk`;
        if given, its predicted line is overlaid.
    ax : matplotlib.axes.Axes, optional
        Axes to draw on; a new figure is created if omitted.
    **kwargs
        Forwarded to the data ``ax.scatter`` call.

    Returns
    -------
    matplotlib.axes.Axes
    """
    if ax is None:
        _, ax = plt.subplots()
    S = np.asarray(S, dtype=np.float64)
    v = np.asarray(v, dtype=np.float64)
    inv_S = 1.0 / S
    ax.scatter(inv_S, 1.0 / v, label="data", **kwargs)
    if fit is not None:
        S_line = np.linspace(np.min(S), np.max(S), 100)
        ax.plot(1.0 / S_line, 1.0 / fit.predict(S_line), color="crimson", label="fit")
        ax.legend()
    ax.set_xlabel("1/[S]")
    ax.set_ylabel("1/v")
    ax.set_title("Lineweaver-Burk plot")
    return ax
