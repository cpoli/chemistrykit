"""Polymer-chemistry plotting: chain-size scaling, molecular-weight distributions, Carothers curves, and free-radical kinetics.

matplotlib-based, following the shared chemistrykit visualizer
convention: thin functions that accept an optional ``ax`` (creating a new
figure if omitted) and return the ``Axes`` they drew on. Free-radical
polymerization trajectories are
:class:`chemistrykit.kinetics.core.base_system.KineticsResult` instances
(since :mod:`chemistrykit.polymer.systems.chain_growth` reuses
:mod:`chemistrykit.kinetics`'s network engine), so
:func:`chemistrykit.kinetics.visualizers.kinetics_plots.plot_concentration_vs_time`
already plots them directly -- :func:`plot_free_radical_kinetics` below
is a thin, polymer-labeled wrapper around it, mirroring
:func:`chemistrykit.photochem.visualizers.photochem_plots.plot_state_populations`.
"""

from __future__ import annotations

import matplotlib.pyplot as plt
import numpy as np

from chemistrykit.kinetics.visualizers.kinetics_plots import plot_concentration_vs_time
from chemistrykit.polymer.core.base_system import PolymerChainModel

__all__ = [
    "plot_chain_scaling",
    "plot_molecular_weight_distribution",
    "plot_carothers_curve",
    "plot_free_radical_kinetics",
]


def plot_chain_scaling(models: dict[str, PolymerChainModel], n_values, b: float, ax=None):
    """Plot end-to-end distance vs. chain length (log-log) for several chain models side by side.

    Parameters
    ----------
    models : dict of str to chemistrykit.polymer.core.base_system.PolymerChainModel
        Label -> model (e.g. ``{"ideal": IdealChain(), "good solvent": RealChain.good_solvent()}``).
    n_values : array-like of float
        Chain lengths (number of segments) to evaluate at.
    b : float
        Segment length, shared across all models.
    ax : matplotlib.axes.Axes, optional

    Returns
    -------
    matplotlib.axes.Axes
    """
    if ax is None:
        _, ax = plt.subplots()
    n_values = np.asarray(n_values, dtype=np.float64)
    for label, model in models.items():
        ax.loglog(n_values, model.end_to_end_distance(n_values, b), label=label)
    ax.set_xlabel("n (segments)")
    ax.set_ylabel("R")
    ax.set_title("Chain-size scaling")
    ax.legend()
    return ax


def plot_molecular_weight_distribution(x, fraction, ax=None, label: str = "", **kwargs):
    """Plot a chain-length (or molar-mass) distribution, e.g. the Flory-Schulz number or weight fraction.

    Parameters
    ----------
    x, fraction : array-like of float
        Chain length (or molar mass) and the corresponding fraction.
    ax : matplotlib.axes.Axes, optional
    label : str, optional
    **kwargs
        Forwarded to ``ax.plot``.

    Returns
    -------
    matplotlib.axes.Axes
    """
    if ax is None:
        _, ax = plt.subplots()
    ax.plot(np.asarray(x, dtype=np.float64), np.asarray(fraction, dtype=np.float64), label=label or None, **kwargs)
    ax.set_xlabel("chain length x")
    ax.set_ylabel("fraction")
    ax.set_title("Molecular-weight distribution")
    if label:
        ax.legend()
    return ax


def plot_carothers_curve(p, Xn, ax=None, **kwargs):
    """Plot the Carothers-equation degree of polymerization against extent of reaction.

    Parameters
    ----------
    p, Xn : array-like of float
    ax : matplotlib.axes.Axes, optional
    **kwargs
        Forwarded to ``ax.plot`` (e.g. ``label=...``).

    Returns
    -------
    matplotlib.axes.Axes
    """
    if ax is None:
        _, ax = plt.subplots()
    ax.plot(np.asarray(p, dtype=np.float64), np.asarray(Xn, dtype=np.float64), **kwargs)
    ax.set_xlabel("extent of reaction p")
    ax.set_ylabel(r"$\bar{X}_n$")
    ax.set_title("Carothers equation")
    return ax


def plot_free_radical_kinetics(result, species=None, ax=None, **kwargs):
    """Plot initiator/monomer/radical/dead-polymer concentration trajectories against time.

    A thin, polymer-labeled wrapper around
    :func:`chemistrykit.kinetics.visualizers.kinetics_plots.plot_concentration_vs_time`.

    Parameters
    ----------
    result : chemistrykit.kinetics.core.base_system.KineticsResult
        The result of integrating
        :func:`chemistrykit.polymer.systems.chain_growth.free_radical_network`.
    species : sequence of str, optional
    ax : matplotlib.axes.Axes, optional
    **kwargs
        Forwarded to ``ax.plot``.

    Returns
    -------
    matplotlib.axes.Axes
    """
    ax = plot_concentration_vs_time(result, species=species, ax=ax, **kwargs)
    ax.set_ylabel("concentration")
    ax.set_title("Free-radical polymerization kinetics")
    return ax
