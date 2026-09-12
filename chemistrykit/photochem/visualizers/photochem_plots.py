"""Photochemistry plotting: Jablonski-diagram population decay, Stern-Volmer plots, and photostationary-state approach.

matplotlib-based, following the shared chemistrykit visualizer
convention: thin functions that accept an optional ``ax`` (creating a new
figure if omitted) and return the ``Axes`` they drew on. Jablonski and
photoswitch population trajectories are
:class:`chemistrykit.kinetics.core.base_system.KineticsResult` instances
(since both reuse :mod:`chemistrykit.kinetics`'s network engine), so
:func:`chemistrykit.kinetics.visualizers.kinetics_plots.plot_concentration_vs_time`
already plots them directly -- the wrapper below just supplies
photochemistry-appropriate axis labels.
"""

from __future__ import annotations

import matplotlib.pyplot as plt
import numpy as np

from chemistrykit.kinetics.visualizers.kinetics_plots import plot_concentration_vs_time

__all__ = ["plot_state_populations", "plot_stern_volmer", "plot_photostationary_approach"]


def plot_state_populations(result, species=None, ax=None, **kwargs):
    """Plot excited-state (or photoswitch) population trajectories against time.

    A thin, photochemistry-labeled wrapper around
    :func:`chemistrykit.kinetics.visualizers.kinetics_plots.plot_concentration_vs_time`.

    Parameters
    ----------
    result : chemistrykit.kinetics.core.base_system.KineticsResult
        The result of integrating a
        :func:`chemistrykit.photochem.systems.jablonski.jablonski_network`
        or :func:`chemistrykit.photochem.systems.photostationary_state.photoswitch_network`.
    species : sequence of str, optional
    ax : matplotlib.axes.Axes, optional
    **kwargs
        Forwarded to ``ax.plot``.

    Returns
    -------
    matplotlib.axes.Axes
    """
    ax = plot_concentration_vs_time(result, species=species, ax=ax, **kwargs)
    ax.set_ylabel("population")
    ax.set_title("Excited-state population vs. time")
    return ax


def plot_stern_volmer(Q, intensity_ratio, fit=None, ax=None, **kwargs):
    """Plot a Stern-Volmer plot: :math:`I_0/I` against quencher concentration.

    Parameters
    ----------
    Q, intensity_ratio : array-like of float
    fit : chemistrykit.photochem.systems.stern_volmer.SternVolmerFit, optional
        If given, its predicted line is overlaid.
    ax : matplotlib.axes.Axes, optional
    **kwargs
        Forwarded to the data ``ax.scatter`` call.

    Returns
    -------
    matplotlib.axes.Axes
    """
    if ax is None:
        _, ax = plt.subplots()
    Q = np.asarray(Q, dtype=np.float64)
    intensity_ratio = np.asarray(intensity_ratio, dtype=np.float64)
    ax.scatter(Q, intensity_ratio, label="data", **kwargs)
    if fit is not None:
        Q_line = np.linspace(0.0, np.max(Q), 100)
        ax.plot(Q_line, fit.predict(Q_line), color="crimson", label="fit")
        ax.legend()
    ax.set_xlabel("[Q]")
    ax.set_ylabel(r"$I_0/I$")
    ax.set_title("Stern-Volmer plot")
    return ax


def plot_photostationary_approach(result, pss_result, ax=None):
    """Plot a photoswitch's approach to its photostationary state.

    Overlays the [B]/[A] population ratio's time evolution against the
    algebraic photostationary-state ratio it should approach.

    Parameters
    ----------
    result : chemistrykit.kinetics.core.base_system.KineticsResult
        Result of integrating
        :func:`chemistrykit.photochem.systems.photostationary_state.photoswitch_network`.
    pss_result : chemistrykit.photochem.core.base_system.PhotostationaryStateResult
    ax : matplotlib.axes.Axes, optional

    Returns
    -------
    matplotlib.axes.Axes
    """
    if ax is None:
        _, ax = plt.subplots()
    ratio = result.concentration("B") / result.concentration("A")
    ax.plot(result.t, ratio, label="[B]/[A]")
    ax.axhline(pss_result.ratio_B_over_A, color="gray", linestyle="--", linewidth=0.8, label="PSS ratio")
    ax.set_xlabel("t")
    ax.set_ylabel("[B]/[A]")
    ax.set_title("Approach to the photostationary state")
    ax.legend()
    return ax
