"""Surface-chemistry plotting: isotherms, linearizations, and Langmuir-Hinshelwood rate curves.

matplotlib-based, following the shared chemistrykit visualizer
convention: thin functions that accept an optional ``ax`` (creating a new
figure if omitted) and return the ``Axes`` they drew on.
"""

from __future__ import annotations

import matplotlib.pyplot as plt
import numpy as np

from chemistrykit.surface.core.base_system import AdsorptionIsotherm

__all__ = ["plot_isotherm", "plot_linearization", "plot_lh_rate_vs_pressure"]


def plot_isotherm(isotherm: AdsorptionIsotherm, P_max: float, P_data=None, q_data=None, ax=None, n_points: int = 200, **kwargs):
    """Plot an adsorption isotherm's loading curve, optionally overlaid with data.

    Parameters
    ----------
    isotherm : chemistrykit.surface.core.base_system.AdsorptionIsotherm
        A fitted or hand-built isotherm instance.
    P_max : float
        Upper end of the pressure range to plot (from 0).
    P_data, q_data : array-like of float, optional
        Measured data to overlay as a scatter.
    ax : matplotlib.axes.Axes, optional
    n_points : int, default 200
        Number of points in the smooth isotherm curve.
    **kwargs
        Forwarded to the isotherm-curve ``ax.plot`` call.

    Returns
    -------
    matplotlib.axes.Axes
    """
    if ax is None:
        _, ax = plt.subplots()
    P = np.linspace(0.0, P_max, n_points)
    ax.plot(P, isotherm.loading(P), label="isotherm", **kwargs)
    if P_data is not None and q_data is not None:
        ax.scatter(P_data, q_data, color="crimson", zorder=3, label="data")
        ax.legend()
    ax.set_xlabel("P")
    ax.set_ylabel("loading")
    ax.set_title(f"{type(isotherm).__name__} isotherm")
    return ax


def plot_linearization(x, y, fit=None, ax=None, xlabel: str = "x", ylabel: str = "y"):
    """Plot a linearized isotherm fit: transformed data with the fitted line overlaid.

    Parameters
    ----------
    x, y : array-like of float
        Linearized coordinates (e.g. ``1/P`` and ``1/q`` for Langmuir).
    fit : tuple of (float, float), optional
        ``(slope, intercept)`` of the fitted line, overlaid if given.
    ax : matplotlib.axes.Axes, optional
    xlabel, ylabel : str

    Returns
    -------
    matplotlib.axes.Axes
    """
    if ax is None:
        _, ax = plt.subplots()
    x = np.asarray(x, dtype=np.float64)
    y = np.asarray(y, dtype=np.float64)
    ax.scatter(x, y, label="data")
    if fit is not None:
        slope, intercept = fit
        x_line = np.linspace(np.min(x), np.max(x), 100)
        ax.plot(x_line, slope * x_line + intercept, color="crimson", label="linear fit")
        ax.legend()
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.set_title("Isotherm linearization")
    return ax


def plot_lh_rate_vs_pressure(P_A, rate, ax=None, **kwargs):
    """Plot a Langmuir-Hinshelwood surface-reaction rate against reactant pressure.

    Parameters
    ----------
    P_A, rate : array-like of float
    ax : matplotlib.axes.Axes, optional
    **kwargs
        Forwarded to ``ax.plot``.

    Returns
    -------
    matplotlib.axes.Axes
    """
    if ax is None:
        _, ax = plt.subplots()
    ax.plot(np.asarray(P_A, dtype=np.float64), np.asarray(rate, dtype=np.float64), **kwargs)
    ax.set_xlabel(r"$P_A$")
    ax.set_ylabel("rate")
    ax.set_title("Langmuir-Hinshelwood rate vs. pressure")
    return ax
