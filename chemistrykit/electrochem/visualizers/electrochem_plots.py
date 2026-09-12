"""Electrochemistry plotting: Tafel plots, Nernst concentration dependence, and battery discharge curves.

matplotlib-based, following the shared chemistrykit visualizer
convention: thin functions that accept an optional ``ax`` (creating a new
figure if omitted) and return the ``Axes`` they drew on.
"""

from __future__ import annotations

import matplotlib.pyplot as plt
import numpy as np

from chemistrykit.electrochem.systems.butler_volmer import butler_volmer_current_density
from chemistrykit.electrochem.systems.nernst import nernst_potential

__all__ = ["plot_tafel", "plot_nernst_concentration_dependence", "plot_discharge_curve"]


def plot_tafel(i0: float, eta_range=(-0.4, 0.4), alpha: float = 0.5, n: int = 1, fit=None, ax=None, n_points: int = 200, **kwargs):
    """Plot a Tafel plot: ``log10|i|`` against overpotential, from the full Butler-Volmer equation.

    Parameters
    ----------
    i0 : float
        Exchange current density.
    eta_range : tuple of float, default (-0.4, 0.4)
        Overpotential range to plot, in V.
    alpha : float, default 0.5
    n : int, default 1
    fit : chemistrykit.electrochem.systems.butler_volmer.TafelFit, optional
        If given, its linear (high-overpotential) prediction is overlaid
        for comparison against the full nonlinear curve.
    ax : matplotlib.axes.Axes, optional
        Axes to draw on; a new figure is created if omitted.
    n_points : int, default 200
    **kwargs
        Forwarded to the main ``ax.plot`` call.

    Returns
    -------
    matplotlib.axes.Axes
    """
    if ax is None:
        _, ax = plt.subplots()
    eta = np.linspace(eta_range[0], eta_range[1], n_points)
    i = butler_volmer_current_density(i0, eta, alpha=alpha, n=n)
    ax.plot(np.log10(np.abs(i)), eta, label="Butler-Volmer", **kwargs)
    if fit is not None:
        i_line = i[np.abs(eta) > 0.15 * max(abs(eta_range[0]), abs(eta_range[1]))]
        if i_line.size:
            ax.plot(np.log10(np.abs(i_line)), fit.predict(i_line), "k--", linewidth=1, label="Tafel fit")
            ax.legend()
    ax.set_xlabel(r"$\log_{10}|i|$")
    ax.set_ylabel(r"$\eta$ (V)")
    ax.set_title("Tafel plot")
    return ax


def plot_nernst_concentration_dependence(E_standard: float, n: int, Q_range=(1e-3, 1e3), ax=None, n_points: int = 200, **kwargs):
    """Plot cell potential (Nernst equation) against reaction quotient `Q`, on a log-Q axis.

    Parameters
    ----------
    E_standard : float
        Standard cell potential, in V.
    n : int
        Electrons transferred.
    Q_range : tuple of float, default (1e-3, 1e3)
        Range of `Q` to plot (log-spaced).
    ax : matplotlib.axes.Axes, optional
    n_points : int, default 200
    **kwargs
        Forwarded to ``ax.plot``.

    Returns
    -------
    matplotlib.axes.Axes
    """
    if ax is None:
        _, ax = plt.subplots()
    Q = np.logspace(np.log10(Q_range[0]), np.log10(Q_range[1]), n_points)
    E = nernst_potential(E_standard, n, Q)
    ax.plot(Q, E, **kwargs)
    ax.axhline(E_standard, color="gray", linestyle="--", linewidth=0.8, label=r"$E^\circ$")
    ax.set_xscale("log")
    ax.set_xlabel("Q")
    ax.set_ylabel("E (V)")
    ax.set_title("Nernst equation: cell potential vs. reaction quotient")
    ax.legend()
    return ax


def plot_discharge_curve(model, t=None, n_points: int = 200, ax=None, **kwargs):
    """Plot a battery's terminal voltage against elapsed discharge time.

    Parameters
    ----------
    model : chemistrykit.electrochem.core.base_system.BatteryDischargeModel
    t : array-like of float, optional
        Times to evaluate at; defaults to a grid spanning the model's
        full discharge time (for a :class:`~chemistrykit.electrochem.systems.battery.ConstantCurrentBattery`).
    n_points : int, default 200
    ax : matplotlib.axes.Axes, optional
    **kwargs
        Forwarded to ``ax.plot``.

    Returns
    -------
    matplotlib.axes.Axes
    """
    if ax is None:
        _, ax = plt.subplots()
    if t is None:
        t_total = model.discharge_time()
        t = np.linspace(0.0, 1.1 * t_total, n_points)
    result = model.discharge_curve(t)
    ax.plot(result.t, result.voltage, **kwargs)
    ax.set_xlabel("t (h)")
    ax.set_ylabel("Terminal voltage (V)")
    ax.set_title("Constant-current discharge curve")
    return ax
