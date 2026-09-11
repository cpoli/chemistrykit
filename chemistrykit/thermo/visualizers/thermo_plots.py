"""Thermo plotting: P-V isotherms, phase boundaries, van't Hoff, and P-x-y diagrams.

matplotlib-based, following the established chemistrykit visualizer
convention (:mod:`chemistrykit.kinetics.visualizers.kinetics_plots`):
thin functions that accept an optional ``ax`` (creating a new figure if
omitted) and return the ``Axes`` they drew on.
"""

from __future__ import annotations

import matplotlib.pyplot as plt
import numpy as np

__all__ = [
    "plot_isotherms",
    "plot_phase_boundary",
    "plot_van_t_hoff",
    "plot_pxy_diagram",
]


def plot_isotherms(eos_list, T: float, Vm_range, labels=None, ax=None, **kwargs):
    """Plot P vs. Vm isotherms for one or more equations of state at a fixed temperature.

    Parameters
    ----------
    eos_list : sequence of chemistrykit.thermo.core.base_system.EquationOfState
        The EOS instances to compare (e.g. an :class:`~chemistrykit.thermo.systems.equations_of_state.IdealGas`,
        :class:`~chemistrykit.thermo.systems.equations_of_state.VanDerWaals`,
        and :class:`~chemistrykit.thermo.systems.equations_of_state.RedlichKwong`
        built with the same substance's parameters).
    T : float
        Temperature of the isotherm, in K.
    Vm_range : array-like of float
        Molar volumes at which to evaluate each EOS's pressure, in m^3/mol.
    labels : sequence of str, optional
        Legend label for each entry in `eos_list`; defaults to each
        instance's class name.
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
    Vm_range = np.asarray(Vm_range, dtype=np.float64)
    labels = [type(eos).__name__ for eos in eos_list] if labels is None else labels
    for eos, label in zip(eos_list, labels):
        ax.plot(Vm_range, eos.pressure(Vm_range, T), label=label, **kwargs)
    ax.set_xlabel("Vm (m^3/mol)")
    ax.set_ylabel("P (Pa)")
    ax.set_title(f"P-V isotherm at T={T:g} K")
    ax.legend()
    return ax


def plot_phase_boundary(clausius_clapeyron, T_range, ax=None, **kwargs):
    """Plot a Clausius-Clapeyron liquid-vapor (or solid-vapor) phase boundary.

    Parameters
    ----------
    clausius_clapeyron : chemistrykit.thermo.systems.phase_equilibria.ClausiusClapeyron
    T_range : array-like of float
        Temperatures at which to evaluate the vapor pressure, in K.
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
    ax.plot(T_range, clausius_clapeyron.pressure(T_range), **kwargs)
    ax.scatter([clausius_clapeyron.T_ref], [clausius_clapeyron.P_ref], color="crimson", zorder=3, label="reference point")
    ax.set_xlabel("T (K)")
    ax.set_ylabel("P (Pa)")
    ax.set_title("Clausius-Clapeyron phase boundary")
    ax.legend()
    return ax


def plot_van_t_hoff(T, K, fit=None, ax=None, **kwargs):
    """Plot a van't Hoff plot: ``ln K`` against ``1/T``, with an optional fitted line.

    Parameters
    ----------
    T, K : array-like of float
        Temperatures and measured equilibrium constants.
    fit : chemistrykit.thermo.systems.equilibrium.VantHoffFit, optional
        A fit from :func:`chemistrykit.thermo.systems.equilibrium.fit_van_t_hoff`;
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
    K = np.asarray(K, dtype=np.float64)
    ax.scatter(1.0 / T, np.log(K), label="data", **kwargs)
    if fit is not None:
        T_line = np.linspace(np.min(T), np.max(T), 100)
        ax.plot(1.0 / T_line, np.log(fit.predict(T_line)), color="crimson", label="fit")
        ax.legend()
    ax.set_xlabel("1/T")
    ax.set_ylabel("ln K")
    ax.set_title("van't Hoff plot")
    return ax


def plot_pxy_diagram(solution, ax=None, n_points: int = 100, **kwargs):
    """Plot a P-x-y diagram for a binary ideal solution.

    Shows total vapor pressure against liquid composition (the "P-x"
    curve) and against the corresponding vapor composition (the "P-y"
    curve), the standard way of visualizing vapor-liquid equilibrium for
    an ideal binary mixture.

    Parameters
    ----------
    solution : chemistrykit.thermo.systems.mixtures.BinaryIdealSolution
    ax : matplotlib.axes.Axes, optional
        Axes to draw on; a new figure is created if omitted.
    n_points : int, default 100
        Number of composition points to evaluate.
    **kwargs
        Forwarded to both ``ax.plot`` calls.

    Returns
    -------
    matplotlib.axes.Axes
    """
    if ax is None:
        _, ax = plt.subplots()
    x_A = np.linspace(0.0, 1.0, n_points)
    P = solution.total_pressure(x_A)
    y_A = solution.vapor_composition(x_A)
    ax.plot(x_A, P, label="P-x (liquid)", **kwargs)
    ax.plot(y_A, P, label="P-y (vapor)", **kwargs)
    ax.set_xlabel("mole fraction of A")
    ax.set_ylabel("P")
    ax.set_title("P-x-y diagram")
    ax.legend()
    return ax
