"""Solutions plotting: titration curves and Debye-Huckel activity-coefficient curves.

matplotlib-based, following the established chemistrykit visualizer
convention (:mod:`chemistrykit.kinetics.visualizers.kinetics_plots`):
thin functions that accept an optional ``ax`` (creating a new figure if
omitted) and return the ``Axes`` they drew on.
"""

from __future__ import annotations

import matplotlib.pyplot as plt
import numpy as np

__all__ = ["plot_titration_curve", "plot_activity_coefficients"]


def plot_titration_curve(titration, Vb, ax=None, mark_equivalence: bool = True, **kwargs):
    """Plot a titration curve (pH vs. titrant volume), optionally marking the equivalence point.

    Parameters
    ----------
    titration : chemistrykit.solutions.core.base_system.Titration
    Vb : array-like of float
        Titrant volumes to evaluate, in L.
    ax : matplotlib.axes.Axes, optional
        Axes to draw on; a new figure is created if omitted.
    mark_equivalence : bool, default True
        If True, mark the numerically detected equivalence point.
    **kwargs
        Forwarded to ``ax.plot``.

    Returns
    -------
    matplotlib.axes.Axes
    """
    if ax is None:
        _, ax = plt.subplots()
    result = titration.curve(Vb)
    ax.plot(result.Vb, result.pH, **kwargs)
    if mark_equivalence:
        Vb_eq = titration.find_equivalence_point(Vb)
        idx = int(np.argmin(np.abs(result.Vb - Vb_eq)))
        ax.scatter([Vb_eq], [result.pH[idx]], color="crimson", zorder=3, label="equivalence point")
        ax.legend()
    ax.set_xlabel("titrant volume")
    ax.set_ylabel("pH")
    ax.set_title("Titration curve")
    return ax


def plot_activity_coefficients(I_range, z_values, extended: bool = False, Ba: float = 1.0, ax=None, **kwargs):
    """Plot Debye-Huckel activity coefficients vs. ionic strength for one or more ion charges.

    Parameters
    ----------
    I_range : array-like of float
        Ionic strengths to evaluate, in mol/L.
    z_values : sequence of float
        Charge numbers to plot one curve per.
    extended : bool, default False
        If True, use the extended Debye-Huckel law; otherwise the
        limiting law.
    Ba : float, default 1.0
        Dimensionless ion-size product, forwarded to the extended law.
    ax : matplotlib.axes.Axes, optional
        Axes to draw on; a new figure is created if omitted.
    **kwargs
        Forwarded to ``ax.plot``.

    Returns
    -------
    matplotlib.axes.Axes
    """
    from chemistrykit.solutions.systems.activity import (
        activity_coefficient_debye_huckel_extended,
        activity_coefficient_debye_huckel_limiting,
    )

    if ax is None:
        _, ax = plt.subplots()
    I_range = np.asarray(I_range, dtype=np.float64)
    for z in z_values:
        if extended:
            gamma = [activity_coefficient_debye_huckel_extended(z, I, Ba=Ba) for I in I_range]
        else:
            gamma = [activity_coefficient_debye_huckel_limiting(z, I) for I in I_range]
        ax.plot(I_range, gamma, label=f"z={z}", **kwargs)
    ax.set_xlabel("ionic strength I (mol/L)")
    ax.set_ylabel("activity coefficient")
    law = "extended" if extended else "limiting"
    ax.set_title(f"Debye-Huckel {law} law")
    ax.legend()
    return ax
