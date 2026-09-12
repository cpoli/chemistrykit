"""Crystal-chemistry plotting: packing fractions, Madelung-sum convergence, XRD patterns, and defect concentrations.

matplotlib-based, following the shared chemistrykit visualizer
convention: thin functions that accept an optional ``ax`` (creating a new
figure if omitted) and return the ``Axes`` they drew on.
"""

from __future__ import annotations

import matplotlib.pyplot as plt
import numpy as np

__all__ = [
    "plot_packing_fractions",
    "plot_madelung_convergence",
    "plot_xrd_pattern",
    "plot_defect_concentration_vs_temperature",
]


def plot_packing_fractions(names, fractions, ax=None):
    """Bar chart comparing several lattices' atomic packing factors.

    Parameters
    ----------
    names : sequence of str
    fractions : sequence of float
    ax : matplotlib.axes.Axes, optional

    Returns
    -------
    matplotlib.axes.Axes
    """
    if ax is None:
        _, ax = plt.subplots()
    ax.bar(names, fractions, color="steelblue")
    ax.set_ylabel("packing fraction")
    ax.set_title("Atomic packing factor by lattice type")
    ax.set_ylim(0, 1)
    return ax


def plot_madelung_convergence(n_shells, values, literature_value=None, ax=None):
    """Plot a Madelung-sum estimate against summation cutoff, showing convergence.

    Parameters
    ----------
    n_shells : array-like of int
    values : array-like of float
    literature_value : float, optional
        Reference value to overlay as a horizontal line.
    ax : matplotlib.axes.Axes, optional

    Returns
    -------
    matplotlib.axes.Axes
    """
    if ax is None:
        _, ax = plt.subplots()
    ax.plot(np.asarray(n_shells), np.asarray(values, dtype=np.float64), marker="o")
    if literature_value is not None:
        ax.axhline(literature_value, color="crimson", linestyle="--", label=f"literature = {literature_value}")
        ax.legend()
    ax.set_xlabel("n_shells")
    ax.set_ylabel("Madelung constant estimate")
    ax.set_title("Evjen-method convergence")
    return ax


def plot_xrd_pattern(peaks, ax=None):
    """Plot a simulated powder-XRD stick pattern from a list of XRDPeak.

    Parameters
    ----------
    peaks : sequence of chemistrykit.crystal.systems.xrd.XRDPeak
    ax : matplotlib.axes.Axes, optional

    Returns
    -------
    matplotlib.axes.Axes
    """
    if ax is None:
        _, ax = plt.subplots()
    two_theta = [p.two_theta for p in peaks]
    intensity = np.array([p.relative_intensity for p in peaks], dtype=np.float64)
    intensity = intensity / intensity.max() * 100.0 if len(intensity) else intensity
    ax.vlines(two_theta, 0.0, intensity, color="darkslateblue")
    for tt, inten, peak in zip(two_theta, intensity, peaks, strict=True):
        ax.annotate(str(peak.hkl), (tt, inten), textcoords="offset points", xytext=(0, 3), fontsize=8, ha="center")
    ax.set_xlabel(r"$2\theta$ (degrees)")
    ax.set_ylabel("relative intensity")
    ax.set_title("Simulated powder XRD pattern")
    return ax


def plot_defect_concentration_vs_temperature(T, n_defects, ax=None, label=None):
    """Plot defect concentration against temperature (typically on a log-y scale).

    Parameters
    ----------
    T : array-like of float
    n_defects : array-like of float
    ax : matplotlib.axes.Axes, optional
    label : str, optional

    Returns
    -------
    matplotlib.axes.Axes
    """
    if ax is None:
        _, ax = plt.subplots()
    ax.semilogy(np.asarray(T, dtype=np.float64), np.asarray(n_defects, dtype=np.float64), marker="o", label=label)
    if label is not None:
        ax.legend()
    ax.set_xlabel("T (K)")
    ax.set_ylabel("defect concentration")
    ax.set_title("Point-defect concentration vs. temperature")
    return ax
