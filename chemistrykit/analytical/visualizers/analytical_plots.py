"""Analytical-chemistry plotting: titration curves, van Deemter curves, chromatograms, and calibration curves.

matplotlib-based, following the shared chemistrykit visualizer
convention: thin functions that accept an optional ``ax`` (creating a new
figure if omitted) and return the ``Axes`` they drew on.
"""

from __future__ import annotations

import matplotlib.pyplot as plt
import numpy as np

__all__ = [
    "plot_titration_curve",
    "plot_van_deemter",
    "plot_chromatogram",
    "plot_calibration_curve",
]


def plot_titration_curve(V, response, V_equiv=None, ax=None, ylabel: str = "response", **kwargs):
    """Plot a generic titration curve (pH, E, or pM vs. titrant volume).

    Parameters
    ----------
    V, response : array-like of float
    V_equiv : float, optional
        Equivalence volume to mark with a vertical dashed line.
    ax : matplotlib.axes.Axes, optional
    ylabel : str, default "response"
    **kwargs
        Forwarded to ``ax.plot``.

    Returns
    -------
    matplotlib.axes.Axes
    """
    if ax is None:
        _, ax = plt.subplots()
    ax.plot(np.asarray(V, dtype=np.float64), np.asarray(response, dtype=np.float64), **kwargs)
    if V_equiv is not None:
        ax.axvline(V_equiv, color="gray", linestyle="--", linewidth=0.8)
    ax.set_xlabel("titrant volume")
    ax.set_ylabel(ylabel)
    ax.set_title("Titration curve")
    return ax


def plot_van_deemter(u, H, u_opt=None, H_min=None, ax=None):
    """Plot the van Deemter curve H(u), optionally marking the optimum.

    Parameters
    ----------
    u, H : array-like of float
    u_opt, H_min : float, optional
        Optimum flow velocity/minimum plate height to mark.
    ax : matplotlib.axes.Axes, optional

    Returns
    -------
    matplotlib.axes.Axes
    """
    if ax is None:
        _, ax = plt.subplots()
    ax.plot(np.asarray(u, dtype=np.float64), np.asarray(H, dtype=np.float64))
    if u_opt is not None and H_min is not None:
        ax.plot([u_opt], [H_min], marker="o", color="crimson")
        ax.axvline(u_opt, color="gray", linestyle="--", linewidth=0.8)
    ax.set_xlabel("linear velocity u")
    ax.set_ylabel("plate height H")
    ax.set_title("van Deemter curve")
    return ax


def plot_chromatogram(t, chrom, ax=None, **kwargs):
    """Plot a simulated chromatogram.

    Parameters
    ----------
    t, chrom : array-like of float
    ax : matplotlib.axes.Axes, optional
    **kwargs
        Forwarded to ``ax.plot``.

    Returns
    -------
    matplotlib.axes.Axes
    """
    if ax is None:
        _, ax = plt.subplots()
    ax.plot(np.asarray(t, dtype=np.float64), np.asarray(chrom, dtype=np.float64), **kwargs)
    ax.set_xlabel("time")
    ax.set_ylabel("detector response")
    ax.set_title("Simulated chromatogram")
    return ax


def plot_calibration_curve(concentration, signal, calibration=None, lod=None, ax=None):
    """Plot calibration data with the fitted line, optionally marking the LOD.

    Parameters
    ----------
    concentration, signal : array-like of float
        Measured calibration standards.
    calibration : chemistrykit.analytical.systems.calibration.LinearCalibration, optional
        Fitted calibration to overlay.
    lod : float, optional
        Limit of detection to mark with a vertical dashed line.
    ax : matplotlib.axes.Axes, optional

    Returns
    -------
    matplotlib.axes.Axes
    """
    if ax is None:
        _, ax = plt.subplots()
    concentration = np.asarray(concentration, dtype=np.float64)
    signal = np.asarray(signal, dtype=np.float64)
    ax.scatter(concentration, signal, label="standards", zorder=3)
    if calibration is not None:
        x_line = np.linspace(0.0, concentration.max(), 100)
        ax.plot(x_line, calibration.predict_signal(x_line), color="crimson", label="fit")
        ax.legend()
    if lod is not None:
        ax.axvline(lod, color="gray", linestyle="--", linewidth=0.8, label="LOD")
    ax.set_xlabel("concentration")
    ax.set_ylabel("signal")
    ax.set_title("Calibration curve")
    return ax
