"""Spectroscopy plotting: stick and broadened spectra, and the Beer-Lambert linearity-deviation curve.

matplotlib-based, following the established chemistrykit visualizer
convention (:mod:`chemistrykit.quantum.visualizers.quantum_plots`): thin
functions that accept an optional ``ax`` (creating a new figure if
omitted) and return the ``Axes`` they drew on.
"""

from __future__ import annotations

import matplotlib.pyplot as plt
import numpy as np

__all__ = ["plot_stick_spectrum", "plot_broadened_spectrum", "plot_beer_lambert_deviation"]


def plot_stick_spectrum(spectrum, ax=None, **kwargs):
    """Draw a spectrum as vertical "stick" lines at each peak position.

    Parameters
    ----------
    spectrum : chemistrykit.spectro.core.base_system.Spectrum
    ax : matplotlib.axes.Axes, optional
    **kwargs
        Forwarded to ``ax.vlines``.

    Returns
    -------
    matplotlib.axes.Axes
    """
    if ax is None:
        _, ax = plt.subplots()
    ax.vlines(spectrum.positions, 0.0, spectrum.intensities, **kwargs)
    ax.set_ylabel("relative intensity")
    ax.set_title("Stick spectrum")
    return ax


def plot_broadened_spectrum(spectrum, x, ax=None, shape: str = "gaussian", fwhm: float = 1.0, fwhm_lorentzian=None, **kwargs):
    """Plot a spectrum broadened into a continuous lineshape.

    Parameters
    ----------
    spectrum : chemistrykit.spectro.core.base_system.Spectrum
    x : array-like of float
        Grid to evaluate the broadened spectrum on.
    ax : matplotlib.axes.Axes, optional
    shape : {"gaussian", "lorentzian", "voigt"}, default "gaussian"
    fwhm : float, default 1.0
    fwhm_lorentzian : float, optional
        Required for `shape="voigt"`.
    **kwargs
        Forwarded to ``ax.plot``.

    Returns
    -------
    matplotlib.axes.Axes
    """
    if ax is None:
        _, ax = plt.subplots()
    x = np.asarray(x, dtype=np.float64)
    y = spectrum.broaden(x, shape=shape, fwhm=fwhm, fwhm_lorentzian=fwhm_lorentzian)
    ax.plot(x, y, **kwargs)
    ax.set_ylabel("intensity")
    ax.set_title(f"Broadened spectrum ({shape})")
    return ax


def plot_beer_lambert_deviation(molar_absorptivity: float, path_length: float, concentrations, stray_light_fraction: float, ax=None):
    """Plot true vs. apparent (stray-light-limited) Beer-Lambert absorbance across a concentration range.

    Parameters
    ----------
    molar_absorptivity : float
        In L mol^-1 cm^-1.
    path_length : float
        In cm.
    concentrations : array-like of float
        Concentrations to evaluate, in mol/L.
    stray_light_fraction : float
        Passed to :func:`chemistrykit.spectro.systems.beer_lambert.apparent_absorbance_with_stray_light`.
    ax : matplotlib.axes.Axes, optional

    Returns
    -------
    matplotlib.axes.Axes
    """
    from chemistrykit.spectro.systems.beer_lambert import absorbance, apparent_absorbance_with_stray_light

    if ax is None:
        _, ax = plt.subplots()
    concentrations = np.asarray(concentrations, dtype=np.float64)
    true_absorbance = [absorbance(molar_absorptivity, c, path_length) for c in concentrations]
    apparent_absorbance = [apparent_absorbance_with_stray_light(molar_absorptivity, c, path_length, stray_light_fraction) for c in concentrations]
    ax.plot(concentrations, true_absorbance, "--", label="ideal Beer-Lambert")
    ax.plot(concentrations, apparent_absorbance, "-", label="apparent (stray light)")
    ax.set_xlabel("concentration (mol/L)")
    ax.set_ylabel("absorbance")
    ax.set_title("Beer-Lambert law: deviation from linearity")
    ax.legend()
    return ax
