r"""The Beer-Lambert absorbance law, and its instrumental (stray-light) deviation from linearity.

See Skoog, Holler & Crouch, *Principles of Instrumental Analysis*, 7th
ed. (2018), Ch. 13.2, or Harris, *Quantitative Chemical Analysis*, 9th
ed., Ch. 18, throughout.
"""

from __future__ import annotations

import numpy as np

__all__ = ["absorbance", "transmittance", "concentration_from_absorbance", "apparent_absorbance_with_stray_light"]


def absorbance(molar_absorptivity: float, concentration, path_length: float):
    r"""The Beer-Lambert law: :math:`A=\varepsilon c l`.

    Linear in concentration by construction -- the exact relationship
    assumed to hold for a dilute, non-interacting, monochromatic-light
    measurement (Skoog, Holler & Crouch, *Principles of Instrumental
    Analysis*, 7th ed., Ch. 13.2). See
    :func:`apparent_absorbance_with_stray_light` for one well-known real
    instrumental effect that breaks this linearity at high absorbance.

    Parameters
    ----------
    molar_absorptivity : float
        Molar absorptivity :math:`\varepsilon`, in L mol^-1 cm^-1.
    concentration : float or array-like of float
        Concentration `c`, in mol/L.
    path_length : float
        Path length `l`, in cm.

    Returns
    -------
    float or ndarray
        Absorbance (dimensionless).

    Examples
    --------
    >>> round(absorbance(molar_absorptivity=5000.0, concentration=2.0e-5, path_length=1.0), 4)
    0.1
    >>> absorbance(5000.0, 4.0e-5, 1.0) == 2.0 * absorbance(5000.0, 2.0e-5, 1.0)
    True
    """
    result = molar_absorptivity * np.asarray(concentration, dtype=np.float64) * path_length
    return float(result) if result.ndim == 0 else result


def transmittance(absorbance_value):
    r"""Transmittance :math:`T=10^{-A}`, the fraction of incident light transmitted.

    Parameters
    ----------
    absorbance_value : float or array-like of float

    Returns
    -------
    float or ndarray
        In :math:`[0, 1]`.

    Examples
    --------
    >>> round(transmittance(0.0), 6)
    1.0
    >>> round(transmittance(1.0), 6)
    0.1
    """
    result = 10.0 ** (-np.asarray(absorbance_value, dtype=np.float64))
    return float(result) if result.ndim == 0 else result


def concentration_from_absorbance(absorbance_value, molar_absorptivity: float, path_length: float):
    r"""Invert the Beer-Lambert law for concentration: :math:`c=A/(\varepsilon l)`.

    Parameters
    ----------
    absorbance_value : float or array-like of float
    molar_absorptivity : float
        In L mol^-1 cm^-1.
    path_length : float
        In cm.

    Returns
    -------
    float or ndarray
        Concentration, in mol/L.

    Examples
    --------
    >>> c = concentration_from_absorbance(absorbance_value=0.1, molar_absorptivity=5000.0, path_length=1.0)
    >>> round(c, 8)
    2e-05
    """
    result = np.asarray(absorbance_value, dtype=np.float64) / (molar_absorptivity * path_length)
    return float(result) if result.ndim == 0 else result


def apparent_absorbance_with_stray_light(molar_absorptivity: float, concentration: float, path_length: float, stray_light_fraction: float) -> float:
    r"""The *apparent* absorbance measured in the presence of instrumental stray light.

    **Approximation flagged explicitly**: real Beer-Lambert-law deviations
    at high concentration have several distinct physical origins (a
    concentration-dependent chemical equilibrium shifting the true
    absorbing species' concentration, use of polychromatic rather than
    strictly monochromatic radiation, and instrumental stray light);
    this function models only the last of those, one specific,
    well-characterized, purely instrumental effect (Skoog, Holler &
    Crouch, *Principles of Instrumental Analysis*, 7th ed., Ch. 13.4).

    A small fraction :math:`p_s` of the light reaching the detector is
    "stray" radiation that never actually passed through the fully
    absorbing sample path (e.g. light scattered around the sample cell),
    so it is *not* attenuated by the sample the way the main beam is.
    The detector reads the combined transmitted-plus-stray intensity, so
    the apparent absorbance is

    .. math::

        A_{\text{apparent}} = -\log_{10}\!\left(10^{-A_{\text{true}}} + p_s\right)

    (normalized so that at `stray_light_fraction=0`, or in the dilute
    limit :math:`A_{\text{true}}\to0`, this reduces exactly to the true
    Beer-Lambert absorbance). Because the stray-light term does not
    shrink as the sample absorbs more strongly, it dominates the
    detector reading at high true absorbance, causing the classic
    *negative* (sublinear, "rolling over") deviation from Beer-Lambert
    linearity seen in real instruments at high concentration.

    Parameters
    ----------
    molar_absorptivity : float
        In L mol^-1 cm^-1.
    concentration : float
        In mol/L.
    path_length : float
        In cm.
    stray_light_fraction : float
        Fraction of incident intensity that reaches the detector as
        stray light, :math:`p_s\in[0,1)`.

    Returns
    -------
    float
        Apparent absorbance, always :math:`\le A_{\text{true}}`.

    Raises
    ------
    ValueError
        If `stray_light_fraction` is not in :math:`[0, 1)`.

    Examples
    --------
    With no stray light, the apparent absorbance is exactly the true one:

    >>> round(apparent_absorbance_with_stray_light(5000.0, 2.0e-5, 1.0, stray_light_fraction=0.0), 6)
    0.1

    At low concentration (low true absorbance) the deviation is
    negligible, but it grows sharply at high concentration -- Beer-
    Lambert linearity breaks down exactly where the textbook warns it
    does:

    >>> low = apparent_absorbance_with_stray_light(5000.0, 1.0e-6, 1.0, stray_light_fraction=0.001)
    >>> low_true = absorbance(5000.0, 1.0e-6, 1.0)
    >>> high = apparent_absorbance_with_stray_light(5000.0, 1.0e-3, 1.0, stray_light_fraction=0.001)
    >>> high_true = absorbance(5000.0, 1.0e-3, 1.0)
    >>> bool(abs(low - low_true) < abs(high - high_true))
    True
    """
    if not (0.0 <= stray_light_fraction < 1.0):
        raise ValueError("stray_light_fraction must be in [0, 1)")
    true_absorbance = absorbance(molar_absorptivity, concentration, path_length)
    return float(-np.log10(10.0**-true_absorbance + stray_light_fraction))
