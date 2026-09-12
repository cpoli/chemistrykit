r"""Lorentzian, Gaussian, and Voigt lineshape functions, and a stick-spectrum broadener.

Shared machinery for every spectroscopic model in this subpackage
(:mod:`chemistrykit.spectro.systems.rotational`,
:mod:`chemistrykit.spectro.systems.vibrational`,
:mod:`chemistrykit.spectro.systems.electronic`) -- each produces a
discrete "stick spectrum" (peak positions and intensities), and this
module turns that into the continuous, broadened spectrum an actual
instrument would record. See Hollas, *Modern Spectroscopy*, 4th ed.
(2004), Ch. 2.3, or Demtroder, *Laser Spectroscopy*, 5th ed., Ch. 3, for
the physical origin of each shape: a Lorentzian from homogeneous
(lifetime/collisional) broadening, a Gaussian from inhomogeneous
(Doppler/instrumental) broadening, and a Voigt profile -- their
convolution -- when both mechanisms contribute together.
"""

from __future__ import annotations

import numpy as np
from scipy.special import wofz

__all__ = ["gaussian", "lorentzian", "voigt", "broaden_stick_spectrum"]


def gaussian(x, x0: float, fwhm: float):
    r"""A normalized Gaussian lineshape.

    .. math::

        g(x) = \frac{2}{\text{FWHM}}\sqrt{\frac{\ln 2}{\pi}}
               \exp\!\left[-4\ln 2\left(\frac{x-x_0}{\text{FWHM}}\right)^2\right]

    Normalized so :math:`\int_{-\infty}^{\infty}g(x)\,dx=1` (Hollas,
    *Modern Spectroscopy*, 4th ed., Ch. 2.3).

    Parameters
    ----------
    x : float or array-like of float
        Evaluation point(s).
    x0 : float
        Peak center.
    fwhm : float
        Full width at half maximum, same units as `x`.

    Returns
    -------
    float or ndarray

    Examples
    --------
    The peak height at `x0` and the normalization both follow the
    closed-form prefactor:

    >>> import numpy as np
    >>> expected_peak = 2.0 / 1.0 * np.sqrt(np.log(2) / np.pi)
    >>> round(float(gaussian(0.0, 0.0, fwhm=1.0)), 6) == round(float(expected_peak), 6)
    True
    >>> x = np.linspace(-20, 20, 200001)
    >>> round(float(np.trapezoid(gaussian(x, 0.0, fwhm=2.0), x)), 3)
    1.0
    """
    x = np.asarray(x, dtype=np.float64)
    prefactor = (2.0 / fwhm) * np.sqrt(np.log(2.0) / np.pi)
    return prefactor * np.exp(-4.0 * np.log(2.0) * ((x - x0) / fwhm) ** 2)


def lorentzian(x, x0: float, fwhm: float):
    r"""A normalized Lorentzian lineshape.

    .. math::

        L(x) = \frac{1}{\pi}\frac{\text{FWHM}/2}{(x-x_0)^2+(\text{FWHM}/2)^2}

    Normalized so :math:`\int_{-\infty}^{\infty}L(x)\,dx=1` (Hollas,
    *Modern Spectroscopy*, 4th ed., Ch. 2.3).

    Parameters
    ----------
    x : float or array-like of float
    x0 : float
        Peak center.
    fwhm : float
        Full width at half maximum, same units as `x`.

    Returns
    -------
    float or ndarray

    Examples
    --------
    >>> import numpy as np
    >>> x = np.linspace(-2000, 2000, 4_000_001)
    >>> round(float(np.trapezoid(lorentzian(x, 0.0, fwhm=2.0), x)), 3)
    1.0

    By definition, the value has dropped to half its peak at
    :math:`x_0\pm\text{FWHM}/2`:

    >>> peak = lorentzian(0.0, 0.0, fwhm=2.0)
    >>> half_max_point = lorentzian(1.0, 0.0, fwhm=2.0)
    >>> round(float(half_max_point), 9) == round(float(peak) / 2.0, 9)
    True
    """
    x = np.asarray(x, dtype=np.float64)
    half_width = fwhm / 2.0
    return (1.0 / np.pi) * half_width / ((x - x0) ** 2 + half_width**2)


def voigt(x, x0: float, fwhm_gaussian: float, fwhm_lorentzian: float):
    r"""A normalized Voigt profile: the convolution of a Gaussian and a Lorentzian.

    .. math::

        V(x) = \mathrm{Re}\!\left[\frac{w(z)}{\sigma\sqrt{2\pi}}\right],
        \qquad z=\frac{(x-x_0)+i\gamma}{\sigma\sqrt2}

    where `w` is the Faddeeva function (:func:`scipy.special.wofz`),
    :math:`\sigma=\text{FWHM}_G/(2\sqrt{2\ln2})` is the Gaussian standard
    deviation, and :math:`\gamma=\text{FWHM}_L/2` is the Lorentzian
    half-width -- the standard numerically stable evaluation of the Voigt
    profile (Hollas, *Modern Spectroscopy*, 4th ed., Ch. 2.3; J.
    Humlicek, *J. Quant. Spectrosc. Radiat. Transfer* 21, 309 (1979)).

    Parameters
    ----------
    x : float or array-like of float
    x0 : float
        Peak center.
    fwhm_gaussian : float
        Gaussian-component FWHM.
    fwhm_lorentzian : float
        Lorentzian-component FWHM.

    Returns
    -------
    float or ndarray

    Examples
    --------
    Normalized, and reduces to a pure Gaussian/Lorentzian in the
    appropriate limit:

    >>> import numpy as np
    >>> x = np.linspace(-1000, 1000, 200001)
    >>> round(float(np.trapezoid(voigt(x, 0.0, fwhm_gaussian=1.0, fwhm_lorentzian=1.0), x)), 2)
    1.0
    >>> x_narrow = np.linspace(-50, 50, 200001)
    >>> pure_lorentzian_limit = voigt(x_narrow, 0.0, fwhm_gaussian=1e-8, fwhm_lorentzian=2.0)
    >>> bool(np.allclose(pure_lorentzian_limit, lorentzian(x_narrow, 0.0, fwhm=2.0), atol=1e-3))
    True
    """
    x = np.asarray(x, dtype=np.float64)
    sigma = fwhm_gaussian / (2.0 * np.sqrt(2.0 * np.log(2.0)))
    gamma = fwhm_lorentzian / 2.0
    z = ((x - x0) + 1j * gamma) / (sigma * np.sqrt(2.0))
    return np.real(wofz(z)) / (sigma * np.sqrt(2.0 * np.pi))


def broaden_stick_spectrum(centers, intensities, x, shape: str = "gaussian", fwhm: float = 1.0, fwhm_lorentzian=None):
    """Sum a lineshape over every (center, intensity) stick to build a continuous spectrum.

    Parameters
    ----------
    centers : array-like of float
        Peak positions.
    intensities : array-like of float
        Relative intensity (area) of each peak, same length as `centers`.
    x : array-like of float
        Grid to evaluate the broadened spectrum on.
    shape : {"gaussian", "lorentzian", "voigt"}, default "gaussian"
        Lineshape to use for every peak.
    fwhm : float, default 1.0
        FWHM (or the Gaussian-component FWHM, for `shape="voigt"`).
    fwhm_lorentzian : float, optional
        Lorentzian-component FWHM, required (and only used) for
        `shape="voigt"`.

    Returns
    -------
    ndarray, shape matching `x`

    Raises
    ------
    ValueError
        For an unrecognized `shape`, or a missing `fwhm_lorentzian` when
        `shape="voigt"`.

    Examples
    --------
    A single stick recovers the plain lineshape function, scaled by its
    intensity:

    >>> import numpy as np
    >>> x = np.linspace(-10, 10, 501)
    >>> spectrum = broaden_stick_spectrum([0.0], [2.0], x, shape="gaussian", fwhm=1.0)
    >>> bool(np.allclose(spectrum, 2.0 * gaussian(x, 0.0, 1.0)))
    True
    """
    x = np.asarray(x, dtype=np.float64)
    centers = np.asarray(centers, dtype=np.float64)
    intensities = np.asarray(intensities, dtype=np.float64)
    spectrum = np.zeros_like(x)
    for center, intensity in zip(centers, intensities, strict=True):
        if shape == "gaussian":
            spectrum += intensity * gaussian(x, center, fwhm)
        elif shape == "lorentzian":
            spectrum += intensity * lorentzian(x, center, fwhm)
        elif shape == "voigt":
            if fwhm_lorentzian is None:
                raise ValueError('fwhm_lorentzian is required for shape="voigt"')
            spectrum += intensity * voigt(x, center, fwhm, fwhm_lorentzian)
        else:
            raise ValueError(f'shape must be "gaussian", "lorentzian", or "voigt", got {shape!r}')
    return spectrum
