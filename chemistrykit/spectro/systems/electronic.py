r"""Electronic (UV-Vis) spectra: Franck-Condon vibronic progressions.

See Atkins & de Paula, *Physical Chemistry*, 11th ed., Ch. 12.7(b), or
Hollas, *Modern Spectroscopy*, 4th ed., Ch. 6.1(c), throughout.

The Franck-Condon principle says an electronic transition is so much
faster than nuclear motion that the nuclei are effectively frozen during
it (J. Franck, *Trans. Faraday Soc.* 21, 536 (1925); E. U. Condon, *Phys.
Rev.* 32, 858 (1928)), so the intensity of the :math:`v''=0\to v'`
vibronic transition (ground vibrational level of the lower electronic
state to level :math:`v'` of the upper one) is governed by the overlap
of the two states' vibrational wavefunctions, the Franck-Condon factor
:math:`|\langle\chi_{v'}|\chi_{v''=0}\rangle|^2`. For two *displaced
harmonic oscillators of equal frequency* (the standard first
approximation -- the excited-state potential has the same curvature as
the ground state's, just shifted by :math:`\Delta Q` along the normal
coordinate), this overlap integral has the closed-form Poisson
distribution

.. math::

    \text{FC}(0\to v') = \frac{e^{-S}S^{v'}}{v'!}

where :math:`S=\frac12\frac{m\omega}{\hbar}\Delta Q^2` is the
dimensionless Huang-Rhys displacement parameter (K. Huang & A. Rhys,
*Proc. R. Soc. Lond. A* 204, 406 (1950)) -- flagged explicitly as an
approximation: unequal ground/excited-state frequencies (a real
excited-state potential is rarely exactly as "stiff" as the ground
state's) require the more general Duschinsky-rotation treatment, not
implemented here.
"""

from __future__ import annotations

import numpy as np
from scipy.special import factorial

from chemistrykit.constants import HBAR
from chemistrykit.spectro.core.base_system import Spectrum

__all__ = ["huang_rhys_factor", "franck_condon_factor", "franck_condon_progression", "franck_condon_spectrum"]


def huang_rhys_factor(displacement: float, mass: float, angular_frequency: float) -> float:
    r"""The dimensionless Huang-Rhys displacement parameter :math:`S=\frac12\frac{m\omega}{\hbar}\Delta Q^2`.

    Parameters
    ----------
    displacement : float
        Excited-minus-ground-state equilibrium displacement along the
        normal coordinate, :math:`\Delta Q`, in m (for a mass-weighted
        coordinate) or the natural length unit of the coordinate used.
    mass : float
        Oscillator (reduced) mass, in kg.
    angular_frequency : float
        Shared vibrational angular frequency :math:`\omega`, in rad/s.

    Returns
    -------
    float

    Examples
    --------
    >>> round(huang_rhys_factor(displacement=0.0, mass=1.6e-27, angular_frequency=5.0e13), 6)
    0.0
    """
    return float(0.5 * mass * angular_frequency / HBAR * displacement**2)


def franck_condon_factor(v: int, S: float) -> float:
    r"""The :math:`0\to v` Franck-Condon factor for two displaced, equal-frequency harmonic oscillators.

    .. math::

        \text{FC}(0\to v) = \frac{e^{-S}S^v}{v!}

    (a Poisson distribution in `v` with mean `S`; see the module
    docstring for the underlying approximation).

    Parameters
    ----------
    v : int
        Final-state vibrational quantum number, >= 0.
    S : float
        Huang-Rhys parameter, >= 0.

    Returns
    -------
    float
        In :math:`[0, 1]`.

    Examples
    --------
    Zero displacement (`S=0`) puts all intensity in the `0->0` transition
    (no vibronic structure at all):

    >>> round(franck_condon_factor(0, S=0.0), 6)
    1.0
    >>> round(franck_condon_factor(1, S=0.0), 6)
    0.0
    """
    if v < 0:
        raise ValueError("v must be >= 0")
    if S < 0:
        raise ValueError("S must be >= 0")
    return float(np.exp(-S) * S**v / factorial(v))


def franck_condon_progression(v_max: int, S: float) -> np.ndarray:
    r"""Tabulate the Franck-Condon factors for :math:`v=0,1,\dots,v_{max}`.

    Parameters
    ----------
    v_max : int
    S : float
        Huang-Rhys parameter.

    Returns
    -------
    ndarray, shape (v_max + 1,)

    Examples
    --------
    The Franck-Condon factors are (approximately) a complete probability
    distribution over the final vibrational level, so they sum to
    (approximately) 1 once `v_max` comfortably exceeds the Poisson
    distribution's mean `S` -- the required normalization check for any
    genuine Franck-Condon calculation:

    >>> round(float(np.sum(franck_condon_progression(v_max=40, S=3.0))), 9)
    1.0

    The most probable (highest-intensity) vibronic transition lands near
    :math:`v'\approx S` -- the qualitative rule of thumb for reading off
    a vibronic progression's displacement from its band shape:

    >>> progression = franck_condon_progression(v_max=20, S=5.0)
    >>> bool(abs(int(np.argmax(progression)) - 5) <= 1)
    True
    """
    v_values = np.arange(0, v_max + 1)
    return np.array([franck_condon_factor(int(v), S) for v in v_values])


def franck_condon_spectrum(origin_wavenumber: float, vibrational_wavenumber: float, S: float, v_max: int) -> Spectrum:
    r"""Build a vibronic-progression :class:`~chemistrykit.spectro.core.base_system.Spectrum`.

    Stick positions are the electronic origin (the :math:`0\to0`
    transition) plus `v'` quanta of the *excited-state* vibrational
    wavenumber; stick intensities are the corresponding Franck-Condon
    factors (Atkins & de Paula, *Physical Chemistry*, 11th ed., Ch.
    12.7(b)).

    Parameters
    ----------
    origin_wavenumber : float
        The :math:`0\to0` (electronic origin) transition wavenumber, in
        cm^-1.
    vibrational_wavenumber : float
        Excited-state vibrational spacing, in cm^-1.
    S : float
        Huang-Rhys parameter.
    v_max : int
        Highest excited-state vibrational level to include.

    Returns
    -------
    Spectrum
        `positions` in cm^-1, `intensities` the Franck-Condon factors.

    Examples
    --------
    >>> spectrum = franck_condon_spectrum(origin_wavenumber=25000.0, vibrational_wavenumber=1500.0, S=1.5, v_max=10)
    >>> round(float(spectrum.positions[0]), 1)
    25000.0
    >>> round(float(spectrum.positions[2] - spectrum.positions[1]), 1)
    1500.0
    """
    v_values = np.arange(0, v_max + 1)
    positions = origin_wavenumber + v_values * vibrational_wavenumber
    intensities = franck_condon_progression(v_max, S)
    labels = [f"0->{v}" for v in v_values]
    return Spectrum(positions=positions, intensities=intensities, labels=labels)
