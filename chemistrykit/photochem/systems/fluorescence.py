r"""Steady-state fluorescence observables: the Stokes shift and the Perrin anisotropy equation.

* :func:`stokes_shift` -- the energy gap between absorption and emission
  maxima (G. G. Stokes, *Phil. Trans. R. Soc. Lond.* 142, 463 (1852)).
* :func:`perrin_anisotropy` and :func:`rotational_correlation_time` --
  the steady-state fluorescence anisotropy of a rotating fluorophore
  (F. Perrin, *J. Phys. Radium* 7, 390 (1926)) with the
  Stokes-Einstein-Debye rotational correlation time. See Lakowicz,
  *Principles of Fluorescence Spectroscopy*, 3rd ed., Chs. 1 and 10.
"""

from __future__ import annotations

import numpy as np

from chemistrykit.constants import K_B

__all__ = ["stokes_shift", "perrin_anisotropy", "rotational_correlation_time"]


def stokes_shift(absorption_max_nm, emission_max_nm):
    r"""Stokes shift in wavenumbers, :math:`\Delta\tilde\nu = 10^7/\lambda_{abs} - 10^7/\lambda_{em}` (cm\ :sup:`-1`).

    Stokes (1852) observed that fluorescence is emitted at longer
    wavelength (lower energy) than the light that excites it; the energy
    difference is lost to vibrational relaxation and solvent
    reorganization in the excited state before emission (Lakowicz, Ch.
    1.4).

    Parameters
    ----------
    absorption_max_nm : float or array-like of float
        Wavelength of the absorption maximum, in nm.
    emission_max_nm : float or array-like of float
        Wavelength of the emission maximum, in nm.

    Returns
    -------
    float or ndarray
        Stokes shift in cm\ :sup:`-1` (positive for red-shifted emission).

    Examples
    --------
    Absorption at 400 nm and emission at 500 nm:

    >>> round(stokes_shift(400.0, 500.0), 6)
    5000.0
    """
    lam_a = np.asarray(absorption_max_nm, dtype=np.float64)
    lam_e = np.asarray(emission_max_nm, dtype=np.float64)
    result = 1.0e7 / lam_a - 1.0e7 / lam_e
    return float(result) if result.ndim == 0 else result


def rotational_correlation_time(viscosity, volume, T):
    r"""Stokes-Einstein-Debye rotational correlation time :math:`\theta = \eta V/(k_B T)`.

    Parameters
    ----------
    viscosity : float or array-like of float
        Solvent viscosity :math:`\eta`, in Pa*s.
    volume : float
        Hydrodynamic volume of the rotating molecule, in m\ :sup:`3`.
    T : float or array-like of float
        Temperature, in K.

    Returns
    -------
    float or ndarray
        :math:`\theta`, in s.

    Examples
    --------
    A ~1 nm\ :sup:`3` fluorophore in water (1 mPa*s) at 298 K rotates in
    about a quarter of a nanosecond:

    >>> round(rotational_correlation_time(1.0e-3, 1.0e-27, 298.15) * 1e9, 3)
    0.243
    """
    result = np.asarray(viscosity, dtype=np.float64) * volume / (K_B * np.asarray(T, dtype=np.float64))
    return float(result) if result.ndim == 0 else result


def perrin_anisotropy(r0: float, tau, theta):
    r"""Perrin equation for steady-state fluorescence anisotropy, :math:`r = r_0/(1+\tau/\theta)`.

    A fluorophore excited by polarized light emits polarized light, but
    rotational diffusion during the excited-state lifetime :math:`\tau`
    scrambles the orientation; for a spherical rotor with rotational
    correlation time :math:`\theta` the time-averaged anisotropy is
    :math:`r_0/r = 1+\tau/\theta` (Perrin 1926; Lakowicz, Ch. 10, eq.
    10.45). :math:`r_0` (at most 0.4 for one-photon excitation) is the
    anisotropy in the absence of rotation.

    Parameters
    ----------
    r0 : float
        Fundamental (rotation-free) anisotropy.
    tau : float or array-like of float
        Fluorescence lifetime.
    theta : float or array-like of float
        Rotational correlation time (same unit as `tau`).

    Returns
    -------
    float or ndarray

    Examples
    --------
    When the lifetime equals the correlation time, the anisotropy is
    halved:

    >>> round(perrin_anisotropy(0.4, tau=4.0, theta=4.0), 6)
    0.2
    """
    result = r0 / (1.0 + np.asarray(tau, dtype=np.float64) / np.asarray(theta, dtype=np.float64))
    return float(result) if result.ndim == 0 else result
