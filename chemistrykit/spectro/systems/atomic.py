r"""Hydrogen-like atomic line positions: the Balmer-Rydberg formula.

For a one-electron atom of nuclear charge `Z`, Bohr's 1913 model (and,
identically for the level energies, the nonrelativistic Schrodinger
equation) gives the wavenumber of the :math:`n_\text{upper}\to
n_\text{lower}` emission line as

.. math::

    \tilde\nu = Z^2 R_M\left(\frac{1}{n_\text{lower}^2}-\frac{1}{n_\text{upper}^2}\right),
    \qquad R_M=\frac{R_\infty}{1+m_e/M},

the empirical formula J. R. Rydberg generalized (1890) from J. J.
Balmer's 1885 fit to the four visible hydrogen lines
(:math:`n_\text{lower}=2`). :math:`R_\infty` is the infinite-nuclear-mass
Rydberg constant and :math:`R_M` its reduced-mass-corrected value for a
nucleus of mass `M`. See Atkins & de Paula, *Physical Chemistry*, 11th
ed., Ch. 9A (hydrogenic atoms) throughout.
"""

from __future__ import annotations

import numpy as np
import scipy.constants as sc

from chemistrykit.constants import ELECTRON_MASS

__all__ = ["rydberg_wavenumber"]


def rydberg_wavenumber(n_lower, n_upper, nuclear_charge: int = 1, nuclear_mass: float | None = None):
    r"""Wavenumber of a hydrogen-like atom's :math:`n_\text{upper}\to n_\text{lower}` line, in cm^-1.

    Parameters
    ----------
    n_lower : int or array-like of int
        Principal quantum number of the lower level (1 = Lyman, 2 =
        Balmer, 3 = Paschen, ...).
    n_upper : int or array-like of int
        Principal quantum number of the upper level (must exceed
        `n_lower`).
    nuclear_charge : int, default 1
        Nuclear charge `Z`.
    nuclear_mass : float, optional
        Nuclear mass `M`, in kg, for the reduced-mass correction
        :math:`R_M=R_\infty/(1+m_e/M)`. ``None`` uses the
        infinite-nuclear-mass constant :math:`R_\infty`.

    Returns
    -------
    float or ndarray
        Vacuum wavenumber, in cm^-1 (positive).

    Examples
    --------
    Balmer's H-alpha line (:math:`3\to2`) of ordinary hydrogen, at a
    vacuum wavelength of 656.47 nm:

    >>> import scipy.constants as sc
    >>> nu = rydberg_wavenumber(2, 3, nuclear_mass=sc.m_p)
    >>> round(1e7 / nu, 2)
    656.47

    The series limit (:math:`n_\text{upper}\to\infty`) of the Lyman
    series is the ionization energy, :math:`R_\infty` itself for an
    infinitely heavy nucleus:

    >>> round(rydberg_wavenumber(1, 10**6) / (sc.Rydberg / 100.0), 9)
    1.0
    """
    n_lo = np.asarray(n_lower, dtype=np.float64)
    n_up = np.asarray(n_upper, dtype=np.float64)
    if np.any(n_lo < 1) or np.any(n_up <= n_lo):
        raise ValueError("require 1 <= n_lower < n_upper")
    rydberg = sc.Rydberg / 100.0  # cm^-1
    if nuclear_mass is not None:
        if nuclear_mass <= 0:
            raise ValueError("nuclear_mass must be positive")
        rydberg = rydberg / (1.0 + ELECTRON_MASS / nuclear_mass)
    result = nuclear_charge**2 * rydberg * (1.0 / n_lo**2 - 1.0 / n_up**2)
    return float(result) if result.ndim == 0 else result
