r"""Schottky and Frenkel point-defect equilibria: Boltzmann-factor defect concentration vs. temperature.

Point defects (vacancies, interstitials) are present in any ionic
crystal in thermodynamic equilibrium at :math:`T>0`, because although
forming a defect costs enthalpy :math:`\Delta H`, it also increases the
crystal's configurational entropy (there are many equivalent ways to
place a handful of defects among many lattice sites). Minimizing the
Gibbs energy :math:`\Delta G=\Delta H-T\Delta S_{config}` with
:math:`\Delta S_{config}=k_B\ln\binom{N}{n}` (Stirling's approximation,
:math:`n\ll N`) over the number of defects `n` gives, in both cases
below, a defect population governed by a Boltzmann factor in
:math:`\Delta H/2` -- the factor of 2 arising because forming *one*
Schottky or Frenkel defect actually requires two independent
"choices" (which cation site is vacated, and independently which
interstitial/anion site is involved), so the degeneracy (and hence the
entropy, and hence the equilibrium condition) is effectively squared
relative to a single-site defect. See West, *Solid State Chemistry and
its Applications*, 2nd ed. (2014), Ch. 1.4, or Kittel, *Introduction to
Solid State Physics*, 8th ed., Ch. 18, for the full derivation.
"""

from __future__ import annotations

import numpy as np

from chemistrykit.constants import K_B

__all__ = ["schottky_defect_concentration", "frenkel_defect_concentration"]


def schottky_defect_concentration(N: float, delta_h: float, T):
    r"""Equilibrium number of Schottky defects (paired cation+anion vacancies).

    .. math::

        n_S = N\exp\!\left(-\frac{\Delta H_S}{2k_BT}\right)

    valid for :math:`n_S\ll N` (West, *Solid State Chemistry and its
    Applications*, 2nd ed., Ch. 1.4, eq. 1.10).

    Parameters
    ----------
    N : float
        Number of cation (equivalently, formula-unit) sites in the crystal.
    delta_h : float
        Enthalpy of formation of one Schottky defect (one cation vacancy
        + one anion vacancy), in J.
    T : float or array-like of float
        Absolute temperature(s), in K.

    Returns
    -------
    float or ndarray
        Number of Schottky defects, :math:`n_S \le N`.

    Examples
    --------
    Defect concentration increases with temperature (more thermal energy
    to pay the formation enthalpy):

    >>> import numpy as np
    >>> T = np.array([300.0, 600.0, 900.0])
    >>> n = schottky_defect_concentration(N=1e20, delta_h=2.0e-19, T=T)
    >>> bool(np.all(np.diff(n) > 0))
    True

    At very low temperature, essentially no defects form:

    >>> round(float(schottky_defect_concentration(N=1e20, delta_h=2.0e-19, T=1.0)), 6)
    0.0
    """
    T = np.asarray(T, dtype=np.float64)
    n = N * np.exp(-delta_h / (2.0 * K_B * T))
    return float(n) if n.ndim == 0 else n


def frenkel_defect_concentration(N: float, N_interstitial: float, delta_h: float, T):
    r"""Equilibrium number of Frenkel defects (a cation displaced to an interstitial site).

    .. math::

        n_F = \sqrt{NN_i}\exp\!\left(-\frac{\Delta H_F}{2k_BT}\right)

    valid for :math:`n_F\ll N,N_i` (West, *Solid State Chemistry and its
    Applications*, 2nd ed., Ch. 1.4, eq. 1.14); reduces to the Schottky
    form with :math:`N_i\to N` up to the geometric- vs. arithmetic-mean
    site count.

    Parameters
    ----------
    N : float
        Number of normal cation lattice sites.
    N_interstitial : float
        Number of available interstitial sites.
    delta_h : float
        Enthalpy of formation of one Frenkel defect (one vacancy + one
        interstitial), in J.
    T : float or array-like of float
        Absolute temperature(s), in K.

    Returns
    -------
    float or ndarray

    Examples
    --------
    >>> import numpy as np
    >>> T = np.array([300.0, 600.0, 900.0])
    >>> n = frenkel_defect_concentration(N=1e20, N_interstitial=1e20, delta_h=2.5e-19, T=T)
    >>> bool(np.all(np.diff(n) > 0))
    True

    With equal site counts N=N_i, this reduces to the Schottky formula's
    functional form (same exponential factor, prefactor N instead of
    sqrt(N*N)=N -- they coincide):

    >>> from chemistrykit.crystal.systems.defects import schottky_defect_concentration
    >>> n_frenkel = frenkel_defect_concentration(N=1e20, N_interstitial=1e20, delta_h=2.5e-19, T=500.0)
    >>> n_schottky = schottky_defect_concentration(N=1e20, delta_h=2.5e-19, T=500.0)
    >>> bool(abs(n_frenkel - n_schottky) / n_schottky < 1e-12)
    True
    """
    T = np.asarray(T, dtype=np.float64)
    n = np.sqrt(N * N_interstitial) * np.exp(-delta_h / (2.0 * K_B * T))
    return float(n) if n.ndim == 0 else n
