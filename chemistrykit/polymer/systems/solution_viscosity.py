r"""Dilute-solution viscosity and molar mass: Staudinger's rule and the Mark-Houwink equation.

Staudinger (H. Staudinger & W. Heuer, *Ber. Dtsch. Chem. Ges.* 63, 222
(1930)) proposed that the specific viscosity of a dilute solution of long
chain molecules grows in direct proportion to their chain length,
:math:`\eta_\text{sp}/c = K_mM` -- viscosity as a molecular-weight probe,
and evidence that the chains are real molecules of definite length.
Mark (1938) and Houwink (*J. Prakt. Chem.* 157, 15 (1940)) generalized
this to the power law

.. math::

    [\eta] = KM^a

with an exponent :math:`0.5\le a\lesssim0.8` for flexible coils
(Staudinger's rule is the special case :math:`a=1`). The Flory-Fox
relation :math:`[\eta]\propto R^3/M` (P. J. Flory & T. G Fox, *J. Am.
Chem. Soc.* 73, 1904 (1951)) with :math:`R\propto M^\nu` links the exponent
to chain statistics, :math:`a=3\nu-1`: :math:`a=1/2` in a theta solvent,
:math:`a=4/5` in a good solvent (Flory's :math:`\nu=3/5`).
"""

from __future__ import annotations

import numpy as np

__all__ = [
    "staudinger_specific_viscosity",
    "mark_houwink_intrinsic_viscosity",
    "mark_houwink_exponent_from_flory",
    "fit_mark_houwink",
]


def staudinger_specific_viscosity(c, M, Km: float):
    r"""Staudinger's viscosity rule, :math:`\eta_\text{sp}=K_mcM`.

    Parameters
    ----------
    c : float or array-like of float
        Polymer concentration.
    M : float or array-like of float
        Molar mass (or degree of polymerization).
    Km : float
        Staudinger constant.

    Returns
    -------
    float or ndarray

    Examples
    --------
    Doubling the chain length doubles the specific viscosity:

    >>> float(staudinger_specific_viscosity(0.01, 2e5, 1e-4) / staudinger_specific_viscosity(0.01, 1e5, 1e-4))
    2.0
    """
    return Km * np.asarray(c, dtype=float) * np.asarray(M, dtype=float)


def mark_houwink_intrinsic_viscosity(M, K: float, a: float):
    r"""Mark-Houwink equation, :math:`[\eta]=KM^a`.

    Parameters
    ----------
    M : float or array-like of float
        Molar mass.
    K, a : float
        Mark-Houwink constants for a given polymer/solvent/temperature.

    Returns
    -------
    float or ndarray

    Examples
    --------
    >>> round(float(mark_houwink_intrinsic_viscosity(1e4, K=0.01, a=0.5)), 6)
    1.0
    """
    return K * np.asarray(M, dtype=float) ** a


def mark_houwink_exponent_from_flory(nu: float) -> float:
    r"""Mark-Houwink exponent implied by a Flory exponent via Flory-Fox, :math:`a=3\nu-1`.

    Parameters
    ----------
    nu : float
        Flory exponent (:math:`R\propto M^\nu`).

    Returns
    -------
    float

    Examples
    --------
    >>> mark_houwink_exponent_from_flory(0.5)
    0.5
    >>> round(mark_houwink_exponent_from_flory(0.6), 6)
    0.8
    """
    return 3.0 * nu - 1.0


def fit_mark_houwink(M, eta) -> tuple[float, float]:
    r"""Fit :math:`[\eta]=KM^a` by least squares on :math:`\log[\eta]` vs :math:`\log M`.

    Parameters
    ----------
    M, eta : array-like of float
        Molar masses and measured intrinsic viscosities.

    Returns
    -------
    K, a : float

    Examples
    --------
    >>> M = np.array([1e4, 1e5, 1e6])
    >>> K, a = fit_mark_houwink(M, 0.02 * M**0.7)
    >>> round(K, 6), round(a, 6)
    (0.02, 0.7)
    """
    a, logK = np.polyfit(np.log(np.asarray(M, dtype=float)), np.log(np.asarray(eta, dtype=float)), 1)
    return float(np.exp(logK)), float(a)
