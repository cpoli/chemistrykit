r"""Polanyi's adsorption potential and the Dubinin-Radushkevich isotherm.

M. Polanyi, *Verh. Dtsch. Phys. Ges.* 16 (1914), 1012; M. M. Dubinin and
L. V. Radushkevich, *Proc. Acad. Sci. USSR, Phys. Chem. Sect.* 55 (1947),
331. Polanyi described adsorption by the *adsorption potential*

.. math::

    A = RT\ln\frac{P_0}{P},

the isothermal work of compressing vapor from its equilibrium pressure
:math:`P` to saturation :math:`P_0`, and postulated that the filled
adsorption volume :math:`W` depends on :math:`A` alone, through a
temperature-independent *characteristic curve* :math:`W(A)`. Dubinin and
Radushkevich gave that curve an explicit form for micropore filling in
activated carbons,

.. math::

    W = W_0\exp\!\left[-\left(\frac{A}{E}\right)^2\right],

with :math:`W_0` the limiting micropore volume and :math:`E` a
characteristic energy (Gregg & Sing, *Adsorption, Surface Area and
Porosity*, 2nd ed., Ch. 4). Because :math:`\ln W` is linear in
:math:`A^2`, :func:`fit_dubinin_radushkevich` recovers :math:`(W_0, E)`
by the same shared least-squares routine as every other linearization in
this package.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from chemistrykit.constants import R
from chemistrykit.surface.core.base_system import AdsorptionIsotherm
from chemistrykit.surface.utils.regression import linear_fit

__all__ = [
    "polanyi_potential",
    "dubinin_radushkevich_loading",
    "DubininRadushkevichIsotherm",
    "DubininRadushkevichFit",
    "fit_dubinin_radushkevich",
]


def polanyi_potential(P, P0: float, T: float, R_gas: float = R):
    r"""Polanyi adsorption potential :math:`A = RT\ln(P_0/P)`, in J/mol.

    Pressures above `P0` (bulk condensation) are clipped to :math:`A=0`.

    Parameters
    ----------
    P : float or array-like of float
        Equilibrium pressure (same units as `P0`).
    P0 : float
        Saturation vapor pressure at temperature `T`.
    T : float
        Absolute temperature, in K.
    R_gas : float, default :data:`chemistrykit.constants.R`

    Returns
    -------
    float or ndarray

    Examples
    --------
    At saturation the potential vanishes; at :math:`P=P_0/e` it equals
    :math:`RT`:

    >>> polanyi_potential(P=1.0, P0=1.0, T=77.0)
    0.0
    >>> round(polanyi_potential(P=np.exp(-1.0), P0=1.0, T=100.0, R_gas=8.0), 10)
    800.0
    """
    P = np.asarray(P, dtype=np.float64)
    result = R_gas * T * np.log(P0 / P)
    result = np.maximum(result, 0.0)
    return float(result) if result.ndim == 0 else result


def dubinin_radushkevich_loading(P, W0: float, E: float, P0: float, T: float, R_gas: float = R):
    r"""Dubinin-Radushkevich micropore filling :math:`W = W_0\exp[-(A/E)^2]`.

    Parameters
    ----------
    P : float or array-like of float
        Equilibrium pressure.
    W0 : float
        Limiting micropore volume (or loading) at :math:`P=P_0`.
    E : float
        Characteristic energy, in J/mol.
    P0 : float
        Saturation vapor pressure at `T`.
    T : float
        Absolute temperature, in K.
    R_gas : float, default :data:`chemistrykit.constants.R`

    Returns
    -------
    float or ndarray

    Examples
    --------
    When the adsorption potential equals the characteristic energy
    (:math:`A=E`), the pores are filled to exactly :math:`W_0/e`:

    >>> T, E = 300.0, 10.0e3
    >>> P = 1.0 * np.exp(-E / (8.31446261815324 * T))
    >>> round(float(dubinin_radushkevich_loading(P, W0=2.0, E=E, P0=1.0, T=T) * np.e), 8)
    2.0
    """
    A = np.asarray(polanyi_potential(P, P0, T, R_gas), dtype=np.float64)
    result = W0 * np.exp(-((A / E) ** 2))
    return float(result) if result.ndim == 0 else result


class DubininRadushkevichIsotherm(AdsorptionIsotherm):
    r"""The Dubinin-Radushkevich isotherm at a fixed temperature.

    Parameters
    ----------
    W0 : float
        Limiting micropore volume (saturation loading).
    E : float
        Characteristic energy, in J/mol.
    P0 : float
        Saturation vapor pressure at `T`.
    T : float
        Absolute temperature, in K.

    Examples
    --------
    >>> iso = DubininRadushkevichIsotherm(W0=0.5, E=12e3, P0=1.0, T=300.0)
    >>> round(float(iso.fractional_coverage(1.0)), 10)
    1.0
    """

    def __init__(self, W0: float, E: float, P0: float, T: float):
        if W0 <= 0 or E <= 0 or P0 <= 0 or T <= 0:
            raise ValueError("W0, E, P0 and T must all be positive")
        self.W0 = float(W0)
        self.E = float(E)
        self.P0 = float(P0)
        self.T = float(T)
        self.saturation_loading = self.W0

    def loading(self, P):
        return dubinin_radushkevich_loading(P, self.W0, self.E, self.P0, self.T)


@dataclass
class DubininRadushkevichFit:
    """Result of fitting loading-vs-pressure data to the Dubinin-Radushkevich isotherm."""

    W0: float
    """float: Fitted limiting micropore volume."""

    E: float
    """float: Fitted characteristic energy, in J/mol."""

    r_squared: float
    """float: Coefficient of determination of the linearized (``ln W`` vs ``A^2``) fit."""


def fit_dubinin_radushkevich(P, W, P0: float, T: float, R_gas: float = R) -> DubininRadushkevichFit:
    r"""Fit data to the Dubinin-Radushkevich isotherm via :math:`\ln W = \ln W_0 - A^2/E^2`.

    Parameters
    ----------
    P : array-like of float
        Equilibrium pressures, all below `P0`.
    W : array-like of float
        Corresponding adsorbed volumes (all positive).
    P0 : float
        Saturation vapor pressure at `T`.
    T : float
        Absolute temperature, in K.
    R_gas : float, default :data:`chemistrykit.constants.R`

    Returns
    -------
    DubininRadushkevichFit

    Examples
    --------
    >>> P = np.array([1e-4, 1e-3, 1e-2, 0.05, 0.1, 0.3])
    >>> W = dubinin_radushkevich_loading(P, W0=0.4, E=9.0e3, P0=1.0, T=77.0)
    >>> fit = fit_dubinin_radushkevich(P, W, P0=1.0, T=77.0)
    >>> round(fit.W0, 6), round(fit.E, 3)
    (0.4, 9000.0)
    """
    A = np.asarray(polanyi_potential(P, P0, T, R_gas), dtype=np.float64)
    W = np.asarray(W, dtype=np.float64)
    fit = linear_fit(A**2, np.log(W))
    return DubininRadushkevichFit(W0=float(np.exp(fit.intercept)), E=float(1.0 / np.sqrt(-fit.slope)), r_squared=fit.r_squared)
