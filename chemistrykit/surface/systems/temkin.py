r"""The Temkin adsorption isotherm, and its origin in a spread of site energies.

M. I. Temkin and V. Pyzhev, "Kinetics of Ammonia Synthesis on Promoted
Iron Catalysts," *Acta Physicochim. URSS* 12 (1940), 327. Temkin replaced
Langmuir's single site energy with a heat of adsorption that *falls
linearly* with coverage. Equivalently, the surface carries a uniform
spread of Langmuir sites whose adsorption energies span a width
:math:`fRT`, so that :math:`\ln K` is uniformly distributed between
:math:`\ln K_{max} - f` and :math:`\ln K_{max}`. Averaging the Langmuir
coverage over that spread gives exactly

.. math::

    \theta = \frac{1}{f}\ln\frac{1+K_{max}P}{1+K_{max}e^{-f}P}

(:func:`uniform_energy_coverage`), which in the broad middle-coverage
range :math:`K_{max}e^{-f}P \ll 1 \ll K_{max}P` collapses to the
logarithmic Temkin isotherm :math:`\theta \approx (1/f)\ln(K_{max}P)`.
In the loading form commonly fitted to data (Atkins & de Paula, *Physical
Chemistry*, 11th ed.),

.. math::

    q = \frac{RT}{b_T}\ln(A_TP),

a straight line in :math:`\ln P` recovered by :func:`fit_temkin`.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from chemistrykit.constants import R
from chemistrykit.surface.core.base_system import AdsorptionIsotherm
from chemistrykit.surface.utils.regression import linear_fit

__all__ = ["uniform_energy_coverage", "temkin_loading", "TemkinIsotherm", "TemkinFit", "fit_temkin"]


def uniform_energy_coverage(K_max: float, f: float, P):
    r"""Langmuir coverage averaged over a uniform spread of site energies of width :math:`fRT`.

    Parameters
    ----------
    K_max : float
        Langmuir constant of the most strongly binding sites.
    f : float
        Dimensionless width of the site-energy spread, :math:`\Delta Q/RT` (> 0).
    P : float or array-like of float
        Equilibrium pressure.

    Returns
    -------
    float or ndarray
        Mean fractional coverage in :math:`[0, 1)`.

    Examples
    --------
    As the spread vanishes the surface is uniform again and the Langmuir
    coverage :math:`KP/(1+KP)` is recovered:

    >>> round(float(uniform_energy_coverage(K_max=2.0, f=1e-8, P=0.5)), 6)
    0.5

    In the middle-coverage range the result is Temkin's logarithm:

    >>> K_max, f, P = 1.0e6, 20.0, 1.0e-2
    >>> theta = uniform_energy_coverage(K_max, f, P)
    >>> bool(abs(theta - np.log(K_max * P) / f) < 1e-3)
    True
    """
    if f <= 0:
        raise ValueError("f must be positive")
    P = np.asarray(P, dtype=np.float64)
    result = (np.log1p(K_max * P) - np.log1p(K_max * np.exp(-f) * P)) / f
    return float(result) if result.ndim == 0 else result


def temkin_loading(A_T: float, b_T: float, P, T: float, R_gas: float = R):
    r"""The Temkin isotherm :math:`q = (RT/b_T)\ln(A_TP)`.

    Valid only in the middle-coverage range; it turns negative for
    :math:`A_TP<1`, where the full expression
    :func:`uniform_energy_coverage` should be used instead.

    Parameters
    ----------
    A_T : float
        Temkin equilibrium binding constant, inverse to `P`'s units.
    b_T : float
        Temkin heat-of-adsorption constant, in J/mol per unit loading.
    P : float or array-like of float
    T : float
        Absolute temperature, in K.
    R_gas : float, default :data:`chemistrykit.constants.R`

    Returns
    -------
    float or ndarray

    Examples
    --------
    The loading vanishes at :math:`P=1/A_T`, and each factor of *e* in
    pressure adds exactly :math:`RT/b_T`:

    >>> round(float(temkin_loading(A_T=4.0, b_T=1.0e3, P=0.25, T=300.0)), 10)
    0.0
    >>> q1 = temkin_loading(4.0, 1.0e3, 1.0, 300.0)
    >>> q2 = temkin_loading(4.0, 1.0e3, np.e, 300.0)
    >>> round(float((q2 - q1) / (8.31446261815324 * 300.0 / 1.0e3)), 10)
    1.0
    """
    P = np.asarray(P, dtype=np.float64)
    result = R_gas * T / b_T * np.log(A_T * P)
    return float(result) if result.ndim == 0 else result


class TemkinIsotherm(AdsorptionIsotherm):
    r"""The Temkin isotherm :math:`q = (RT/b_T)\ln(A_TP)` at fixed temperature.

    Parameters
    ----------
    A_T : float
        Temkin binding constant.
    b_T : float
        Temkin heat-of-adsorption constant, in J/mol.
    T : float
        Absolute temperature, in K.

    Examples
    --------
    >>> round(float(TemkinIsotherm(A_T=2.0, b_T=500.0, T=300.0).loading(0.5)), 10)
    0.0
    """

    def __init__(self, A_T: float, b_T: float, T: float):
        if A_T <= 0 or b_T <= 0 or T <= 0:
            raise ValueError("A_T, b_T and T must all be positive")
        self.A_T = float(A_T)
        self.b_T = float(b_T)
        self.T = float(T)

    def loading(self, P):
        return temkin_loading(self.A_T, self.b_T, P, self.T)

    def fractional_coverage(self, P):
        raise NotImplementedError("the logarithmic Temkin isotherm has no saturation loading to normalize by")


@dataclass
class TemkinFit:
    """Result of fitting loading-vs-pressure data to the Temkin isotherm."""

    A_T: float
    """float: Fitted Temkin binding constant."""

    b_T: float
    """float: Fitted Temkin heat-of-adsorption constant, in J/mol."""

    r_squared: float
    """float: Coefficient of determination of the ``q`` vs ``ln P`` fit."""


def fit_temkin(P, q, T: float, R_gas: float = R) -> TemkinFit:
    r"""Fit data to the Temkin isotherm via :math:`q = (RT/b_T)\ln A_T + (RT/b_T)\ln P`.

    Parameters
    ----------
    P, q : array-like of float
        Equilibrium pressures and loadings in the middle-coverage range.
    T : float
        Absolute temperature, in K.
    R_gas : float, default :data:`chemistrykit.constants.R`

    Returns
    -------
    TemkinFit

    Examples
    --------
    >>> P = np.array([0.5, 1.0, 2.0, 5.0, 10.0])
    >>> q = temkin_loading(A_T=3.0, b_T=800.0, P=P, T=300.0)
    >>> fit = fit_temkin(P, q, T=300.0)
    >>> round(fit.A_T, 6), round(fit.b_T, 6)
    (3.0, 800.0)
    """
    P = np.asarray(P, dtype=np.float64)
    q = np.asarray(q, dtype=np.float64)
    fit = linear_fit(np.log(P), q)
    b_T = R_gas * T / fit.slope
    A_T = np.exp(fit.intercept / fit.slope)
    return TemkinFit(A_T=float(A_T), b_T=float(b_T), r_squared=fit.r_squared)
