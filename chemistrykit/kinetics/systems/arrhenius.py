"""The Arrhenius equation and activation-energy fitting.

See Atkins & de Paula, *Physical Chemistry*, 11th ed., Ch. 20.3 ("The
temperature dependence of reaction rates"), and the original: S.
Arrhenius, *Z. Phys. Chem.* 4, 226 (1889).
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from chemistrykit.constants import R
from chemistrykit.kinetics.utils.linear_regression import linear_fit

__all__ = ["arrhenius_rate_constant", "ArrheniusFit", "fit_arrhenius"]


def arrhenius_rate_constant(A: float, Ea: float, T, R_gas: float = R):
    r"""Evaluate the Arrhenius equation :math:`k = A e^{-E_a/(RT)}`.

    Parameters
    ----------
    A : float
        Pre-exponential ("frequency") factor, in the same units as `k`.
    Ea : float
        Activation energy, in J/mol.
    T : float or array-like of float
        Absolute temperature(s), in K.
    R_gas : float, default :data:`chemistrykit.constants.R`
        Gas constant, in J mol^-1 K^-1.

    Returns
    -------
    float or ndarray
        Rate constant(s) `k`, in the same units as `A`.

    Examples
    --------
    >>> round(float(arrhenius_rate_constant(A=1e13, Ea=50e3, T=300.0)), 2)
    19696.84
    """
    T = np.asarray(T, dtype=np.float64)
    return A * np.exp(-Ea / (R_gas * T))


@dataclass
class ArrheniusFit:
    """Result of fitting rate-constant-vs-temperature data to the Arrhenius equation."""

    Ea: float
    """float: Fitted activation energy, in J/mol."""

    A: float
    """float: Fitted pre-exponential factor."""

    r_squared: float
    """float: Coefficient of determination of the linear (ln k vs 1/T) fit."""

    R_gas: float = R
    """float: Gas constant used in the fit, in J mol^-1 K^-1."""

    def predict(self, T):
        """Evaluate the fitted Arrhenius equation at temperature(s) `T`.

        Parameters
        ----------
        T : float or array-like of float
            Absolute temperature(s), in K.

        Returns
        -------
        float or ndarray
        """
        return arrhenius_rate_constant(self.A, self.Ea, T, R_gas=self.R_gas)


def fit_arrhenius(T, k, R_gas: float = R) -> ArrheniusFit:
    r"""Fit rate constant vs. temperature data to the Arrhenius equation.

    Linearizes :math:`\ln k = \ln A - E_a/R \cdot (1/T)` and fits by
    ordinary least squares -- the standard "Arrhenius plot" method (Atkins
    & de Paula, *Physical Chemistry*, 11th ed., Ch. 20.3): plotting
    :math:`\ln k` against :math:`1/T` gives a straight line of slope
    :math:`-E_a/R` and intercept :math:`\ln A`.

    Parameters
    ----------
    T : array-like of float
        Absolute temperatures, in K (at least 2 distinct values).
    k : array-like of float
        Rate constants measured at each temperature in `T`, same units
        throughout.
    R_gas : float, default :data:`chemistrykit.constants.R`
        Gas constant, in J mol^-1 K^-1.

    Returns
    -------
    ArrheniusFit

    Examples
    --------
    Generate exact data from a known (Ea, A) and recover them:

    >>> import numpy as np
    >>> T = np.array([280.0, 300.0, 320.0, 340.0, 360.0])
    >>> k = arrhenius_rate_constant(A=5e12, Ea=60e3, T=T)
    >>> fit = fit_arrhenius(T, k)
    >>> round(fit.Ea, 2)
    60000.0
    >>> round(fit.r_squared, 6)
    1.0
    """
    T = np.asarray(T, dtype=np.float64)
    k = np.asarray(k, dtype=np.float64)
    fit = linear_fit(1.0 / T, np.log(k))
    Ea = -fit.slope * R_gas
    A = np.exp(fit.intercept)
    return ArrheniusFit(Ea=float(Ea), A=float(A), r_squared=fit.r_squared, R_gas=R_gas)
