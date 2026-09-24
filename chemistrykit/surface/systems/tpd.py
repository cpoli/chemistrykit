r"""Temperature-programmed desorption (TPD) and Redhead's peak analysis.

A monolayer heated at a constant rate :math:`\beta = dT/dt` desorbs
according to the Polanyi-Wigner rate equation

.. math::

    -\frac{d\theta}{dT} = \frac{\nu}{\beta}\,\theta^n\,e^{-E_d/RT},

with :math:`E_d` the desorption activation energy, :math:`\nu` the
pre-exponential factor and :math:`n` the desorption order (Masel,
*Principles of Adsorption and Reaction on Solid Surfaces*, Wiley, 1996,
Ch. 7). For :math:`n=1,2` it integrates in closed form given
:math:`I(T)=\int_{T_0}^{T}e^{-E_d/RT'}\,dT'`:
:math:`\theta=\theta_0e^{-\nu I/\beta}` (first order) and
:math:`\theta=\theta_0/(1+\theta_0\nu I/\beta)` (second order), which
:func:`simulate_tpd` evaluates by cumulative quadrature.

Setting the derivative of the first-order desorption rate to zero gives
the exact peak condition :math:`E_d/(RT_p^2) = (\nu/\beta)e^{-E_d/RT_p}`
(:func:`first_order_peak_temperature`). P. A. Redhead, "Thermal
Desorption of Gases," *Vacuum* 12 (1962), 203, approximated its solution
as

.. math::

    E_d = RT_p\left[\ln\frac{\nu T_p}{\beta} - 3.64\right],

accurate to within a few percent for :math:`10^8 < \nu/\beta < 10^{13}\ \mathrm{K^{-1}}`
and typical binding energies (:func:`redhead_desorption_energy`) -- turning a single measured peak
temperature into a binding energy.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from scipy.integrate import cumulative_trapezoid
from scipy.optimize import brentq

from chemistrykit.constants import R

__all__ = ["TPDResult", "simulate_tpd", "first_order_peak_temperature", "redhead_desorption_energy"]


@dataclass
class TPDResult:
    """Result of a simulated temperature-programmed desorption run."""

    T: np.ndarray
    """ndarray: Temperature grid, in K."""

    coverage: np.ndarray
    """ndarray: Fractional coverage :math:`\\theta(T)`."""

    desorption_rate: np.ndarray
    """ndarray: Desorption rate per kelvin, :math:`-d\\theta/dT`."""

    peak_temperature: float
    """float: Temperature of maximum desorption rate, in K."""


def simulate_tpd(
    Ed: float,
    nu: float,
    beta: float,
    theta0: float = 1.0,
    order: int = 1,
    T_start: float = 100.0,
    T_end: float = 800.0,
    n_points: int = 20001,
    R_gas: float = R,
) -> TPDResult:
    r"""Simulate a TPD spectrum from the Polanyi-Wigner equation (first or second order).

    Parameters
    ----------
    Ed : float
        Desorption activation energy, in J/mol.
    nu : float
        Pre-exponential factor, in 1/s (first order) or 1/(s * coverage) (second order).
    beta : float
        Linear heating rate, in K/s.
    theta0 : float, default 1.0
        Initial coverage.
    order : {1, 2}, default 1
        Desorption order.
    T_start, T_end : float
        Temperature range of the ramp, in K.
    n_points : int, default 20001
        Number of grid points.
    R_gas : float, default :data:`chemistrykit.constants.R`

    Returns
    -------
    TPDResult

    Examples
    --------
    The simulated first-order peak sits where the exact peak condition
    puts it:

    >>> res = simulate_tpd(Ed=100e3, nu=1e13, beta=10.0)
    >>> bool(abs(res.peak_temperature - first_order_peak_temperature(100e3, 1e13, 10.0)) < 0.05)
    True
    """
    if order not in (1, 2):
        raise ValueError("order must be 1 or 2")
    T = np.linspace(T_start, T_end, n_points)
    boltzmann = np.exp(-Ed / (R_gas * T))
    integral = cumulative_trapezoid(boltzmann, T, initial=0.0)
    x = nu * integral / beta
    theta = theta0 * np.exp(-x) if order == 1 else theta0 / (1.0 + theta0 * x)
    rate = nu / beta * theta**order * boltzmann
    i = int(np.argmax(rate))
    T_peak = float(T[i])
    if 0 < i < n_points - 1:
        y0, y1, y2 = rate[i - 1], rate[i], rate[i + 1]
        denom = y0 - 2.0 * y1 + y2
        if denom != 0.0:
            T_peak += 0.5 * (y0 - y2) / denom * (T[1] - T[0])
    return TPDResult(T=T, coverage=theta, desorption_rate=rate, peak_temperature=T_peak)


def first_order_peak_temperature(Ed: float, nu: float, beta: float, R_gas: float = R) -> float:
    r"""Exact first-order TPD peak temperature from :math:`E_d/(RT_p^2)=(\nu/\beta)e^{-E_d/RT_p}`.

    Parameters
    ----------
    Ed : float
        Desorption activation energy, in J/mol.
    nu : float
        Pre-exponential factor, in 1/s.
    beta : float
        Heating rate, in K/s.
    R_gas : float, default :data:`chemistrykit.constants.R`

    Returns
    -------
    float
        Peak temperature, in K.

    Examples
    --------
    >>> Tp = first_order_peak_temperature(Ed=100e3, nu=1e13, beta=10.0)
    >>> lhs = 100e3 / (8.31446261815324 * Tp**2)
    >>> rhs = 1e13 / 10.0 * np.exp(-100e3 / (8.31446261815324 * Tp))
    >>> round(float(lhs / rhs), 8)
    1.0
    """

    def g(T):
        return np.log(Ed / (R_gas * T**2)) - np.log(nu / beta) + Ed / (R_gas * T)

    return float(brentq(g, 1.0, 1.0e6, xtol=1e-12))


def redhead_desorption_energy(T_peak, nu: float, beta: float, R_gas: float = R):
    r"""Redhead's estimate :math:`E_d = RT_p[\ln(\nu T_p/\beta) - 3.64]`, in J/mol.

    Parameters
    ----------
    T_peak : float or array-like of float
        Measured first-order TPD peak temperature, in K.
    nu : float
        Assumed pre-exponential factor, in 1/s (commonly :math:`10^{13}`).
    beta : float
        Heating rate, in K/s.
    R_gas : float, default :data:`chemistrykit.constants.R`

    Returns
    -------
    float or ndarray

    Examples
    --------
    Redhead's formula recovers the true desorption energy to within about
    one percent at the common choice :math:`\nu=10^{13}` s^-1:

    >>> Tp = first_order_peak_temperature(Ed=100e3, nu=1e13, beta=10.0)
    >>> bool(abs(redhead_desorption_energy(Tp, nu=1e13, beta=10.0) / 100e3 - 1.0) < 0.015)
    True
    """
    T_peak = np.asarray(T_peak, dtype=np.float64)
    result = R_gas * T_peak * (np.log(nu * T_peak / beta) - 3.64)
    return float(result) if result.ndim == 0 else result
