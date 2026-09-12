r"""The Langmuir adsorption isotherm, and its standard linearization.

See Atkins & de Paula, *Physical Chemistry*, 11th ed., Ch. 22 ("Processes
at solid surfaces"), and the original: I. Langmuir, *J. Am. Chem. Soc.*
40, 1361 (1918).

The Langmuir model assumes a fixed number of identical, independent
surface sites, each holding at most one adsorbate molecule, in dynamic
equilibrium between adsorption (rate :math:`k_a P (1-\theta)`,
proportional to the pressure and the fraction of *empty* sites) and
desorption (rate :math:`k_d \theta`). Setting these equal and defining
:math:`K = k_a/k_d` gives the fractional surface coverage

.. math::

    \theta(P) = \frac{KP}{1+KP}

which saturates at :math:`\theta \to 1` as :math:`P \to \infty` (a true
monolayer limit) -- the key qualitative difference from the Freundlich
isotherm (:mod:`chemistrykit.surface.systems.freundlich`), which has no
saturation limit.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from chemistrykit.surface.core.base_system import AdsorptionIsotherm
from chemistrykit.surface.utils.regression import LinearFit, linear_fit

__all__ = ["langmuir_coverage", "LangmuirIsotherm", "LangmuirFit", "fit_langmuir"]


def langmuir_coverage(K: float, P):
    r"""The Langmuir isotherm's fractional surface coverage :math:`\theta = KP/(1+KP)`.

    Parameters
    ----------
    K : float
        Langmuir equilibrium (adsorption) constant, inverse to `P`'s
        units (e.g. 1/Pa for `P` in Pa).
    P : float or array-like of float
        Equilibrium (partial) pressure of the adsorbate.

    Returns
    -------
    float or ndarray
        Fractional coverage :math:`\theta \in [0, 1)`.

    Examples
    --------
    Coverage is exactly one-half at the half-saturation pressure
    :math:`P = 1/K` -- the defining, exactly-solvable feature of the
    Langmuir isotherm (Atkins & de Paula, *Physical Chemistry*, 11th ed.,
    Ch. 22.1):

    >>> round(float(langmuir_coverage(K=2.0, P=0.5)), 10)
    0.5

    Coverage saturates toward 1 at high pressure and vanishes at zero
    pressure:

    >>> round(float(langmuir_coverage(K=2.0, P=0.0)), 10)
    0.0
    >>> round(float(langmuir_coverage(K=2.0, P=1.0e6)), 6)
    1.0
    """
    P = np.asarray(P, dtype=np.float64)
    result = K * P / (1.0 + K * P)
    return float(result) if result.ndim == 0 else result


class LangmuirIsotherm(AdsorptionIsotherm):
    r"""The Langmuir monolayer adsorption isotherm :math:`q(P) = q_{max}\,KP/(1+KP)`.

    Parameters
    ----------
    K : float
        Langmuir equilibrium constant, inverse to `P`'s units.
    qmax : float, default 1.0
        Monolayer (saturation) capacity, in whatever amount unit `q` is
        measured in (e.g. mol/g).

    Examples
    --------
    >>> iso = LangmuirIsotherm(K=2.0, qmax=5.0)
    >>> round(float(iso.loading(P=0.5)), 6)
    2.5
    >>> round(float(iso.fractional_coverage(P=0.5)), 6)
    0.5
    """

    def __init__(self, K: float, qmax: float = 1.0):
        if K <= 0:
            raise ValueError("K must be positive")
        if qmax <= 0:
            raise ValueError("qmax must be positive")
        self.K = float(K)
        self.qmax = float(qmax)
        self.saturation_loading = self.qmax

    def loading(self, P):
        return self.qmax * langmuir_coverage(self.K, P)

    def half_saturation_pressure(self) -> float:
        """Return the pressure :math:`P=1/K` at which coverage is exactly one-half.

        Returns
        -------
        float

        Examples
        --------
        >>> LangmuirIsotherm(K=4.0).half_saturation_pressure()
        0.25
        """
        return 1.0 / self.K


@dataclass
class LangmuirFit:
    """Result of fitting loading-vs-pressure data to the Langmuir isotherm."""

    K: float
    """float: Fitted Langmuir equilibrium constant."""

    qmax: float
    """float: Fitted monolayer capacity."""

    r_squared: float
    """float: Coefficient of determination of the linearized (``1/q`` vs ``1/P``) fit."""

    def to_isotherm(self) -> LangmuirIsotherm:
        """Return a :class:`LangmuirIsotherm` built from the fitted parameters."""
        return LangmuirIsotherm(K=self.K, qmax=self.qmax)

    def predict(self, P):
        """Evaluate the fitted Langmuir isotherm at pressure(s) `P`.

        Parameters
        ----------
        P : float or array-like of float

        Returns
        -------
        float or ndarray
        """
        return self.to_isotherm().loading(P)


def fit_langmuir(P, q) -> LangmuirFit:
    r"""Fit loading-vs-pressure data to the Langmuir isotherm via its linearization.

    The Langmuir equation :math:`q = q_{max}KP/(1+KP)` is inverted to a
    straight line (Atkins & de Paula, *Physical Chemistry*, 11th ed., Ch.
    22.1):

    .. math::

        \frac{1}{q} = \frac{1}{q_{max}} + \frac{1}{q_{max}K}\cdot\frac{1}{P}

    so plotting :math:`1/q` against :math:`1/P` gives a straight line of
    intercept :math:`1/q_{max}` and slope :math:`1/(q_{max}K)`.

    Parameters
    ----------
    P : array-like of float
        Equilibrium pressures (at least 2 distinct positive values).
    q : array-like of float
        Corresponding measured loadings (all positive).

    Returns
    -------
    LangmuirFit

    Examples
    --------
    Generate exact data from a known :math:`(K, q_{max})` and recover them:

    >>> P = np.array([0.1, 0.5, 1.0, 2.0, 5.0, 10.0])
    >>> q = LangmuirIsotherm(K=3.0, qmax=8.0).loading(P)
    >>> fit = fit_langmuir(P, q)
    >>> round(fit.K, 6), round(fit.qmax, 6)
    (3.0, 8.0)
    >>> round(fit.r_squared, 6)
    1.0
    """
    P = np.asarray(P, dtype=np.float64)
    q = np.asarray(q, dtype=np.float64)
    fit: LinearFit = linear_fit(1.0 / P, 1.0 / q)
    qmax = 1.0 / fit.intercept
    K = fit.intercept / fit.slope
    return LangmuirFit(K=float(K), qmax=float(qmax), r_squared=fit.r_squared)
