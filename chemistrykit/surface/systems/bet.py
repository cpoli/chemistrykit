r"""The BET multilayer adsorption isotherm, and its standard linearization.

See Atkins & de Paula, *Physical Chemistry*, 11th ed., Ch. 22.1, and the
original: S. Brunauer, P. H. Emmett & E. Teller, *J. Am. Chem. Soc.* 60,
309 (1938). BET generalizes Langmuir
(:mod:`chemistrykit.surface.systems.langmuir`) by allowing *multilayer*
adsorption: once a site is occupied, a second adsorbate molecule can
condense on top of the first (with the same energetics as ordinary
liquefaction), a third on top of that, and so on, up to the saturation
vapor pressure :math:`P_0` at which bulk condensation occurs and the
adsorbed volume formally diverges.

.. math::

    \frac{V}{V_m} = \frac{Cx}{(1-x)(1-x+Cx)}, \qquad x \equiv \frac{P}{P_0}

where :math:`V_m` is the monolayer capacity and :math:`C` (roughly,
:math:`e^{(E_1-E_L)/RT}`, the ratio of first-layer to liquefaction
adsorption energies) controls how sharply the isotherm's initial
(Langmuir-like) rise is separated from its later upswing toward bulk
condensation.

**BET reduces exactly to Langmuir as** :math:`P_0 \to \infty` **at fixed**
:math:`K \equiv C/P_0`: in that limit :math:`x=P/P_0\to 0` while
:math:`Cx = KP` stays finite, so :math:`(1-x)\to1` and
:math:`(1-x+Cx)\to 1+KP`, giving :math:`V/V_m \to KP/(1+KP)` --
:func:`chemistrykit.surface.systems.langmuir.langmuir_coverage`. Physically
this is the limit in which the vapor is always far below its saturation
pressure, so multilayer condensation never has a chance to set in and
only the monolayer term survives. This reduction is checked numerically
in :mod:`chemistrykit.surface.tests.test_bet` (and was verified
numerically -- to a relative error below :math:`10^{-4}` at
:math:`P_0/K \sim 10^7` -- before this docstring was written).
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from chemistrykit.surface.core.base_system import AdsorptionIsotherm
from chemistrykit.surface.utils.regression import LinearFit, linear_fit

__all__ = ["bet_loading", "BETIsotherm", "BETFit", "fit_bet"]


def bet_loading(Vm: float, C: float, P, P0: float):
    r"""The BET isotherm :math:`V(P) = V_m Cx/[(1-x)(1-x+Cx)]`, :math:`x=P/P_0`.

    Parameters
    ----------
    Vm : float
        Monolayer capacity (e.g. cm^3(STP)/g or mol/g).
    C : float
        BET energy constant (:math:`C \gg 1` gives a sharp "knee" at
        monolayer completion; :math:`C \to 1` washes it out).
    P : float or array-like of float
        Equilibrium pressure, in the same units as `P0`.
    P0 : float
        Saturation vapor pressure of the adsorbate at the experiment's
        temperature.

    Returns
    -------
    float or ndarray
        Volume (or amount) adsorbed. Diverges as :math:`P \to P_0`
        (formal bulk condensation); only valid for :math:`P < P_0`.

    Examples
    --------
    At zero pressure nothing is adsorbed:

    >>> round(float(bet_loading(Vm=1.0, C=100.0, P=0.0, P0=10.0)), 10)
    0.0

    BET reduces to the Langmuir isotherm when the saturation pressure is
    far above the working pressure range (multilayer condensation
    suppressed) -- see the module docstring for the derivation:

    >>> from chemistrykit.surface.systems.langmuir import langmuir_coverage
    >>> K = 2.0
    >>> P0 = 1.0e7
    >>> Vm = 1.0
    >>> C = K * P0
    >>> P = np.array([0.01, 0.1, 0.5, 1.0, 2.0])
    >>> V_bet = bet_loading(Vm, C, P, P0)
    >>> theta_langmuir = langmuir_coverage(K, P)
    >>> bool(np.allclose(V_bet, Vm * theta_langmuir, rtol=1e-4))
    True
    """
    P = np.asarray(P, dtype=np.float64)
    x = P / P0
    result = Vm * C * x / ((1.0 - x) * (1.0 - x + C * x))
    return float(result) if result.ndim == 0 else result


class BETIsotherm(AdsorptionIsotherm):
    r"""The BET multilayer adsorption isotherm.

    Parameters
    ----------
    Vm : float
        Monolayer capacity.
    C : float
        BET energy constant.
    P0 : float
        Saturation vapor pressure of the adsorbate.

    Examples
    --------
    >>> iso = BETIsotherm(Vm=2.0, C=100.0, P0=10.0)
    >>> round(float(iso.loading(P=1.0)), 4)
    2.0387
    """

    def __init__(self, Vm: float, C: float, P0: float):
        if Vm <= 0:
            raise ValueError("Vm must be positive")
        if C <= 0:
            raise ValueError("C must be positive")
        if P0 <= 0:
            raise ValueError("P0 must be positive")
        self.Vm = float(Vm)
        self.C = float(C)
        self.P0 = float(P0)
        self.saturation_loading = self.Vm

    def loading(self, P):
        return bet_loading(self.Vm, self.C, P, self.P0)


@dataclass
class BETFit:
    """Result of fitting loading-vs-pressure data to the BET isotherm."""

    Vm: float
    """float: Fitted monolayer capacity."""

    C: float
    """float: Fitted BET energy constant."""

    P0: float
    """float: Saturation vapor pressure used (an input to the fit, not itself fitted)."""

    r_squared: float
    """float: Coefficient of determination of the linearized fit."""

    def to_isotherm(self) -> BETIsotherm:
        """Return a :class:`BETIsotherm` built from the fitted parameters."""
        return BETIsotherm(Vm=self.Vm, C=self.C, P0=self.P0)

    def predict(self, P):
        """Evaluate the fitted BET isotherm at pressure(s) `P`.

        Parameters
        ----------
        P : float or array-like of float

        Returns
        -------
        float or ndarray
        """
        return self.to_isotherm().loading(P)


def fit_bet(P, V, P0: float) -> BETFit:
    r"""Fit loading-vs-pressure data to the BET isotherm via its linearization.

    The BET equation rearranges to a straight line in :math:`x=P/P_0`
    (Atkins & de Paula, *Physical Chemistry*, 11th ed., Ch. 22.1;
    Brunauer, Emmett & Teller 1938):

    .. math::

        \frac{x}{V(1-x)} = \frac{1}{V_m C} + \frac{C-1}{V_m C}\,x

    so plotting :math:`x/[V(1-x)]` against `x` gives a straight line of
    intercept :math:`1/(V_mC)` and slope :math:`(C-1)/(V_mC)`; solving
    the two simultaneously gives :math:`V_m = 1/(\text{slope}+\text{intercept})`
    and :math:`C = 1 + \text{slope}/\text{intercept}`. Fitting is normally
    restricted to the range :math:`0.05 \lesssim x \lesssim 0.35`, where
    the BET assumptions are most reliable (Atkins & de Paula, *loc. cit.*)
    -- not enforced here, since this function only performs the
    regression on whatever data it is given.

    Parameters
    ----------
    P : array-like of float
        Equilibrium pressures (at least 2 distinct values, all ``< P0``).
    V : array-like of float
        Corresponding measured adsorbed volumes/amounts.
    P0 : float
        Saturation vapor pressure of the adsorbate.

    Returns
    -------
    BETFit

    Examples
    --------
    Generate exact data from a known :math:`(V_m, C)` and recover them:

    >>> P0 = 10.0
    >>> P = np.array([0.5, 1.0, 1.5, 2.0, 2.5, 3.0])
    >>> V = BETIsotherm(Vm=5.0, C=80.0, P0=P0).loading(P)
    >>> fit = fit_bet(P, V, P0=P0)
    >>> round(fit.Vm, 4), round(fit.C, 2)
    (5.0, 80.0)
    >>> round(fit.r_squared, 6)
    1.0
    """
    P = np.asarray(P, dtype=np.float64)
    V = np.asarray(V, dtype=np.float64)
    x = P / P0
    y = x / (V * (1.0 - x))
    fit: LinearFit = linear_fit(x, y)
    Vm = 1.0 / (fit.slope + fit.intercept)
    C = 1.0 + fit.slope / fit.intercept
    return BETFit(Vm=float(Vm), C=float(C), P0=float(P0), r_squared=fit.r_squared)
