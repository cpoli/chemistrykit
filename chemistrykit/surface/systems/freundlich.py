r"""The Freundlich adsorption isotherm, and its standard linearization.

See Atkins & de Paula, *Physical Chemistry*, 11th ed., Ch. 22.1, and the
original: H. Freundlich, *Z. Phys. Chem.* 57, 385 (1907). Unlike Langmuir
(:mod:`chemistrykit.surface.systems.langmuir`), the Freundlich isotherm is
a purely empirical power law with no saturation limit -- it is a good
approximate description of adsorption on a *heterogeneous* surface (a
distribution of site binding energies) over the pressure range it was
fitted to, but should not be extrapolated to very high pressure, where
real adsorption always saturates.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from chemistrykit.surface.core.base_system import AdsorptionIsotherm
from chemistrykit.surface.utils.regression import LinearFit, linear_fit

__all__ = ["freundlich_loading", "FreundlichIsotherm", "FreundlichFit", "fit_freundlich"]


def freundlich_loading(Kf: float, n: float, P):
    r"""The Freundlich isotherm :math:`q = K_f P^{1/n}`.

    Parameters
    ----------
    Kf : float
        Freundlich capacity constant (the loading at :math:`P=1`, in the
        pressure units used).
    n : float
        Freundlich heterogeneity exponent (:math:`n>1` is the usual,
        "favorable" case; :math:`n=1` recovers a linear (Henry's-law-like)
        isotherm).
    P : float or array-like of float
        Equilibrium (partial) pressure of the adsorbate.

    Returns
    -------
    float or ndarray

    Examples
    --------
    At unit pressure, the loading equals :math:`K_f` exactly, by
    construction:

    >>> round(float(freundlich_loading(Kf=2.5, n=3.0, P=1.0)), 10)
    2.5

    With ``n=1`` the isotherm is exactly linear in `P`:

    >>> P = np.array([0.5, 1.0, 2.0, 4.0])
    >>> q = freundlich_loading(Kf=3.0, n=1.0, P=P)
    >>> bool(np.allclose(q, 3.0 * P))
    True
    """
    P = np.asarray(P, dtype=np.float64)
    result = Kf * P ** (1.0 / n)
    return float(result) if result.ndim == 0 else result


class FreundlichIsotherm(AdsorptionIsotherm):
    r"""The Freundlich adsorption isotherm :math:`q(P) = K_f P^{1/n}`.

    Parameters
    ----------
    Kf : float
        Freundlich capacity constant.
    n : float
        Freundlich heterogeneity exponent.

    Notes
    -----
    Unlike :class:`~chemistrykit.surface.systems.langmuir.LangmuirIsotherm`,
    the Freundlich isotherm has no saturation (monolayer) loading, so
    :meth:`~chemistrykit.surface.core.base_system.AdsorptionIsotherm.fractional_coverage`
    is not physically meaningful here and raises :class:`NotImplementedError`.

    Examples
    --------
    >>> iso = FreundlichIsotherm(Kf=2.5, n=3.0)
    >>> round(float(iso.loading(P=1.0)), 6)
    2.5
    """

    def __init__(self, Kf: float, n: float):
        if Kf <= 0:
            raise ValueError("Kf must be positive")
        if n <= 0:
            raise ValueError("n must be positive")
        self.Kf = float(Kf)
        self.n = float(n)

    def loading(self, P):
        return freundlich_loading(self.Kf, self.n, P)

    def fractional_coverage(self, P):
        raise NotImplementedError("the Freundlich isotherm has no saturation (monolayer) loading to normalize by")


@dataclass
class FreundlichFit:
    """Result of fitting loading-vs-pressure data to the Freundlich isotherm."""

    Kf: float
    """float: Fitted Freundlich capacity constant."""

    n: float
    """float: Fitted Freundlich heterogeneity exponent."""

    r_squared: float
    """float: Coefficient of determination of the linearized (``log q`` vs ``log P``) fit."""

    def to_isotherm(self) -> FreundlichIsotherm:
        """Return a :class:`FreundlichIsotherm` built from the fitted parameters."""
        return FreundlichIsotherm(Kf=self.Kf, n=self.n)

    def predict(self, P):
        """Evaluate the fitted Freundlich isotherm at pressure(s) `P`.

        Parameters
        ----------
        P : float or array-like of float

        Returns
        -------
        float or ndarray
        """
        return self.to_isotherm().loading(P)


def fit_freundlich(P, q) -> FreundlichFit:
    r"""Fit loading-vs-pressure data to the Freundlich isotherm via its linearization.

    Taking logarithms of :math:`q = K_f P^{1/n}` gives a straight line
    (Atkins & de Paula, *Physical Chemistry*, 11th ed., Ch. 22.1):

    .. math::

        \ln q = \ln K_f + \frac{1}{n}\ln P

    so plotting :math:`\ln q` against :math:`\ln P` gives a straight line
    of intercept :math:`\ln K_f` and slope :math:`1/n`.

    Parameters
    ----------
    P : array-like of float
        Equilibrium pressures (at least 2 distinct positive values).
    q : array-like of float
        Corresponding measured loadings (all positive).

    Returns
    -------
    FreundlichFit

    Examples
    --------
    Generate exact data from a known :math:`(K_f, n)` and recover them:

    >>> P = np.array([0.1, 0.5, 1.0, 2.0, 5.0, 10.0])
    >>> q = FreundlichIsotherm(Kf=4.0, n=2.5).loading(P)
    >>> fit = fit_freundlich(P, q)
    >>> round(fit.Kf, 6), round(fit.n, 6)
    (4.0, 2.5)
    >>> round(fit.r_squared, 6)
    1.0
    """
    P = np.asarray(P, dtype=np.float64)
    q = np.asarray(q, dtype=np.float64)
    fit: LinearFit = linear_fit(np.log(P), np.log(q))
    n = 1.0 / fit.slope
    Kf = np.exp(fit.intercept)
    return FreundlichFit(Kf=float(Kf), n=float(n), r_squared=fit.r_squared)
