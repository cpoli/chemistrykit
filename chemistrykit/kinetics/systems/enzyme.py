"""Michaelis-Menten enzyme kinetics: Lineweaver-Burk linearization and inhibition.

See Atkins & de Paula, *Physical Chemistry*, 11th ed., Ch. 20.6 ("Enzyme-
catalysed reactions"); the original: L. Michaelis & M. L. Menten, *Biochem.
Z.* 49, 333 (1913); the linearization: H. Lineweaver & D. Burk, *J. Am.
Chem. Soc.* 56, 658 (1934).
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from numba import njit
from numpy.typing import NDArray

from chemistrykit.kinetics.core.base_system import ReactionNetwork
from chemistrykit.kinetics.utils.linear_regression import linear_fit

__all__ = [
    "michaelis_menten_rate",
    "competitive_inhibition_rate",
    "noncompetitive_inhibition_rate",
    "MichaelisMentenFit",
    "fit_lineweaver_burk",
    "MichaelisMentenProgress",
]


def michaelis_menten_rate(S, Vmax: float, Km: float):
    r"""Michaelis-Menten rate law :math:`v = V_{max}[S] / (K_m + [S])`.

    Parameters
    ----------
    S : float or array-like of float
        Substrate concentration(s).
    Vmax : float
        Maximum reaction rate (attained as :math:`[S] \to \infty`).
    Km : float
        Michaelis constant: the substrate concentration at which
        :math:`v = V_{max}/2`.

    Returns
    -------
    float or ndarray

    Examples
    --------
    At ``S = Km`` the rate is exactly half of ``Vmax``, by construction:

    >>> round(float(michaelis_menten_rate(S=2.0, Vmax=10.0, Km=2.0)), 6)
    5.0
    """
    S = np.asarray(S, dtype=np.float64)
    return Vmax * S / (Km + S)


def competitive_inhibition_rate(S, I, Vmax: float, Km: float, Ki: float):
    r"""Competitive-inhibition Michaelis-Menten rate law.

    The inhibitor competes with substrate for the active site, which is
    equivalent to inflating the apparent :math:`K_m` by a factor
    :math:`\alpha = 1 + [I]/K_i` while leaving :math:`V_{max}` unchanged
    (Atkins & de Paula, *Physical Chemistry*, 11th ed., Ch. 20.6):

    .. math::

        v = \frac{V_{max}[S]}{K_m (1 + [I]/K_i) + [S]}

    Parameters
    ----------
    S : float or array-like of float
        Substrate concentration(s).
    I : float
        Inhibitor concentration.
    Vmax, Km : float
        Uninhibited Michaelis-Menten parameters.
    Ki : float
        Inhibitor dissociation constant.

    Returns
    -------
    float or ndarray

    Examples
    --------
    With no inhibitor (``I=0``) this reduces exactly to the uninhibited
    rate law:

    >>> float(competitive_inhibition_rate(S=2.0, I=0.0, Vmax=10.0, Km=2.0, Ki=1.0))
    5.0
    """
    S = np.asarray(S, dtype=np.float64)
    alpha = 1.0 + I / Ki
    return Vmax * S / (Km * alpha + S)


def noncompetitive_inhibition_rate(S, I, Vmax: float, Km: float, Ki: float):
    r"""Noncompetitive-inhibition Michaelis-Menten rate law.

    The inhibitor binds a separate site with equal affinity for the free
    enzyme and the enzyme-substrate complex, which is equivalent to
    deflating the apparent :math:`V_{max}` by :math:`\alpha = 1 + [I]/K_i`
    while leaving :math:`K_m` unchanged (Atkins & de Paula, *Physical
    Chemistry*, 11th ed., Ch. 20.6):

    .. math::

        v = \frac{V_{max}[S]}{\alpha (K_m + [S])}

    Parameters
    ----------
    S : float or array-like of float
        Substrate concentration(s).
    I : float
        Inhibitor concentration.
    Vmax, Km : float
        Uninhibited Michaelis-Menten parameters.
    Ki : float
        Inhibitor dissociation constant.

    Returns
    -------
    float or ndarray

    Examples
    --------
    With no inhibitor (``I=0``) this reduces exactly to the uninhibited
    rate law:

    >>> float(noncompetitive_inhibition_rate(S=2.0, I=0.0, Vmax=10.0, Km=2.0, Ki=1.0))
    5.0
    """
    S = np.asarray(S, dtype=np.float64)
    alpha = 1.0 + I / Ki
    return Vmax * S / (alpha * (Km + S))


@dataclass
class MichaelisMentenFit:
    """Result of a Lineweaver-Burk (double-reciprocal) linear fit."""

    Vmax: float
    """float: Fitted maximum rate."""

    Km: float
    """float: Fitted Michaelis constant."""

    r_squared: float
    """float: Coefficient of determination of the linear (1/v vs 1/S) fit."""

    def predict(self, S):
        """Evaluate the fitted Michaelis-Menten rate law at substrate concentration(s) `S`.

        Parameters
        ----------
        S : float or array-like of float

        Returns
        -------
        float or ndarray
        """
        return michaelis_menten_rate(S, self.Vmax, self.Km)


def fit_lineweaver_burk(S, v) -> MichaelisMentenFit:
    r"""Fit substrate/rate data via the Lineweaver-Burk linearization.

    Inverting the Michaelis-Menten equation gives the "double-reciprocal"
    linear form (Lineweaver & Burk, 1934; Atkins & de Paula, *Physical
    Chemistry*, 11th ed., Ch. 20.6):

    .. math::

        \frac{1}{v} = \frac{K_m}{V_{max}} \cdot \frac{1}{[S]} + \frac{1}{V_{max}}

    so plotting :math:`1/v` against :math:`1/[S]` gives a straight line
    of slope :math:`K_m/V_{max}` and intercept :math:`1/V_{max}`, fit here
    by ordinary least squares.

    Parameters
    ----------
    S : array-like of float
        Substrate concentrations (at least 2 distinct, positive values).
    v : array-like of float
        Initial rates measured at each concentration in `S`.

    Returns
    -------
    MichaelisMentenFit

    Examples
    --------
    Generate exact data from known (Vmax, Km) and recover them:

    >>> import numpy as np
    >>> S = np.array([0.5, 1.0, 2.0, 4.0, 8.0])
    >>> v = michaelis_menten_rate(S, Vmax=8.0, Km=1.5)
    >>> fit = fit_lineweaver_burk(S, v)
    >>> round(fit.Vmax, 6)
    8.0
    >>> round(fit.Km, 6)
    1.5
    """
    S = np.asarray(S, dtype=np.float64)
    v = np.asarray(v, dtype=np.float64)
    fit = linear_fit(1.0 / S, 1.0 / v)
    Vmax = 1.0 / fit.intercept
    Km = fit.slope * Vmax
    return MichaelisMentenFit(Vmax=float(Vmax), Km=float(Km), r_squared=fit.r_squared)


@njit(cache=True)
def _mm_progress_rhs(state: NDArray[np.float64], t: float, params: NDArray[np.float64]) -> NDArray[np.float64]:
    """Substrate-depletion vector field ``d[S]/dt = -Vmax*[S]/(Km+[S])``.

    Parameters
    ----------
    state : ndarray of float, shape (1,)
        State vector ``([S],)``.
    t : float
        Current time (unused; the system is autonomous).
    params : ndarray of float, shape (2,)
        Parameters ``(Vmax, Km)``.

    Returns
    -------
    ndarray of float, shape (1,)
    """
    Vmax, Km = params[0], params[1]
    S = state[0]
    out = np.empty(1)
    out[0] = -Vmax * S / (Km + S)
    return out


class MichaelisMentenProgress(ReactionNetwork):
    """Substrate-depletion progress curve under Michaelis-Menten kinetics.

    Integrates ``d[S]/dt = -Vmax*[S]/(Km+[S])`` forward in time via
    :mod:`chemistrykit.integrators` -- the "progress curve" an enzyme assay
    actually measures, as distinct from the *initial*-rate
    :func:`michaelis_menten_rate` used for a Lineweaver-Burk fit (Atkins &
    de Paula, *Physical Chemistry*, 11th ed., Ch. 20.6). In the two
    asymptotic limits this reduces to the closed-form rate laws in
    :mod:`chemistrykit.kinetics.systems.rate_laws`: zeroth order in `S`
    when ``S0 >> Km`` (:class:`~chemistrykit.kinetics.systems.rate_laws.ZeroOrder`)
    and first order when ``S0 << Km``
    (:class:`~chemistrykit.kinetics.systems.rate_laws.FirstOrder`).

    Parameters
    ----------
    S0 : float
        Initial substrate concentration.
    Vmax, Km : float
        Michaelis-Menten parameters.
    """

    species = ("S",)

    def __init__(self, S0: float, Vmax: float, Km: float):
        self.Vmax = float(Vmax)
        self.Km = float(Km)
        self.params = np.array([self.Vmax, self.Km])
        self._rhs_njit = _mm_progress_rhs
        super().__init__([S0])

    def rhs(self, state, t):
        return np.asarray(_mm_progress_rhs(np.asarray(state, dtype=np.float64), t, self.params))
