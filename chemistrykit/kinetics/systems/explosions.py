r"""Branched-chain explosions: Semenov's net-branching criterion and the H2/O2 explosion limits.

In a branched chain reaction each propagation cycle can make more than one
chain carrier. For the carrier concentration :math:`n`,

.. math::

    \frac{dn}{dt} = w_0 + \varphi\,n, \qquad \varphi = f - g

with :math:`w_0` the initiation rate, :math:`f` the branching and
:math:`g` the termination rate coefficient. If :math:`\varphi < 0`, :math:`n`
settles at :math:`w_0/|\varphi|` and the reaction proceeds steadily. If
:math:`\varphi > 0`, :math:`n` grows exponentially and the mixture explodes
(N. N. Semenov, *Chemical Kinetics and Chain Reactions*, Oxford, 1935;
Nobel Prize 1956, shared with Hinshelwood).

For hydrogen-oxygen the rate-limiting branching step is
:math:`H + O_2 \to OH + O`, which (followed by the fast
:math:`O + H_2 \to OH + H` and :math:`OH + H_2 \to H_2O + H`) turns one H
atom into three, so :math:`f = 2k_b[O_2]`. H atoms are removed at the
vessel wall (first order, :math:`k_w`) and by the termolecular
:math:`H + O_2 + M \to HO_2 + M` (:math:`k_t[O_2][M]`). With
:math:`[O_2] = x_{O_2}c` and :math:`[M] = c = P/RT`,

.. math::

    \varphi(c) = 2k_bx_{O_2}c - k_w - k_tx_{O_2}c^2

is a downward parabola in :math:`c`. It is positive between two roots:
the *first* explosion limit, where branching overtakes wall loss, and the
*second*, where three-body termination catches up (as :math:`k_w\to 0`,
:math:`c_2\to 2k_b/k_t`). Their temperature dependence traces the
explosion peninsula (Atkins & de Paula, *Physical Chemistry*, 11th ed.,
Ch. 19.3; Lewis & von Elbe, *Combustion, Flames and Explosions of Gases*,
3rd ed., Ch. 2). The third, thermal limit at still higher pressure arises
from self-heating, which this isothermal model does not include. The wall
term is taken as kinetically controlled (constant :math:`k_w`) rather
than diffusion-limited.
"""

from __future__ import annotations

from collections.abc import Callable
from dataclasses import dataclass
from typing import Optional

import numpy as np

from chemistrykit.constants import R

__all__ = ["ChainBranchingExplosion"]

RateConstant = float | Callable[[float], float]


def _at(k: RateConstant, T: float) -> float:
    return float(k(T)) if callable(k) else float(k)


@dataclass
class ChainBranchingExplosion:
    r"""Semenov branching-versus-termination model of the H2/O2 first and second explosion limits.

    Each rate constant may be a number or a callable ``k(T)`` (e.g. a
    lambda wrapping :func:`chemistrykit.kinetics.arrhenius_rate_constant`)
    so the limits can be traced across temperature.

    Parameters
    ----------
    k_branch : float or callable
        :math:`k_b` of :math:`H + O_2 \to OH + O`, in m^3 mol^-1 s^-1.
    k_wall : float or callable
        First-order wall-termination rate constant :math:`k_w`, in s^-1.
    k_termination : float or callable
        :math:`k_t` of :math:`H + O_2 + M \to HO_2 + M`, in m^6 mol^-2 s^-1.
    x_O2 : float, default 1/3
        Mole fraction of O2 (stoichiometric 2H2 + O2 by default).

    Examples
    --------
    Explosion happens only between the two limits:

    >>> model = ChainBranchingExplosion(k_branch=1.0e3, k_wall=50.0, k_termination=1.0e3)
    >>> P1, P2 = model.explosion_limits(T=800.0)
    >>> model.explodes(0.5 * P1, 800.0), model.explodes((P1 * P2) ** 0.5, 800.0), model.explodes(2 * P2, 800.0)
    (False, True, False)
    """

    k_branch: RateConstant
    k_wall: RateConstant
    k_termination: RateConstant
    x_O2: float = 1.0 / 3.0

    def branching_factor(self, P, T: float, R_gas: float = R):
        r"""Net branching factor :math:`\varphi = 2k_bx_{O_2}c - k_w - k_tx_{O_2}c^2`, in s^-1.

        Parameters
        ----------
        P : float or array-like of float
            Total pressure, in Pa.
        T : float
            Temperature, in K.
        R_gas : float, default :data:`chemistrykit.constants.R`

        Returns
        -------
        float or ndarray
        """
        c = np.asarray(P, dtype=np.float64) / (R_gas * T)
        kb, kw, kt = _at(self.k_branch, T), _at(self.k_wall, T), _at(self.k_termination, T)
        return 2.0 * kb * self.x_O2 * c - kw - kt * self.x_O2 * c**2

    def explodes(self, P: float, T: float) -> bool:
        """Whether carriers multiply without bound (:math:`\\varphi > 0`) at `P`, `T`.

        Parameters
        ----------
        P : float
            Pressure, in Pa.
        T : float
            Temperature, in K.

        Returns
        -------
        bool
        """
        return bool(self.branching_factor(P, T) > 0.0)

    def explosion_limits(self, T: float, R_gas: float = R) -> Optional[tuple[float, float]]:
        r"""First and second explosion-limit pressures, the roots of :math:`\varphi=0`.

        .. math::

            c_{1,2} = \frac{k_b}{k_t} \mp \sqrt{\frac{k_b^2}{k_t^2} - \frac{k_w}{k_tx_{O_2}}}

        Parameters
        ----------
        T : float
            Temperature, in K.
        R_gas : float, default :data:`chemistrykit.constants.R`

        Returns
        -------
        (P1, P2) or None
            Limit pressures in Pa, or ``None`` if :math:`\varphi<0` at every
            pressure (below the peninsula's tip temperature).
        """
        kb, kw, kt = _at(self.k_branch, T), _at(self.k_wall, T), _at(self.k_termination, T)
        center = kb / kt
        disc = center**2 - kw / (kt * self.x_O2)
        if disc <= 0.0:
            return None
        root = np.sqrt(disc)
        return float((center - root) * R_gas * T), float((center + root) * R_gas * T)

    def carrier_concentration(self, t, P: float, T: float, w0: float, n0: float = 0.0):
        r"""Carrier concentration :math:`n(t)=(n_0+w_0/\varphi)e^{\varphi t}-w_0/\varphi`.

        The exact solution of :math:`dn/dt=w_0+\varphi n`: it levels off
        when :math:`\varphi<0` and grows exponentially when :math:`\varphi>0`.

        Parameters
        ----------
        t : float or array-like of float
            Time(s), in s.
        P, T : float
            Pressure (Pa) and temperature (K).
        w0 : float
            Initiation rate, in mol m^-3 s^-1.
        n0 : float, default 0.0
            Initial carrier concentration.

        Returns
        -------
        float or ndarray
        """
        t = np.asarray(t, dtype=np.float64)
        phi = float(self.branching_factor(P, T))
        if phi == 0.0:
            return n0 + w0 * t
        return (n0 + w0 / phi) * np.exp(phi * t) - w0 / phi
