r"""Lindemann-Hinshelwood theory of unimolecular reactions and their pressure falloff.

A "unimolecular" gas reaction A -> P still needs collisions to supply its
activation energy. F. A. Lindemann (*Trans. Faraday Soc.* 17, 598 (1922))
and C. N. Hinshelwood (*Proc. R. Soc. Lond. A* 113, 230 (1927)) proposed

.. math::

    A + M \underset{k_{-1}}{\overset{k_1}{\rightleftharpoons}} A^* + M,
    \qquad A^* \xrightarrow{k_2} P

where M is any collision partner. The steady-state approximation for the
energized molecule :math:`A^*` gives an effective first-order rate
constant that depends on the bath-gas concentration :math:`[M]`:

.. math::

    k_{uni} = \frac{k_1k_2[M]}{k_{-1}[M]+k_2}

which is first order in :math:`[M]` at low pressure (:math:`k_0 = k_1`,
activation is rate-limiting) and saturates at :math:`k_\infty =
k_1k_2/k_{-1}` at high pressure (:math:`A^*` decay is rate-limiting); the
crossover, the *falloff* region, is centered at :math:`[M]_{1/2} =
k_2/k_{-1}` (Atkins & de Paula, *Physical Chemistry*, 11th ed., Ch. 18.5;
Gilbert & Smith, *Theory of Unimolecular and Recombination Reactions*,
Blackwell, 1990, Ch. 1).
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from chemistrykit.kinetics.systems.networks import StoichiometricNetwork

__all__ = ["LindemannHinshelwood"]


@dataclass
class LindemannHinshelwood:
    r"""The Lindemann-Hinshelwood mechanism for a unimolecular reaction A -> P.

    Parameters
    ----------
    k1 : float
        Collisional activation rate constant (second order, A + M).
    k_minus1 : float
        Collisional deactivation rate constant (second order, A* + M).
    k2 : float
        Unimolecular decay rate constant of A* (first order).

    Examples
    --------
    :math:`k_{uni}` rises linearly at low :math:`[M]`, is exactly half of
    :math:`k_\infty` at :math:`[M]_{1/2}`, and saturates at high :math:`[M]`:

    >>> lh = LindemannHinshelwood(k1=1.0, k_minus1=10.0, k2=100.0)
    >>> lh.k_inf, lh.M_half
    (10.0, 10.0)
    >>> float(lh.rate_constant(lh.M_half)) == lh.k_inf / 2
    True
    """

    k1: float
    k_minus1: float
    k2: float

    @property
    def k0(self) -> float:
        r"""float: Low-pressure limiting (second-order) rate constant, :math:`k_1`."""
        return self.k1

    @property
    def k_inf(self) -> float:
        r"""float: High-pressure limiting (first-order) rate constant, :math:`k_1k_2/k_{-1}`."""
        return self.k1 * self.k2 / self.k_minus1

    @property
    def M_half(self) -> float:
        r"""float: Falloff center :math:`[M]_{1/2}=k_2/k_{-1}`, where :math:`k_{uni}=k_\infty/2`."""
        return self.k2 / self.k_minus1

    def rate_constant(self, M):
        r"""Effective first-order rate constant :math:`k_{uni}([M])`.

        Parameters
        ----------
        M : float or array-like of float
            Bath-gas concentration (e.g. :math:`P/RT`).

        Returns
        -------
        float or ndarray
        """
        M = np.asarray(M, dtype=np.float64)
        return self.k1 * self.k2 * M / (self.k_minus1 * M + self.k2)

    def reduced_pressure(self, M):
        r"""Reduced pressure :math:`P_r=k_0[M]/k_\infty`, so that :math:`k_{uni}/k_\infty = P_r/(1+P_r)`.

        Parameters
        ----------
        M : float or array-like of float

        Returns
        -------
        float or ndarray
        """
        return self.k0 * np.asarray(M, dtype=np.float64) / self.k_inf

    def inverse_rate_constant(self, M):
        r"""The Lindemann linearization :math:`1/k_{uni} = 1/k_\infty + 1/(k_1[M])`.

        Plotting :math:`1/k_{uni}` against :math:`1/[M]` gives a straight
        line whose deviations in real data expose the theory's neglect of
        energy-dependent :math:`k_2` (the RRK/RRKM refinement).

        Parameters
        ----------
        M : float or array-like of float

        Returns
        -------
        float or ndarray
        """
        M = np.asarray(M, dtype=np.float64)
        return 1.0 / self.k_inf + 1.0 / (self.k1 * M)

    def network(self, M: float, A0: float = 1.0) -> StoichiometricNetwork:
        """The full three-step mechanism at fixed :math:`[M]`, for numerical integration.

        `M` is folded into pseudo-first-order constants :math:`k_1[M]`
        and :math:`k_{-1}[M]`, so the network needs no steady-state
        approximation; its long-time decay rate of A approaches
        :meth:`rate_constant` when :math:`A^*` stays scarce.

        Parameters
        ----------
        M : float
            Bath-gas concentration.
        A0 : float, default 1.0
            Initial concentration of A.

        Returns
        -------
        StoichiometricNetwork
            Species ``("A", "A*", "P")``.
        """
        stoich = [[-1.0, 1.0, 0.0], [1.0, -1.0, -1.0], [0.0, 0.0, 1.0]]
        orders = [[1.0, 0.0, 0.0], [0.0, 1.0, 1.0], [0.0, 0.0, 0.0]]
        return StoichiometricNetwork(("A", "A*", "P"), stoich, [self.k1 * M, self.k_minus1 * M, self.k2], orders, [A0, 0.0, 0.0])
