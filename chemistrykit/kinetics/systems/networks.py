"""A general stoichiometry-matrix reaction-network engine.

:class:`StoichiometricNetwork` models an arbitrary set of elementary
mass-action reactions from a stoichiometry matrix and per-reaction rate
constants, integrated forward in time via :mod:`chemistrykit.integrators`.
For a network of ``R`` reactions among ``N`` species, reaction ``j``'s
rate is the mass-action law

.. math::

    \\text{rate}_j = k_j \\prod_i [X_i]^{\\nu_{ij}^{react}}

(product over each species' *reactant*-side stoichiometric coefficient
only) and the net rate of change of species `i` is

.. math::

    \\frac{d[X_i]}{dt} = \\sum_j \\nu_{ij} \\, \\text{rate}_j

where :math:`\\nu_{ij}` is the (signed, net) stoichiometric coefficient of
species `i` in reaction `j` -- see Atkins & de Paula, *Physical
Chemistry*, 11th ed., Ch. 20 for the elementary-reaction rate laws this
composes, or Espenson, *Chemical Kinetics and Reaction Mechanisms*, 2nd
ed., Ch. 4-5 for parallel/consecutive/reversible mechanisms specifically.

The classmethods below (:meth:`~StoichiometricNetwork.parallel`,
:meth:`~StoichiometricNetwork.consecutive`,
:meth:`~StoichiometricNetwork.reversible`) build the specific named
topologies the mechanism has a closed-form solution for (see the
``*_analytic`` functions), so numerical results can be checked against
them; the general constructor accepts any stoichiometry matrix.
"""

from __future__ import annotations

from collections.abc import Sequence

import numpy as np
from numba import njit
from numpy.typing import NDArray

from chemistrykit.kinetics.core.base_system import ReactionNetwork

__all__ = [
    "StoichiometricNetwork",
    "consecutive_analytic",
    "reversible_analytic",
    "ssa_intermediate_concentration",
]


def _make_mass_action_rhs(stoich_matrix: np.ndarray, rate_constants: np.ndarray, reactant_orders: np.ndarray):
    """Build a standalone ``@njit`` mass-action right-hand side.

    Factory that closes over the network's structural arrays (its
    stoichiometry matrix, rate constants, and reactant orders) and returns
    a genuine ``@njit`` dispatcher usable as
    :attr:`chemistrykit.kinetics.core.base_system.ReactionNetwork._rhs_njit`
    -- the same factory-closure pattern physicskit uses for e.g.
    ``physicskit.classical.systems.newtonian._make_force_njit``, needed
    because Numba's nopython mode can only call another genuine njit
    dispatcher from inside a compiled function, not a bound Python method.

    Parameters
    ----------
    stoich_matrix : ndarray, shape (n_species, n_reactions)
        Net stoichiometric coefficient of each species in each reaction.
    rate_constants : ndarray, shape (n_reactions,)
        Rate constant of each reaction.
    reactant_orders : ndarray, shape (n_species, n_reactions)
        Reactant-side order of each species in each reaction's mass-action
        rate law (0 for a species that is not a reactant in that
        reaction).

    Returns
    -------
    callable
        ``@njit`` dispatcher ``rhs(state, t, params) -> dstate`` (the
        ``params`` argument is accepted, per
        :data:`chemistrykit.integrators.RHSFunc`, but unused: every
        numeric parameter is already baked into this closure).
    """
    stoich_matrix = np.ascontiguousarray(stoich_matrix, dtype=np.float64)
    rate_constants = np.ascontiguousarray(rate_constants, dtype=np.float64)
    reactant_orders = np.ascontiguousarray(reactant_orders, dtype=np.float64)
    n_species, n_reactions = stoich_matrix.shape

    @njit(cache=False)
    def rhs(state, t, params):
        rates = np.empty(n_reactions)
        for j in range(n_reactions):
            r = rate_constants[j]
            for i in range(n_species):
                order = reactant_orders[i, j]
                if order != 0.0:
                    r *= state[i] ** order
            rates[j] = r
        dstate = np.zeros(n_species)
        for i in range(n_species):
            acc = 0.0
            for j in range(n_reactions):
                acc += stoich_matrix[i, j] * rates[j]
            dstate[i] = acc
        return dstate

    return rhs


class StoichiometricNetwork(ReactionNetwork):
    """A network of elementary mass-action reactions, integrated numerically.

    Parameters
    ----------
    species : sequence of str
        Ordered species names.
    stoich_matrix : array-like, shape (n_species, n_reactions)
        Net stoichiometric coefficient of each species in each reaction
        (negative for a net-consumed species, positive for net-produced).
    rate_constants : array-like, shape (n_reactions,)
        Rate constant of each reaction.
    reactant_orders : array-like, shape (n_species, n_reactions)
        Reactant-side order of each species in each reaction's mass-action
        rate law (0 for a species that is not a reactant in that
        reaction; for an elementary reaction this equals the reactant's
        stoichiometric coefficient, but is given independently so
        non-elementary/pseudo-order rate laws can also be modeled).
    state0 : array-like, shape (n_species,)
        Initial concentrations, in `species` order.
    """

    def __init__(
        self,
        species: Sequence[str],
        stoich_matrix,
        rate_constants,
        reactant_orders,
        state0,
    ):
        self.species = tuple(species)
        self.stoich_matrix = np.asarray(stoich_matrix, dtype=np.float64)
        self.rate_constants = np.asarray(rate_constants, dtype=np.float64)
        self.reactant_orders = np.asarray(reactant_orders, dtype=np.float64)
        self._rhs_njit = _make_mass_action_rhs(self.stoich_matrix, self.rate_constants, self.reactant_orders)
        self.params = np.empty(0)
        super().__init__(state0)

    def rhs(self, state: NDArray[np.float64], t: float) -> NDArray[np.float64]:
        return np.asarray(self._rhs_njit(np.asarray(state, dtype=np.float64), t, self.params))

    @classmethod
    def parallel(cls, k1: float, k2: float, A0: float = 1.0) -> StoichiometricNetwork:
        """Build the parallel (competing) mechanism ``A -> B`` (k1), ``A -> C`` (k2).

        Parameters
        ----------
        k1, k2 : float
            First-order rate constants of the two competing channels.
        A0 : float
            Initial concentration of A.

        Returns
        -------
        StoichiometricNetwork

        Examples
        --------
        The product ratio [B]/[C] is exactly k1/k2 at every time, the
        classic parallel-reaction result (Espenson, *Chemical Kinetics
        and Reaction Mechanisms*, 2nd ed., Ch. 4):

        >>> import numpy as np
        >>> net = StoichiometricNetwork.parallel(k1=2.0, k2=1.0, A0=1.0)
        >>> result = net.integrate((0.0, 5.0), dt=1e-3, method="rk4")
        >>> round(float(result.concentration("B")[-1] / result.concentration("C")[-1]), 4)
        2.0
        """
        stoich = [[-1.0, -1.0], [1.0, 0.0], [0.0, 1.0]]
        orders = [[1.0, 1.0], [0.0, 0.0], [0.0, 0.0]]
        return cls(("A", "B", "C"), stoich, [k1, k2], orders, [A0, 0.0, 0.0])

    @classmethod
    def consecutive(cls, k1: float, k2: float, A0: float = 1.0) -> StoichiometricNetwork:
        """Build the consecutive (chain) mechanism ``A -> B`` (k1), ``B -> C`` (k2).

        Parameters
        ----------
        k1, k2 : float
            First-order rate constants of the two steps.
        A0 : float
            Initial concentration of A.

        Returns
        -------
        StoichiometricNetwork

        Examples
        --------
        >>> net = StoichiometricNetwork.consecutive(k1=1.0, k2=0.3, A0=1.0)
        >>> result = net.integrate((0.0, 20.0), dt=1e-3, method="rk4")
        >>> round(float(result.concentration("A")[-1] + result.concentration("B")[-1] + result.concentration("C")[-1]), 6)
        1.0
        """
        stoich = [[-1.0, 0.0], [1.0, -1.0], [0.0, 1.0]]
        orders = [[1.0, 0.0], [0.0, 1.0], [0.0, 0.0]]
        return cls(("A", "B", "C"), stoich, [k1, k2], orders, [A0, 0.0, 0.0])

    @classmethod
    def reversible(cls, kf: float, kr: float, A0: float = 1.0, B0: float = 0.0) -> StoichiometricNetwork:
        """Build the reversible mechanism ``A <-> B`` (forward kf, reverse kr).

        Parameters
        ----------
        kf, kr : float
            Forward and reverse first-order rate constants.
        A0, B0 : float
            Initial concentrations.

        Returns
        -------
        StoichiometricNetwork

        Examples
        --------
        >>> net = StoichiometricNetwork.reversible(kf=2.0, kr=1.0, A0=1.0)
        >>> result = net.integrate((0.0, 20.0), dt=1e-3, method="rk4")
        >>> round(float(result.concentration("A")[-1]), 3)
        0.333
        """
        stoich = [[-1.0, 1.0], [1.0, -1.0]]
        orders = [[1.0, 0.0], [0.0, 1.0]]
        return cls(("A", "B"), stoich, [kf, kr], orders, [A0, B0])


def consecutive_analytic(A0: float, k1: float, k2: float, t):
    """Closed-form ``A -> B -> C`` concentrations (Bateman's equations).

    See Espenson, *Chemical Kinetics and Reaction Mechanisms*, 2nd ed.,
    Ch. 4.2, or the original: H. Bateman, *Proc. Cambridge Philos. Soc.*
    15, 423 (1910) (radioactive decay chains have the identical
    mathematical structure). The ``k1 != k2`` formula below has a
    removable singularity at ``k1 == k2``; the degenerate limit
    :math:`[B](t) = A_0 k_1 t\\,e^{-k_1 t}` is used automatically in that
    case (obtained by L'Hopital's rule / directly integrating
    :math:`d[B]/dt = k_1[A] - k_1[B]` with equal rate constants).

    Parameters
    ----------
    A0 : float
        Initial concentration of A.
    k1, k2 : float
        Rate constants of ``A -> B`` and ``B -> C`` respectively.
    t : float or array-like of float
        Time(s) at which to evaluate the concentrations.

    Returns
    -------
    A, B, C : float or ndarray

    Examples
    --------
    Mass balance holds exactly at every time:

    >>> import numpy as np
    >>> A, B, C = consecutive_analytic(A0=1.0, k1=1.0, k2=0.3, t=np.array([0.0, 1.0, 5.0, 20.0]))
    >>> np.allclose(A + B + C, 1.0)
    True
    """
    t = np.asarray(t, dtype=np.float64)
    A = A0 * np.exp(-k1 * t)
    if k1 == k2:
        B = A0 * k1 * t * np.exp(-k1 * t)
    else:
        B = A0 * k1 / (k2 - k1) * (np.exp(-k1 * t) - np.exp(-k2 * t))
    C = A0 - A - B
    return A, B, C


def ssa_intermediate_concentration(A0: float, k1: float, k2: float, t):
    r"""Steady-state-approximation (SSA) estimate of the intermediate [B].

    When the second step is much faster than the first (:math:`k_2 \gg
    k_1`), the intermediate B is consumed almost as fast as it forms, so
    :math:`d[B]/dt \approx 0` at all times after a brief induction period
    -- the steady-state approximation (Espenson, *Chemical Kinetics and
    Reaction Mechanisms*, 2nd ed., Ch. 5.1). Applying it to
    :math:`d[B]/dt = k_1[A] - k_2[B] = 0` gives
    :math:`[B]_{ssa} = (k_1/k_2)[A](t) = (k_1/k_2) A_0 e^{-k_1 t}`, which
    should agree closely with :func:`consecutive_analytic`'s exact
    :math:`[B](t)` whenever :math:`k_2/k_1` is large.

    Parameters
    ----------
    A0 : float
        Initial concentration of A.
    k1, k2 : float
        Rate constants of ``A -> B`` and ``B -> C`` respectively.
    t : float or array-like of float
        Time(s) at which to evaluate the approximation.

    Returns
    -------
    float or ndarray

    Examples
    --------
    The SSA error shrinks as k2/k1 grows:

    >>> import numpy as np
    >>> t = np.linspace(0.5, 10.0, 50)
    >>> errs = []
    >>> for k2 in (2.0, 10.0, 100.0):
    ...     _, B_exact, _ = consecutive_analytic(A0=1.0, k1=1.0, k2=k2, t=t)
    ...     B_ssa = ssa_intermediate_concentration(A0=1.0, k1=1.0, k2=k2, t=t)
    ...     errs.append(np.max(np.abs(B_exact - B_ssa)))
    >>> bool(errs[0] > errs[1] > errs[2])
    True
    """
    t = np.asarray(t, dtype=np.float64)
    return (k1 / k2) * A0 * np.exp(-k1 * t)


def reversible_analytic(A0: float, kf: float, kr: float, t, B0: float = 0.0):
    r"""Closed-form ``A <-> B`` relaxation to equilibrium.

    For :math:`d[A]/dt = -k_f[A] + k_r[B]` with :math:`[A]+[B] = A_0+B_0`
    conserved, the solution relaxes exponentially to equilibrium with
    relaxation rate :math:`k_f + k_r` (Atkins & de Paula, *Physical
    Chemistry*, 11th ed., Ch. 20.4):

    .. math::

        [A](t) = [A]_{eq} + ([A]_0 - [A]_{eq}) e^{-(k_f+k_r)t}, \quad
        [A]_{eq} = \frac{k_r ([A]_0+[B]_0)}{k_f+k_r}

    Parameters
    ----------
    A0 : float
        Initial concentration of A.
    kf, kr : float
        Forward and reverse rate constants.
    t : float or array-like of float
        Time(s) at which to evaluate the concentrations.
    B0 : float
        Initial concentration of B.

    Returns
    -------
    A, B : float or ndarray

    Examples
    --------
    >>> import numpy as np
    >>> A, B = reversible_analytic(A0=1.0, kf=2.0, kr=1.0, t=1e6)
    >>> round(float(A), 3), round(float(B), 3)
    (0.333, 0.667)
    """
    t = np.asarray(t, dtype=np.float64)
    total = A0 + B0
    A_eq = kr * total / (kf + kr)
    A = A_eq + (A0 - A_eq) * np.exp(-(kf + kr) * t)
    B = total - A
    return A, B
