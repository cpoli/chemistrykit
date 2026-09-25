r"""Gillespie's stochastic simulation algorithm (SSA) for reaction networks.

At low copy numbers a reaction network is a continuous-time Markov jump
process on integer molecule counts, not a set of ODEs for smooth
concentrations. Gillespie's "direct method" samples exact trajectories of
that process: with propensities :math:`a_j(\mathbf{n})` and total
:math:`a_0 = \sum_j a_j`, the waiting time to the next reaction is
exponential with rate :math:`a_0`, and reaction :math:`j` is the one that
fires with probability :math:`a_j/a_0` (D. T. Gillespie, *J. Comput.
Phys.* 22, 403 (1976); *J. Phys. Chem.* 81, 2340 (1977)).

Mass-action propensities use the combinatorial count of distinct reactant
combinations, :math:`a_j = c_j \prod_i \binom{n_i}{\nu_{ij}}`, where
:math:`c_j` is the stochastic rate constant of reaction :math:`j`.
"""

from __future__ import annotations

from collections.abc import Sequence
from dataclasses import dataclass, field
from math import factorial

import numpy as np

__all__ = ["StochasticTrajectory", "gillespie_ssa"]


@dataclass
class StochasticTrajectory:
    """One exact sample path of a stochastic reaction network."""

    t: np.ndarray
    """ndarray, shape (n_events + 1,): Reaction-event times (starting at 0)."""

    counts: np.ndarray
    """ndarray of int, shape (n_events + 1, n_species): Molecule counts
    immediately after each event (row 0 is the initial state)."""

    species: Sequence[str] = field(default_factory=tuple)
    """tuple of str: Species names, in column order matching ``counts``."""

    def count(self, name: str) -> np.ndarray:
        """Return the count trajectory of a single named species.

        Parameters
        ----------
        name : str

        Returns
        -------
        ndarray of int
        """
        return self.counts[:, list(self.species).index(name)]

    def sample(self, times) -> np.ndarray:
        """Evaluate the piecewise-constant trajectory at arbitrary times.

        Parameters
        ----------
        times : array-like of float
            Times in ``[0, t_max]``.

        Returns
        -------
        ndarray of int, shape (len(times), n_species)
        """
        idx = np.searchsorted(self.t, np.asarray(times, dtype=np.float64), side="right") - 1
        return self.counts[np.clip(idx, 0, len(self.t) - 1)]


def gillespie_ssa(
    stoich_matrix,
    rate_constants,
    reactant_orders,
    n0,
    t_max: float,
    species: Sequence[str] | None = None,
    seed=None,
    max_events: int = 10_000_000,
) -> StochasticTrajectory:
    r"""Simulate a mass-action network exactly with Gillespie's direct method.

    Parameters
    ----------
    stoich_matrix : array-like of int, shape (n_species, n_reactions)
        Net change in each species' count when each reaction fires.
    rate_constants : array-like of float, shape (n_reactions,)
        Stochastic rate constants :math:`c_j`.
    reactant_orders : array-like of int, shape (n_species, n_reactions)
        Number of molecules of each species consumed as reactants by each
        reaction (0 if not a reactant).
    n0 : array-like of int, shape (n_species,)
        Initial molecule counts.
    t_max : float
        Simulate until this time (or until no reaction can fire).
    species : sequence of str, optional
        Species names; defaults to ``("S0", "S1", ...)``.
    seed : int or numpy.random.Generator, optional
        Seed or generator for reproducible trajectories.
    max_events : int, default 10_000_000
        Safety cap on the number of reaction events.

    Returns
    -------
    StochasticTrajectory

    Examples
    --------
    Pure decay ``A -> 0`` of 50 molecules: counts only ever drop by one,
    and a single trajectory ends with every molecule gone:

    >>> import numpy as np
    >>> traj = gillespie_ssa([[-1]], [1.0], [[1]], n0=[50], t_max=100.0, seed=0)
    >>> int(traj.count("S0")[0]), int(traj.count("S0")[-1]), len(traj.t)
    (50, 0, 51)
    >>> bool(np.all(np.diff(traj.count("S0")) == -1))
    True
    """
    stoich = np.asarray(stoich_matrix, dtype=np.int64)
    c = np.asarray(rate_constants, dtype=np.float64)
    orders = np.asarray(reactant_orders, dtype=np.int64)
    n = np.asarray(n0, dtype=np.int64).copy()
    n_species, n_reactions = stoich.shape
    if species is None:
        species = tuple(f"S{i}" for i in range(n_species))
    rng = seed if isinstance(seed, np.random.Generator) else np.random.default_rng(seed)

    # (species, reaction, order, 1/order!) for every reactant entry
    terms = [(i, j, int(orders[i, j]), 1.0 / factorial(int(orders[i, j]))) for i, j in zip(*np.nonzero(orders), strict=True)]

    ts = [0.0]
    states = [n.copy()]
    t = 0.0
    for _ in range(max_events):
        a = c.copy()
        for i, j, order, inv_fact in terms:
            falling = 1.0
            for m in range(order):
                falling *= n[i] - m
            a[j] *= max(falling, 0.0) * inv_fact
        a0 = a.sum()
        if a0 <= 0.0:
            break
        t += rng.exponential(1.0 / a0)
        if t > t_max:
            break
        j = int(np.searchsorted(np.cumsum(a), rng.random() * a0, side="right"))
        n += stoich[:, min(j, n_reactions - 1)]
        ts.append(t)
        states.append(n.copy())
    return StochasticTrajectory(t=np.array(ts), counts=np.array(states), species=tuple(species))
