r"""Periodic-boundary neighbor lists for pairwise force evaluation.

See Allen & Tildesley, *Computer Simulation of Liquids*, 2nd ed., Ch. 5.3
(the Verlet neighbor list; L. Verlet, *Phys. Rev.* 159, 98 (1967)), or
Frenkel & Smit, *Understanding Molecular Simulation*, 2nd ed., Ch. 5.3.
"""

from __future__ import annotations

from typing import Optional

import numpy as np

from chemistrykit.md.utils.pbc import minimum_image_displacement

__all__ = ["build_neighbor_list", "VerletNeighborList"]


def build_neighbor_list(positions, box_length: Optional[float], cutoff: float):
    """Build the list of all particle pairs within `cutoff` of each other.

    A direct, vectorized :math:`O(N^2)` all-pairs distance computation
    under the minimum-image convention -- exact, and fast enough for the
    particle counts this package targets (up to a few thousand). A
    production MD code instead uses a cell (linked-list) decomposition to
    reduce this to :math:`O(N)`; that further optimization is out of
    scope here. See :class:`VerletNeighborList` for the one optimization
    this module *does* provide: reusing a list built with a distance
    "skin" across several integration steps instead of rebuilding every
    step.

    Parameters
    ----------
    positions : ndarray, shape (N, d)
    box_length : float or None
        Cubic periodic box side length; ``None`` or ``<= 0`` for an
        unbounded (non-periodic) system.
    cutoff : float
        Pairs farther apart than this are excluded.

    Returns
    -------
    pairs_i, pairs_j : ndarray of int64, shape (n_pairs,)
        Indices such that ``pairs_i < pairs_j`` for every returned pair.

    Examples
    --------
    >>> import numpy as np
    >>> positions = np.array([[0.0, 0.0, 0.0], [0.5, 0.0, 0.0], [5.0, 0.0, 0.0]])
    >>> pairs_i, pairs_j = build_neighbor_list(positions, box_length=None, cutoff=1.0)
    >>> list(zip(pairs_i.tolist(), pairs_j.tolist()))
    [(0, 1)]
    """
    positions = np.asarray(positions, dtype=np.float64)
    n = positions.shape[0]
    iu, ju = np.triu_indices(n, k=1)
    if iu.size == 0:
        return iu.astype(np.int64), ju.astype(np.int64)
    diff = minimum_image_displacement(positions[iu], positions[ju], box_length)
    r2 = np.sum(diff * diff, axis=-1)
    mask = r2 < cutoff * cutoff
    return iu[mask].astype(np.int64), ju[mask].astype(np.int64)


class VerletNeighborList:
    """A Verlet ("skin") neighbor list, rebuilt only every few calls.

    Building the pair list with a cutoff enlarged by a `skin` distance
    means the list stays a superset of the true `cutoff`-pairs for
    several integration steps, since a pair can only enter the true
    cutoff sphere after its separation has changed by up to `skin` --
    the standard MD performance trick (Allen & Tildesley, *Computer
    Simulation of Liquids*, 2nd ed., Ch. 5.3). This implementation
    rebuilds unconditionally every `rebuild_every` calls rather than
    tracking accumulated displacement, which is simpler but requires the
    caller to choose `rebuild_every`/`skin` conservatively for their
    timestep and temperature -- a stricter implementation would instead
    track the two largest particle displacements since the last rebuild
    and force an early rebuild once their sum exceeds `skin`.

    Parameters
    ----------
    cutoff : float
        The physical interaction cutoff.
    skin : float, default 0.3
        Extra distance added to `cutoff` when building the list.
    rebuild_every : int, default 20
        Rebuild after this many calls to :meth:`pairs`.
    """

    def __init__(self, cutoff: float, skin: float = 0.3, rebuild_every: int = 20):
        if cutoff <= 0:
            raise ValueError("cutoff must be positive")
        if skin < 0:
            raise ValueError("skin must be non-negative")
        self.cutoff = float(cutoff)
        self.skin = float(skin)
        self.rebuild_every = int(rebuild_every)
        self._pairs = None
        self._calls_since_rebuild = 0

    def pairs(self, positions, box_length: Optional[float]):
        """Return (possibly cached) ``(pairs_i, pairs_j)`` for `positions`.

        Parameters
        ----------
        positions : ndarray, shape (N, d)
        box_length : float or None

        Returns
        -------
        pairs_i, pairs_j : ndarray of int64
        """
        if self._pairs is None or self._calls_since_rebuild >= self.rebuild_every:
            self._pairs = build_neighbor_list(positions, box_length, self.cutoff + self.skin)
            self._calls_since_rebuild = 0
        self._calls_since_rebuild += 1
        return self._pairs

    def reset(self) -> None:
        """Force the next call to :meth:`pairs` to rebuild."""
        self._pairs = None
        self._calls_since_rebuild = 0
