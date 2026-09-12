r"""Evjen's method: a genuinely converging summation for conditionally-convergent ionic lattice sums.

The electrostatic (Madelung) lattice sum for an ionic crystal,

.. math::

    S = \sum_{j\neq0} \frac{q_j/q_0}{r_j/r_0}

is only *conditionally* convergent: it does not have a single well-defined
value as an unordered sum, so its numerical value depends on the order/
shape in which terms are added, and a naive sum over an expanding cube or
sphere of lattice points does not converge as the cutoff grows -- it
oscillates, because the outermost shell added at each step is never
electrically neutral by itself. This module implements Evjen's method
(H. M. Evjen, *Phys. Rev.* 39, 675 (1932)), the standard fix: terms whose
lattice indices lie exactly on the truncation cube's boundary are given
fractional weight (1/2 for a face, 1/4 for an edge, 1/8 for a corner of
the cube) so that every finite partial sum is built from electrically
neutral "shells" (each truncation cube plus its fractionally-weighted
boundary carries zero net charge), which restores genuine, rapidly
converging convergence to the *correct* conditionally-convergent limit.
See also West, *Solid State Chemistry and its Applications*, 2nd ed.
(2014), Ch. 1.4, for the same method applied to NaCl.

This is verified against the literature NaCl value in
:mod:`chemistrykit.crystal.systems.madelung`; see that module's docstring
for a demonstration that the naive (unweighted) truncated sum this module
avoids does *not* converge.
"""

from __future__ import annotations

import numpy as np

__all__ = ["evjen_lattice_sum_cubic_alternating"]


def evjen_lattice_sum_cubic_alternating(n_shells: int) -> float:
    r"""Evjen-weighted sum, over a simple cubic lattice of alternating unit charges, of :math:`\pm1/r`.

    Sums :math:`\sum (-1)^{i+j+k}/\sqrt{i^2+j^2+k^2}` over all integer
    triples :math:`(i,j,k)\neq(0,0,0)` with :math:`|i|,|j|,|k|\le n` (`n`
    = `n_shells`), Evjen-weighting any term with one or more indices
    exactly at the cutoff :math:`\pm n` by :math:`1/2^{(\text{number of
    boundary indices})}` (so a face term counts 1/2, an edge term 1/4, a
    corner term 1/8) so that the finite region summed is charge-neutral.

    This particular alternating-sign-by-parity charge pattern on a simple
    cubic lattice is exactly the combined Na+/Cl- ion arrangement of the
    NaCl (rock-salt) structure once both interpenetrating FCC sublattices
    are re-expressed on the finer simple-cubic mesh of half the
    conventional cell edge (Ashcroft & Mermin, *Solid State Physics*,
    1976, problem 20.1) -- see :mod:`chemistrykit.crystal.systems.madelung`
    for that identification and the resulting Madelung constant.

    Parameters
    ----------
    n_shells : int
        Half-width of the cubic summation region, in lattice units (must
        be >= 1). Larger values converge closer to the true (infinite-
        lattice) conditionally-convergent limit.

    Returns
    -------
    float

    Examples
    --------
    The magnitude grows and stabilizes quickly with `n_shells` (contrast
    with an unweighted/naive truncated sum, which oscillates indefinitely
    -- see :mod:`chemistrykit.crystal.systems.madelung`):

    >>> vals = [evjen_lattice_sum_cubic_alternating(n) for n in (4, 8, 12)]
    >>> bool(abs(vals[2] - vals[1]) < abs(vals[1] - vals[0]))
    True
    """
    if n_shells < 1:
        raise ValueError("n_shells must be >= 1")
    n = n_shells
    idx = np.arange(-n, n + 1)
    i, j, k = np.meshgrid(idx, idx, idx, indexing="ij")
    r = np.sqrt(i.astype(np.float64) ** 2 + j.astype(np.float64) ** 2 + k.astype(np.float64) ** 2)
    origin = (i == 0) & (j == 0) & (k == 0)
    r_safe = np.where(origin, 1.0, r)
    sign = np.where((i + j + k) % 2 == 0, 1.0, -1.0)
    n_boundary = (np.abs(i) == n).astype(np.int64) + (np.abs(j) == n).astype(np.int64) + (np.abs(k) == n).astype(np.int64)
    weight = 1.0 / (2.0**n_boundary)
    term = np.where(origin, 0.0, weight * sign / r_safe)
    return float(np.sum(term))
