r"""Stereoregularity (tacticity) of vinyl polymers: Bernoullian dyad/triad statistics.

In a vinyl polymer such as polypropylene, consecutive stereocenters form
either *meso* (m, same configuration) or *racemo* (r) dyads. Natta's
stereospecific Ziegler-type catalysts (G. Natta et al., *J. Am. Chem. Soc.*
77, 1708 (1955)) make nearly all dyads meso (isotactic), whereas a
free-radical mechanism places each monomer with nearly random
orientation (atactic). If each placement is an independent event with
meso probability :math:`P_m` (Bernoullian statistics; F. A. Bovey & G. V.
D. Tiers, *J. Polym. Sci.* 44, 173 (1960)), the triad fractions
observable by NMR are

.. math::

    [mm]=P_m^2, \qquad [mr]=2P_m(1-P_m), \qquad [rr]=(1-P_m)^2

and the mean isotactic run length (consecutive m dyads) is
:math:`1/(1-P_m)`.
"""

from __future__ import annotations

import numpy as np

__all__ = [
    "bernoullian_triad_fractions",
    "mean_isotactic_run_length",
    "sample_dyad_sequence",
]


def bernoullian_triad_fractions(Pm) -> tuple:
    r"""Triad fractions :math:`([mm],[mr],[rr])` for Bernoullian placement with meso probability `Pm`.

    Parameters
    ----------
    Pm : float or array-like of float
        Probability of a meso placement.

    Returns
    -------
    mm, mr, rr : float or ndarray
        They sum to 1.

    Examples
    --------
    Atactic (:math:`P_m=1/2`) gives the 1:2:1 ratio:

    >>> tuple(round(float(v), 3) for v in bernoullian_triad_fractions(0.5))
    (0.25, 0.5, 0.25)
    """
    Pm = np.asarray(Pm, dtype=float)
    return Pm**2, 2.0 * Pm * (1.0 - Pm), (1.0 - Pm) ** 2


def mean_isotactic_run_length(Pm: float) -> float:
    r"""Mean length of a run of consecutive meso dyads, :math:`1/(1-P_m)` (Bernoullian).

    Parameters
    ----------
    Pm : float
        Meso placement probability, :math:`<1`.

    Returns
    -------
    float

    Examples
    --------
    >>> round(mean_isotactic_run_length(0.99), 6)
    100.0
    """
    return 1.0 / (1.0 - Pm)


def sample_dyad_sequence(n_dyads: int, Pm: float, rng=None) -> np.ndarray:
    r"""Sample a Bernoullian dyad sequence: ``True`` for meso, ``False`` for racemo.

    Parameters
    ----------
    n_dyads : int
    Pm : float
        Meso placement probability.
    rng : int, numpy.random.Generator, or None, optional

    Returns
    -------
    ndarray of bool, shape (n_dyads,)

    Examples
    --------
    >>> s = sample_dyad_sequence(100000, 0.9, rng=0)
    >>> bool(abs(s.mean() - 0.9) < 0.01)
    True
    """
    rng = np.random.default_rng(rng)
    return rng.random(n_dyads) < Pm
