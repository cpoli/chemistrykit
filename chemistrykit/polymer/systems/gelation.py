r"""Gelation of branching step-growth polymerizations: Carothers and Flory-Stockmayer gel points.

When monomers carry more than two functional groups, step-growth
polymerization produces branched molecules and, at a critical extent of
reaction, an "infinite" network -- the gel point. For the
self-condensation of an :math:`f`-functional monomer :math:`A_f`, Flory's
tree (branching-process) statistics (P. J. Flory, *J. Am. Chem. Soc.*
63, 3083 (1941); W. H. Stockmayer, *J. Chem. Phys.* 11, 45 (1943)) give

.. math::

    \bar X_n = \frac{1}{1-fp/2}, \qquad
    \bar X_w = \frac{1+p}{1-(f-1)p}, \qquad
    p_c = \frac{1}{f-1}

so the *weight*-average degree of polymerization diverges at
:math:`p_c` while the number average stays finite. Carothers's cruder
estimate (*Trans. Faraday Soc.* 32, 39 (1936)) sets :math:`\bar X_n=\infty`
instead, :math:`p_c=2/f_\text{avg}`, overestimating the gel point. With
:math:`f=2` both reduce to linear step growth (Flory-Schulz, :math:`p_c=1`).
"""

from __future__ import annotations

import numpy as np

__all__ = [
    "carothers_gel_point",
    "flory_stockmayer_gel_point",
    "branching_number_average_DP",
    "branching_weight_average_DP",
]


def carothers_gel_point(f_avg: float) -> float:
    r"""Carothers gel point, :math:`p_c=2/f_\text{avg}` (where :math:`\bar X_n\to\infty`).

    Parameters
    ----------
    f_avg : float
        Average monomer functionality.

    Returns
    -------
    float

    Examples
    --------
    >>> round(carothers_gel_point(3.0), 6)
    0.666667
    """
    return 2.0 / f_avg


def flory_stockmayer_gel_point(f: float) -> float:
    r"""Flory-Stockmayer gel point for :math:`A_f` self-condensation, :math:`p_c=1/(f-1)`.

    Parameters
    ----------
    f : float
        Monomer functionality (> 2 for gelation at :math:`p<1`).

    Returns
    -------
    float

    Examples
    --------
    >>> flory_stockmayer_gel_point(3)
    0.5
    """
    return 1.0 / (f - 1.0)


def branching_number_average_DP(p, f: float):
    r"""Number-average degree of polymerization of :math:`A_f` self-condensation, :math:`1/(1-fp/2)`.

    Parameters
    ----------
    p : float or array-like of float
        Extent of reaction (fraction of groups reacted), below the gel point.
    f : float
        Functionality.

    Returns
    -------
    float or ndarray

    Examples
    --------
    Finite at the Flory-Stockmayer gel point, :math:`2(f-1)/(f-2)`:

    >>> round(float(branching_number_average_DP(0.5, f=3)), 6)
    4.0
    """
    return 1.0 / (1.0 - f * np.asarray(p, dtype=float) / 2.0)


def branching_weight_average_DP(p, f: float):
    r"""Weight-average degree of polymerization of :math:`A_f` self-condensation, :math:`(1+p)/(1-(f-1)p)`.

    Parameters
    ----------
    p : float or array-like of float
        Extent of reaction, :math:`p<p_c`.
    f : float
        Functionality.

    Returns
    -------
    float or ndarray

    Examples
    --------
    With ``f=2`` this is the Flory-Schulz :math:`(1+p)/(1-p)`:

    >>> round(float(branching_weight_average_DP(0.9, f=2)), 6)
    19.0
    """
    p = np.asarray(p, dtype=float)
    return (1.0 + p) / (1.0 - (f - 1.0) * p)
