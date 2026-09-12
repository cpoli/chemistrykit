r"""Generic number- and weight-average moment helpers for a chain-length/molar-mass distribution.

Both a measured molar-mass distribution
(:func:`chemistrykit.polymer.systems.molecular_weight_distribution.number_average_molar_mass` /
:func:`~chemistrykit.polymer.systems.molecular_weight_distribution.weight_average_molar_mass`)
and a numerical check of the analytic Flory-Schulz distribution's
number-/weight-average degree of polymerization (see
:mod:`chemistrykit.polymer.tests.test_molecular_weight_distribution`)
reduce to the same two moment ratios of a distribution over chain length
(or molar mass) -- this is that one shared numerical routine, following
the house convention of factoring numerics used by more than one
``systems/`` computation into ``utils/`` (cf.
``chemistrykit.kinetics.utils.linear_regression``,
``chemistrykit.thermo.utils.cubic_roots``).
"""

from __future__ import annotations

import numpy as np

__all__ = ["number_average", "weight_average"]


def number_average(x, N_x) -> float:
    r"""Number-average of `x` weighted by population counts `N_x`: :math:`\bar{x}_n=\sum N_ix_i/\sum N_i`.

    Parameters
    ----------
    x : array-like of float
        Values of the quantity being averaged (e.g. chain length or molar
        mass of each species/fraction `i`).
    N_x : array-like of float
        Number of chains (or moles) at each value of `x`.

    Returns
    -------
    float

    Examples
    --------
    >>> number_average(x=[1.0, 2.0, 3.0], N_x=[1.0, 1.0, 1.0])
    2.0
    """
    x = np.asarray(x, dtype=np.float64)
    N_x = np.asarray(N_x, dtype=np.float64)
    return float(np.sum(N_x * x) / np.sum(N_x))


def weight_average(x, N_x) -> float:
    r"""Weight-average of `x` weighted by population counts `N_x`: :math:`\bar{x}_w=\sum N_ix_i^2/\sum N_ix_i`.

    Weighting by :math:`N_ix_i` (the total mass/length in fraction `i`,
    not just its number of chains) rather than by :math:`N_i` alone is
    what distinguishes a weight average from a number average -- longer
    chains contribute more to the weight average because there is simply
    more material at that chain length (Odian, *Principles of
    Polymerization*, 4th ed., Ch. 2.3).

    Parameters
    ----------
    x, N_x : array-like of float
        As in :func:`number_average`.

    Returns
    -------
    float

    Examples
    --------
    A monodisperse population has equal number- and weight-averages:

    >>> weight_average(x=[5.0, 5.0, 5.0], N_x=[1.0, 1.0, 1.0])
    5.0

    A polydisperse population's weight-average exceeds its number-average
    (longer chains are over-represented by mass):

    >>> x, N_x = [1.0, 2.0, 3.0], [3.0, 2.0, 1.0]
    >>> weight_average(x, N_x) > number_average(x, N_x)
    True
    """
    x = np.asarray(x, dtype=np.float64)
    N_x = np.asarray(N_x, dtype=np.float64)
    return float(np.sum(N_x * x**2) / np.sum(N_x * x))
