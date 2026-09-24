r"""Living (termination-free) polymerization and the Poisson chain-length distribution.

M. Szwarc, "'Living' Polymers," *Nature* 178, 1168 (1956); the
distribution itself is P. J. Flory, *J. Am. Chem. Soc.* 62, 1561 (1940).

If every chain is started at once by a fast initiation and then grows by
adding monomer with no termination or transfer, each chain accumulates
monomers as independent random events, so the number :math:`k` of
monomers added to a chain is Poisson distributed with mean :math:`\nu`
(monomers consumed per initiator). Counting the initiator-bearing first
unit, a chain has :math:`x=k+1` units with number fraction

.. math::

    N_x = \frac{e^{-\nu}\nu^{x-1}}{(x-1)!}, \qquad
    \bar X_n = 1+\nu, \qquad
    \frac{\bar X_w}{\bar X_n} = 1+\frac{\nu}{(1+\nu)^2}

so the polydispersity index tends to 1 (roughly :math:`1+1/\nu`) for
long chains -- far narrower than the Flory-Schulz limit of 2.
"""

from __future__ import annotations

import numpy as np
from scipy.special import gammaln

__all__ = [
    "poisson_number_fraction",
    "poisson_number_average_DP",
    "poisson_weight_average_DP",
    "poisson_pdi",
    "simulate_living_polymerization",
]


def poisson_number_fraction(x, nu: float):
    r"""Poisson number fraction :math:`N_x=e^{-\nu}\nu^{x-1}/(x-1)!` of chains with `x` units.

    Parameters
    ----------
    x : int or array-like of int
        Chain length (:math:`\ge1`).
    nu : float
        Mean number of monomers added per chain.

    Returns
    -------
    float or ndarray

    Examples
    --------
    >>> import numpy as np
    >>> round(float(np.sum(poisson_number_fraction(np.arange(1, 400), 50.0))), 10)
    1.0
    """
    k = np.asarray(x, dtype=float) - 1.0
    return np.exp(-nu + k * np.log(nu) - gammaln(k + 1.0))


def poisson_number_average_DP(nu: float) -> float:
    r"""Number-average degree of polymerization, :math:`1+\nu`.

    Parameters
    ----------
    nu : float

    Returns
    -------
    float

    Examples
    --------
    >>> poisson_number_average_DP(99.0)
    100.0
    """
    return 1.0 + nu


def poisson_weight_average_DP(nu: float) -> float:
    r"""Weight-average degree of polymerization, :math:`(\nu^2+3\nu+1)/(\nu+1)`.

    Parameters
    ----------
    nu : float

    Returns
    -------
    float

    Examples
    --------
    >>> round(poisson_weight_average_DP(1.0), 6)
    2.5
    """
    return (nu**2 + 3.0 * nu + 1.0) / (nu + 1.0)


def poisson_pdi(nu: float) -> float:
    r"""Polydispersity index of the Poisson distribution, :math:`1+\nu/(1+\nu)^2`.

    Parameters
    ----------
    nu : float

    Returns
    -------
    float

    Examples
    --------
    >>> round(poisson_pdi(99.0), 6)
    1.0099
    """
    return 1.0 + nu / (1.0 + nu) ** 2


def simulate_living_polymerization(n_chains: int, n_monomers: int, rng=None) -> np.ndarray:
    r"""Stochastically grow `n_chains` living chains by adding `n_monomers` one at a time to random chains.

    Each monomer attaches to a uniformly chosen active chain (equal
    reactivity, no termination), so chain lengths follow a multinomial
    distribution that tends to the Poisson distribution
    (:func:`poisson_number_fraction`) with :math:`\nu=` `n_monomers` /
    `n_chains`.

    Parameters
    ----------
    n_chains : int
        Number of initiated chains (each starts with one unit).
    n_monomers : int
        Number of monomers added in total.
    rng : int, numpy.random.Generator, or None, optional
        Seed or generator.

    Returns
    -------
    ndarray of int, shape (n_chains,)
        Final chain lengths :math:`x` (units per chain, including the first).

    Examples
    --------
    >>> x = simulate_living_polymerization(1000, 50000, rng=0)
    >>> int(x.sum())
    51000
    """
    rng = np.random.default_rng(rng)
    counts = np.bincount(rng.integers(0, n_chains, size=n_monomers), minlength=n_chains)
    return counts + 1
