r"""Molecular-weight-distribution statistics, and the Flory-Schulz (most-probable) distribution.

See Odian, *Principles of Polymerization*, 4th ed., Ch. 2.3 ("Molecular
weight distribution in linear step polymerization"), and the original:
P. J. Flory, *J. Am. Chem. Soc.* 58, 1877 (1936).

A real polymer sample is a mixture of chains of many different lengths;
:func:`number_average_molar_mass` (:math:`M_n`) and
:func:`weight_average_molar_mass` (:math:`M_w`) are the two standard
summary moments of that distribution, and their ratio
:func:`polydispersity_index` (PDI, :math:`\text{\dj}=M_w/M_n\ge1`)
quantifies how broad it is (PDI = 1 for a perfectly monodisperse sample).

For an ideal linear step-growth polymerization carried to extent of
reaction :math:`p` (see :mod:`chemistrykit.polymer.systems.step_growth`),
elementary statistics of the condensation process give the chain-length
distribution in closed form -- the **Flory-Schulz** (or "most probable")
distribution: the probability that a randomly chosen *chain* (not
monomer unit) has exactly `x` repeat units is
:math:`(1-p)p^{x-1}` (a geometric distribution: `x-1` successful
condensation steps followed by one chain end that didn't react further),
and the corresponding mass/weight fraction is
:math:`w_x=x(1-p)^2p^{x-1}`. Their moments,

.. math::

    \bar{X}_n = \sum_x xN_x\Big/\sum_x N_x = \frac{1}{1-p}, \qquad
    \bar{X}_w = \sum_x xw_x = \frac{1+p}{1-p}

reproduce the Carothers equation for :math:`\bar X_n` (see
:mod:`chemistrykit.polymer.systems.step_growth`) and give
:math:`\text{PDI}=\bar X_w/\bar X_n=1+p\to2` as :math:`p\to1` -- the
textbook result that an ideal step-growth polymerization's PDI approaches
(but never exceeds) 2, however far the reaction is driven.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from chemistrykit.polymer.utils.moments import number_average, weight_average

__all__ = [
    "number_average_molar_mass",
    "weight_average_molar_mass",
    "polydispersity_index",
    "MolecularWeightDistribution",
    "flory_schulz_number_fraction",
    "flory_schulz_weight_fraction",
    "flory_schulz_number_average_DP",
    "flory_schulz_weight_average_DP",
    "flory_schulz_pdi",
]


def number_average_molar_mass(N_i, M_i) -> float:
    r"""Number-average molar mass :math:`M_n=\sum N_iM_i/\sum N_i`.

    Parameters
    ----------
    N_i : array-like of float
        Number (or mole count) of chains in each molar-mass fraction.
    M_i : array-like of float
        Molar mass of each fraction.

    Returns
    -------
    float

    Examples
    --------
    >>> number_average_molar_mass(N_i=[10.0, 5.0], M_i=[1000.0, 2000.0])
    1333.3333333333333
    """
    return number_average(M_i, N_i)


def weight_average_molar_mass(N_i, M_i) -> float:
    r"""Weight-average molar mass :math:`M_w=\sum N_iM_i^2/\sum N_iM_i`.

    Parameters
    ----------
    N_i, M_i : array-like of float
        As in :func:`number_average_molar_mass`.

    Returns
    -------
    float

    Examples
    --------
    :math:`M_w \ge M_n` always (Cauchy-Schwarz), with equality only for a
    monodisperse sample:

    >>> N_i, M_i = [10.0, 5.0], [1000.0, 2000.0]
    >>> weight_average_molar_mass(N_i, M_i) >= number_average_molar_mass(N_i, M_i)
    True
    """
    return weight_average(M_i, N_i)


def polydispersity_index(Mn: float, Mw: float) -> float:
    r"""Polydispersity index :math:`\text{\dj} = M_w/M_n`.

    Parameters
    ----------
    Mn, Mw : float
        Number- and weight-average molar mass.

    Returns
    -------
    float
        :math:`\ge 1`, with 1 for a perfectly monodisperse sample.

    Examples
    --------
    >>> polydispersity_index(Mn=1000.0, Mw=1000.0)
    1.0
    """
    return Mw / Mn


@dataclass
class MolecularWeightDistribution:
    """Summary statistics (Mn, Mw, PDI) of a molar-mass distribution."""

    Mn: float
    """float: Number-average molar mass."""

    Mw: float
    """float: Weight-average molar mass."""

    PDI: float
    """float: Polydispersity index, Mw/Mn."""

    @classmethod
    def from_counts(cls, N_i, M_i) -> MolecularWeightDistribution:
        """Build from raw (count, molar mass) fraction data.

        Parameters
        ----------
        N_i, M_i : array-like of float
            As in :func:`number_average_molar_mass`.

        Returns
        -------
        MolecularWeightDistribution

        Examples
        --------
        >>> mwd = MolecularWeightDistribution.from_counts(N_i=[10.0, 5.0], M_i=[1000.0, 2000.0])
        >>> round(mwd.Mn, 2), round(mwd.Mw, 2)
        (1333.33, 1500.0)
        >>> round(mwd.PDI, 4)
        1.125
        """
        Mn = number_average_molar_mass(N_i, M_i)
        Mw = weight_average_molar_mass(N_i, M_i)
        return cls(Mn=Mn, Mw=Mw, PDI=polydispersity_index(Mn, Mw))


def flory_schulz_number_fraction(x, p):
    r"""Flory-Schulz number fraction of chains of length `x`: :math:`N_x=(1-p)p^{x-1}`.

    Parameters
    ----------
    x : int or array-like of int
        Chain length (degree of polymerization), :math:`x\ge1`.
    p : float
        Extent of reaction, :math:`0\le p<1`.

    Returns
    -------
    float or ndarray

    Examples
    --------
    The number fractions over all chain lengths sum to exactly 1 (a
    proper probability distribution):

    >>> import numpy as np
    >>> x = np.arange(1, 5000)
    >>> round(float(np.sum(flory_schulz_number_fraction(x, p=0.98))), 6)
    1.0
    """
    x = np.asarray(x, dtype=np.float64)
    result = (1.0 - p) * p ** (x - 1.0)
    return float(result) if result.ndim == 0 else result


def flory_schulz_weight_fraction(x, p):
    r"""Flory-Schulz weight fraction of chains of length `x`: :math:`w_x=x(1-p)^2p^{x-1}`.

    Parameters
    ----------
    x : int or array-like of int
    p : float

    Returns
    -------
    float or ndarray

    Examples
    --------
    The weight fractions over all chain lengths sum to exactly 1:

    >>> import numpy as np
    >>> x = np.arange(1, 5000)
    >>> round(float(np.sum(flory_schulz_weight_fraction(x, p=0.98))), 6)
    1.0
    """
    x = np.asarray(x, dtype=np.float64)
    result = x * (1.0 - p) ** 2 * p ** (x - 1.0)
    return float(result) if result.ndim == 0 else result


def flory_schulz_number_average_DP(p: float) -> float:
    r"""Flory-Schulz number-average degree of polymerization :math:`\bar X_n=1/(1-p)`.

    Identical to the Carothers equation,
    :func:`chemistrykit.polymer.systems.step_growth.degree_of_polymerization`
    -- this is the same quantity derived from the chain-length
    distribution's first moment rather than from a mass-balance argument.

    Parameters
    ----------
    p : float
        Extent of reaction, :math:`0\le p<1`.

    Returns
    -------
    float

    Examples
    --------
    >>> round(flory_schulz_number_average_DP(0.9), 6)
    10.0
    """
    return 1.0 / (1.0 - p)


def flory_schulz_weight_average_DP(p: float) -> float:
    r"""Flory-Schulz weight-average degree of polymerization :math:`\bar X_w=(1+p)/(1-p)`.

    Parameters
    ----------
    p : float
        Extent of reaction, :math:`0\le p<1`.

    Returns
    -------
    float

    Examples
    --------
    >>> round(flory_schulz_weight_average_DP(0.9), 6)
    19.0
    """
    return (1.0 + p) / (1.0 - p)


def flory_schulz_pdi(p: float) -> float:
    r"""Flory-Schulz polydispersity index :math:`\text{\dj}=\bar X_w/\bar X_n=1+p`.

    Parameters
    ----------
    p : float
        Extent of reaction, :math:`0\le p\le1`.

    Returns
    -------
    float
        In :math:`[1, 2]`.

    Examples
    --------
    The PDI is exactly 1 with no reaction (every "chain" is a lone
    monomer -- trivially monodisperse):

    >>> flory_schulz_pdi(0.0)
    1.0

    The PDI approaches exactly 2 in the high-conversion limit
    (:math:`p\to1`) -- the textbook result that an ideal step-growth
    polymerization's molecular-weight distribution can never exceed a
    polydispersity of 2, however far the reaction is driven:

    >>> flory_schulz_pdi(1.0)
    2.0
    """
    return 1.0 + p
