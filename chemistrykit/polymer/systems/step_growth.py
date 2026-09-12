r"""Step-growth (condensation) polymerization kinetics: the Carothers equation.

See Odian, *Principles of Polymerization*, 4th ed., Ch. 2.2 ("Kinetics of
step polymerization"), and the original: W. H. Carothers, *Trans.
Faraday Soc.* 32, 39 (1936).

In an ideal step-growth polymerization (e.g. a diol + diacid
polyesterification), any two molecules with compatible functional groups
-- monomers, dimers, or growing chains of any length -- can react with
each other, so the average chain length grows only slowly, as more and
more of the *functional groups* (not molecules) react. If :math:`N_0` is
the initial number of monomer molecules (equivalently, functional
groups) and :math:`p` is the fraction of functional groups that have
reacted (the **extent of reaction**), then after reaction there are
:math:`N=N_0(1-p)` molecules left (every reaction event fuses two
molecules into one, so :math:`N_0-N` reaction events have occurred, each
consuming one pair of functional groups, i.e. :math:`N_0p/2` events --
consistent with :math:`N=N_0-N_0p/2\cdot2/1`... more directly: each
reacted functional group belongs to a molecule that merged with another,
so the molecule count drops by exactly one per reacted group-*pair*,
giving :math:`N=N_0(1-p)` directly). The number-average degree of
polymerization is then the **Carothers equation**:

.. math::

    \bar{X}_n = \frac{N_0}{N} = \frac{1}{1-p}

which diverges as :math:`p\to1` -- reaching a useful high molecular
weight requires driving step-growth reactions to very high (often
>99%) conversion, a key qualitative difference from chain-growth
polymerization (:mod:`chemistrykit.polymer.systems.chain_growth`), where
high molecular weight chains form immediately even at low monomer
conversion.
"""

from __future__ import annotations

__all__ = [
    "degree_of_polymerization",
    "extent_of_reaction_for_DP",
    "degree_of_polymerization_stoichiometric_imbalance",
]


def degree_of_polymerization(p: float) -> float:
    r"""The Carothers equation: :math:`\bar{X}_n=1/(1-p)`.

    Parameters
    ----------
    p : float
        Extent of reaction (fraction of functional groups reacted),
        :math:`0\le p<1`.

    Returns
    -------
    float

    Examples
    --------
    At 90% conversion the average chain is 10 repeat units long:

    >>> round(degree_of_polymerization(0.9), 6)
    10.0

    A "high" molecular weight polymer genuinely requires extreme
    conversion -- reaching :math:`\bar X_n=100` needs 99% conversion, not
    90%:

    >>> round(degree_of_polymerization(0.99), 1)
    100.0
    """
    return 1.0 / (1.0 - p)


def extent_of_reaction_for_DP(Xn: float) -> float:
    r"""Invert the Carothers equation: the extent of reaction needed for a target :math:`\bar X_n`.

    Parameters
    ----------
    Xn : float
        Target number-average degree of polymerization, :math:`\ge1`.

    Returns
    -------
    float
        Required extent of reaction `p`.

    Examples
    --------
    Exactly inverts :func:`degree_of_polymerization`:

    >>> round(extent_of_reaction_for_DP(degree_of_polymerization(0.95)), 10)
    0.95
    """
    return 1.0 - 1.0 / Xn


def degree_of_polymerization_stoichiometric_imbalance(p: float, r: float) -> float:
    r"""Carothers equation generalized to a stoichiometric imbalance of the two functional groups.

    When the two reactive functional groups (e.g. -OH and -COOH) are not
    present in exactly equal number -- either a deliberate stoichiometric
    imbalance, or the deliberate addition of a monofunctional "chain
    stopper" -- the maximum attainable degree of polymerization is capped
    even at :math:`p=1` (Odian, *Principles of Polymerization*, 4th ed.,
    Ch. 2.2, eq. 2-70):

    .. math::

        \bar{X}_n = \frac{1+r}{1+r-2rp}

    where :math:`r\le1` is the ratio of the minority to majority
    functional-group count (or, with a monofunctional chain stopper of
    concentration :math:`N_B'`, :math:`r=N_A/(N_B+2N_B')`; see Odian for
    the full derivation). Setting :math:`r=1` (perfect stoichiometric
    balance) recovers the ordinary :func:`degree_of_polymerization`
    exactly.

    Parameters
    ----------
    p : float
        Extent of reaction of the limiting (minority) functional group,
        :math:`0\le p\le1`.
    r : float
        Stoichiometric ratio, :math:`0<r\le1`.

    Returns
    -------
    float

    Examples
    --------
    With ``r=1`` this reduces exactly to the ordinary Carothers equation:

    >>> round(degree_of_polymerization_stoichiometric_imbalance(p=0.9, r=1.0), 10)
    10.0

    A stoichiometric imbalance caps the attainable degree of
    polymerization even at complete conversion of the limiting group
    (:math:`p=1`):

    >>> round(degree_of_polymerization_stoichiometric_imbalance(p=1.0, r=0.98), 4)
    99.0
    """
    return (1.0 + r) / (1.0 + r - 2.0 * r * p)
