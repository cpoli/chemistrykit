r"""A curated table of standard reduction potentials, and cell-potential combination.

See Bard & Faulkner, *Electrochemical Methods: Fundamentals and
Applications*, 2nd ed., Appendix C.3, or Atkins & de Paula, *Physical
Chemistry*, 11th ed., Table 6.4, for the tabulated values below (all
referenced to the standard hydrogen electrode, SHE, defined as exactly 0
V at all temperatures) and Ch. 6.9-6.10 (Atkins & de Paula) for how two
half-reactions combine into a cell potential.

This is a small, curated reference table in the same spirit as
:mod:`chemistrykit.periodic_table` -- plain data, not a comprehensive
electrochemical-series database. Values are standard (298.15 K, unit
activity) reduction potentials in volts.
"""

from __future__ import annotations

import math
from dataclasses import dataclass

__all__ = [
    "HalfReaction",
    "STANDARD_REDUCTION_POTENTIALS",
    "cell_potential",
    "standard_cell_potential",
    "balance_redox_reaction",
    "is_spontaneous",
]


@dataclass(frozen=True)
class HalfReaction:
    """A single tabulated standard reduction half-reaction.

    Parameters
    ----------
    name : str
        Short label, e.g. ``"Cu2+/Cu"`` for :math:`Cu^{2+} + 2e^- \to Cu`.
    n : int
        Number of electrons transferred in the half-reaction as written
        (in the reduction direction).
    E_standard : float
        Standard reduction potential, in V vs. SHE.
    """

    name: str
    n: int
    E_standard: float


#: Standard reduction potentials at 298.15 K vs. the standard hydrogen
#: electrode (SHE = 0 V exactly), in V, ordered from most negative
#: (strongest reducing agent in its reduced form) to most positive
#: (strongest oxidizing agent in its oxidized form). Bard & Faulkner,
#: *Electrochemical Methods*, 2nd ed., Appendix C.3; Atkins & de Paula,
#: *Physical Chemistry*, 11th ed., Table 6.4.
STANDARD_REDUCTION_POTENTIALS: dict[str, HalfReaction] = {
    entry.name: entry
    for entry in (
        HalfReaction("Li+/Li", 1, -3.04),
        HalfReaction("K+/K", 1, -2.93),
        HalfReaction("Ca2+/Ca", 2, -2.87),
        HalfReaction("Na+/Na", 1, -2.71),
        HalfReaction("Mg2+/Mg", 2, -2.36),
        HalfReaction("Al3+/Al", 3, -1.66),
        HalfReaction("Zn2+/Zn", 2, -0.76),
        HalfReaction("Fe2+/Fe", 2, -0.44),
        HalfReaction("Ni2+/Ni", 2, -0.26),
        HalfReaction("Sn2+/Sn", 2, -0.14),
        HalfReaction("Pb2+/Pb", 2, -0.13),
        HalfReaction("H+/H2", 2, 0.00),
        HalfReaction("Sn4+/Sn2+", 2, 0.15),
        HalfReaction("Cu2+/Cu", 2, 0.34),
        HalfReaction("I2/I-", 2, 0.54),
        HalfReaction("Fe3+/Fe2+", 1, 0.77),
        HalfReaction("Ag+/Ag", 1, 0.80),
        HalfReaction("Br2/Br-", 2, 1.09),
        HalfReaction("O2/H2O", 4, 1.23),
        HalfReaction("Cr2O7^2-/Cr3+", 6, 1.33),
        HalfReaction("Cl2/Cl-", 2, 1.36),
        HalfReaction("MnO4-/Mn2+", 5, 1.51),
        HalfReaction("F2/F-", 2, 2.87),
    )
}


def cell_potential(cathode: HalfReaction, anode: HalfReaction) -> float:
    r"""Standard cell potential from two half-reactions: :math:`E^\circ_{cell} = E^\circ_{cathode} - E^\circ_{anode}`.

    Both potentials are used exactly as tabulated *reduction* potentials
    -- the anode's is not negated by hand and neither is scaled by the
    number of electrons it transfers, because electrode potential is an
    intensive quantity (a per-electron driving force), unaffected by how
    many electrons the balanced overall reaction happens to require
    (Atkins & de Paula, *Physical Chemistry*, 11th ed., Ch. 6.9). See
    :func:`balance_redox_reaction` for the electron-balancing arithmetic
    this potential calculation deliberately does *not* need.

    Parameters
    ----------
    cathode : HalfReaction
        The half-reaction that runs as a reduction (gains electrons).
    anode : HalfReaction
        The half-reaction that runs as an oxidation (loses electrons);
        its tabulated *reduction* potential is passed in unchanged.

    Returns
    -------
    float
        Standard cell potential, in V. Positive means the reaction as
        paired (cathode reduction + anode oxidation) is spontaneous.

    Examples
    --------
    The classic Daniell cell, Cu cathode / Zn anode:

    >>> from chemistrykit.electrochem.systems.standard_potentials import STANDARD_REDUCTION_POTENTIALS as T
    >>> round(cell_potential(T["Cu2+/Cu"], T["Zn2+/Zn"]), 2)
    1.1
    """
    return cathode.E_standard - anode.E_standard


def standard_cell_potential(cathode_name: str, anode_name: str) -> float:
    """Look up two half-reactions by name in :data:`STANDARD_REDUCTION_POTENTIALS` and combine them.

    Parameters
    ----------
    cathode_name, anode_name : str
        Keys into :data:`STANDARD_REDUCTION_POTENTIALS`.

    Returns
    -------
    float
        Standard cell potential, in V.

    Examples
    --------
    >>> round(standard_cell_potential("Cu2+/Cu", "Zn2+/Zn"), 2)
    1.1
    """
    return cell_potential(STANDARD_REDUCTION_POTENTIALS[cathode_name], STANDARD_REDUCTION_POTENTIALS[anode_name])


def balance_redox_reaction(cathode: HalfReaction, anode: HalfReaction) -> tuple[int, int, int]:
    r"""Electron-balancing multiples for combining two half-reactions into one overall redox reaction.

    Half-reaction potentials combine as-is (see :func:`cell_potential`),
    but *mass-balancing* the overall reaction -- and therefore computing
    e.g. how much of each species is consumed per mole of overall
    reaction via Faraday's laws
    (:mod:`chemistrykit.electrochem.systems.electrolysis`) -- requires
    each half-reaction to transfer the same number of electrons. The
    standard recipe (Atkins & de Paula, *Physical Chemistry*, 11th ed.,
    Ch. 6.9) multiplies the cathode half-reaction by
    :math:`n_{anode}/\gcd(n_{cathode},n_{anode})` and the anode
    half-reaction by :math:`n_{cathode}/\gcd(n_{cathode},n_{anode})`, so
    both sides transfer :math:`\mathrm{lcm}(n_{cathode},n_{anode})`
    electrons.

    Parameters
    ----------
    cathode, anode : HalfReaction

    Returns
    -------
    cathode_multiple, anode_multiple, n_total : int
        Multiples applied to the cathode and anode half-reactions, and
        the resulting total number of electrons transferred by the
        balanced overall reaction.

    Examples
    --------
    :math:`Zn + Cu^{2+} \to Zn^{2+} + Cu` needs no rebalancing (both
    transfer 2 electrons already):

    >>> from chemistrykit.electrochem.systems.standard_potentials import STANDARD_REDUCTION_POTENTIALS as T
    >>> balance_redox_reaction(T["Cu2+/Cu"], T["Zn2+/Zn"])
    (1, 1, 2)

    :math:`MnO_4^- + 5Fe^{2+} \to \dots` needs Fe balanced 5-fold against
    Mn's single 5-electron step:

    >>> balance_redox_reaction(T["MnO4-/Mn2+"], T["Fe3+/Fe2+"])
    (1, 5, 5)
    """
    g = math.gcd(cathode.n, anode.n)
    cathode_multiple = anode.n // g
    anode_multiple = cathode.n // g
    n_total = cathode.n * cathode_multiple
    return cathode_multiple, anode_multiple, n_total


def is_spontaneous(E_cell: float) -> bool:
    r"""Whether a cell reaction is spontaneous as written, from its cell potential.

    :math:`E_{cell} > 0 \iff \Delta G = -nFE_{cell} < 0`, the spontaneity
    criterion (Atkins & de Paula, *Physical Chemistry*, 11th ed., Ch.
    6.9). A spontaneous cell reaction discharges *galvanically*
    (delivers electrical work); a non-spontaneous one must be driven
    *electrolytically* by an externally applied voltage exceeding
    :math:`|E_{cell}|` in the opposing sense -- see
    :mod:`chemistrykit.electrochem.systems.electrolysis`.

    Parameters
    ----------
    E_cell : float
        Cell potential (standard or Nernst-corrected), in V.

    Returns
    -------
    bool

    Examples
    --------
    >>> is_spontaneous(1.10)
    True
    >>> is_spontaneous(-0.76)
    False
    """
    return E_cell > 0.0
