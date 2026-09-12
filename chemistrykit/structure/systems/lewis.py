r"""Formal charge and oxidation-state assignment from a Lewis structure.

See Atkins & de Paula, *Physical Chemistry*, 11th ed., Ch. 9.1(a) and
Ch. 2.9, or Housecroft & Sharpe, *Inorganic Chemistry*, 5th ed., Ch. 2,
throughout. Both quantities are computed *algorithmically* from a
specified Lewis structure (connectivity + bond orders + lone-pair
counts) and, for oxidation state, the atoms' electronegativities from
:mod:`chemistrykit.periodic_table` -- never looked up by molecular
formula.

**Formal charge** answers "how would this atom's electron count compare
to the free atom if every bonding pair were split evenly, regardless of
which atom is more electronegative?" -- a bookkeeping device for
comparing alternative Lewis structures of the *same* molecule (G. N.
Lewis, *J. Am. Chem. Soc.* 38, 762 (1916)):

.. math::

    FC = V - N - \frac{B}{2}

where `V` is the free atom's valence-electron count, `N` its nonbonding
(lone-pair) electrons in this structure, and `B` its total bonding
electrons (twice the sum of bond orders at that atom).

**Oxidation state** instead answers "how would this atom's electron
count compare to the free atom in the fully ionic limit, where every
bonding pair is assigned entirely to the more electronegative partner?"
(a homonuclear bond, no electronegativity difference, splits its
electrons evenly, contributing no net shift to either atom):

.. math::

    OS = V - \left(2N + \sum_{\text{bonds at this atom}} w \cdot(\text{bond order}\times2)\right)

where `w` is `1` if this atom is the more electronegative partner of
that bond, `0` if the less electronegative, and `1/2` for a homonuclear
bond.
"""

from __future__ import annotations

from dataclasses import dataclass, field

from chemistrykit.periodic_table import electronegativity, valence_electrons

__all__ = ["LewisStructure"]


@dataclass
class LewisStructure:
    """A Lewis structure: element symbols, bond orders, and lone-pair counts.

    Parameters
    ----------
    symbols : sequence of str
        Element symbol of each atom (main-group elements only -- see
        :func:`chemistrykit.periodic_table.valence_electrons`), 0-indexed.
    bond_orders : dict of tuple(int, int) to int
        Bond order for each bonded atom pair, e.g. ``{(0, 1): 1, (0, 2): 2}``
        for one single and one double bond from atom 0. Each pair should
        appear once (order doesn't matter, `(i, j)` and `(j, i)` are
        treated identically).
    lone_pairs : dict of int to int, optional
        Number of lone *pairs* (not electrons) on each atom, defaulting
        to 0 for any atom not given explicitly.

    Examples
    --------
    Water: O has 2 single bonds (to atoms 1, 2) and 2 lone pairs.

    >>> water = LewisStructure(symbols=["O", "H", "H"], bond_orders={(0, 1): 1, (0, 2): 1}, lone_pairs={0: 2})
    >>> water.formal_charges()
    {0: 0.0, 1: 0.0, 2: 0.0}
    >>> water.oxidation_states()
    {0: -2.0, 1: 1.0, 2: 1.0}
    """

    symbols: list
    bond_orders: dict
    lone_pairs: dict = field(default_factory=dict)

    def __post_init__(self):
        self.symbols = list(self.symbols)
        self.bond_orders = {tuple(sorted(pair)): order for pair, order in self.bond_orders.items()}
        n = len(self.symbols)
        for i, j in self.bond_orders:
            if not (0 <= i < n and 0 <= j < n):
                raise ValueError(f"bond ({i}, {j}) references an atom index out of range for {n} atoms")
        for i in self.lone_pairs:
            if not (0 <= i < n):
                raise ValueError(f"lone_pairs references atom index {i} out of range for {n} atoms")

    def _bonds_at(self, atom: int):
        for (i, j), order in self.bond_orders.items():
            if i == atom:
                yield j, order
            elif j == atom:
                yield i, order

    def formal_charges(self) -> dict:
        r"""Compute the formal charge :math:`FC=V-N-B/2` for every atom.

        Returns
        -------
        dict of int to float

        Examples
        --------
        Ammonium, :math:`NH_4^+`: nitrogen has 4 single bonds and no lone
        pairs, giving the textbook +1 formal charge:

        >>> nh4_plus = LewisStructure(symbols=["N", "H", "H", "H", "H"], bond_orders={(0, k): 1 for k in (1, 2, 3, 4)})
        >>> nh4_plus.formal_charges()[0]
        1.0
        """
        charges = {}
        for i, symbol in enumerate(self.symbols):
            V = valence_electrons(symbol)
            N = 2 * self.lone_pairs.get(i, 0)
            B = sum(order for _, order in self._bonds_at(i))
            charges[i] = float(V - N - B)
        return charges

    def total_formal_charge(self) -> float:
        """Sum of all atoms' formal charges -- must equal the molecule's net charge.

        A useful self-consistency check on a proposed Lewis structure: any
        valid structure's formal charges must sum to the actual net
        charge of the species (0 for a neutral molecule).

        Returns
        -------
        float
        """
        return sum(self.formal_charges().values())

    def oxidation_states(self) -> dict:
        r"""Compute the oxidation state of every atom via the electronegativity-based algorithm.

        For each bond, both bonding electrons (per unit of bond order)
        are assigned entirely to the more electronegative atom (split
        evenly for a homonuclear bond); nonbonding electrons always stay
        with their own atom. See the module docstring for the full
        formula.

        Returns
        -------
        dict of int to float

        Raises
        ------
        KeyError
            If any atom's element has no tabulated Pauling
            electronegativity (see
            :func:`chemistrykit.periodic_table.electronegativity`).

        Examples
        --------
        Hydrogen peroxide, :math:`H_2O_2` (H-O-O-H): the O-O bond is
        homonuclear (no net electron shift), each O keeps its lone pairs
        and wins its one O-H bond, giving the well-known -1 oxidation
        state (intermediate between water's -2 and O2's 0):

        >>> h2o2 = LewisStructure(symbols=["H", "O", "O", "H"], bond_orders={(0, 1): 1, (1, 2): 1, (2, 3): 1}, lone_pairs={1: 2, 2: 2})
        >>> h2o2.oxidation_states()
        {0: 1.0, 1: -1.0, 2: -1.0, 3: 1.0}
        """
        owned = {i: 2.0 * self.lone_pairs.get(i, 0) for i in range(len(self.symbols))}
        for (i, j), order in self.bond_orders.items():
            electrons = 2.0 * order
            en_i = electronegativity(self.symbols[i])
            en_j = electronegativity(self.symbols[j])
            if en_i > en_j:
                owned[i] += electrons
            elif en_j > en_i:
                owned[j] += electrons
            else:
                owned[i] += electrons / 2.0
                owned[j] += electrons / 2.0
        return {i: float(valence_electrons(symbol) - owned[i]) for i, symbol in enumerate(self.symbols)}
