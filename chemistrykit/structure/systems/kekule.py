r"""Kekulé structures: enumerating the alternating single/double-bond patterns of a conjugated molecule.

Kekulé proposed in 1865 that benzene is a ring of six carbon atoms joined
by alternating single and double bonds (A. Kekulé, *Bull. Soc. Chim.
Paris* 3, 98 (1865)). In graph terms, a Kekulé structure of a conjugated
hydrocarbon's carbon skeleton is a *perfect matching*: a set of double
bonds such that every carbon takes part in exactly one of them. Benzene
has two (Kekulé's own "oscillating" pair), naphthalene three, anthracene
four and phenanthrene five (see Gutman & Cyvin, *Introduction to the
Theory of Benzenoid Hydrocarbons* (Springer, 1989), Ch. 2). The count,
the Kekulé structure count :math:`K`, is a classic index of aromatic
stability.

The enumeration here is a plain depth-first search: pick the lowest-
numbered atom that has no double bond yet, try each unmatched neighbor as
its double-bond partner, and recurse.
"""

from __future__ import annotations

__all__ = ["kekule_structures", "count_kekule_structures"]


def kekule_structures(n_atoms: int, bonds) -> list:
    """Enumerate every Kekulé structure (perfect matching) of a conjugated skeleton.

    Parameters
    ----------
    n_atoms : int
        Number of conjugated atoms (e.g. the carbons of a benzenoid
        hydrocarbon), labelled ``0 .. n_atoms - 1``.
    bonds : sequence of tuple(int, int)
        The sigma-bond connectivity of the skeleton.

    Returns
    -------
    list of list of tuple(int, int)
        Each entry is one Kekulé structure: the sorted list of bonds that
        are double bonds in it. Every atom appears in exactly one of them.
        An odd-membered or otherwise unmatchable skeleton returns ``[]``.

    Examples
    --------
    Benzene has exactly Kekulé's two structures:

    >>> ring = [(i, (i + 1) % 6) for i in range(6)]
    >>> for s in kekule_structures(6, ring):
    ...     print(s)
    [(0, 1), (2, 3), (4, 5)]
    [(0, 5), (1, 2), (3, 4)]
    """
    neighbors: dict = {i: set() for i in range(n_atoms)}
    for i, j in bonds:
        if i == j or not (0 <= i < n_atoms and 0 <= j < n_atoms):
            raise ValueError(f"invalid bond ({i}, {j}) for {n_atoms} atoms")
        neighbors[i].add(j)
        neighbors[j].add(i)
    if n_atoms % 2:
        return []

    structures: list = []
    matched = [False] * n_atoms

    def search(current: list) -> None:
        try:
            first = matched.index(False)
        except ValueError:
            structures.append(sorted(current))
            return
        matched[first] = True
        for partner in sorted(neighbors[first]):
            if not matched[partner]:
                matched[partner] = True
                current.append((min(first, partner), max(first, partner)))
                search(current)
                current.pop()
                matched[partner] = False
        matched[first] = False

    search([])
    return structures


def count_kekule_structures(n_atoms: int, bonds) -> int:
    """Number of Kekulé structures :math:`K` of a conjugated skeleton.

    Parameters
    ----------
    n_atoms : int
    bonds : sequence of tuple(int, int)
        See :func:`kekule_structures`.

    Returns
    -------
    int

    Examples
    --------
    Naphthalene (two fused six-membered rings, 10 carbons) has three:

    >>> naphthalene = [(0, 1), (1, 2), (2, 3), (3, 4), (4, 9), (9, 0),
    ...                (4, 5), (5, 6), (6, 7), (7, 8), (8, 9)]
    >>> count_kekule_structures(10, naphthalene)
    3
    """
    return len(kekule_structures(n_atoms, bonds))
