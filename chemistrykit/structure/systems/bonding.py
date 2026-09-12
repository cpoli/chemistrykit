r"""Bond order: the Pauling bond-order/bond-length correlation, and Huckel-theory pi bond order.

Two independent routes to a bond order are implemented, matching two
independent chemistry traditions:

* :func:`bond_order_from_length` / :func:`bond_length_from_order` --
  Pauling's empirical correlation between an experimentally measured
  bond length and the (otherwise not directly observable) bond order,
  calibrated against known-integer-bond-order reference bonds (L.
  Pauling, *J. Am. Chem. Soc.* 69, 542 (1947); *The Nature of the
  Chemical Bond*, 3rd ed. (1960), Ch. 7). **Flagged as an approximation**:
  the correlation constant `c` is empirical and bond-type-specific (this
  module's default, 0.71 Å, is Pauling's own value for carbon-carbon
  bonds specifically); a different bond type (C-N, C-O, N-N, ...) needs
  its own fitted `c`, which this module does not tabulate.
* :func:`coulson_pi_bond_order` -- the Coulson bond order from Huckel
  molecular-orbital theory (C. A. Coulson, *Proc. R. Soc. Lond. A* 169,
  413 (1939)), computed *from the MO coefficients* of a
  :class:`chemistrykit.quantum.systems.huckel.HuckelSystem` solve, not
  from bond length -- reusing chemistrykit.quantum's Huckel machinery
  rather than re-deriving it, per the spec's explicit guidance.
"""

from __future__ import annotations

import numpy as np

__all__ = ["bond_order_from_length", "bond_length_from_order", "coulson_pi_bond_order"]

#: Pauling's original carbon-carbon bond-order/length correlation
#: constant, in angstrom (L. Pauling, *J. Am. Chem. Soc.* 69, 542 (1947)).
#: The default for `c` in :func:`bond_order_from_length` /
#: :func:`bond_length_from_order`; pass a different value for other bond
#: types (this module does not tabulate them).
PAULING_C_C_CONSTANT = 0.71


def bond_order_from_length(single_bond_length: float, observed_length: float, c: float = PAULING_C_C_CONSTANT) -> float:
    r"""Estimate a bond order from an observed bond length, via Pauling's empirical correlation.

    .. math::

        D(n) = D(1) - c\ln n \quad\Longrightarrow\quad n = \exp\!\left(\frac{D(1)-D(n)}{c}\right)

    where :math:`D(1)` is the reference single-bond length and `c` is an
    empirical, bond-type-specific constant (L. Pauling, *J. Am. Chem.
    Soc.* 69, 542 (1947)). Shorter bonds correspond to a higher (order
    > 1) bond order; this is an empirical correlation, not a first-
    principles bonding calculation -- see the module docstring.

    Parameters
    ----------
    single_bond_length : float
        Reference single-bond (`n=1`) length :math:`D(1)`, in any
        consistent length unit (e.g. angstrom).
    observed_length : float
        The bond length whose order is being estimated, :math:`D(n)`,
        same units as `single_bond_length`.
    c : float, default 0.71 (angstrom, Pauling's C-C value)
        The empirical correlation constant; must be re-fit for bond types
        other than carbon-carbon.

    Returns
    -------
    float
        Estimated (generally non-integer) bond order.

    Examples
    --------
    Benzene's C-C bond (1.397 Å) is intermediate between a single
    (1.54 Å) and double (1.34 Å) bond, giving an estimated order between
    1 and 2 (the classic evidence for delocalized aromatic bonding):

    >>> n = bond_order_from_length(single_bond_length=1.54, observed_length=1.397)
    >>> bool(1.0 < n < 2.0)
    True

    A bond exactly at the reference single-bond length has order 1:

    >>> round(bond_order_from_length(1.54, 1.54), 6)
    1.0
    """
    return float(np.exp((single_bond_length - observed_length) / c))


def bond_length_from_order(single_bond_length: float, bond_order: float, c: float = PAULING_C_C_CONSTANT) -> float:
    r"""Invert :func:`bond_order_from_length`: predict a bond length from a bond order.

    .. math::

        D(n) = D(1) - c\ln n

    Parameters
    ----------
    single_bond_length : float
        Reference single-bond length :math:`D(1)`.
    bond_order : float
        Bond order `n` (> 0).
    c : float, default 0.71 (angstrom, Pauling's C-C value)

    Returns
    -------
    float
        Predicted bond length, same units as `single_bond_length`.

    Raises
    ------
    ValueError
        If `bond_order` is not positive.

    Examples
    --------
    A round trip through :func:`bond_order_from_length` recovers the
    original length exactly (the two functions are exact inverses):

    >>> D1, Dn = 1.54, 1.34
    >>> n = bond_order_from_length(D1, Dn)
    >>> round(bond_length_from_order(D1, n), 9) == round(Dn, 9)
    True

    Higher bond order predicts a shorter bond -- triple shorter than
    double shorter than single:

    >>> lengths = [bond_length_from_order(1.54, n) for n in (1, 2, 3)]
    >>> bool(lengths[0] > lengths[1] > lengths[2])
    True
    """
    if bond_order <= 0:
        raise ValueError("bond_order must be positive")
    return float(single_bond_length - c * np.log(bond_order))


def coulson_pi_bond_order(coefficients: np.ndarray, occupations, i: int, j: int) -> float:
    r"""The Coulson pi bond order between atoms `i` and `j` from Huckel molecular-orbital coefficients.

    .. math::

        p_{ij} = \sum_k n_k c_{ik} c_{jk}

    summed over molecular orbitals `k` with occupation number
    :math:`n_k\in\{0,1,2\}` (C. A. Coulson, *Proc. R. Soc. Lond. A* 169,
    413 (1939); Streitwieser, *Molecular Orbital Theory for Organic
    Chemists*, Ch. 2). This measures the pi-bonding contribution only
    (Huckel theory does not model the sigma framework at all, per
    :mod:`chemistrykit.quantum.systems.huckel`'s module docstring), so a
    formally "single" sigma-bonded pair with a fractional
    :math:`p_{ij}` (e.g. benzene's 2/3) has a *total* bond order of
    :math:`1+p_{ij}`.

    Parameters
    ----------
    coefficients : ndarray, shape (n_atoms, n_mo)
        MO coefficient matrix, columns are individual MOs -- e.g.
        :attr:`chemistrykit.quantum.core.base_system.EigenstateResult.coefficients`
        from a :meth:`chemistrykit.quantum.systems.huckel.HuckelSystem.solve` call.
    occupations : array-like of float, length n_mo
        Occupation number of each MO (0, 1, or 2 electrons), same column
        order as `coefficients`.
    i, j : int
        0-indexed atom (basis-function) indices.

    Returns
    -------
    float

    Examples
    --------
    Ethene: the single pi bond is fully formed, :math:`p_{01}=1`:

    >>> from chemistrykit.quantum.systems.huckel import HuckelSystem
    >>> ethene = HuckelSystem(n_atoms=2, bonds=[(0, 1)])
    >>> result = ethene.solve()
    >>> order = np.argsort(result.energies)
    >>> occ = np.zeros(2)
    >>> occ[order[0]] = 2.0
    >>> round(coulson_pi_bond_order(result.coefficients, occ, 0, 1), 6)
    1.0

    Benzene: the textbook Coulson bond order of 2/3 for every adjacent
    carbon pair (delocalization spreads the pi bonding evenly around the
    ring, weaker than a localized double bond):

    >>> benzene = HuckelSystem.cyclic_polyene(n_atoms=6)
    >>> result = benzene.solve()
    >>> order = np.argsort(result.energies)
    >>> occ = np.zeros(6)
    >>> occ[order[:3]] = 2.0
    >>> round(coulson_pi_bond_order(result.coefficients, occ, 0, 1), 4)
    0.6667
    """
    coefficients = np.asarray(coefficients, dtype=np.float64)
    occupations = np.asarray(occupations, dtype=np.float64)
    return float(np.sum(occupations * coefficients[i, :] * coefficients[j, :]))
