r"""Simple Huckel molecular-orbital (HMO) theory for conjugated pi systems.

See Atkins & de Paula, *Physical Chemistry*, 11th ed., Ch. 9.4-9.5, or
Levine, *Quantum Chemistry*, 7th ed., Ch. 16.1-16.3, for Huckel theory,
and Atkins & de Paula, *Physical Chemistry*, 11th ed., Ch. 9.5(b) for
Huckel's :math:`4n+2` aromaticity rule (E. Huckel, *Z. Phys.* 70, 204
(1931)).

Huckel theory treats each sp2 carbon's unhybridized 2p_z atomic orbital as
one basis function and makes three simplifying approximations, all
flagged here explicitly rather than silently baked in:

1. **Zero differential overlap**: :math:`S_{ij}=\delta_{ij}` -- the pi
   atomic orbitals are treated as orthonormal, so the secular equation
   :math:`HC=SCE` collapses to the ordinary eigenvalue problem
   :math:`HC=CE` (:func:`chemistrykit.quantum.utils.secular_equation.solve_secular_equation`
   with ``S=None``).
2. **Nearest-neighbor-only resonance integrals**: :math:`H_{ii}=\alpha`
   (the same Coulomb integral for every carbon 2p_z orbital) and
   :math:`H_{ij}=\beta` only for `i`, `j` directly pi-bonded (`0`
   otherwise) -- no through-space or long-range interaction between
   non-bonded p_z orbitals.
3. **sigma-pi separability**: the sigma-bond framework is not modeled at
   all; only the pi system's own energy is computed.
"""

from __future__ import annotations

import numpy as np

from chemistrykit.quantum.core.base_system import EigenstateResult, VariationalSolver
from chemistrykit.quantum.utils.secular_equation import solve_secular_equation

__all__ = ["HuckelSystem", "linear_polyene_eigenvalues", "cyclic_polyene_eigenvalues", "is_aromatic_by_huckel_rule"]


class HuckelSystem(VariationalSolver):
    r"""A conjugated pi system's Huckel Hamiltonian, built from an explicit bond list.

    Parameters
    ----------
    n_atoms : int
        Number of conjugated (sp2, pi-contributing) atoms.
    bonds : sequence of tuple(int, int)
        Pi bonds, as 0-indexed atom pairs ``(i, j)``.
    alpha : float, default 0.0
        The Coulomb integral :math:`\alpha`, in energy units (every atom
        is assumed identical -- carbon 2p_z -- so a single scalar
        suffices; heteroatom Huckel parameters are not modeled here).
    beta : float, default -1.0
        The (negative, by convention) resonance integral :math:`\beta`
        for a pi bond.
    labels : sequence of str, optional
        Human-readable atom labels (defaults to ``C1``, ``C2``, ...).

    Examples
    --------
    Ethene (a single pi bond, 2 atoms) has the trivial 2x2 Huckel result
    :math:`\alpha\pm\beta`:

    >>> ethene = HuckelSystem(n_atoms=2, bonds=[(0, 1)])
    >>> result = ethene.solve()
    >>> np.round(result.energies, 6)
    array([-1.,  1.])
    """

    def __init__(self, n_atoms: int, bonds, alpha: float = 0.0, beta: float = -1.0, labels=None):
        if n_atoms < 1:
            raise ValueError("n_atoms must be >= 1")
        if beta >= 0:
            raise ValueError("beta must be negative (a bonding resonance integral)")
        self.n_atoms = int(n_atoms)
        self.bonds = [(int(i), int(j)) for i, j in bonds]
        self.alpha = float(alpha)
        self.beta = float(beta)
        self.labels = tuple(labels) if labels is not None else tuple(f"C{i + 1}" for i in range(n_atoms))
        if len(self.labels) != n_atoms:
            raise ValueError("labels must have length n_atoms")

    @classmethod
    def linear_polyene(cls, n_atoms: int, alpha: float = 0.0, beta: float = -1.0) -> HuckelSystem:
        """Build an open (acyclic) conjugated chain of `n_atoms` sp2 carbons (e.g. butadiene, ``n_atoms=4``).

        Parameters
        ----------
        n_atoms : int
        alpha : float, default 0.0
        beta : float, default -1.0

        Returns
        -------
        HuckelSystem
        """
        bonds = [(i, i + 1) for i in range(n_atoms - 1)]
        return cls(n_atoms=n_atoms, bonds=bonds, alpha=alpha, beta=beta)

    @classmethod
    def cyclic_polyene(cls, n_atoms: int, alpha: float = 0.0, beta: float = -1.0) -> HuckelSystem:
        """Build a cyclic conjugated ring of `n_atoms` sp2 carbons (e.g. benzene, ``n_atoms=6``).

        Parameters
        ----------
        n_atoms : int
            Must be >= 3.
        alpha : float, default 0.0
        beta : float, default -1.0

        Returns
        -------
        HuckelSystem
        """
        if n_atoms < 3:
            raise ValueError("a cyclic polyene needs n_atoms >= 3")
        bonds = [(i, (i + 1) % n_atoms) for i in range(n_atoms)]
        return cls(n_atoms=n_atoms, bonds=bonds, alpha=alpha, beta=beta)

    def hamiltonian(self) -> np.ndarray:
        """Build the Huckel Hamiltonian matrix: :math:`\\alpha` on the diagonal, :math:`\\beta` for bonded pairs.

        Returns
        -------
        ndarray, shape (n_atoms, n_atoms)
        """
        H = np.full((self.n_atoms, self.n_atoms), 0.0)
        np.fill_diagonal(H, self.alpha)
        for i, j in self.bonds:
            H[i, j] = self.beta
            H[j, i] = self.beta
        return H

    def solve(self) -> EigenstateResult:
        """Diagonalize the Huckel Hamiltonian (orthonormal basis: ``S=I``).

        Returns
        -------
        EigenstateResult
            ``energies`` ascending (most bonding, i.e. most negative for
            ``beta<0``, first); ``coefficients`` columns are the
            molecular-orbital pi-electron-density expansion coefficients.
        """
        energies, coefficients = solve_secular_equation(self.hamiltonian(), S=None)
        return EigenstateResult(energies=energies, coefficients=coefficients, basis_labels=self.labels, overlap=None)

    def pi_electron_energy(self, n_pi_electrons: int) -> float:
        r"""Total pi-electron energy, filling MOs lowest-energy-first, two electrons per orbital.

        :math:`E_\pi=\sum_kn_k\varepsilon_k`, with occupation numbers
        :math:`n_k\in\{0,1,2\}` (Aufbau + Pauli, the same filling rule
        used throughout this subpackage, e.g.
        :func:`chemistrykit.quantum.systems.particle_in_box.conjugated_dye_absorption_wavelength`).

        Parameters
        ----------
        n_pi_electrons : int
            Total number of pi electrons to place, ``<= 2 * n_atoms``.

        Returns
        -------
        float
            Total pi-electron energy, in the same units as `alpha`/`beta`.
        """
        if not (0 <= n_pi_electrons <= 2 * self.n_atoms):
            raise ValueError("n_pi_electrons out of range")
        energies = np.sort(self.solve().energies)
        n_full = n_pi_electrons // 2
        remainder = n_pi_electrons % 2
        total = 2.0 * np.sum(energies[:n_full])
        if remainder:
            total += energies[n_full]
        return float(total)

    def delocalization_energy(self, n_pi_electrons: int) -> float:
        r"""Huckel delocalization (resonance) energy relative to isolated (localized) pi bonds.

        :math:`E_{deloc}=E_\pi(\text{Huckel})-n_{bonds}\times2(\alpha+\beta)`,
        comparing the delocalized Huckel pi-energy to the energy of the
        same number of pi electrons confined to isolated, non-interacting
        double bonds (each an ethene-like 2-orbital system contributing
        :math:`2(\alpha+\beta)` when filled with 2 electrons) -- the
        standard measure of aromatic/conjugative stabilization (Atkins &
        de Paula, *Physical Chemistry*, 11th ed., Ch. 9.5(a)).

        Parameters
        ----------
        n_pi_electrons : int
            Must be even (whole pi bonds) and equal to the number of pi
            bonds in `bonds` times 2.

        Returns
        -------
        float
            More negative means more delocalization-stabilized.
        """
        if n_pi_electrons % 2 != 0:
            raise ValueError("n_pi_electrons must be even to compare against localized double bonds")
        n_bonds = n_pi_electrons // 2
        localized = n_bonds * 2.0 * (self.alpha + self.beta)
        return self.pi_electron_energy(n_pi_electrons) - localized

    def frontier_electron_density(self, n_pi_electrons: int, orbital: str = "homo", tol: float = 1e-6) -> np.ndarray:
        r"""Fukui's frontier electron density on each atom, :math:`f_r=2c_{r,\mathrm{F}}^2`.

        Fukui, Yonezawa, and Shingu's frontier-orbital reactivity index
        (K. Fukui, T. Yonezawa, H. Shingu, *J. Chem. Phys.* 20, 722
        (1952)): an electrophile attacks the atom with the largest
        HOMO density and a nucleophile the atom with the largest LUMO
        density. The factor 2 makes the densities sum to 2 over all atoms
        (two frontier electrons).

        Parameters
        ----------
        n_pi_electrons : int
            Total pi electrons (must be even, so the HOMO is doubly occupied).
        orbital : {"homo", "lumo"}, default "homo"
            Which frontier orbital to use.
        tol : float, default 1e-6
            A frontier orbital degenerate with its neighbor within `tol`
            is ambiguous and raises ``ValueError``.

        Returns
        -------
        ndarray, shape (n_atoms,)

        Examples
        --------
        Butadiene's HOMO is concentrated on the terminal carbons, the
        sites where electrophiles add:

        >>> density = HuckelSystem.linear_polyene(4).frontier_electron_density(4)
        >>> np.round(density, 3)
        array([0.724, 0.276, 0.276, 0.724])
        """
        if n_pi_electrons <= 0 or n_pi_electrons % 2 != 0 or n_pi_electrons >= 2 * self.n_atoms:
            raise ValueError("n_pi_electrons must be even, positive, and leave at least one empty orbital")
        if orbital not in ("homo", "lumo"):
            raise ValueError("orbital must be 'homo' or 'lumo'")
        result = self.solve()
        order = np.argsort(result.energies)
        energies = result.energies[order]
        index = n_pi_electrons // 2 - 1 if orbital == "homo" else n_pi_electrons // 2
        neighbors = [k for k in (index - 1, index + 1) if 0 <= k < self.n_atoms]
        if any(abs(energies[k] - energies[index]) < tol for k in neighbors):
            raise ValueError("the frontier orbital is degenerate; its density is not uniquely defined")
        c = result.coefficients[:, order[index]]
        return 2.0 * c**2


def linear_polyene_eigenvalues(n_atoms: int, alpha: float = 0.0, beta: float = -1.0) -> np.ndarray:
    r"""Closed-form Huckel eigenvalues of a linear (acyclic) conjugated chain.

    .. math::

        E_k = \alpha+2\beta\cos\!\left(\frac{k\pi}{n+1}\right), \qquad k=1,\dots,n

    (Coulson's formula; C. A. Coulson, *Proc. R. Soc. Lond. A* 169, 413
    (1939); Levine, *Quantum Chemistry*, 7th ed., eq. 16.36.) Provided as
    an independent closed-form cross-check of
    :meth:`HuckelSystem.solve`'s numerical diagonalization for
    :meth:`HuckelSystem.linear_polyene`.

    Parameters
    ----------
    n_atoms : int
    alpha : float, default 0.0
    beta : float, default -1.0

    Returns
    -------
    ndarray, shape (n_atoms,)
        Ascending.

    Examples
    --------
    Butadiene's Huckel eigenvalues (:math:`n=4`):

    >>> np.round(linear_polyene_eigenvalues(4), 4)
    array([-1.618, -0.618,  0.618,  1.618])
    """
    k = np.arange(1, n_atoms + 1)
    return np.sort(alpha + 2.0 * beta * np.cos(k * np.pi / (n_atoms + 1)))


def cyclic_polyene_eigenvalues(n_atoms: int, alpha: float = 0.0, beta: float = -1.0) -> np.ndarray:
    r"""Closed-form Huckel eigenvalues of a cyclic conjugated ring.

    .. math::

        E_k = \alpha+2\beta\cos\!\left(\frac{2\pi k}{n}\right), \qquad k=0,1,\dots,n-1

    (the "Frost circle" mnemonic; A. A. Frost & B. Musulin, *J. Chem.
    Phys.* 21, 572 (1953); Levine, *Quantum Chemistry*, 7th ed., eq.
    16.62.) Provided as an independent closed-form cross-check of
    :meth:`HuckelSystem.solve`'s numerical diagonalization for
    :meth:`HuckelSystem.cyclic_polyene`.

    Parameters
    ----------
    n_atoms : int
    alpha : float, default 0.0
    beta : float, default -1.0

    Returns
    -------
    ndarray, shape (n_atoms,)
        Ascending.

    Examples
    --------
    Benzene's famous Huckel pattern -- :math:`\alpha+2\beta`,
    :math:`\alpha+\beta` (doubly degenerate), :math:`\alpha-\beta`
    (doubly degenerate), :math:`\alpha-2\beta`:

    >>> np.round(cyclic_polyene_eigenvalues(6), 4)
    array([-2., -1., -1.,  1.,  1.,  2.])
    """
    k = np.arange(0, n_atoms)
    return np.sort(alpha + 2.0 * beta * np.cos(2.0 * np.pi * k / n_atoms))


def is_aromatic_by_huckel_rule(n_pi_electrons: int, energies, tol: float = 1e-6) -> bool:
    r"""Check Huckel's :math:`4n+2` rule against an *actual computed* Huckel spectrum.

    A cyclic, fully conjugated, planar system is predicted aromatic when
    (Atkins & de Paula, *Physical Chemistry*, 11th ed., Ch. 9.5(b); E.
    Huckel, *Z. Phys.* 70, 204 (1931)):

    1. it has :math:`4n+2` pi electrons for some non-negative integer `n`, *and*
    2. that electron count exactly fills a set of Huckel MOs with no
       unpaired electron in a partially-filled degenerate level (the
       closed-shell condition) -- checked here directly against the
       degeneracies of the *computed* `energies`, not assumed.

    Condition 2 is what makes this a genuine check against the
    diagonalized spectrum rather than a bare electron-counting rule:
    a system could satisfy :math:`N=4n+2` yet still land mid-degenerate-level
    for an unusual (non-uniform-ring) Huckel Hamiltonian, which this
    function would correctly flag as not closed-shell.

    Parameters
    ----------
    n_pi_electrons : int
        Number of pi electrons.
    energies : array-like of float
        The full computed Huckel spectrum (e.g. from
        :meth:`HuckelSystem.solve`'s ``.energies``), any order.
    tol : float, default 1e-6
        Absolute energy tolerance for treating two levels as degenerate.

    Returns
    -------
    bool

    Examples
    --------
    Benzene: 6 pi electrons, :math:`4(1)+2`, and the computed spectrum's
    filling is closed-shell (both the doubly degenerate HOMO orbitals are
    fully occupied):

    >>> benzene_energies = cyclic_polyene_eigenvalues(6)
    >>> is_aromatic_by_huckel_rule(6, benzene_energies)
    True

    Cyclobutadiene: 4 pi electrons is not :math:`4n+2` for any integer
    `n` (it is :math:`4n`, the antiaromatic count), so this correctly
    returns False regardless of the degeneracy structure:

    >>> cbd_energies = cyclic_polyene_eigenvalues(4)
    >>> is_aromatic_by_huckel_rule(4, cbd_energies)
    False
    """
    if n_pi_electrons <= 0 or n_pi_electrons % 2 != 0:
        return False
    n = (n_pi_electrons - 2) / 4.0
    if abs(n - round(n)) > 1e-9 or n < 0:
        return False

    energies = np.sort(np.asarray(energies, dtype=np.float64))
    n_orbitals_needed = n_pi_electrons // 2
    if n_orbitals_needed > len(energies):
        return False
    if n_orbitals_needed == len(energies):
        return True
    homo_energy = energies[n_orbitals_needed - 1]
    lumo_energy = energies[n_orbitals_needed]
    # Closed shell iff the last occupied orbital is NOT degenerate with
    # the first unoccupied one (otherwise electrons would have to split
    # unpaired across a degenerate HOMO/LUMO set -- Hund's rule territory,
    # not a closed shell).
    return bool(abs(homo_energy - lumo_energy) > tol)
