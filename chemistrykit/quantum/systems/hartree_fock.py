r"""Restricted Hartree-Fock SCF for closed-shell molecules, and a variational H2+ warm-up.

:class:`RestrictedHartreeFock` is the self-consistent-field (SCF)
Roothaan-Hall procedure of Szabo & Ostlund, *Modern Quantum Chemistry*,
1st ed. rev., Ch. 3.4-3.5, for closed-shell molecules in a basis of
contracted s-type Gaussians (e.g. STO-3G 1s on H2 and HeH+,
:meth:`RestrictedHartreeFock.h2` and :meth:`RestrictedHartreeFock.heh_plus`).
:class:`H2PlusVariational`, below, is the one-electron special case with
no electron-electron repulsion and hence no SCF loop:

The simplest possible molecule with a real two-center chemical bond,
H2+ (one electron, two protons), is exactly solvable in confocal
elliptic coordinates -- but the point of this module is the *variational
method itself*, the way real molecular electronic-structure calculations
work: pick a finite basis, build the Hamiltonian and overlap matrices,
solve the secular equation :math:`HC=SCE`
(:func:`chemistrykit.quantum.utils.secular_equation.solve_secular_equation`),
and then *improve* the result by varying a real basis parameter, exactly
as the Rayleigh-Ritz variational theorem promises: any trial wavefunction
gives an energy that is an upper bound on the true ground-state energy,
so lowering the computed energy by varying a parameter can only move it
closer to the truth (Szabo & Ostlund, *Modern Quantum Chemistry*, 1st ed.
rev., Ch. 1.3 and Ch. 3.4).

Basis: one un-contracted, normalized s-type Gaussian primitive
(:class:`chemistrykit.quantum.utils.basis_sets.GaussianPrimitive`) on
each proton, sharing a single variational orbital exponent
:math:`\alpha` -- a "GTO-1G" minimal basis (a single Gaussian standing in
for the true 1s Slater-type orbital each proton would carry in a real
LCAO-MO treatment; Levine, *Quantum Chemistry*, 7th ed., Ch. 13.1-13.2
covers the (Slater-orbital) LCAO-MO treatment of H2+ this mirrors). The
molecular orbitals are then

.. math::

    \psi_\pm = c_A\chi_A\pm c_B\chi_B

the bonding (+) and antibonding (-) combinations, obtained directly as
the two eigenvectors of the 2x2 secular equation.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from scipy.optimize import minimize_scalar

from chemistrykit.constants import ELEMENTARY_CHARGE, VACUUM_PERMITTIVITY
from chemistrykit.quantum.core.base_system import EigenstateResult, VariationalSolver
from chemistrykit.quantum.utils.basis_sets import (
    BOHR_RADIUS,
    ContractedGaussian,
    GaussianPrimitive,
    contracted_electron_repulsion,
    contracted_one_electron,
    kinetic_integral,
    nuclear_attraction_integral,
    overlap_integral,
    sto3g_1s,
)
from chemistrykit.quantum.utils.secular_equation import solve_secular_equation

__all__ = ["H2PlusVariational", "ExponentOptimizationResult", "RestrictedHartreeFock", "SCFResult"]

#: Coulomb's-law energy scale e^2/(4*pi*epsilon_0), in J*m (proton-proton
#: nuclear repulsion is Z_A*Z_B*this / R, with Z_A=Z_B=1).
_COULOMB_CONSTANT = ELEMENTARY_CHARGE**2 / (4.0 * np.pi * VACUUM_PERMITTIVITY)


class H2PlusVariational(VariationalSolver):
    r"""H2+ (one electron, two protons separated by `bond_length`) in a minimal 2-Gaussian LCAO basis.

    Parameters
    ----------
    bond_length : float
        Proton-proton separation `R`, in m.

    Examples
    --------
    At a plausible bond length, the bonding MO is lower in energy than
    the antibonding one -- the basic LCAO-MO picture of a covalent bond:

    >>> h2plus = H2PlusVariational(bond_length=106e-12)
    >>> result = h2plus.solve(alpha=1.0 / (5.29e-11 ** 2))
    >>> bool(result.energies[0] < result.energies[1])
    True
    """

    def __init__(self, bond_length: float):
        if bond_length <= 0:
            raise ValueError("bond_length must be positive")
        self.bond_length = float(bond_length)
        self._center_a = np.array([0.0, 0.0, -bond_length / 2.0])
        self._center_b = np.array([0.0, 0.0, bond_length / 2.0])

    @property
    def nuclear_repulsion(self) -> float:
        r"""float: The proton-proton repulsion energy, :math:`e^2/(4\pi\varepsilon_0R)`, in J."""
        return _COULOMB_CONSTANT / self.bond_length

    def _basis(self, alpha: float):
        return GaussianPrimitive(alpha, self._center_a), GaussianPrimitive(alpha, self._center_b)

    def hamiltonian(self, alpha: float = 1.0) -> np.ndarray:
        r"""Build the 2x2 core Hamiltonian (kinetic + both nuclear attractions) at exponent `alpha`.

        Parameters
        ----------
        alpha : float, default 1.0
            Shared Gaussian orbital exponent, in m^-2.

        Returns
        -------
        ndarray, shape (2, 2)
            Energy, in J.
        """
        chi_a, chi_b = self._basis(alpha)
        H = np.empty((2, 2))
        for i, gi in enumerate((chi_a, chi_b)):
            for j, gj in enumerate((chi_a, chi_b)):
                kinetic = kinetic_integral(gi, gj)
                attraction = nuclear_attraction_integral(gi, gj, Z=1.0, nucleus_center=self._center_a) + nuclear_attraction_integral(
                    gi, gj, Z=1.0, nucleus_center=self._center_b
                )
                H[i, j] = kinetic + attraction
        return H

    def overlap_matrix(self, alpha: float = 1.0) -> np.ndarray:
        """Build the 2x2 overlap matrix at exponent `alpha`.

        Parameters
        ----------
        alpha : float, default 1.0

        Returns
        -------
        ndarray, shape (2, 2)
        """
        chi_a, chi_b = self._basis(alpha)
        S = np.eye(2)
        S[0, 1] = S[1, 0] = overlap_integral(chi_a, chi_b)
        return S

    def solve(self, alpha: float = 1.0) -> EigenstateResult:
        r"""Solve the secular equation :math:`HC=SCE` at a given (fixed) orbital exponent.

        Parameters
        ----------
        alpha : float, default 1.0
            Shared Gaussian orbital exponent, in m^-2.

        Returns
        -------
        EigenstateResult
            Two states (bonding, antibonding); ``extra["nuclear_repulsion"]``
            and ``extra["total_energy"]`` (electronic ground state +
            nuclear repulsion) are also attached.
        """
        H = self.hamiltonian(alpha)
        S = self.overlap_matrix(alpha)
        energies, coefficients = solve_secular_equation(H, S)
        total_energy = float(energies[0] + self.nuclear_repulsion)
        return EigenstateResult(
            energies=energies,
            coefficients=coefficients,
            basis_labels=("H_A", "H_B"),
            overlap=S,
            extra={"nuclear_repulsion": self.nuclear_repulsion, "total_energy": total_energy, "alpha": float(alpha)},
        )

    def total_energy(self, alpha: float) -> float:
        r"""Total energy (bonding-orbital electronic energy + nuclear repulsion) at exponent `alpha`.

        Parameters
        ----------
        alpha : float
            Shared Gaussian orbital exponent, in m^-2 (must be positive).

        Returns
        -------
        float
            Energy, in J.
        """
        if alpha <= 0:
            raise ValueError("alpha must be positive")
        return self.solve(alpha).extra["total_energy"]

    def optimize_exponent(self, alpha_guess: float = 1.0 / (5.29e-11**2)) -> ExponentOptimizationResult:
        r"""Variationally optimize the shared Gaussian exponent to minimize the total energy.

        By the variational theorem, the computed ground-state energy at
        *any* value of `alpha` is an upper bound on the true H2+
        ground-state energy, so minimizing over `alpha` with
        :func:`scipy.optimize.minimize_scalar` produces a strictly better
        (lower or equal) energy than any single fixed guess -- a genuine,
        if minimal, variational calculation (Szabo & Ostlund, *Modern
        Quantum Chemistry*, 1st ed. rev., Ch. 1.3).

        Parameters
        ----------
        alpha_guess : float, default the exponent matching a hydrogen-atom-like 1s size
            Initial bracket center for the 1D minimization, in m^-2.

        Returns
        -------
        ExponentOptimizationResult

        Examples
        --------
        Optimizing the exponent never makes the energy worse than the
        naive starting guess -- the variational principle in action:

        >>> h2plus = H2PlusVariational(bond_length=106e-12)
        >>> naive_energy = h2plus.total_energy(alpha_guess := 1.0 / (5.29e-11 ** 2))
        >>> result = h2plus.optimize_exponent(alpha_guess)
        >>> bool(result.optimized_energy <= naive_energy)
        True
        """
        naive_energy = self.total_energy(alpha_guess)

        def objective(log_alpha):
            return self.total_energy(np.exp(log_alpha))

        log_alpha_guess = np.log(alpha_guess)
        res = minimize_scalar(objective, bounds=(log_alpha_guess - 3.0, log_alpha_guess + 3.0), method="bounded")
        optimized_alpha = float(np.exp(res.x))
        optimized_energy = float(res.fun)
        return ExponentOptimizationResult(
            optimized_alpha=optimized_alpha,
            optimized_energy=optimized_energy,
            naive_alpha=float(alpha_guess),
            naive_energy=float(naive_energy),
            converged=bool(res.success),
        )


@dataclass
class ExponentOptimizationResult:
    """Result of :meth:`H2PlusVariational.optimize_exponent`."""

    optimized_alpha: float
    """float: The variationally optimized Gaussian exponent, in m^-2."""

    optimized_energy: float
    """float: Total energy at ``optimized_alpha``, in J."""

    naive_alpha: float
    """float: The starting-guess exponent, in m^-2."""

    naive_energy: float
    """float: Total energy at ``naive_alpha``, in J."""

    converged: bool
    """bool: Whether the underlying 1D optimizer reported success."""

    @property
    def improvement(self) -> float:
        """float: How much lower the optimized energy is than the naive guess's, in J (>= 0 by the variational theorem)."""
        return self.naive_energy - self.optimized_energy


@dataclass
class SCFResult:
    """Result of :meth:`RestrictedHartreeFock.scf`."""

    total_energy: float
    """float: Electronic energy + nuclear repulsion, in J."""

    electronic_energy: float
    """float: :math:`\\frac12\\sum_{\\mu\\nu}P_{\\nu\\mu}(H_{\\mu\\nu}+F_{\\mu\\nu})`, in J."""

    nuclear_repulsion: float
    """float: Nucleus-nucleus Coulomb repulsion, in J."""

    orbital_energies: np.ndarray
    """ndarray, shape (n_basis,): Canonical MO energies, ascending, in J."""

    coefficients: np.ndarray
    """ndarray, shape (n_basis, n_basis): MO coefficients as columns (:math:`C^TSC=I`)."""

    density: np.ndarray
    """ndarray, shape (n_basis, n_basis): Converged density matrix :math:`P=2\\sum_a^{occ}C_{\\mu a}C_{\\nu a}`."""

    n_iterations: int
    """int: SCF iterations performed."""

    converged: bool
    """bool: Whether the density converged within the tolerance."""

    energy_history: np.ndarray
    """ndarray, shape (n_iterations,): Total energy after each iteration, in J."""


class RestrictedHartreeFock:
    r"""Closed-shell restricted Hartree-Fock (RHF) by the Roothaan-Hall SCF procedure.

    Each doubly occupied spatial orbital is a linear combination of basis
    functions, :math:`\psi_a=\sum_\mu C_{\mu a}\phi_\mu`, with coefficients
    solving the Roothaan-Hall equations :math:`FC=SC\varepsilon` (C. C. J.
    Roothaan, *Rev. Mod. Phys.* 23, 69 (1951); G. G. Hall, *Proc. R. Soc.
    Lond. A* 205, 541 (1951)). The Fock matrix

    .. math::

        F_{\mu\nu} = H^{core}_{\mu\nu} + \sum_{\lambda\sigma}P_{\lambda\sigma}
            \left[(\mu\nu|\sigma\lambda) - \tfrac12(\mu\lambda|\sigma\nu)\right]

    depends on the density :math:`P` built from its own solution, so the
    equations are iterated to self-consistency from a core-Hamiltonian
    (:math:`P=0`) guess (Szabo & Ostlund, *Modern Quantum Chemistry*, 1st
    ed. rev., Ch. 3.4.6). Each step reuses
    :func:`chemistrykit.quantum.utils.secular_equation.solve_secular_equation`
    for :math:`FC=SC\varepsilon`.

    Parameters
    ----------
    basis : sequence of ContractedGaussian
        s-type basis functions.
    nuclei : sequence of (float, array-like)
        ``(Z, center)`` for each nucleus; centers in m.
    n_electrons : int
        Total electron count; must be even (closed shell).

    Examples
    --------
    H2 at 1.4 bohr in STO-3G reproduces Szabo & Ostlund's -1.117 hartree:

    >>> from chemistrykit.quantum.systems.helium import HARTREE_ENERGY
    >>> result = RestrictedHartreeFock.h2().scf()
    >>> round(result.total_energy / HARTREE_ENERGY, 4)
    -1.1167
    """

    def __init__(self, basis, nuclei, n_electrons: int):
        if n_electrons <= 0 or n_electrons % 2:
            raise ValueError("RHF needs a positive, even number of electrons")
        if n_electrons // 2 > len(basis):
            raise ValueError("basis too small for the number of occupied orbitals")
        self.basis: list[ContractedGaussian] = list(basis)
        self.nuclei = [(float(Z), np.asarray(center, dtype=np.float64)) for Z, center in nuclei]
        self.n_electrons = int(n_electrons)

    @classmethod
    def h2(cls, bond_length: float = 1.4 * BOHR_RADIUS, zeta: float = 1.24) -> RestrictedHartreeFock:
        """H2 in the STO-3G minimal basis (Szabo & Ostlund, Ch. 3.5.2).

        Parameters
        ----------
        bond_length : float, default 1.4 bohr
            H-H distance, in m.
        zeta : float, default 1.24
            Slater exponent of each H 1s function.

        Returns
        -------
        RestrictedHartreeFock
        """
        a, b = np.array([0.0, 0.0, 0.0]), np.array([0.0, 0.0, bond_length])
        return cls([sto3g_1s(zeta, a), sto3g_1s(zeta, b)], [(1.0, a), (1.0, b)], n_electrons=2)

    @classmethod
    def heh_plus(cls, bond_length: float = 1.4632 * BOHR_RADIUS, zeta_he: float = 2.0925, zeta_h: float = 1.24) -> RestrictedHartreeFock:
        """HeH+ in the STO-3G minimal basis (Szabo & Ostlund, Ch. 3.5.3 and Appendix B).

        Parameters
        ----------
        bond_length : float, default 1.4632 bohr
            He-H distance, in m.
        zeta_he, zeta_h : float
            Slater exponents of the He and H 1s functions.

        Returns
        -------
        RestrictedHartreeFock
        """
        he, h = np.array([0.0, 0.0, 0.0]), np.array([0.0, 0.0, bond_length])
        return cls([sto3g_1s(zeta_he, he), sto3g_1s(zeta_h, h)], [(2.0, he), (1.0, h)], n_electrons=2)

    @property
    def nuclear_repulsion(self) -> float:
        """float: :math:`\\sum_{A<B}Z_AZ_Be^2/(4\\pi\\varepsilon_0R_{AB})`, in J."""
        total = 0.0
        for i, (Za, Ra) in enumerate(self.nuclei):
            for Zb, Rb in self.nuclei[i + 1 :]:
                total += Za * Zb * _COULOMB_CONSTANT / float(np.linalg.norm(Ra - Rb))
        return total

    def overlap_matrix(self) -> np.ndarray:
        """Overlap matrix :math:`S_{\\mu\\nu}`.

        Returns
        -------
        ndarray, shape (n_basis, n_basis)
        """
        return self._one_electron_matrix(lambda a, b: contracted_one_electron(overlap_integral, a, b))

    def core_hamiltonian(self) -> np.ndarray:
        """Core Hamiltonian :math:`H^{core}=T+\\sum_AV_A` (kinetic + nuclear attraction), in J.

        Returns
        -------
        ndarray, shape (n_basis, n_basis)
        """

        def element(a, b):
            h = contracted_one_electron(kinetic_integral, a, b)
            for Z, center in self.nuclei:
                h += contracted_one_electron(nuclear_attraction_integral, a, b, Z, center)
            return h

        return self._one_electron_matrix(element)

    def two_electron_integrals(self) -> np.ndarray:
        """All two-electron integrals :math:`(\\mu\\nu|\\lambda\\sigma)` (chemists' notation), in J.

        Returns
        -------
        ndarray, shape (n_basis, n_basis, n_basis, n_basis)
        """
        n = len(self.basis)
        eri = np.empty((n, n, n, n))
        for i in range(n):
            for j in range(i + 1):
                for k in range(n):
                    for l in range(k + 1):
                        if i * (i + 1) // 2 + j < k * (k + 1) // 2 + l:
                            continue
                        value = contracted_electron_repulsion(self.basis[i], self.basis[j], self.basis[k], self.basis[l])
                        for p, q, r, t in ((i, j, k, l), (j, i, k, l), (i, j, l, k), (j, i, l, k)):
                            eri[p, q, r, t] = eri[r, t, p, q] = value
        return eri

    def _one_electron_matrix(self, element) -> np.ndarray:
        n = len(self.basis)
        M = np.empty((n, n))
        for i in range(n):
            for j in range(i + 1):
                M[i, j] = M[j, i] = element(self.basis[i], self.basis[j])
        return M

    def scf(self, max_iter: int = 100, tol: float = 1e-10) -> SCFResult:
        """Iterate the Roothaan-Hall equations to self-consistency.

        Parameters
        ----------
        max_iter : int, default 100
        tol : float, default 1e-10
            Convergence threshold on the largest change of any density
            matrix element between iterations (dimensionless).

        Returns
        -------
        SCFResult
        """
        if max_iter < 1:
            raise ValueError("max_iter must be at least 1")
        S = self.overlap_matrix()
        H = self.core_hamiltonian()
        eri = self.two_electron_integrals()
        E_nuc = self.nuclear_repulsion
        n_occ = self.n_electrons // 2

        def fock(P):
            return H + np.einsum("mnls,sl->mn", eri, P) - 0.5 * np.einsum("mlsn,sl->mn", eri, P)

        P = np.zeros_like(S)
        history = []
        converged = False
        for _ in range(max_iter):
            eps, C = solve_secular_equation(fock(P), S)
            P_new = 2.0 * C[:, :n_occ] @ C[:, :n_occ].T
            history.append(0.5 * float(np.sum(P_new * (H + fock(P_new)))) + E_nuc)
            delta = float(np.max(np.abs(P_new - P)))
            P = P_new
            if delta < tol:
                converged = True
                break
        return SCFResult(
            total_energy=history[-1],
            electronic_energy=history[-1] - E_nuc,
            nuclear_repulsion=E_nuc,
            orbital_energies=eps,
            coefficients=C,
            density=P,
            n_iterations=len(history),
            converged=converged,
            energy_history=np.array(history),
        )
