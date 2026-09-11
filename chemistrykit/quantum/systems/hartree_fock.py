r"""A minimal variational (Hartree-Fock-style) LCAO treatment of H2+ in a 2-Gaussian basis.

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
from chemistrykit.quantum.utils.basis_sets import GaussianPrimitive, kinetic_integral, nuclear_attraction_integral, overlap_integral
from chemistrykit.quantum.utils.secular_equation import solve_secular_equation

__all__ = ["H2PlusVariational", "ExponentOptimizationResult"]

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
