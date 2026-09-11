r"""Abstract base classes for quantum-chemistry models, and the eigenstate result container.

:mod:`chemistrykit.quantum` mixes two genuinely different numerical
shapes, unlike the single-shape domains built so far
(:mod:`chemistrykit.kinetics`'s ODE right-hand-sides,
:mod:`chemistrykit.md`'s velocity-Verlet trajectories):

* **Exactly solvable models** (:mod:`~chemistrykit.quantum.systems.particle_in_box`,
  :mod:`~chemistrykit.quantum.systems.harmonic_oscillator`,
  :mod:`~chemistrykit.quantum.systems.rigid_rotor`,
  :mod:`~chemistrykit.quantum.systems.hydrogenlike`) have a closed-form
  energy eigenvalue for each quantum number, the way
  :class:`chemistrykit.kinetics.core.base_system.RateLaw` has a
  closed-form integrated concentration -- no numerical linear algebra is
  needed at all. :class:`QuantumSystem` declares this common interface.
* **Variational/secular-equation models**
  (:mod:`~chemistrykit.quantum.systems.huckel`,
  :mod:`~chemistrykit.quantum.systems.hartree_fock`) build an explicit
  matrix Hamiltonian (and, where the basis is non-orthogonal, an overlap
  matrix) and obtain their spectrum by genuinely diagonalizing it --
  :func:`chemistrykit.quantum.utils.secular_equation.solve_secular_equation`
  solves :math:`HC=SCE` via :func:`numpy.linalg.eigh` (orthonormal basis)
  or :func:`scipy.linalg.eigh` (generalized problem), exactly the
  Hamiltonian-diagonalization pattern used throughout condensed-matter
  tight-binding theory (physicskit's
  ``physicskit.condensed.tight_binding.Hamiltonian.bands``) -- Huckel
  molecular-orbital theory *is* a tight-binding model of the pi system.
  :class:`VariationalSolver` declares this common interface, and
  :class:`EigenstateResult` is its stable return type (mirroring
  :class:`chemistrykit.kinetics.core.base_system.KineticsResult`).

See Atkins & de Paula, *Physical Chemistry*, 11th ed., Ch. 7-9, and
Levine, *Quantum Chemistry*, 7th ed., Ch. 2-4 and 11, for the general
framing of both families of model.
"""

from __future__ import annotations

from abc import ABC, abstractmethod
from collections.abc import Sequence
from dataclasses import dataclass, field
from typing import Optional

import numpy as np

__all__ = ["QuantumSystem", "VariationalSolver", "EigenstateResult"]


class QuantumSystem(ABC):
    """Common interface for an exactly solvable bound-state model with a discrete spectrum.

    Concrete subclasses (:class:`~chemistrykit.quantum.systems.particle_in_box.ParticleInBox1D`,
    :class:`~chemistrykit.quantum.systems.harmonic_oscillator.QuantumHarmonicOscillator`,
    :class:`~chemistrykit.quantum.systems.rigid_rotor.RigidRotor`,
    :class:`~chemistrykit.quantum.systems.hydrogenlike.HydrogenLikeAtom`, ...)
    implement :meth:`energy` with whatever quantum number(s) are natural
    for that model (a single principal quantum number `n`, or `(n, l)`,
    etc.) -- mirroring
    :class:`chemistrykit.kinetics.core.base_system.RateLaw`'s closed-form
    ``concentration(t)``, there is nothing to numerically integrate or
    diagonalize here.
    """

    @abstractmethod
    def energy(self, *quantum_numbers):
        """Return the energy eigenvalue for the given quantum number(s).

        Parameters
        ----------
        *quantum_numbers
            One or more quantum numbers, model-specific (e.g. a single
            `n` for a 1D box, `(n, l)` for a hydrogen-like atom).

        Returns
        -------
        float or ndarray
            Energy, in J.
        """


@dataclass
class EigenstateResult:
    """Container for the output of a :meth:`VariationalSolver.solve` call.

    Mirrors :class:`chemistrykit.kinetics.core.base_system.KineticsResult`
    and :class:`chemistrykit.md.core.base_system.MDResult` (a stable
    dataclass return type consumed by visualizers/tests), specialized to
    a diagonalized-Hamiltonian spectrum: ordered eigenvalues plus the
    matrix of eigenvectors (molecular-orbital coefficients), each column
    one eigenstate.
    """

    energies: np.ndarray
    """ndarray, shape (n_states,): Eigenvalues (orbital/state energies),
    ascending."""

    coefficients: np.ndarray
    """ndarray, shape (n_basis, n_states): Eigenvector matrix; column `i`
    is the basis-function expansion coefficients of the state with energy
    ``energies[i]``."""

    basis_labels: Sequence[str] = field(default_factory=tuple)
    """tuple of str: Human-readable label for each basis function (e.g.
    atom labels for Huckel pi orbitals, center labels for an LCAO
    basis), same order as ``coefficients``' rows."""

    overlap: Optional[np.ndarray] = None
    """ndarray, shape (n_basis, n_basis), or None: The overlap matrix `S`
    used to obtain this result, or ``None`` when the basis was assumed
    orthonormal (``S = I``), as in standard Huckel theory."""

    extra: dict = field(default_factory=dict)
    """dict: Free-form slot for additional diagnostics a solver chooses
    to attach (e.g. total electronic energy, nuclear repulsion)."""

    def orbital_energy(self, i: int) -> float:
        """Return the energy of the `i`-th state (0-indexed, ascending).

        Parameters
        ----------
        i : int

        Returns
        -------
        float
        """
        return float(self.energies[i])

    def orbital_coefficients(self, i: int) -> np.ndarray:
        """Return the basis-expansion coefficients of the `i`-th state.

        Parameters
        ----------
        i : int

        Returns
        -------
        ndarray, shape (n_basis,)
        """
        return self.coefficients[:, i]

    def is_degenerate(self, i: int, j: int, tol: float = 1e-9) -> bool:
        """Return whether states `i` and `j` are degenerate to within `tol`.

        Parameters
        ----------
        i, j : int
        tol : float, default 1e-9
            Absolute energy tolerance, in the same units as ``energies``.

        Returns
        -------
        bool

        Examples
        --------
        >>> import numpy as np
        >>> result = EigenstateResult(energies=np.array([-2.0, -1.0, -1.0, 1.0]), coefficients=np.eye(4))
        >>> result.is_degenerate(1, 2)
        True
        >>> result.is_degenerate(0, 1)
        False
        """
        return bool(abs(self.energies[i] - self.energies[j]) <= tol)


class VariationalSolver(ABC):
    """Common interface for a model whose spectrum comes from diagonalizing a matrix Hamiltonian.

    Concrete subclasses (:class:`~chemistrykit.quantum.systems.huckel.HuckelSystem`,
    :class:`~chemistrykit.quantum.systems.hartree_fock.H2PlusVariational`)
    build a Hamiltonian matrix `H` (and, for a non-orthogonal basis, an
    overlap matrix `S`) from their physical parameters and delegate the
    actual eigenvalue problem to
    :func:`chemistrykit.quantum.utils.secular_equation.solve_secular_equation`.
    """

    @abstractmethod
    def hamiltonian(self) -> np.ndarray:
        """Build and return this model's Hamiltonian matrix.

        Returns
        -------
        ndarray, shape (n_basis, n_basis)
        """

    @abstractmethod
    def solve(self) -> EigenstateResult:
        """Diagonalize the Hamiltonian and return the resulting spectrum.

        Returns
        -------
        EigenstateResult
        """
