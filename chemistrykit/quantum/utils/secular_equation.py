r"""Shared generalized-eigenvalue ("secular equation") solver.

Every variational model in this subpackage --
:class:`chemistrykit.quantum.systems.huckel.HuckelSystem` (an orthonormal
atomic-orbital basis) and
:class:`chemistrykit.quantum.systems.hartree_fock.H2PlusVariational` (a
non-orthogonal Gaussian basis) -- reduces to the same linear-algebra
problem: expand a trial wavefunction in a finite basis,
:math:`|\psi\rangle=\sum_i c_i|\chi_i\rangle`, and require the expansion
coefficients to make the energy stationary. This gives the matrix
generalized eigenvalue problem

.. math::

    HC = SCE

(the "secular equation" or "Roothaan equation"; Szabo & Ostlund, *Modern
Quantum Chemistry*, 1st ed. rev., Ch. 3.4; Levine, *Quantum Chemistry*,
7th ed., Ch. 11.1), where :math:`H_{ij}=\langle\chi_i|\hat H|\chi_j\rangle`
is the Hamiltonian matrix, :math:`S_{ij}=\langle\chi_i|\chi_j\rangle` is
the overlap matrix, `C` collects the eigenvectors as columns, and `E` is
the diagonal matrix of orbital energies. When the basis is already
orthonormal (:math:`S=I`, standard Huckel theory's approximation), this
reduces to the ordinary eigenvalue problem :math:`HC=CE`.

Both `H` and `S` are Hermitian by construction (real-symmetric for the
real bases used throughout this subpackage), so
:func:`numpy.linalg.eigh`/:func:`scipy.linalg.eigh` -- which exploit that
structure for speed and numerical stability, and guarantee real
eigenvalues -- are the correct and only tool needed; this is exactly the
diagonalization pattern physicskit's tight-binding Bloch Hamiltonians use
(``physicskit.condensed.tight_binding.Hamiltonian.bands``), since Huckel
theory *is* a (molecular, real-space) tight-binding model of the pi
system.
"""

from __future__ import annotations

import numpy as np
import scipy.linalg as sla

__all__ = ["solve_secular_equation"]


def solve_secular_equation(H, S=None):
    r"""Solve the secular equation :math:`HC=SCE` for a Hermitian `H` and (optional) overlap `S`.

    Parameters
    ----------
    H : array-like, shape (n_basis, n_basis)
        Hamiltonian matrix in the chosen basis. Must be Hermitian (only
        the values are checked for approximate symmetry, not enforced).
    S : array-like, shape (n_basis, n_basis), optional
        Overlap matrix. If omitted (or ``None``), the basis is assumed
        orthonormal (:math:`S=I`) and the ordinary eigenvalue problem is
        solved via :func:`numpy.linalg.eigh` -- the approximation
        standard Huckel theory makes (neglect of differential overlap
        between atomic pi orbitals). If given, the true generalized
        problem is solved via :func:`scipy.linalg.eigh`, which requires
        `S` to be symmetric positive-definite; a non-positive-definite
        `S` means the basis is linearly dependent (e.g. two basis
        functions placed on top of each other) and raises `ValueError`.

    Returns
    -------
    energies : ndarray, shape (n_basis,)
        Eigenvalues, ascending.
    coefficients : ndarray, shape (n_basis, n_basis)
        Eigenvectors as columns, `S`-orthonormal (:math:`C^TSC=I`) when
        `S` is given, orthonormal (:math:`C^TC=I`) otherwise.

    Raises
    ------
    ValueError
        If `S` is given but is not symmetric positive-definite (a
        linearly dependent basis).

    Examples
    --------
    A trivial 2x2 orthonormal-basis case (`S` omitted) recovers the plain
    eigenvalues of a symmetric matrix:

    >>> import numpy as np
    >>> H = np.array([[0.0, -1.0], [-1.0, 0.0]])
    >>> energies, C = solve_secular_equation(H)
    >>> np.round(energies, 6)
    array([-1.,  1.])

    A non-orthogonal basis (`S != I`) genuinely changes the eigenvalues
    relative to the naive (ordinary) eigenvalue problem:

    >>> S = np.array([[1.0, 0.5], [0.5, 1.0]])
    >>> energies_S, _ = solve_secular_equation(H, S)
    >>> bool(np.any(np.abs(np.sort(energies_S) - np.sort(energies)) > 1e-9))
    True
    """
    H = np.asarray(H, dtype=np.float64)
    if H.ndim != 2 or H.shape[0] != H.shape[1]:
        raise ValueError("H must be a square matrix")
    if not np.allclose(H, H.T, atol=1e-8):
        raise ValueError("H must be (numerically) symmetric/Hermitian")

    if S is None:
        energies, coefficients = np.linalg.eigh(H)
        return energies, coefficients

    S = np.asarray(S, dtype=np.float64)
    if S.shape != H.shape:
        raise ValueError("S must have the same shape as H")
    if not np.allclose(S, S.T, atol=1e-8):
        raise ValueError("S must be (numerically) symmetric")
    try:
        energies, coefficients = sla.eigh(H, S)
    except sla.LinAlgError as exc:
        raise ValueError("S is not symmetric positive-definite -- the basis is linearly dependent") from exc
    return energies, coefficients
