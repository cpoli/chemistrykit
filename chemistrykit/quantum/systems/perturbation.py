r"""Rayleigh-Schrodinger perturbation theory applied to the anharmonic oscillator.

See Atkins & de Paula, *Physical Chemistry*, 11th ed., Ch. 8 ("The
chemical bond") for anharmonicity, and Levine, *Quantum Chemistry*, 7th
ed., Ch. 9.1-9.3 for Rayleigh-Schrodinger perturbation theory generally
(Griffiths, *Introduction to Quantum Mechanics*, 2nd ed., Ch. 6, works
this exact cubic/quartic-oscillator example in detail).

For :math:`\hat H=\hat H_0+\lambda\hat H'`, the first-order energy
correction is :math:`E_n^{(1)}=\langle n|\hat H'|n\rangle`, evaluated
here in the unperturbed harmonic-oscillator eigenbasis using the ladder
operators :math:`x=\sqrt{\hbar/(2m\omega)}(a+a^\dagger)` (:func:`position_operator_matrix`),
represented as a truncated matrix so the same
:func:`chemistrykit.quantum.utils.secular_equation.solve_secular_equation`
used by Huckel/Hartree-Fock also gives the *exact* (within the
truncation) numerically diagonalized spectrum of the full anharmonic
Hamiltonian, for comparison against the perturbative estimate.
"""

from __future__ import annotations

import numpy as np

from chemistrykit.constants import HBAR
from chemistrykit.quantum.utils.secular_equation import solve_secular_equation

__all__ = [
    "position_operator_matrix",
    "cubic_perturbation_first_order_correction",
    "quartic_perturbation_first_order_correction",
    "anharmonic_hamiltonian_matrix",
    "anharmonic_energy_levels",
]


def position_operator_matrix(n_basis: int, mass: float, omega: float) -> np.ndarray:
    r"""Matrix of the position operator :math:`\hat x` in the truncated harmonic-oscillator number basis.

    .. math::

        \hat x = \sqrt{\frac{\hbar}{2m\omega}}\left(\hat a+\hat a^\dagger\right)

    with :math:`\langle n-1|\hat a|n\rangle=\langle n|\hat a^\dagger|n-1\rangle=\sqrt n`
    the standard ladder-operator matrix elements (Levine, *Quantum
    Chemistry*, 7th ed., eq. 2.44), giving a tridiagonal matrix in the
    :math:`\{|0\rangle,\dots,|n_{basis}-1\rangle\}` basis. Powers of this
    matrix (via ordinary matrix multiplication) give the matrix elements
    of :math:`x^2`, :math:`x^3`, :math:`x^4`, etc. -- exact within the
    truncated basis, and increasingly accurate for the lowest states as
    `n_basis` grows.

    Parameters
    ----------
    n_basis : int
        Truncated basis size (>= 2).
    mass : float
        Oscillator mass, in kg.
    omega : float
        Angular frequency, in rad/s.

    Returns
    -------
    ndarray, shape (n_basis, n_basis)
        Symmetric tridiagonal, in m.

    Examples
    --------
    >>> X = position_operator_matrix(3, mass=1.6e-27, omega=1.0e14)
    >>> bool(np.allclose(X, X.T))
    True
    >>> bool(np.allclose(np.diag(X), 0.0))  # <n|x|n> = 0 for every n (parity)
    True
    """
    if n_basis < 2:
        raise ValueError("n_basis must be >= 2")
    beta = np.sqrt(HBAR / (2.0 * mass * omega))
    X = np.zeros((n_basis, n_basis))
    for i in range(1, n_basis):
        val = beta * np.sqrt(i)
        X[i, i - 1] = val
        X[i - 1, i] = val
    return X


def cubic_perturbation_first_order_correction(n: int) -> float:
    r"""First-order energy correction from a cubic perturbation :math:`\hat H'=c\hat x^3`: always exactly zero.

    :math:`\langle n|\hat x^3|n\rangle=0` for every harmonic-oscillator
    eigenstate `n`, by parity: :math:`\hat x^3` is an odd function of
    :math:`\hat x`, while :math:`|n\rangle` has definite parity
    :math:`(-1)^n`, so the integrand :math:`\psi_n\hat x^3\psi_n` is
    always odd and integrates to zero over the symmetric domain
    :math:`(-\infty,\infty)` (Griffiths, *Introduction to Quantum
    Mechanics*, 2nd ed., Ch. 6.3). A real cubic term (e.g. the leading
    anharmonic correction to a Morse-like bond potential) therefore only
    shifts the energy levels at *second* order in perturbation theory --
    not implemented here, but this is why :func:`quartic_perturbation_first_order_correction`
    is the one that matters at first order.

    Parameters
    ----------
    n : int
        Vibrational quantum number, >= 0 (unused -- the result is zero
        for every `n`; kept as a parameter for a uniform calling
        convention with :func:`quartic_perturbation_first_order_correction`).

    Returns
    -------
    float
        Exactly ``0.0``.

    Examples
    --------
    >>> cubic_perturbation_first_order_correction(5)
    0.0
    """
    if n < 0:
        raise ValueError("n must be >= 0")
    return 0.0


def quartic_perturbation_first_order_correction(n: int, mass: float, omega: float, b: float) -> float:
    r"""First-order energy correction from a quartic perturbation :math:`\hat H'=b\hat x^4`.

    .. math::

        E_n^{(1)} = b\langle n|\hat x^4|n\rangle
                  = \frac{3b\hbar^2}{4m^2\omega^2}\left(2n^2+2n+1\right)

    using :math:`\langle n|\hat x^4|n\rangle=3\left(\frac{\hbar}{2m\omega}\right)^2(2n^2+2n+1)`
    (Griffiths, *Introduction to Quantum Mechanics*, 2nd ed., Ch. 6
    problem on the quartic perturbation; at :math:`n=0` this reduces to
    the Gaussian ground state's fourth moment :math:`\langle x^4\rangle=3\langle x^2\rangle^2`).
    Unlike the cubic term, this survives at first order because
    :math:`\hat x^4` is even, matching :math:`|n\rangle`'s squared
    (always-even) probability density.

    Parameters
    ----------
    n : int
        Vibrational quantum number, >= 0.
    mass : float
        Oscillator mass, in kg.
    omega : float
        Angular frequency, in rad/s.
    b : float
        Quartic perturbation strength, in J/m^4.

    Returns
    -------
    float
        Energy correction, in J.

    Examples
    --------
    The correction grows with `n` (higher states are more delocalized,
    probing the quartic term's steep walls more):

    >>> correction_0 = quartic_perturbation_first_order_correction(0, mass=1.6e-27, omega=1.0e14, b=1.0e18)
    >>> correction_5 = quartic_perturbation_first_order_correction(5, mass=1.6e-27, omega=1.0e14, b=1.0e18)
    >>> bool(correction_5 > correction_0 > 0)
    True
    """
    if n < 0:
        raise ValueError("n must be >= 0")
    return b * 3.0 * HBAR**2 / (4.0 * mass**2 * omega**2) * (2.0 * n**2 + 2.0 * n + 1.0)


def anharmonic_hamiltonian_matrix(n_basis: int, mass: float, omega: float, b: float, c: float = 0.0) -> np.ndarray:
    r"""Build the matrix of :math:`\hat H=\hat H_0+c\hat x^3+b\hat x^4` in the truncated harmonic-oscillator basis.

    :math:`\hat H_0` is diagonal (:math:`\hbar\omega(n+1/2)`); the
    perturbation terms are built from powers of
    :func:`position_operator_matrix`. This is the *unperturbed-basis
    matrix-mechanics* approach to the anharmonic oscillator: truncating
    at finite `n_basis` is itself an approximation (states near the top
    of the truncated basis are contaminated by the missing higher
    states), so only the lowest handful of eigenvalues (well below
    `n_basis`) should be trusted -- :func:`anharmonic_energy_levels`
    exists specifically to hand back only those.

    Parameters
    ----------
    n_basis : int
        Truncated basis size.
    mass : float
        Oscillator mass, in kg.
    omega : float
        Angular frequency, in rad/s.
    b : float
        Quartic perturbation strength, in J/m^4.
    c : float, default 0.0
        Cubic perturbation strength, in J/m^3.

    Returns
    -------
    ndarray, shape (n_basis, n_basis)
        Symmetric, in J.
    """
    n = np.arange(n_basis)
    H0 = np.diag(HBAR * omega * (n + 0.5))
    X = position_operator_matrix(n_basis, mass, omega)
    X3 = X @ X @ X
    X4 = X3 @ X
    return H0 + c * X3 + b * X4


def anharmonic_energy_levels(n_basis: int, mass: float, omega: float, b: float, c: float = 0.0, n_levels: int = 3) -> np.ndarray:
    r"""Exact (within the truncated basis) anharmonic energy levels, by direct numerical diagonalization.

    Diagonalizes :func:`anharmonic_hamiltonian_matrix` via
    :func:`chemistrykit.quantum.utils.secular_equation.solve_secular_equation`
    (an ordinary eigenvalue problem -- the harmonic-oscillator basis is
    orthonormal) and returns the lowest `n_levels` eigenvalues, the
    numerically "exact" benchmark against which
    :func:`quartic_perturbation_first_order_correction`'s perturbative
    estimate is checked.

    Parameters
    ----------
    n_basis : int
        Truncated basis size; should be well above `n_levels` for the
        returned levels to be trustworthy (see
        :func:`anharmonic_hamiltonian_matrix`'s truncation caveat).
    mass : float
        Oscillator mass, in kg.
    omega : float
        Angular frequency, in rad/s.
    b : float
        Quartic perturbation strength, in J/m^4.
    c : float, default 0.0
        Cubic perturbation strength, in J/m^3.
    n_levels : int, default 3
        Number of lowest levels to return.

    Returns
    -------
    ndarray, shape (n_levels,)
        Ascending.

    Examples
    --------
    For a weak quartic perturbation, exact diagonalization matches the
    first-order perturbative estimate to within a small relative error:

    >>> mass, omega, b = 1.6e-27, 1.0e14, 1.0e18  # a weak perturbation
    >>> exact = anharmonic_energy_levels(n_basis=40, mass=mass, omega=omega, b=b, n_levels=1)
    >>> E0_harmonic = 0.5 * HBAR * omega
    >>> E0_perturbative = E0_harmonic + quartic_perturbation_first_order_correction(0, mass, omega, b)
    >>> relative_error = abs(exact[0] - E0_perturbative) / abs(E0_perturbative)
    >>> bool(relative_error < 1.0e-3)
    True
    """
    if n_levels > n_basis:
        raise ValueError("n_levels must be <= n_basis")
    H = anharmonic_hamiltonian_matrix(n_basis, mass, omega, b, c)
    energies, _ = solve_secular_equation(H, S=None)
    return np.sort(energies)[:n_levels]
