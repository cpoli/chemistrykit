"""Tests for chemistrykit.quantum.systems.perturbation: RS perturbation theory vs. exact diagonalization."""

import numpy as np
import pytest

from chemistrykit.constants import HBAR
from chemistrykit.quantum.systems.perturbation import (
    anharmonic_energy_levels,
    anharmonic_hamiltonian_matrix,
    cubic_perturbation_first_order_correction,
    position_operator_matrix,
    quartic_perturbation_first_order_correction,
)

MASS = 1.6e-27
OMEGA = 1.0e14


def test_cubic_first_order_correction_always_zero_by_parity():
    for n in range(10):
        assert cubic_perturbation_first_order_correction(n) == 0.0


def test_cubic_correction_rejects_negative_n():
    with pytest.raises(ValueError):
        cubic_perturbation_first_order_correction(-1)


def test_quartic_correction_matches_ground_state_gaussian_fourth_moment():
    # For the HO ground state, <x^4> = 3<x^2>^2 (Gaussian fourth moment),
    # with <x^2>_0 = hbar/(2*m*omega).
    b = 1.0
    correction = quartic_perturbation_first_order_correction(0, MASS, OMEGA, b)
    x2 = HBAR / (2.0 * MASS * OMEGA)
    expected = b * 3.0 * x2**2
    assert correction == pytest.approx(expected)


def test_quartic_correction_increases_with_n():
    b = 1.0e18
    corrections = [quartic_perturbation_first_order_correction(n, MASS, OMEGA, b) for n in range(6)]
    for i in range(len(corrections) - 1):
        assert corrections[i + 1] > corrections[i]


def test_quartic_correction_rejects_negative_n():
    with pytest.raises(ValueError):
        quartic_perturbation_first_order_correction(-1, MASS, OMEGA, 1.0)


def test_position_operator_matrix_is_symmetric_tridiagonal_zero_diagonal():
    X = position_operator_matrix(5, MASS, OMEGA)
    assert np.allclose(X, X.T)
    assert np.allclose(np.diag(X), 0.0)
    # Off-tridiagonal entries are zero.
    assert np.allclose(X[0, 2:], 0.0)


def test_position_operator_matrix_matches_ladder_operator_formula():
    X = position_operator_matrix(4, MASS, OMEGA)
    beta = np.sqrt(HBAR / (2.0 * MASS * OMEGA))
    for i in range(1, 4):
        assert X[i, i - 1] == pytest.approx(beta * np.sqrt(i))


def test_position_operator_matrix_rejects_too_small_basis():
    with pytest.raises(ValueError):
        position_operator_matrix(1, MASS, OMEGA)


def test_anharmonic_hamiltonian_reduces_to_harmonic_when_b_and_c_zero():
    H = anharmonic_hamiltonian_matrix(6, MASS, OMEGA, b=0.0, c=0.0)
    n = np.arange(6)
    expected_diag = HBAR * OMEGA * (n + 0.5)
    np.testing.assert_allclose(np.diag(H), expected_diag)
    off_diag = H - np.diag(np.diag(H))
    assert np.allclose(off_diag, 0.0)


def test_exact_diagonalization_matches_first_order_perturbation_theory_for_weak_coupling():
    b = 1.0e18  # weak perturbation
    exact = anharmonic_energy_levels(n_basis=40, mass=MASS, omega=OMEGA, b=b, n_levels=1)
    E0_perturbative = 0.5 * HBAR * OMEGA + quartic_perturbation_first_order_correction(0, MASS, OMEGA, b)
    relative_error = abs(exact[0] - E0_perturbative) / abs(E0_perturbative)
    assert relative_error < 1.0e-3


def test_perturbation_theory_degrades_for_strong_coupling():
    b_weak, b_strong = 1.0e15, 1.0e21
    errors = {}
    for b in (b_weak, b_strong):
        exact = anharmonic_energy_levels(n_basis=40, mass=MASS, omega=OMEGA, b=b, n_levels=1)
        E0_perturbative = 0.5 * HBAR * OMEGA + quartic_perturbation_first_order_correction(0, MASS, OMEGA, b)
        errors[b] = abs(exact[0] - E0_perturbative) / abs(E0_perturbative)
    assert errors[b_strong] > errors[b_weak]


def test_anharmonic_energy_levels_are_ascending():
    levels = anharmonic_energy_levels(n_basis=20, mass=MASS, omega=OMEGA, b=1.0e18, n_levels=5)
    assert np.all(np.diff(levels) > 0)


def test_anharmonic_energy_levels_rejects_too_many_levels():
    with pytest.raises(ValueError):
        anharmonic_energy_levels(n_basis=5, mass=MASS, omega=OMEGA, b=1.0, n_levels=10)
