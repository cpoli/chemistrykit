"""Tests for chemistrykit.quantum.utils.secular_equation."""

import numpy as np
import pytest

from chemistrykit.quantum.utils.secular_equation import solve_secular_equation


def test_orthonormal_basis_matches_plain_eigh():
    H = np.array([[2.0, -1.0], [-1.0, 2.0]])
    energies, coefficients = solve_secular_equation(H)
    expected_energies, expected_vectors = np.linalg.eigh(H)
    np.testing.assert_allclose(energies, expected_energies)
    np.testing.assert_allclose(np.abs(coefficients), np.abs(expected_vectors))


def test_generalized_problem_differs_from_ordinary_when_s_not_identity():
    H = np.array([[0.0, -1.0], [-1.0, 0.0]])
    S = np.array([[1.0, 0.3], [0.3, 1.0]])
    energies_ordinary, _ = solve_secular_equation(H)
    energies_generalized, _ = solve_secular_equation(H, S)
    assert not np.allclose(np.sort(energies_ordinary), np.sort(energies_generalized))


def test_generalized_eigenvectors_are_s_orthonormal():
    H = np.array([[1.0, 0.5], [0.5, -1.0]])
    S = np.array([[1.0, 0.2], [0.2, 1.0]])
    _, C = solve_secular_equation(H, S)
    np.testing.assert_allclose(C.T @ S @ C, np.eye(2), atol=1.0e-8)


def test_secular_equation_reproduces_hc_equals_sce():
    H = np.array([[1.0, 0.3], [0.3, -0.5]])
    S = np.array([[1.0, 0.1], [0.1, 1.0]])
    energies, C = solve_secular_equation(H, S)
    for i in range(2):
        lhs = H @ C[:, i]
        rhs = energies[i] * (S @ C[:, i])
        np.testing.assert_allclose(lhs, rhs, atol=1.0e-8)


def test_rejects_nonsquare_h():
    with pytest.raises(ValueError):
        solve_secular_equation(np.array([[1.0, 2.0, 3.0], [4.0, 5.0, 6.0]]))


def test_rejects_nonsymmetric_h():
    with pytest.raises(ValueError):
        solve_secular_equation(np.array([[1.0, 2.0], [0.0, 1.0]]))


def test_rejects_mismatched_s_shape():
    H = np.eye(2)
    S = np.eye(3)
    with pytest.raises(ValueError):
        solve_secular_equation(H, S)


def test_rejects_non_positive_definite_overlap():
    H = np.eye(2)
    S = np.array([[1.0, 2.0], [2.0, 1.0]])  # eigenvalues -1, 3: not PD
    with pytest.raises(ValueError):
        solve_secular_equation(H, S)
