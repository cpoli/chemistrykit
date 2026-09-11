"""Tests for chemistrykit.quantum.systems.hartree_fock (the minimal H2+ variational solver)."""

import numpy as np
import pytest

from chemistrykit.constants import ELECTRONVOLT
from chemistrykit.quantum.systems.hartree_fock import H2PlusVariational

BOHR_RADIUS = 5.29177e-11
H2_PLUS_BOND_LENGTH = 106.0e-12  # experimental equilibrium bond length, in m


def test_bonding_orbital_is_lower_than_antibonding():
    h2plus = H2PlusVariational(bond_length=H2_PLUS_BOND_LENGTH)
    energies = h2plus.solve(alpha=1.0 / BOHR_RADIUS**2).energies
    assert energies[0] < energies[1]


def test_overlap_matrix_is_symmetric_positive_definite_for_reasonable_alpha():
    h2plus = H2PlusVariational(bond_length=H2_PLUS_BOND_LENGTH)
    S = h2plus.overlap_matrix(alpha=1.0 / BOHR_RADIUS**2)
    assert np.allclose(S, S.T)
    eigenvalues = np.linalg.eigvalsh(S)
    assert np.all(eigenvalues > 0)


def test_overlap_diagonal_is_one_normalized_basis():
    h2plus = H2PlusVariational(bond_length=H2_PLUS_BOND_LENGTH)
    S = h2plus.overlap_matrix(alpha=1.0 / BOHR_RADIUS**2)
    assert S[0, 0] == pytest.approx(1.0)
    assert S[1, 1] == pytest.approx(1.0)


def test_nuclear_repulsion_matches_coulomb_law():
    h2plus = H2PlusVariational(bond_length=1.0e-10)
    from chemistrykit.constants import ELEMENTARY_CHARGE, VACUUM_PERMITTIVITY

    expected = ELEMENTARY_CHARGE**2 / (4.0 * np.pi * VACUUM_PERMITTIVITY * 1.0e-10)
    assert h2plus.nuclear_repulsion == pytest.approx(expected)


def test_optimized_exponent_improves_on_naive_guess():
    h2plus = H2PlusVariational(bond_length=H2_PLUS_BOND_LENGTH)
    alpha_naive = 1.0 / BOHR_RADIUS**2
    naive_energy = h2plus.total_energy(alpha_naive)
    result = h2plus.optimize_exponent(alpha_naive)
    assert result.optimized_energy <= naive_energy
    assert result.improvement >= 0.0
    assert result.converged is True


def test_optimized_total_energy_is_a_real_negative_bound_state():
    h2plus = H2PlusVariational(bond_length=H2_PLUS_BOND_LENGTH)
    result = h2plus.optimize_exponent()
    energy_ev = result.optimized_energy / ELECTRONVOLT
    # A real, physically sensible (if crude, single-Gaussian) bound-state
    # total energy: negative, and within an order of magnitude of the
    # experimental H2+ value (~-16.4 eV total electronic+nuclear energy).
    assert -30.0 < energy_ev < -5.0


def test_bond_length_rejects_nonpositive():
    with pytest.raises(ValueError):
        H2PlusVariational(bond_length=0.0)


def test_total_energy_rejects_nonpositive_alpha():
    h2plus = H2PlusVariational(bond_length=H2_PLUS_BOND_LENGTH)
    with pytest.raises(ValueError):
        h2plus.total_energy(0.0)


def test_solve_result_has_expected_labels_and_extra():
    h2plus = H2PlusVariational(bond_length=H2_PLUS_BOND_LENGTH)
    result = h2plus.solve(alpha=1.0 / BOHR_RADIUS**2)
    assert result.basis_labels == ("H_A", "H_B")
    assert "total_energy" in result.extra
    assert "nuclear_repulsion" in result.extra
