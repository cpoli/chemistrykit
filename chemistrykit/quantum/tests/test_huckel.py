"""Tests for chemistrykit.quantum.systems.huckel against known closed-form Huckel spectra."""

import numpy as np
import pytest

from chemistrykit.quantum.systems.huckel import (
    HuckelSystem,
    cyclic_polyene_eigenvalues,
    is_aromatic_by_huckel_rule,
    linear_polyene_eigenvalues,
)


def test_ethene_eigenvalues():
    ethene = HuckelSystem(n_atoms=2, bonds=[(0, 1)])
    energies = ethene.solve().energies
    np.testing.assert_allclose(np.sort(energies), [-1.0, 1.0], atol=1.0e-9)


def test_butadiene_matches_coulson_closed_form():
    butadiene = HuckelSystem.linear_polyene(4)
    numeric = np.sort(butadiene.solve().energies)
    closed_form = linear_polyene_eigenvalues(4)
    np.testing.assert_allclose(numeric, closed_form, atol=1.0e-9)
    # Golden-ratio structure of the linear tetraene Huckel spectrum.
    phi = (1.0 + np.sqrt(5.0)) / 2.0
    np.testing.assert_allclose(sorted(np.abs(numeric)), [phi - 1.0, phi - 1.0, phi, phi], atol=1.0e-6)


def test_benzene_matches_known_degenerate_pattern():
    benzene = HuckelSystem.cyclic_polyene(6)
    numeric = np.sort(benzene.solve().energies)
    expected = np.array([-2.0, -1.0, -1.0, 1.0, 1.0, 2.0])
    np.testing.assert_allclose(numeric, expected, atol=1.0e-9)


def test_benzene_matches_frost_circle_closed_form():
    numeric = np.sort(HuckelSystem.cyclic_polyene(6).solve().energies)
    closed_form = cyclic_polyene_eigenvalues(6)
    np.testing.assert_allclose(numeric, closed_form, atol=1.0e-9)


def test_hamiltonian_is_symmetric_with_correct_diagonal_and_bonds():
    system = HuckelSystem.linear_polyene(4, alpha=-2.0, beta=-0.7)
    H = system.hamiltonian()
    assert np.allclose(H, H.T)
    assert np.allclose(np.diag(H), -2.0)
    assert H[0, 1] == pytest.approx(-0.7)
    assert H[0, 2] == pytest.approx(0.0)


def test_benzene_is_aromatic_by_huckel_rule():
    energies = HuckelSystem.cyclic_polyene(6).solve().energies
    assert is_aromatic_by_huckel_rule(6, energies) is True


def test_cyclobutadiene_is_not_aromatic_4n_electron_count():
    energies = HuckelSystem.cyclic_polyene(4).solve().energies
    assert is_aromatic_by_huckel_rule(4, energies) is False


def test_cyclooctatetraene_8_electrons_not_closed_shell():
    # COT: 8 pi electrons is not 4n+2 for any integer n (it's 4n), so
    # this should fail on the electron-count check alone.
    energies = HuckelSystem.cyclic_polyene(8).solve().energies
    assert is_aromatic_by_huckel_rule(8, energies) is False


def test_cyclopentadienyl_cation_4_pi_electrons_not_aromatic():
    # A 4n electron count on a 5-membered Huckel ring: not aromatic.
    energies = HuckelSystem.cyclic_polyene(5).solve().energies
    assert is_aromatic_by_huckel_rule(4, energies) is False


def test_aromaticity_rule_rejects_odd_electron_count():
    energies = HuckelSystem.cyclic_polyene(6).solve().energies
    assert is_aromatic_by_huckel_rule(5, energies) is False


def test_pi_electron_energy_benzene_matches_known_value():
    # Benzene ground state: 6 electrons fill the 3 lowest MOs
    # (alpha+2beta once, alpha+beta twice), so E_pi = 2(alpha+2beta) + 4(alpha+beta) = 6*alpha + 8*beta.
    benzene = HuckelSystem.cyclic_polyene(6, alpha=0.0, beta=-1.0)
    assert benzene.pi_electron_energy(6) == pytest.approx(8.0 * -1.0)


def test_benzene_delocalization_energy_is_negative_stabilizing():
    benzene = HuckelSystem.cyclic_polyene(6, alpha=0.0, beta=-1.0)
    # 3 formal double bonds -> 6 pi electrons.
    assert benzene.delocalization_energy(6) < 0.0


def test_linear_polyene_requires_at_least_one_atom():
    with pytest.raises(ValueError):
        HuckelSystem(n_atoms=0, bonds=[])


def test_cyclic_polyene_requires_at_least_three_atoms():
    with pytest.raises(ValueError):
        HuckelSystem.cyclic_polyene(2)


def test_rejects_nonnegative_beta():
    with pytest.raises(ValueError):
        HuckelSystem(n_atoms=2, bonds=[(0, 1)], beta=0.5)


def test_labels_default_and_custom():
    system = HuckelSystem.linear_polyene(3)
    assert system.labels == ("C1", "C2", "C3")
    custom = HuckelSystem(n_atoms=2, bonds=[(0, 1)], labels=["A", "B"])
    assert custom.solve().basis_labels == ("A", "B")


NAPHTHALENE_BONDS = [(0, 1), (1, 2), (2, 3), (3, 4), (4, 5), (5, 6), (6, 7), (7, 8), (8, 9), (9, 0), (4, 9)]


def test_naphthalene_homo_density_favors_alpha_position():
    # Fukui, Yonezawa & Shingu (1952): HOMO coefficients 0.425 (alpha), 0.263 (beta), 0 (bridgehead).
    naphthalene = HuckelSystem(n_atoms=10, bonds=NAPHTHALENE_BONDS)
    density = naphthalene.frontier_electron_density(10, "homo")
    assert density.sum() == pytest.approx(2.0)
    alpha_density = (5.0 + np.sqrt(5.0)) / 20.0  # 2 c_alpha^2 = 0.3618
    np.testing.assert_allclose(density[[0, 3, 5, 8]], alpha_density, atol=1e-9)
    np.testing.assert_allclose(density[[1, 2, 6, 7]], 0.5 - alpha_density, atol=1e-9)
    np.testing.assert_allclose(density[[4, 9]], 0.0, atol=1e-12)


def test_butadiene_frontier_densities_closed_form():
    butadiene = HuckelSystem.linear_polyene(4)
    c = np.sqrt(2.0 / 5.0) * np.sin(np.arange(1, 5) * 2.0 * np.pi / 5.0)  # HOMO, k=2
    np.testing.assert_allclose(butadiene.frontier_electron_density(4, "homo"), 2.0 * c**2, atol=1e-9)
    np.testing.assert_allclose(butadiene.frontier_electron_density(4, "lumo").sum(), 2.0)


def test_frontier_density_rejects_degenerate_and_bad_input():
    benzene = HuckelSystem.cyclic_polyene(6)
    with pytest.raises(ValueError):
        benzene.frontier_electron_density(6)
    with pytest.raises(ValueError):
        benzene.frontier_electron_density(5)
    with pytest.raises(ValueError):
        HuckelSystem.linear_polyene(4).frontier_electron_density(4, "somo")
