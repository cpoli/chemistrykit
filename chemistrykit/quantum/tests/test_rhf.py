"""Tests for the restricted Hartree-Fock SCF solver against Szabo & Ostlund's STO-3G results."""

import numpy as np
import pytest

from chemistrykit.quantum.systems.hartree_fock import RestrictedHartreeFock
from chemistrykit.quantum.systems.helium import HARTREE_ENERGY
from chemistrykit.quantum.utils.basis_sets import (
    BOHR_RADIUS,
    GaussianPrimitive,
    contracted_one_electron,
    electron_repulsion_integral,
    overlap_integral,
    sto3g_1s,
)


def test_h2_sto3g_matches_szabo_ostlund():
    result = RestrictedHartreeFock.h2().scf()
    assert result.converged
    assert result.total_energy / HARTREE_ENERGY == pytest.approx(-1.1167, abs=1e-4)
    assert result.electronic_energy / HARTREE_ENERGY == pytest.approx(-1.8310, abs=1e-4)
    np.testing.assert_allclose(result.orbital_energies / HARTREE_ENERGY, [-0.5782, 0.6703], atol=1e-4)


def test_heh_plus_sto3g_matches_szabo_ostlund():
    result = RestrictedHartreeFock.heh_plus().scf()
    assert result.converged
    assert result.total_energy / HARTREE_ENERGY == pytest.approx(-2.860662, abs=1e-5)
    np.testing.assert_allclose(result.orbital_energies / HARTREE_ENERGY, [-1.597448, -0.061652], atol=1e-4)


def test_density_is_idempotent_and_holds_all_electrons():
    rhf = RestrictedHartreeFock.heh_plus()
    result = rhf.scf()
    S = rhf.overlap_matrix()
    assert np.trace(result.density @ S) == pytest.approx(2.0)
    PS = result.density @ S
    np.testing.assert_allclose(PS @ PS, 2.0 * PS, atol=1e-10)


def test_h2_binding_curve_has_minimum_near_1_35_bohr():
    R = np.linspace(1.0, 2.0, 21)
    E = [RestrictedHartreeFock.h2(bond_length=r * BOHR_RADIUS).scf().total_energy for r in R]
    assert R[int(np.argmin(E))] == pytest.approx(1.35, abs=0.05)


def test_eri_self_repulsion_closed_form():
    alpha = 1.0e20
    g = GaussianPrimitive(alpha, [0.0, 0.0, 0.0])
    coulomb = 8.9875517923e9 * 1.602176634e-19**2
    assert electron_repulsion_integral(g, g, g, g) == pytest.approx(coulomb * 2 * np.sqrt(alpha / np.pi), rel=1e-8)


def test_eri_permutational_symmetry():
    a = GaussianPrimitive(1.0e20, [0.0, 0.0, 0.0])
    b = GaussianPrimitive(3.0e20, [0.0, 0.0, 1.0e-10])
    c = GaussianPrimitive(2.0e20, [0.5e-10, 0.0, 0.0])
    ref = electron_repulsion_integral(a, b, c, a)
    for perm in [(b, a, c, a), (a, b, a, c), (c, a, a, b), (a, c, b, a)]:
        assert electron_repulsion_integral(*perm) == pytest.approx(ref, rel=1e-12)


def test_sto3g_function_is_normalized():
    phi = sto3g_1s(zeta=2.0925, center=[0.0, 0.0, 0.0])
    assert contracted_one_electron(overlap_integral, phi, phi) == pytest.approx(1.0, abs=1e-5)


@pytest.mark.parametrize("n_electrons", [0, 3])
def test_invalid_electron_counts_raise(n_electrons):
    phi = sto3g_1s(1.24, [0.0, 0.0, 0.0])
    with pytest.raises(ValueError):
        RestrictedHartreeFock([phi], [(1.0, [0.0, 0.0, 0.0])], n_electrons=n_electrons)
