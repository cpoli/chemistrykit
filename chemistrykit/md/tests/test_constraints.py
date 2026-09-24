"""Tests for chemistrykit.md.systems.constraints (SHAKE)."""

import numpy as np
import pytest

from chemistrykit.md.systems.constraints import ShakeMolecule, shake


def test_shake_restores_constraint_along_old_bond_with_mass_weighting():
    old = np.array([[0.0, 0.0, 0.0], [1.0, 0.0, 0.0]])
    new = np.array([[0.0, 0.0, 0.0], [1.3, 0.0, 0.0]])
    fixed, _ = shake(new, old, [(0, 1, 1.0)], masses=[3.0, 1.0])
    # correction of -0.3 split inversely to mass: atom 0 moves +0.075, atom 1 moves -0.225
    np.testing.assert_allclose(fixed[:, 0], [0.075, 1.075], atol=1e-12)
    np.testing.assert_allclose(3.0 * fixed[0] + fixed[1], 3.0 * new[0] + new[1])  # center of mass unchanged


def test_rigid_rotor_keeps_bond_length_and_angular_velocity():
    pos = np.array([[0.0, 0.0, 0.0], [1.0, 0.0, 0.0]])
    vel = np.array([[0.0, -0.5, 0.0], [0.0, 0.5, 0.0]])
    rotor = ShakeMolecule(pos, vel, [1.0, 1.0], constraints=[(0, 1, 1.0)])
    for _ in range(2000):
        rotor.step(0.01)
    assert rotor.bond_lengths()[0] == pytest.approx(1.0, abs=1e-9)
    # free rigid rotor: kinetic energy (and so angular speed 1.0) conserved exactly
    assert rotor.kinetic_energy() == pytest.approx(0.25, rel=1e-8)
    assert rotor.degrees_of_freedom() == 2


def test_rigid_bonds_flexible_angle_conserve_energy():
    theta0 = np.deg2rad(104.5)
    pos = np.array([[np.sin(theta0 / 2), np.cos(theta0 / 2), 0.0], [0.0, 0.0, 0.0], [-np.sin(theta0 / 2), np.cos(theta0 / 2), 0.0]])
    vel = np.zeros_like(pos)
    vel[0] = [0.4, -0.2, 0.1]
    vel[2] = [-0.1, 0.3, 0.0]
    mol = ShakeMolecule(pos, vel, [1.0, 16.0, 1.0], constraints=[(0, 1, 1.0), (1, 2, 1.0)], angles=[(0, 1, 2, 50.0, theta0)])
    e0 = mol.kinetic_energy() + mol.potential_energy()
    for _ in range(3000):
        mol.step(0.002)
    np.testing.assert_allclose(mol.bond_lengths(), 1.0, atol=1e-9)
    assert mol.kinetic_energy() + mol.potential_energy() == pytest.approx(e0, rel=1e-3)


def test_shake_raises_when_not_converged():
    old = np.array([[0.0, 0.0, 0.0], [1.0, 0.0, 0.0]])
    new = np.array([[0.0, 0.0, 0.0], [5.0, 0.0, 0.0]])
    with pytest.raises(RuntimeError):
        shake(new, old, [(0, 1, 1.0)], masses=[1.0, 1.0], max_iter=1)
