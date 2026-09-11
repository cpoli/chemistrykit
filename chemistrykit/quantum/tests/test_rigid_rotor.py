"""Tests for chemistrykit.quantum.systems.rigid_rotor against closed-form results."""

import pytest
import scipy.constants as sc

from chemistrykit.quantum.systems.rigid_rotor import RigidRotor


def test_ground_state_energy_is_zero():
    rotor = RigidRotor(moment_of_inertia=1.45e-46)
    assert rotor.energy(0) == pytest.approx(0.0)


def test_energy_matches_closed_form():
    rotor = RigidRotor(moment_of_inertia=1.45e-46)
    B = rotor.rotational_constant
    for J in range(6):
        assert rotor.energy(J) == pytest.approx(J * (J + 1) * B)


def test_degeneracy_is_two_j_plus_one():
    rotor = RigidRotor(moment_of_inertia=1.0e-46)
    for J in range(6):
        assert rotor.degeneracy(J) == 2 * J + 1


def test_transition_energies_evenly_spaced_by_2b():
    rotor = RigidRotor(moment_of_inertia=1.45e-46)
    B = rotor.rotational_constant
    for J in range(5):
        assert rotor.transition_energy(J) == pytest.approx(2.0 * B * (J + 1))


def test_from_diatomic_reduced_mass_and_moment_of_inertia():
    m1, m2, r = 1.008 * sc.atomic_mass, 34.97 * sc.atomic_mass, 127.5e-12
    rotor = RigidRotor.from_diatomic(m1, m2, r)
    mu = m1 * m2 / (m1 + m2)
    assert rotor.moment_of_inertia == pytest.approx(mu * r**2)


def test_rejects_nonpositive_moment_of_inertia():
    with pytest.raises(ValueError):
        RigidRotor(moment_of_inertia=0.0)


def test_rejects_negative_quantum_number():
    rotor = RigidRotor(moment_of_inertia=1.0e-46)
    with pytest.raises(ValueError):
        rotor.energy(-1)
    with pytest.raises(ValueError):
        rotor.degeneracy(-1)
    with pytest.raises(ValueError):
        rotor.transition_energy(-1)
