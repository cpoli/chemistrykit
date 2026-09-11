"""Tests for chemistrykit.md.systems.pair_potentials against closed-form/known results."""

import numpy as np
import pytest

from chemistrykit.md.systems.pair_potentials import (
    Buckingham,
    DiatomicOscillator,
    HarmonicAngle,
    HarmonicBond,
    HarmonicMolecule,
    Morse,
)


def test_morse_minimum_at_re_with_depth_de_and_zero_force():
    morse = Morse(De=5.0, a=2.0, re=1.2)
    assert float(morse.energy(morse.re)) == pytest.approx(-5.0)
    assert float(morse.force_scalar(morse.re)) == pytest.approx(0.0, abs=1e-10)


def test_morse_dissociates_to_zero_at_large_r():
    morse = Morse(De=5.0, a=2.0, re=1.2)
    assert float(morse.energy(50.0)) == pytest.approx(0.0, abs=1e-6)


def test_morse_force_constant_matches_curvature():
    """k = 2*De*a^2 should match the numerical second derivative of U at re."""
    morse = Morse(De=5.0, a=2.0, re=1.2)
    h = 1e-5
    r = morse.re
    second_derivative = (float(morse.energy(r + h)) - 2.0 * float(morse.energy(r)) + float(morse.energy(r - h))) / h**2
    assert morse.force_constant == pytest.approx(second_derivative, rel=1e-3)


def test_morse_from_force_constant_round_trips():
    morse = Morse.from_force_constant(De=5.0, k=12.0, re=1.0)
    assert morse.force_constant == pytest.approx(12.0)


def test_morse_rejects_nonpositive_de_or_a():
    with pytest.raises(ValueError):
        Morse(De=0.0, a=1.0, re=1.0)
    with pytest.raises(ValueError):
        Morse(De=1.0, a=-1.0, re=1.0)


def test_buckingham_attractive_tail_at_large_r():
    pot = Buckingham(A=1.0e4, B=3.0, C=1.0)
    assert float(pot.energy(10.0)) < 0.0
    assert float(pot.force_scalar(10.0)) < 0.0


def test_buckingham_repulsive_at_short_range():
    """Repulsive in the physically-relevant short-range window -- not the same
    as r -> 0, where the -C/r^6 term's unphysical divergence eventually
    dominates and makes the potential (incorrectly) attractive again -- see
    the class docstring's note on the Buckingham potential's "catastrophe"."""
    pot = Buckingham(A=1.0e4, B=3.0, C=1.0)
    assert float(pot.force_scalar(0.6)) > 0.0


def test_buckingham_rejects_invalid_parameters():
    with pytest.raises(ValueError):
        Buckingham(A=-1.0, B=1.0, C=1.0)
    with pytest.raises(ValueError):
        Buckingham(A=1.0, B=1.0, C=-1.0)


def test_harmonic_bond_zero_at_equilibrium():
    bond = HarmonicBond(k=50.0, r0=1.0)
    assert float(bond.energy(1.0)) == 0.0
    assert float(bond.force_scalar(1.0)) == 0.0


def test_harmonic_bond_restoring_force_direction():
    bond = HarmonicBond(k=50.0, r0=1.0)
    assert float(bond.force_scalar(1.1)) < 0.0  # stretched -> attractive (pulls back in)
    assert float(bond.force_scalar(0.9)) > 0.0  # compressed -> repulsive (pushes back out)


def test_harmonic_angle_zero_at_equilibrium():
    angle = HarmonicAngle(k_theta=100.0, theta0=1.9106)
    assert float(angle.energy(angle.theta0)) == 0.0
    assert float(angle.torque(angle.theta0)) == pytest.approx(0.0)


def test_diatomic_oscillator_harmonic_period():
    k, m1, m2, r0 = 200.0, 1.0, 1.0, 1.2
    bond = HarmonicBond(k=k, r0=1.0)
    oscillator = DiatomicOscillator(bond, m1=m1, m2=m2, r0=r0)
    period = DiatomicOscillator.harmonic_period(k, m1, m2)
    dt = period / 2000.0
    n_steps = 2000  # one full period
    initial_length = oscillator.bond_length()
    for _ in range(n_steps):
        oscillator.step(dt)
    final_length = oscillator.bond_length()
    assert final_length == pytest.approx(initial_length, abs=0.01)


def test_diatomic_oscillator_conserves_energy():
    bond = HarmonicBond(k=200.0, r0=1.0)
    oscillator = DiatomicOscillator(bond, m1=1.0, m2=2.0, r0=1.3)
    e0 = oscillator.kinetic_energy() + oscillator.potential_energy()
    period = DiatomicOscillator.harmonic_period(200.0, 1.0, 2.0)
    dt = period / 1000.0
    for _ in range(1000):
        oscillator.step(dt)
    e1 = oscillator.kinetic_energy() + oscillator.potential_energy()
    assert e1 == pytest.approx(e0, rel=1e-3)


def test_diatomic_oscillator_morse_dissociates_with_enough_energy():
    morse = Morse(De=1.0, a=2.0, re=1.0)
    oscillator = DiatomicOscillator(morse, m1=1.0, m2=1.0, r0=1.0, v_rel0=3.0)
    for _ in range(2000):
        oscillator.step(0.001)
    assert oscillator.bond_length() > 3.0 * morse.re


def test_diatomic_oscillator_rejects_unsupported_potential():
    class DummyPotential:
        pass

    with pytest.raises(TypeError):
        DiatomicOscillator(DummyPotential(), m1=1.0, m2=1.0, r0=1.0)


def test_harmonic_molecule_bond_only_conserves_energy():
    positions = np.array([[0.0, 0.0, 0.0], [1.2, 0.0, 0.0]])
    velocities = np.zeros((2, 3))
    masses = [1.0, 1.0]
    molecule = HarmonicMolecule(positions, velocities, masses, bonds=[(0, 1, 200.0, 1.0)])
    e0 = molecule.kinetic_energy() + molecule.potential_energy()
    for _ in range(2000):
        molecule.step(0.001)
    e1 = molecule.kinetic_energy() + molecule.potential_energy()
    assert e1 == pytest.approx(e0, rel=1e-3)


def test_harmonic_molecule_bent_triatomic_conserves_energy():
    """A water-like bent triatomic with harmonic bonds + a harmonic angle term."""
    theta0 = 1.8239  # ~104.5 degrees
    r0 = 0.9584
    positions = np.array(
        [
            [r0 * np.cos(theta0 / 2), r0 * np.sin(theta0 / 2), 0.0],
            [0.0, 0.0, 0.0],
            [r0 * np.cos(theta0 / 2), -r0 * np.sin(theta0 / 2), 0.0],
        ]
    )
    positions[0, 0] += 0.05  # perturb away from equilibrium so it actually vibrates
    velocities = np.zeros((3, 3))
    masses = [1.008, 16.0, 1.008]
    molecule = HarmonicMolecule(
        positions,
        velocities,
        masses,
        bonds=[(0, 1, 450.0, r0), (1, 2, 450.0, r0)],
        angles=[(0, 1, 2, 50.0, theta0)],
    )
    e0 = molecule.kinetic_energy() + molecule.potential_energy()
    for _ in range(3000):
        molecule.step(0.001)
    e1 = molecule.kinetic_energy() + molecule.potential_energy()
    assert e1 == pytest.approx(e0, rel=1e-2)


def test_harmonic_molecule_angle_only_bounded_by_energy_conservation():
    """Started from rest, energy conservation bounds |theta(t)-theta0| by its initial value at every t.

    With no bond terms, the angle-bending force is exactly tangential (it
    conserves each vertex-to-atom distance -- see the module docstring's
    force derivation), so the angle potential ``0.5*k_theta*(theta-theta0)**2``
    is the *only* potential energy in this system. Starting at rest, all
    energy is initially potential, so ``PE(t) <= PE(0)`` for all t (since
    ``KE(t) >= 0``), which -- because PE is monotonic in the deviation --
    means the angle can never wander farther from theta0 than its started.
    """
    theta0 = np.pi / 2.0
    positions = np.array([[1.0, 0.3, 0.0], [0.0, 0.0, 0.0], [0.0, 1.0, 0.0]])
    velocities = np.zeros((3, 3))
    masses = [1.0, 1.0, 1.0]
    molecule = HarmonicMolecule(positions, velocities, masses, angles=[(0, 1, 2, 30.0, theta0)])

    def angle_now():
        v1 = molecule.positions[0] - molecule.positions[1]
        v2 = molecule.positions[2] - molecule.positions[1]
        cos_t = np.dot(v1, v2) / (np.linalg.norm(v1) * np.linalg.norm(v2))
        return np.arccos(np.clip(cos_t, -1.0, 1.0))

    initial_deviation = abs(angle_now() - theta0)
    max_deviation = initial_deviation
    for _ in range(200):
        molecule.step(0.005)
        max_deviation = max(max_deviation, abs(angle_now() - theta0))
    assert max_deviation <= 1.05 * initial_deviation + 1e-3
