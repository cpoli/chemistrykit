"""Tests for chemistrykit.quantum.systems.particle_in_box against closed-form results."""

import numpy as np
import pytest
import scipy.constants as sc

from chemistrykit.quantum.systems.particle_in_box import ParticleInBox1D, ParticleInBox3D, conjugated_dye_absorption_wavelength

ELECTRON_MASS = sc.m_e
H = sc.h
C = sc.c


def test_1d_ground_state_energy_matches_closed_form():
    box = ParticleInBox1D(length=1.0e-9)
    expected = 1**2 * H**2 / (8.0 * ELECTRON_MASS * (1.0e-9) ** 2)
    assert box.energy(1) == pytest.approx(expected)


def test_1d_energy_scales_as_n_squared():
    box = ParticleInBox1D(length=1.0e-9)
    assert box.energy(2) == pytest.approx(4.0 * box.energy(1))
    assert box.energy(3) == pytest.approx(9.0 * box.energy(1))


def test_1d_energy_scales_as_inverse_length_squared():
    box1 = ParticleInBox1D(length=1.0e-9)
    box2 = ParticleInBox1D(length=2.0e-9)
    assert box1.energy(1) == pytest.approx(4.0 * box2.energy(1))


def test_1d_rejects_invalid_quantum_number():
    box = ParticleInBox1D(length=1.0e-9)
    with pytest.raises(ValueError):
        box.energy(0)


def test_1d_rejects_nonpositive_parameters():
    with pytest.raises(ValueError):
        ParticleInBox1D(length=0.0)
    with pytest.raises(ValueError):
        ParticleInBox1D(length=1.0, mass=-1.0)


def test_1d_wavefunction_normalized():
    box = ParticleInBox1D(length=1.0e-9)
    x = np.linspace(0.0, 1.0e-9, 200001)
    dx = x[1] - x[0]
    for n in (1, 2, 3):
        norm = np.trapezoid(box.probability_density(x, n), dx=dx)
        assert norm == pytest.approx(1.0, abs=1.0e-3)


def test_1d_wavefunction_has_n_minus_1_nodes():
    box = ParticleInBox1D(length=1.0e-9)
    x = np.linspace(1.0e-12, 1.0e-9 - 1.0e-12, 5000)
    for n in (1, 2, 3, 4):
        psi = box.wavefunction(x, n)
        sign_changes = np.sum(np.diff(np.sign(psi)) != 0)
        assert sign_changes == n - 1


def test_3d_cubic_box_degenerate_states_share_energy():
    box = ParticleInBox3D(Lx=1.0e-9, Ly=1.0e-9, Lz=1.0e-9)
    e211 = box.energy(2, 1, 1)
    e121 = box.energy(1, 2, 1)
    e112 = box.energy(1, 1, 2)
    assert e211 == pytest.approx(e121)
    assert e211 == pytest.approx(e112)


def test_3d_ground_state_energy_matches_closed_form():
    box = ParticleInBox3D(Lx=1.0e-9, Ly=2.0e-9, Lz=3.0e-9)
    expected = (H**2 / (8.0 * ELECTRON_MASS)) * (1.0 / 1.0e-9**2 + 1.0 / 2.0e-9**2 + 1.0 / 3.0e-9**2)
    assert box.energy(1, 1, 1) == pytest.approx(expected)


def test_3d_degeneracy_counts_cubic_box_first_excited_triplet():
    box = ParticleInBox3D(Lx=1.0e-9, Ly=1.0e-9, Lz=1.0e-9)
    degeneracies = box.degeneracy(n_max=2)
    assert 3 in degeneracies.values()


def test_3d_rejects_nonpositive_dimensions():
    with pytest.raises(ValueError):
        ParticleInBox3D(Lx=0.0, Ly=1.0, Lz=1.0)


def test_dye_absorption_wavelength_positive_and_in_reasonable_range():
    wavelength = conjugated_dye_absorption_wavelength(box_length=1.2e-9, n_pi_electrons=8)
    assert 1.0e-8 < wavelength < 1.0e-5


def test_dye_absorption_redshifts_with_more_conjugation():
    short = conjugated_dye_absorption_wavelength(box_length=0.8e-9, n_pi_electrons=6)
    long = conjugated_dye_absorption_wavelength(box_length=1.4e-9, n_pi_electrons=10)
    assert long > short


def test_dye_absorption_matches_direct_formula():
    box_length, n_pi = 1.0e-9, 8
    delta_e = (H**2 / (8.0 * ELECTRON_MASS * box_length**2)) * (n_pi + 1)
    expected = H * C / delta_e
    assert conjugated_dye_absorption_wavelength(box_length, n_pi) == pytest.approx(expected)


def test_dye_absorption_rejects_odd_electron_count():
    with pytest.raises(ValueError):
        conjugated_dye_absorption_wavelength(box_length=1.0e-9, n_pi_electrons=7)
