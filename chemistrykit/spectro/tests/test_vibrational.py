"""Tests for chemistrykit.spectro.systems.vibrational against closed-form results and physical sanity checks."""

import numpy as np
import pytest
import scipy.constants as sc

from chemistrykit.quantum.systems.harmonic_oscillator import MorseOscillator
from chemistrykit.spectro.systems.vibrational import (
    TriatomicNormalModes,
    anharmonicity_from_overtones,
    harmonic_fundamental_wavenumber,
    morse_transition_wavenumbers,
)


def test_harmonic_fundamental_matches_closed_form():
    k, mu = 480.0, (1.008 * 34.97) / (1.008 + 34.97) * sc.atomic_mass
    nu = harmonic_fundamental_wavenumber(k, mu)
    expected = (1.0 / (2.0 * np.pi * sc.c)) * np.sqrt(k / mu) / 100.0
    assert nu == pytest.approx(expected)


def test_morse_overtone_below_twice_fundamental():
    morse = MorseOscillator(mass=1.6e-27, force_constant=500.0, dissociation_energy=7.0e-19)
    nu01, nu02 = morse_transition_wavenumbers(morse, v_max=2)
    assert nu02 < 2.0 * nu01


def test_morse_transitions_decrease_in_spacing():
    morse = MorseOscillator(mass=1.6e-27, force_constant=500.0, dissociation_energy=7.0e-19)
    wavenumbers = morse_transition_wavenumbers(morse, v_max=4)
    spacings = np.diff(np.concatenate([[0.0], wavenumbers]))
    assert np.all(np.diff(spacings) < 0)


def test_anharmonicity_from_overtones_round_trip():
    morse = MorseOscillator(mass=1.6e-27, force_constant=500.0, dissociation_energy=7.0e-19)
    nu01, nu02 = morse_transition_wavenumbers(morse, v_max=2)
    omega_e, omega_e_xe = anharmonicity_from_overtones(nu01, nu02)
    reconstructed_01 = omega_e - 2.0 * omega_e_xe
    reconstructed_02 = 2.0 * omega_e - 6.0 * omega_e_xe
    assert reconstructed_01 == pytest.approx(nu01)
    assert reconstructed_02 == pytest.approx(nu02)


def test_anharmonicity_constant_is_positive_for_real_anharmonic_potential():
    morse = MorseOscillator(mass=1.6e-27, force_constant=500.0, dissociation_energy=7.0e-19)
    nu01, nu02 = morse_transition_wavenumbers(morse, v_max=2)
    _, omega_e_xe = anharmonicity_from_overtones(nu01, nu02)
    assert omega_e_xe > 0.0


def _water_normal_modes(k_r=770.0, k_theta=0.7e-18):
    return TriatomicNormalModes.bent(
        mass_terminal=1.008 * sc.atomic_mass,
        mass_central=15.999 * sc.atomic_mass,
        bond_length=95.8e-12,
        bond_angle_degrees=104.5,
        k_r=k_r,
        k_theta=k_theta,
    )


def _co2_normal_modes(k_r=1600.0, k_theta=0.85e-18):
    return TriatomicNormalModes.linear(
        mass_terminal=15.999 * sc.atomic_mass,
        mass_central=12.011 * sc.atomic_mass,
        bond_length=116.3e-12,
        k_r=k_r,
        k_theta=k_theta,
    )


def test_water_normal_modes_all_positive_and_correctly_ordered():
    result = _water_normal_modes().solve()
    assert result.is_linear is False
    assert np.all(result.wavenumbers > 0.0)
    # bend well below both stretches, in the right ballpark for a real
    # bent triatomic (H2O's real modes are 1595, 3657, 3756 cm^-1).
    assert result.wavenumbers[0] < result.wavenumbers[1]
    assert 500.0 < result.wavenumbers[0] < 3000.0
    assert 2000.0 < result.wavenumbers[1] < 5000.0
    assert 2000.0 < result.wavenumbers[2] < 5000.0


def test_co2_normal_modes_bend_lower_than_stretches():
    result = _co2_normal_modes().solve()
    assert result.is_linear is True
    assert np.all(result.wavenumbers > 0.0)
    bend, sym_stretch, asym_stretch = result.wavenumbers
    assert bend < sym_stretch < asym_stretch
    # Real CO2: bend 667, symmetric stretch 1388, antisymmetric stretch 2349 cm^-1;
    # the simple diagonal valence force field should land in the right ballpark.
    assert 300.0 < bend < 1200.0
    assert 900.0 < sym_stretch < 2000.0
    assert 1800.0 < asym_stretch < 3200.0


def test_stiffer_bend_force_constant_raises_bend_frequency():
    soft = _water_normal_modes(k_theta=0.3e-18).solve()
    stiff = _water_normal_modes(k_theta=1.2e-18).solve()
    assert stiff.wavenumbers[0] > soft.wavenumbers[0]


def test_stiffer_stretch_force_constant_raises_stretch_frequencies():
    soft = _water_normal_modes(k_r=400.0).solve()
    stiff = _water_normal_modes(k_r=1200.0).solve()
    assert stiff.wavenumbers[1] > soft.wavenumbers[1]
    assert stiff.wavenumbers[2] > soft.wavenumbers[2]


def test_triatomic_normal_modes_rejects_bad_shapes():
    with pytest.raises(ValueError):
        TriatomicNormalModes(masses=[1.0, 1.0], equilibrium_coordinates=np.zeros((3, 3)), k_r1=1.0, k_r2=1.0, k_theta=1.0)
    with pytest.raises(ValueError):
        TriatomicNormalModes(masses=[1.0, 1.0, 1.0], equilibrium_coordinates=np.zeros((2, 3)), k_r1=1.0, k_r2=1.0, k_theta=1.0)


def test_triatomic_normal_modes_rejects_nonpositive_force_constants():
    coords = np.array([[0.0, 0.0, -1.0], [0.0, 0.0, 0.0], [0.0, 0.0, 1.0]])
    with pytest.raises(ValueError):
        TriatomicNormalModes(masses=[1.0, 1.0, 1.0], equilibrium_coordinates=coords, k_r1=-1.0, k_r2=1.0, k_theta=1.0)
