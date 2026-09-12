"""Tests for chemistrykit.spectro.systems.rotational against closed-form results."""

import numpy as np
import pytest
import scipy.constants as sc

from chemistrykit.quantum.systems.rigid_rotor import RigidRotor
from chemistrykit.spectro.systems.rotational import (
    energy_to_wavenumber,
    isotope_shift_ratio,
    rotational_line_wavenumbers,
    rotational_spectrum,
)


def _hcl_rotor():
    return RigidRotor.from_diatomic(mass1=1.008 * sc.atomic_mass, mass2=34.97 * sc.atomic_mass, bond_length=127.5e-12)


def test_energy_to_wavenumber_closed_form():
    energy = sc.h * sc.c * 500.0 * 100.0  # 500 cm^-1
    assert energy_to_wavenumber(energy) == pytest.approx(500.0)


def test_rotational_lines_evenly_spaced_by_2b():
    rotor = _hcl_rotor()
    lines = rotational_line_wavenumbers(rotor, J_max=6)
    spacings = np.diff(lines)
    assert spacings == pytest.approx(spacings[0] * np.ones_like(spacings))


def test_rotational_lines_match_direct_energy_conversion():
    rotor = _hcl_rotor()
    lines = rotational_line_wavenumbers(rotor, J_max=3)
    for J in range(4):
        expected = energy_to_wavenumber(rotor.transition_energy(J))
        assert lines[J] == pytest.approx(expected)


def test_rotational_spectrum_intensities_are_nonnegative():
    rotor = _hcl_rotor()
    spectrum = rotational_spectrum(rotor, J_max=10, temperature=300.0)
    assert np.all(spectrum.intensities >= 0.0)
    assert len(spectrum.positions) == 11


def test_rotational_spectrum_higher_temperature_populates_higher_j():
    rotor = _hcl_rotor()
    cold = rotational_spectrum(rotor, J_max=20, temperature=50.0)
    hot = rotational_spectrum(rotor, J_max=20, temperature=500.0)
    peak_cold = int(np.argmax(cold.intensities))
    peak_hot = int(np.argmax(hot.intensities))
    assert peak_hot > peak_cold


def test_isotope_shift_ratio_less_than_one_for_heavier_isotope():
    ratio = isotope_shift_ratio(mass1a=1.008 * sc.atomic_mass, mass2=34.97 * sc.atomic_mass, mass1b=2.014 * sc.atomic_mass)
    assert 0.0 < ratio < 1.0


def test_isotope_shift_ratio_matches_direct_rotor_comparison():
    r = 127.5e-12
    rotor_hcl = RigidRotor.from_diatomic(1.008 * sc.atomic_mass, 34.97 * sc.atomic_mass, r)
    rotor_dcl = RigidRotor.from_diatomic(2.014 * sc.atomic_mass, 34.97 * sc.atomic_mass, r)
    ratio_direct = rotor_dcl.rotational_constant / rotor_hcl.rotational_constant
    ratio_formula = isotope_shift_ratio(1.008 * sc.atomic_mass, 34.97 * sc.atomic_mass, 2.014 * sc.atomic_mass)
    assert ratio_direct == pytest.approx(ratio_formula)


def test_isotope_shift_ratio_is_one_for_identical_masses():
    ratio = isotope_shift_ratio(mass1a=10.0, mass2=20.0, mass1b=10.0)
    assert ratio == pytest.approx(1.0)
