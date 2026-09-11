"""Tests for chemistrykit.quantum.systems.harmonic_oscillator against closed-form results."""

import numpy as np
import pytest

from chemistrykit.constants import HBAR
from chemistrykit.quantum.systems.harmonic_oscillator import MorseOscillator, QuantumHarmonicOscillator, compare_harmonic_vs_morse


def test_zero_point_energy_matches_closed_form():
    ho = QuantumHarmonicOscillator(mass=1.6e-27, force_constant=500.0)
    assert ho.zero_point_energy == pytest.approx(0.5 * HBAR * ho.angular_frequency)
    assert ho.energy(0) == pytest.approx(ho.zero_point_energy)


def test_level_spacing_is_exactly_hbar_omega():
    ho = QuantumHarmonicOscillator(mass=1.6e-27, force_constant=500.0)
    for v in range(10):
        spacing = ho.energy(v + 1) - ho.energy(v)
        assert spacing == pytest.approx(HBAR * ho.angular_frequency)


def test_angular_frequency_matches_sqrt_k_over_m():
    ho = QuantumHarmonicOscillator(mass=2.0, force_constant=8.0)
    assert ho.angular_frequency == pytest.approx(2.0)


def test_rejects_nonpositive_parameters():
    with pytest.raises(ValueError):
        QuantumHarmonicOscillator(mass=0.0, force_constant=1.0)
    with pytest.raises(ValueError):
        QuantumHarmonicOscillator(mass=1.0, force_constant=-1.0)


def test_wavefunction_ground_state_is_gaussian_and_normalized():
    ho = QuantumHarmonicOscillator(mass=1.6e-27, force_constant=500.0)
    beta = np.sqrt(ho.mass * ho.angular_frequency / HBAR)
    x = np.linspace(-10.0 / beta, 10.0 / beta, 200001)
    dx = x[1] - x[0]
    psi0 = ho.wavefunction(x, v=0)
    norm = np.trapezoid(psi0**2, dx=dx)
    assert norm == pytest.approx(1.0, abs=1.0e-3)


def test_wavefunctions_are_orthogonal():
    ho = QuantumHarmonicOscillator(mass=1.6e-27, force_constant=500.0)
    beta = np.sqrt(ho.mass * ho.angular_frequency / HBAR)
    x = np.linspace(-12.0 / beta, 12.0 / beta, 300001)
    dx = x[1] - x[0]
    psi0 = ho.wavefunction(x, v=0)
    psi1 = ho.wavefunction(x, v=1)
    overlap = np.trapezoid(psi0 * psi1, dx=dx)
    assert overlap == pytest.approx(0.0, abs=1.0e-6)


def test_morse_ground_state_close_to_harmonic_ground_state():
    mass, k, De = 1.6e-27, 500.0, 7.0e-19
    ho = QuantumHarmonicOscillator(mass, k)
    morse = MorseOscillator(mass, k, De)
    # Small relative difference at v=0 since anharmonicity is a small correction there.
    assert morse.energy(0) == pytest.approx(ho.energy(0), rel=0.05)


def test_morse_level_spacing_shrinks_with_v():
    morse = MorseOscillator(mass=1.6e-27, force_constant=500.0, dissociation_energy=7.0e-19)
    spacings = [morse.energy(v + 1) - morse.energy(v) for v in range(morse.v_max)]
    for i in range(len(spacings) - 1):
        assert spacings[i] > spacings[i + 1]


def test_morse_v_max_brackets_the_level_spacing_sign_change():
    morse = MorseOscillator(mass=1.6e-27, force_constant=500.0, dissociation_energy=7.0e-19)
    # v_max marks (to within rounding) where the (unphysical, beyond-dissociation)
    # quadratic-in-v formula peaks and spacing turns negative: comfortably
    # below v_max the levels are still climbing, comfortably above they've
    # turned over.
    spacing_below = morse.energy(morse.v_max) - morse.energy(morse.v_max - 1)
    spacing_above = morse.energy(morse.v_max + 2) - morse.energy(morse.v_max + 1)
    assert spacing_below > 0
    assert spacing_above < 0


def test_morse_rejects_nonpositive_parameters():
    with pytest.raises(ValueError):
        MorseOscillator(mass=1.0, force_constant=1.0, dissociation_energy=0.0)


def test_compare_harmonic_vs_morse_agreement_worsens_with_v():
    v, E_harmonic, E_morse = compare_harmonic_vs_morse(mass=1.6e-27, force_constant=500.0, dissociation_energy=7.0e-19, v_max=6)
    relative_diff = np.abs(E_harmonic - E_morse) / E_harmonic
    assert relative_diff[-1] > relative_diff[0]
    assert np.all(np.diff(relative_diff) >= -1.0e-12)  # monotonically non-decreasing


def test_compare_harmonic_vs_morse_shapes():
    v, E_harmonic, E_morse = compare_harmonic_vs_morse(mass=1.6e-27, force_constant=500.0, dissociation_energy=7.0e-19, v_max=4)
    assert v.shape == (5,)
    assert E_harmonic.shape == (5,)
    assert E_morse.shape == (5,)
