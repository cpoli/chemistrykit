"""Tests for chemistrykit.spectro.systems.electronic against closed-form Franck-Condon results."""

import numpy as np
import pytest

from chemistrykit.spectro.systems.electronic import (
    franck_condon_factor,
    franck_condon_progression,
    franck_condon_spectrum,
    huang_rhys_factor,
)


def test_zero_displacement_gives_pure_0_0_transition():
    assert franck_condon_factor(0, S=0.0) == pytest.approx(1.0)
    assert franck_condon_factor(1, S=0.0) == pytest.approx(0.0)
    assert franck_condon_factor(5, S=0.0) == pytest.approx(0.0)


def test_huang_rhys_factor_zero_at_zero_displacement():
    assert huang_rhys_factor(displacement=0.0, mass=1.6e-27, angular_frequency=5.0e13) == pytest.approx(0.0)


def test_huang_rhys_factor_scales_with_displacement_squared():
    s1 = huang_rhys_factor(displacement=1e-11, mass=1.6e-27, angular_frequency=5.0e13)
    s2 = huang_rhys_factor(displacement=2e-11, mass=1.6e-27, angular_frequency=5.0e13)
    assert s2 == pytest.approx(4.0 * s1)


@pytest.mark.parametrize("S", [0.5, 1.5, 3.0, 8.0])
def test_franck_condon_progression_sums_to_one(S):
    progression = franck_condon_progression(v_max=80, S=S)
    assert np.sum(progression) == pytest.approx(1.0, abs=1e-9)


def test_franck_condon_progression_matches_poisson_pmf():
    from scipy.stats import poisson

    S = 2.5
    progression = franck_condon_progression(v_max=15, S=S)
    expected = poisson.pmf(np.arange(0, 16), mu=S)
    assert progression == pytest.approx(expected)


def test_franck_condon_peak_near_huang_rhys_parameter():
    for S in (1.0, 3.0, 6.0):
        progression = franck_condon_progression(v_max=30, S=S)
        peak_v = int(np.argmax(progression))
        assert abs(peak_v - S) <= 1


def test_franck_condon_factor_rejects_negative_v():
    with pytest.raises(ValueError):
        franck_condon_factor(-1, S=1.0)


def test_franck_condon_factor_rejects_negative_s():
    with pytest.raises(ValueError):
        franck_condon_factor(0, S=-1.0)


def test_franck_condon_spectrum_positions_evenly_spaced():
    spectrum = franck_condon_spectrum(origin_wavenumber=20000.0, vibrational_wavenumber=1200.0, S=2.0, v_max=8)
    spacings = np.diff(spectrum.positions)
    assert spacings == pytest.approx(1200.0 * np.ones_like(spacings))


def test_franck_condon_spectrum_intensities_match_progression():
    spectrum = franck_condon_spectrum(origin_wavenumber=20000.0, vibrational_wavenumber=1200.0, S=2.0, v_max=8)
    expected = franck_condon_progression(v_max=8, S=2.0)
    assert spectrum.intensities == pytest.approx(expected)
