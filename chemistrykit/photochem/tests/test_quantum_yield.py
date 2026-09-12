"""Tests for chemistrykit.photochem.systems.quantum_yield against closed-form/known results."""

import numpy as np
import pytest

from chemistrykit.photochem.systems.quantum_yield import (
    fluorescence_quantum_yield,
    intersystem_crossing_yield,
    phosphorescence_quantum_yield,
    photochemical_quantum_yield,
    photons_absorbed,
)


@pytest.mark.parametrize(
    "kf,kic,kisc",
    [(1.0, 0.0, 0.0), (2.0, 1.0, 0.5), (0.1, 5.0, 3.0), (10.0, 0.0, 10.0)],
)
def test_fluorescence_yield_bounded_in_unit_interval(kf, kic, kisc):
    phi = fluorescence_quantum_yield(kf, kic, kisc)
    assert 0.0 <= phi <= 1.0


def test_fluorescence_yield_is_one_with_no_competing_pathway():
    assert fluorescence_quantum_yield(kf=1.0, kic=0.0, kisc=0.0) == pytest.approx(1.0)


def test_decay_channel_yields_sum_to_one():
    kf, kic, kisc = 2.0, 1.0, 0.5
    phi_f = fluorescence_quantum_yield(kf, kic, kisc)
    phi_isc = intersystem_crossing_yield(kisc, kf, kic)
    phi_ic = kic / (kf + kic + kisc)
    assert phi_f + phi_isc + phi_ic == pytest.approx(1.0)


def test_phosphorescence_yield_bounded_and_requires_intersystem_crossing():
    phi_p = phosphorescence_quantum_yield(kisc=0.0, kf=1.0, kic=0.5, kp=1.0, kic_T=0.5)
    assert phi_p == pytest.approx(0.0)

    phi_p_full = phosphorescence_quantum_yield(kisc=1.0, kf=0.0, kic=0.0, kp=1.0, kic_T=0.0)
    assert phi_p_full == pytest.approx(1.0)

    phi_p_partial = phosphorescence_quantum_yield(kisc=1.0, kf=1.0, kic=0.0, kp=1.0, kic_T=1.0)
    assert 0.0 <= phi_p_partial <= 1.0


def test_photons_absorbed_matches_beer_lambert_transmittance():
    from chemistrykit.spectro.systems.beer_lambert import transmittance

    I0, A = 5.0, 0.3
    expected = I0 * (1.0 - transmittance(A))
    assert photons_absorbed(I0, A) == pytest.approx(expected)


def test_photons_absorbed_bounds():
    assert photons_absorbed(2.0, absorbance=0.0) == pytest.approx(0.0)
    assert photons_absorbed(2.0, absorbance=20.0) == pytest.approx(2.0, rel=1e-9)


def test_photons_absorbed_array_input():
    I0 = np.array([1.0, 2.0, 3.0])
    result = photons_absorbed(I0, absorbance=0.5)
    assert result.shape == (3,)
    np.testing.assert_allclose(result / I0, result[0] / I0[0])


def test_photochemical_quantum_yield_simple_ratio():
    assert photochemical_quantum_yield(2.0, 4.0) == pytest.approx(0.5)
    assert photochemical_quantum_yield(1.0, 1.0) == pytest.approx(1.0)
