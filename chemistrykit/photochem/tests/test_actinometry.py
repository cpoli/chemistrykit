"""Tests for chemistrykit.photochem.systems.actinometry."""

import pytest

from chemistrykit.photochem.systems.actinometry import ferrioxalate_fe2_moles, ferrioxalate_photon_flux


def test_fe2_moles_beer_lambert():
    assert ferrioxalate_fe2_moles(0.555, volume_L=0.02, path_length_cm=0.5) == pytest.approx(0.555 * 0.02 / (1.11e4 * 0.5))


def test_photon_flux_roundtrip_and_partial_absorption():
    q_true, t, A, phi = 3.0e-9, 60.0, 0.5, 1.21
    n = q_true * t * phi * (1.0 - 10**-A)
    assert ferrioxalate_photon_flux(n, t, absorbance=A, quantum_yield=phi) == pytest.approx(q_true)
    assert ferrioxalate_photon_flux(1.21e-6, 100.0) == pytest.approx(1e-8)
