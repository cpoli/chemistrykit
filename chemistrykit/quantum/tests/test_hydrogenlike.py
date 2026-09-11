"""Tests for chemistrykit.quantum.systems.hydrogenlike against closed-form/known results."""

import numpy as np
import pytest

from chemistrykit.constants import ELECTRONVOLT
from chemistrykit.quantum.systems.hydrogenlike import HydrogenLikeAtom


def test_hydrogen_ground_state_is_minus_13p6_ev():
    h_atom = HydrogenLikeAtom(Z=1)
    energy_ev = -h_atom.energy(1) / ELECTRONVOLT
    assert energy_ev == pytest.approx(13.6, abs=0.02)


def test_energy_scales_as_z_squared():
    h_atom = HydrogenLikeAtom(Z=1)
    he_plus = HydrogenLikeAtom(Z=2)
    assert he_plus.energy(1) == pytest.approx(4.0 * h_atom.energy(1))


def test_energy_scales_as_inverse_n_squared():
    h_atom = HydrogenLikeAtom(Z=1)
    assert h_atom.energy(2) == pytest.approx(h_atom.energy(1) / 4.0)
    assert h_atom.energy(3) == pytest.approx(h_atom.energy(1) / 9.0)


def test_energy_is_negative_and_increases_toward_zero_with_n():
    h_atom = HydrogenLikeAtom(Z=1)
    e1, e2, e3 = h_atom.energy(1), h_atom.energy(2), h_atom.energy(3)
    assert e1 < e2 < e3 < 0.0


def test_bohr_radius_matches_known_value():
    h_atom = HydrogenLikeAtom(Z=1)
    assert h_atom.bohr_radius == pytest.approx(5.29177e-11, rel=1.0e-4)


def test_reduced_mass_shifts_ground_state_slightly():
    import scipy.constants as sc

    mu = sc.m_e * sc.m_p / (sc.m_e + sc.m_p)
    h_infinite_mass = HydrogenLikeAtom(Z=1, reduced_mass=sc.m_e)
    h_reduced_mass = HydrogenLikeAtom(Z=1, reduced_mass=mu)
    # The reduced-mass correction lowers the binding energy magnitude very slightly.
    assert abs(h_reduced_mass.energy(1)) < abs(h_infinite_mass.energy(1))
    assert abs(h_reduced_mass.energy(1)) == pytest.approx(abs(h_infinite_mass.energy(1)), rel=1.0e-2)


def test_radial_wavefunction_1s_matches_closed_form():
    h_atom = HydrogenLikeAtom(Z=1)
    a0 = h_atom.bohr_radius
    r = np.array([0.0, a0, 2.0 * a0])
    expected = 2.0 * a0**-1.5 * np.exp(-r / a0)
    np.testing.assert_allclose(h_atom.radial_wavefunction(r, n=1, l=0), expected, rtol=1.0e-10)


@pytest.mark.parametrize("n,l", [(1, 0), (2, 0), (2, 1), (3, 0), (3, 1), (3, 2)])
def test_radial_wavefunction_is_normalized(n, l):
    h_atom = HydrogenLikeAtom(Z=1)
    assert h_atom.check_radial_normalization(n, l) == pytest.approx(1.0, abs=1.0e-4)


def test_radial_wavefunction_rejects_invalid_l():
    h_atom = HydrogenLikeAtom(Z=1)
    with pytest.raises(ValueError):
        h_atom.radial_wavefunction(1.0e-10, n=1, l=1)
    with pytest.raises(ValueError):
        h_atom.radial_wavefunction(1.0e-10, n=2, l=-1)


def test_1s_radial_distribution_peaks_at_bohr_radius():
    h_atom = HydrogenLikeAtom(Z=1)
    a0 = h_atom.bohr_radius
    r = np.linspace(1.0e-4 * a0, 6.0 * a0, 200000)
    P = h_atom.radial_distribution_function(r, n=1, l=0)
    r_peak = r[np.argmax(P)]
    assert r_peak == pytest.approx(a0, rel=1.0e-3)


def test_angular_wavefunction_s_orbital_is_isotropic():
    h_atom = HydrogenLikeAtom(Z=1)
    Y00_a = h_atom.angular_wavefunction(theta=0.3, phi=0.1, l=0, m=0)
    Y00_b = h_atom.angular_wavefunction(theta=2.1, phi=4.0, l=0, m=0)
    assert Y00_a == pytest.approx(Y00_b)


def test_rejects_invalid_z_and_mass():
    with pytest.raises(ValueError):
        HydrogenLikeAtom(Z=0)
    with pytest.raises(ValueError):
        HydrogenLikeAtom(Z=1, reduced_mass=-1.0)


def test_energy_rejects_invalid_n():
    h_atom = HydrogenLikeAtom(Z=1)
    with pytest.raises(ValueError):
        h_atom.energy(0)
