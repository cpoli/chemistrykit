"""Tests for chemistrykit.crystal.systems.defects: Boltzmann-factor scaling of defect concentration."""

import numpy as np
import pytest

from chemistrykit.crystal.systems.defects import frenkel_defect_concentration, schottky_defect_concentration


def test_schottky_concentration_increases_with_temperature():
    T = np.linspace(300.0, 1500.0, 20)
    n = schottky_defect_concentration(N=1e22, delta_h=2.0e-19, T=T)
    assert np.all(np.diff(n) > 0)


def test_schottky_concentration_matches_closed_form():
    from chemistrykit.constants import K_B

    N, delta_h, T = 1e22, 2.0e-19, 500.0
    expected = N * np.exp(-delta_h / (2.0 * K_B * T))
    assert schottky_defect_concentration(N, delta_h, T) == pytest.approx(expected)


def test_schottky_concentration_never_exceeds_N():
    T = np.linspace(300.0, 5000.0, 50)
    n = schottky_defect_concentration(N=1e22, delta_h=1.0e-19, T=T)
    assert np.all(n <= 1e22)


def test_frenkel_concentration_increases_with_temperature():
    T = np.linspace(300.0, 1500.0, 20)
    n = frenkel_defect_concentration(N=1e22, N_interstitial=1e21, delta_h=2.5e-19, T=T)
    assert np.all(np.diff(n) > 0)


def test_frenkel_reduces_to_schottky_form_with_equal_site_counts():
    N, delta_h, T = 1e22, 2.5e-19, 500.0
    n_frenkel = frenkel_defect_concentration(N=N, N_interstitial=N, delta_h=delta_h, T=T)
    n_schottky = schottky_defect_concentration(N=N, delta_h=delta_h, T=T)
    assert n_frenkel == pytest.approx(n_schottky)


def test_frenkel_scales_as_geometric_mean_of_site_counts():
    delta_h, T = 2.5e-19, 500.0
    n1 = frenkel_defect_concentration(N=1e22, N_interstitial=1e20, delta_h=delta_h, T=T)
    n2 = frenkel_defect_concentration(N=1e22, N_interstitial=4e20, delta_h=delta_h, T=T)
    # N_interstitial quadruples -> sqrt(N*N_i) doubles.
    assert n2 / n1 == pytest.approx(2.0)


def test_higher_formation_enthalpy_gives_fewer_defects_at_fixed_T():
    T = 800.0
    n_low = schottky_defect_concentration(N=1e22, delta_h=1.0e-19, T=T)
    n_high = schottky_defect_concentration(N=1e22, delta_h=3.0e-19, T=T)
    assert n_high < n_low
