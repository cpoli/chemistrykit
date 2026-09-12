"""Tests for chemistrykit.polymer.systems.molecular_weight_distribution against closed-form/known results.

Includes a numerical cross-check: summing the Flory-Schulz number/weight
fractions directly (via chemistrykit.polymer.utils.moments) must
reproduce the closed-form Xn = 1/(1-p), Xw = (1+p)/(1-p) results.
"""

import numpy as np
import pytest

from chemistrykit.polymer.systems.molecular_weight_distribution import (
    MolecularWeightDistribution,
    flory_schulz_number_average_DP,
    flory_schulz_number_fraction,
    flory_schulz_pdi,
    flory_schulz_weight_average_DP,
    flory_schulz_weight_fraction,
    number_average_molar_mass,
    polydispersity_index,
    weight_average_molar_mass,
)
from chemistrykit.polymer.utils.moments import number_average, weight_average


def test_number_average_molar_mass_known_value():
    assert number_average_molar_mass(N_i=[10.0, 5.0], M_i=[1000.0, 2000.0]) == pytest.approx(1333.333333, rel=1e-6)


def test_weight_average_molar_mass_known_value():
    assert weight_average_molar_mass(N_i=[10.0, 5.0], M_i=[1000.0, 2000.0]) == pytest.approx(1500.0, rel=1e-9)


def test_weight_average_at_least_number_average():
    N_i, M_i = [10.0, 5.0], [1000.0, 2000.0]
    Mn = number_average_molar_mass(N_i, M_i)
    Mw = weight_average_molar_mass(N_i, M_i)
    assert Mw >= Mn


def test_polydispersity_index_is_one_for_monodisperse_sample():
    assert polydispersity_index(Mn=1000.0, Mw=1000.0) == pytest.approx(1.0)


def test_molecular_weight_distribution_from_counts():
    mwd = MolecularWeightDistribution.from_counts(N_i=[10.0, 5.0], M_i=[1000.0, 2000.0])
    assert mwd.Mn == pytest.approx(1333.333333, rel=1e-6)
    assert mwd.Mw == pytest.approx(1500.0, rel=1e-9)
    assert mwd.PDI == pytest.approx(1.125, rel=1e-9)


def test_flory_schulz_number_fraction_sums_to_one():
    x = np.arange(1, 5000)
    total = np.sum(flory_schulz_number_fraction(x, p=0.98))
    assert total == pytest.approx(1.0, abs=1e-6)


def test_flory_schulz_weight_fraction_sums_to_one():
    x = np.arange(1, 5000)
    total = np.sum(flory_schulz_weight_fraction(x, p=0.98))
    assert total == pytest.approx(1.0, abs=1e-6)


@pytest.mark.parametrize("p", [0.5, 0.9, 0.95, 0.98])
def test_flory_schulz_number_average_DP_matches_direct_summation(p):
    x = np.arange(1, 20000)
    N_x = flory_schulz_number_fraction(x, p)
    Xn_numeric = number_average(x, N_x)
    assert Xn_numeric == pytest.approx(flory_schulz_number_average_DP(p), rel=1e-3)


@pytest.mark.parametrize("p", [0.5, 0.9, 0.95, 0.98])
def test_flory_schulz_weight_average_DP_matches_direct_summation(p):
    x = np.arange(1, 20000)
    N_x = flory_schulz_number_fraction(x, p)
    Xw_numeric = weight_average(x, N_x)
    assert Xw_numeric == pytest.approx(flory_schulz_weight_average_DP(p), rel=1e-3)


def test_flory_schulz_pdi_is_ratio_of_averages():
    p = 0.85
    assert flory_schulz_pdi(p) == pytest.approx(flory_schulz_weight_average_DP(p) / flory_schulz_number_average_DP(p))


def test_flory_schulz_pdi_is_one_with_no_reaction():
    assert flory_schulz_pdi(0.0) == pytest.approx(1.0)


def test_flory_schulz_pdi_approaches_exactly_two_at_full_conversion():
    assert flory_schulz_pdi(1.0) == pytest.approx(2.0)
    # And monotonically approaches it as p -> 1.
    p_values = np.array([0.9, 0.99, 0.999, 0.9999])
    pdi_values = flory_schulz_pdi(p_values)
    assert np.all(np.diff(pdi_values) > 0)
    assert np.all(pdi_values < 2.0)
