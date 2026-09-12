"""Tests for chemistrykit.surface.systems.catalysis against closed-form/known results."""

import numpy as np
import pytest

from chemistrykit.constants import R
from chemistrykit.surface.systems.catalysis import compare_catalyzed_rate, turnover_frequency, turnover_number


def test_turnover_frequency_is_rate_over_site_concentration():
    assert turnover_frequency(rate=5.0e-3, active_site_concentration=2.0e-6) == pytest.approx(2500.0)


def test_turnover_number_is_moles_converted_over_moles_catalyst():
    assert turnover_number(moles_converted=0.5, moles_catalyst=1.0e-4) == pytest.approx(5000.0)


def test_compare_catalyzed_rate_enhancement_matches_exponential_of_delta_Ea():
    Ea_uncat, Ea_cat, T = 80e3, 50e3, 298.15
    comparison = compare_catalyzed_rate(Ea_uncatalyzed=Ea_uncat, Ea_catalyzed=Ea_cat, T=T, A_uncatalyzed=1e13)
    expected = np.exp((Ea_uncat - Ea_cat) / (R * T))
    assert comparison.rate_enhancement == pytest.approx(expected, rel=1e-10)
    assert comparison.delta_Ea == pytest.approx(Ea_uncat - Ea_cat)


def test_compare_catalyzed_rate_defaults_A_catalyzed_to_A_uncatalyzed():
    c1 = compare_catalyzed_rate(Ea_uncatalyzed=80e3, Ea_catalyzed=50e3, T=300.0, A_uncatalyzed=1e12)
    c2 = compare_catalyzed_rate(Ea_uncatalyzed=80e3, Ea_catalyzed=50e3, T=300.0, A_uncatalyzed=1e12, A_catalyzed=1e12)
    assert c1.rate_enhancement == pytest.approx(c2.rate_enhancement)


def test_compare_catalyzed_rate_enhancement_exceeds_one_when_Ea_is_reduced():
    comparison = compare_catalyzed_rate(Ea_uncatalyzed=100e3, Ea_catalyzed=60e3, T=298.15, A_uncatalyzed=1e13)
    assert comparison.rate_enhancement > 1.0
    assert comparison.k_catalyzed > comparison.k_uncatalyzed
