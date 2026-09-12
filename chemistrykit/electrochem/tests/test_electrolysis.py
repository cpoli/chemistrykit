"""Tests for chemistrykit.electrochem.systems.electrolysis against closed-form/known results."""

import pytest

from chemistrykit.electrochem.systems.electrolysis import (
    charge_from_current,
    faradays_law_mass,
    mass_from_charge,
    minimum_applied_voltage_electrolytic,
    moles_from_charge,
)


def test_charge_from_current():
    assert charge_from_current(current=5.0, time=100.0) == pytest.approx(500.0)


def test_mass_proportional_to_charge():
    """Faraday's first law: mass deposited is directly proportional to charge passed."""
    m1 = mass_from_charge(charge=1000.0, molar_mass=63.546, n=2)
    m2 = mass_from_charge(charge=2000.0, molar_mass=63.546, n=2)
    assert m2 == pytest.approx(2.0 * m1)


def test_moles_from_charge_proportional_to_charge():
    m1 = moles_from_charge(charge=96485.332, n=1)
    m2 = moles_from_charge(charge=2 * 96485.332, n=1)
    assert m2 == pytest.approx(2.0 * m1)
    assert m1 == pytest.approx(1.0, rel=1e-6)


def test_faradays_law_mass_matches_known_silver_deposition():
    """1 F (96,485 C) deposits exactly 1 mole (107.87 g) of Ag via Ag+ + e- -> Ag."""
    m = faradays_law_mass(current=96485.332, time=1.0, molar_mass=107.87, n=1)
    assert m == pytest.approx(107.87, rel=1e-5)


def test_faradays_law_mass_scales_with_current_and_time():
    base = faradays_law_mass(current=2.0, time=1000.0, molar_mass=58.69, n=2)
    double_current = faradays_law_mass(current=4.0, time=1000.0, molar_mass=58.69, n=2)
    double_time = faradays_law_mass(current=2.0, time=2000.0, molar_mass=58.69, n=2)
    assert double_current == pytest.approx(2.0 * base)
    assert double_time == pytest.approx(2.0 * base)


def test_faradays_law_mass_efficiency_scaling():
    ideal = faradays_law_mass(10.0, 3600.0, 58.69, n=2, efficiency=1.0)
    reduced = faradays_law_mass(10.0, 3600.0, 58.69, n=2, efficiency=0.8)
    assert reduced == pytest.approx(0.8 * ideal)


def test_faradays_law_mass_rejects_bad_efficiency():
    with pytest.raises(ValueError):
        faradays_law_mass(1.0, 1.0, 10.0, n=1, efficiency=0.0)
    with pytest.raises(ValueError):
        faradays_law_mass(1.0, 1.0, 10.0, n=1, efficiency=1.5)


def test_minimum_applied_voltage_is_magnitude_of_negative_cell_potential():
    assert minimum_applied_voltage_electrolytic(-2.5) == pytest.approx(2.5)
    assert minimum_applied_voltage_electrolytic(2.5) == pytest.approx(2.5)
