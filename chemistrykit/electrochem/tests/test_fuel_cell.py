"""Tests for chemistrykit.electrochem.systems.fuel_cell against closed-form/known results."""

import pytest

from chemistrykit.electrochem.systems.fuel_cell import (
    fuel_cell_efficiency_limit,
    reversible_cell_voltage,
    reversible_cell_voltage_at_temperature,
)
from chemistrykit.electrochem.systems.standard_potentials import standard_cell_potential


def test_hydrogen_oxygen_voltage_matches_standard_potentials():
    E = reversible_cell_voltage(-237.13e3, n=2)
    assert E == pytest.approx(standard_cell_potential("O2/H2O", "H+/H2"), abs=5e-3)


def test_temperature_form_reduces_to_gibbs_form():
    dH, dS, T = -285.83e3, -163.3, 298.15
    assert reversible_cell_voltage_at_temperature(dH, dS, 2, T) == pytest.approx(reversible_cell_voltage(dH - T * dS, 2))


def test_efficiency_limit():
    assert fuel_cell_efficiency_limit(-237.13e3, -285.83e3) == pytest.approx(0.8296, abs=1e-4)
