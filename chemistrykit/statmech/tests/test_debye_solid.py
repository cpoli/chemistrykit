"""Tests for chemistrykit.statmech.systems.debye_solid against closed-form limits."""

import numpy as np
import pytest

from chemistrykit.constants import K_B, R
from chemistrykit.statmech.systems.debye_solid import DebyeSolid


def test_high_temperature_limit_is_dulong_petit():
    solid = DebyeSolid(debye_temperature=300.0)
    assert solid.heat_capacity_v(1.0e5) == pytest.approx(3.0 * R, rel=1e-6)


def test_low_temperature_limit_is_t_cubed_law():
    solid = DebyeSolid(debye_temperature=300.0)
    T = np.array([0.5, 1.0, 2.0])
    assert solid.heat_capacity_v(T) == pytest.approx(solid.low_temperature_heat_capacity(T), rel=1e-9)


def test_heat_capacity_is_derivative_of_internal_energy():
    solid = DebyeSolid(debye_temperature=343.0)
    T, dT = 150.0, 1e-3
    numeric = (solid.internal_energy(T + dT, N=1.0) - solid.internal_energy(T - dT, N=1.0)) / (2 * dT)
    assert numeric == pytest.approx(solid.heat_capacity_v(T, N=1.0), rel=1e-6)


def test_high_temperature_energy_is_classical_minus_zero_point():
    # U - U0 -> 3NkT - (9/8) N k Theta_D as T -> infinity
    solid = DebyeSolid(debye_temperature=100.0)
    T = 1.0e5
    assert solid.internal_energy(T, N=1.0) == pytest.approx(3 * K_B * T - 9.0 / 8.0 * K_B * 100.0, rel=1e-8)


def test_rejects_nonpositive_inputs():
    with pytest.raises(ValueError):
        DebyeSolid(0.0)
    with pytest.raises(ValueError):
        DebyeSolid(100.0).heat_capacity_v(0.0)
