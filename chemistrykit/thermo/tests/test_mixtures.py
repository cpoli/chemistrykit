"""Tests for chemistrykit.thermo.systems.mixtures against closed-form results."""

import numpy as np
import pytest

from chemistrykit.thermo.systems.mixtures import (
    BinaryIdealSolution,
    boiling_point_elevation,
    freezing_point_depression,
    henry_law_pressure,
    osmotic_pressure,
    raoult_vapor_pressure,
)


def test_raoult_vapor_pressure_pure_component_limit():
    assert raoult_vapor_pressure(x=1.0, P_pure=500.0) == pytest.approx(500.0)


def test_raoult_vapor_pressure_scales_linearly():
    assert raoult_vapor_pressure(x=0.3, P_pure=200.0) == pytest.approx(60.0)


def test_binary_ideal_solution_equal_volatility_gives_equal_composition():
    solution = BinaryIdealSolution(P_A_star=100.0, P_B_star=100.0)
    x_A = np.linspace(0.0, 1.0, 11)
    np.testing.assert_allclose(solution.vapor_composition(x_A), x_A)


def test_binary_ideal_solution_enriches_vapor_in_more_volatile_component():
    solution = BinaryIdealSolution(P_A_star=300.0, P_B_star=50.0)
    x_A = 0.5
    assert solution.vapor_composition(x_A) > x_A


def test_binary_ideal_solution_total_pressure_at_pure_ends():
    solution = BinaryIdealSolution(P_A_star=300.0, P_B_star=50.0)
    assert solution.total_pressure(1.0) == pytest.approx(300.0)
    assert solution.total_pressure(0.0) == pytest.approx(50.0)


def test_henry_law_pressure_scales_linearly():
    assert henry_law_pressure(x=0.002, K_H=1.0e5) == pytest.approx(200.0)


def test_freezing_point_depression_van_t_hoff_factor_scaling():
    """Doubling the van't Hoff factor doubles the depression, at fixed molality."""
    dT1 = freezing_point_depression(Kf=1.86, b=0.5, i=1.0)
    dT2 = freezing_point_depression(Kf=1.86, b=0.5, i=2.0)
    assert dT2 == pytest.approx(2.0 * dT1)


def test_boiling_point_elevation_matches_formula():
    assert boiling_point_elevation(Kb=0.51, b=2.0, i=1.0) == pytest.approx(1.02)


def test_osmotic_pressure_matches_ideal_gas_law_analogy():
    """Pi = i*M*R*T -- structurally identical to the ideal gas law."""
    from chemistrykit.constants import R

    M, T = 50.0, 300.0
    assert osmotic_pressure(M=M, T=T, i=1.0) == pytest.approx(R * T * M)


def test_osmotic_pressure_scales_linearly_with_concentration():
    p1 = osmotic_pressure(M=10.0, T=298.15)
    p2 = osmotic_pressure(M=20.0, T=298.15)
    assert p2 == pytest.approx(2.0 * p1)
