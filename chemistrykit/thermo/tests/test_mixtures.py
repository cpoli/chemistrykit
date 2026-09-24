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


def test_margules_with_zero_parameter_is_ideal():
    from chemistrykit.thermo.systems.mixtures import MargulesSolution

    x = np.linspace(0.0, 1.0, 11)
    ideal = BinaryIdealSolution(P_A_star=30.0, P_B_star=10.0)
    marg = MargulesSolution(A=0.0, P1_star=30.0, P2_star=10.0)
    np.testing.assert_allclose(marg.total_pressure(x), ideal.total_pressure(x))
    np.testing.assert_allclose(marg.vapor_composition(x[1:]), ideal.vapor_composition(x[1:]))


def test_margules_limiting_laws():
    """Raoult's law as x1 -> 1 and Henry's law K_H = P1* exp(A) as x1 -> 0."""
    from chemistrykit.thermo.systems.mixtures import MargulesSolution

    sol = MargulesSolution(A=1.5, P1_star=30.0, P2_star=10.0)
    P1_dilute, _ = sol.partial_pressures(1e-6)
    assert P1_dilute / 1e-6 == pytest.approx(sol.henry_constant_1, rel=1e-5)
    g1, _ = sol.activity_coefficients(1.0)
    assert g1 == pytest.approx(1.0)


def test_margules_satisfies_gibbs_duhem():
    """x1 dln(g1)/dx1 + x2 dln(g2)/dx1 = 0."""
    from chemistrykit.thermo.systems.mixtures import MargulesSolution

    sol = MargulesSolution(A=-0.8, P1_star=1.0, P2_star=1.0)
    x = np.linspace(0.05, 0.95, 19)
    h = 1e-6
    lg1p, lg2p = (np.log(g) for g in sol.activity_coefficients(x + h))
    lg1m, lg2m = (np.log(g) for g in sol.activity_coefficients(x - h))
    residual = x * (lg1p - lg1m) / (2 * h) + (1 - x) * (lg2p - lg2m) / (2 * h)
    np.testing.assert_allclose(residual, 0.0, atol=1e-6)


def test_margules_excess_gibbs_closed_form():
    from chemistrykit.thermo.systems.mixtures import MargulesSolution

    sol = MargulesSolution(A=2.0, P1_star=1.0, P2_star=1.0)
    assert sol.excess_gibbs(0.5, T=300.0) == pytest.approx(8.314462618 * 300.0 * 2.0 * 0.25, rel=1e-6)
