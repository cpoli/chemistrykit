"""Tests for chemistrykit.thermo.systems.phase_equilibria against closed-form results."""

import numpy as np
import pytest

from chemistrykit.thermo.systems.phase_equilibria import ClausiusClapeyron, gibbs_phase_rule


def test_clausius_clapeyron_reference_point_is_exact():
    model = ClausiusClapeyron(delta_h_vap=40700.0, T_ref=373.15, P_ref=101325.0)
    assert model.pressure(373.15) == pytest.approx(101325.0)


def test_clausius_clapeyron_pressure_decreases_below_reference_temperature():
    model = ClausiusClapeyron(delta_h_vap=40700.0, T_ref=373.15, P_ref=101325.0)
    assert model.pressure(350.0) < model.pressure(373.15)


def test_clausius_clapeyron_boiling_point_inverts_pressure():
    model = ClausiusClapeyron(delta_h_vap=35000.0, T_ref=300.0, P_ref=1.0e4)
    T = 320.0
    P = model.pressure(T)
    assert model.boiling_point(P) == pytest.approx(T, rel=1e-8)


def test_clausius_clapeyron_from_two_points_recovers_delta_h():
    delta_h_true = 42000.0
    model = ClausiusClapeyron(delta_h_vap=delta_h_true, T_ref=290.0, P_ref=2.0e4)
    T2 = 310.0
    P2 = model.pressure(T2)
    fit = ClausiusClapeyron.from_two_points(290.0, 2.0e4, T2, P2)
    assert fit.delta_h_vap == pytest.approx(delta_h_true, rel=1e-8)


def test_clausius_clapeyron_vectorized_temperatures():
    model = ClausiusClapeyron(delta_h_vap=40700.0, T_ref=373.15, P_ref=101325.0)
    T = np.array([350.0, 373.15, 400.0])
    P = model.pressure(T)
    assert np.all(np.diff(P) > 0)  # pressure increases monotonically with T


def test_gibbs_phase_rule_pure_substance_triple_point():
    assert gibbs_phase_rule(n_components=1, n_phases=3) == 0


def test_gibbs_phase_rule_pure_substance_two_phases():
    assert gibbs_phase_rule(n_components=1, n_phases=2) == 1


def test_gibbs_phase_rule_pure_substance_single_phase():
    assert gibbs_phase_rule(n_components=1, n_phases=1) == 2


def test_gibbs_phase_rule_with_reaction_constraint():
    """Each independent reaction equilibrium removes one degree of freedom."""
    assert gibbs_phase_rule(n_components=3, n_phases=1, reactions=1) == 3
