"""Tests for the Wilson/NRTL/UNIQUAC activity models, the shared VLE base, and the Antoine equation."""

import numpy as np
import pytest

from chemistrykit.thermo.systems.mixtures import MargulesSolution, NRTLSolution, UNIQUACSolution, WilsonSolution
from chemistrykit.thermo.systems.phase_equilibria import AntoineEquation, ClausiusClapeyron

MODELS = [
    MargulesSolution(A=1.5, P1_star=30.0, P2_star=20.0),
    WilsonSolution(Lambda12=0.3, Lambda21=0.8, P1_star=30.0, P2_star=20.0),
    NRTLSolution(tau12=0.5, tau21=1.2, P1_star=30.0, P2_star=20.0, alpha=0.3),
    UNIQUACSolution(r1=2.1055, r2=0.92, q1=1.972, q2=1.40, tau12=0.8, tau21=0.6, P1_star=30.0, P2_star=20.0),
]


@pytest.mark.parametrize("model", MODELS, ids=lambda m: type(m).__name__)
def test_gibbs_duhem(model):
    x = np.linspace(0.05, 0.95, 19)
    h = 1e-6
    g1p, g2p = model.activity_coefficients(x + h)
    g1m, g2m = model.activity_coefficients(x - h)
    d1 = (np.log(g1p) - np.log(g1m)) / (2 * h)
    d2 = (np.log(g2p) - np.log(g2m)) / (2 * h)
    np.testing.assert_allclose(x * d1 + (1 - x) * d2, 0.0, atol=1e-7)


@pytest.mark.parametrize("model", MODELS, ids=lambda m: type(m).__name__)
def test_pure_component_limits(model):
    g1, _ = model.activity_coefficients(1.0)
    _, g2 = model.activity_coefficients(0.0)
    assert g1 == pytest.approx(1.0)
    assert g2 == pytest.approx(1.0)
    assert model.total_pressure(1.0) == pytest.approx(model.P1_star)
    assert model.total_pressure(0.0) == pytest.approx(model.P2_star)


@pytest.mark.parametrize("model", MODELS, ids=lambda m: type(m).__name__)
def test_azeotrope_has_equal_liquid_and_vapor_composition(model):
    x_az, P_az = model.azeotrope()
    assert model.vapor_composition(x_az) == pytest.approx(x_az, abs=1e-9)
    x = np.linspace(0.01, 0.99, 99)
    assert P_az >= model.total_pressure(x).max() - 1e-6  # positive deviation -> maximum-pressure azeotrope


def test_symmetric_margules_azeotrope_is_at_half():
    x_az, _ = MargulesSolution(A=2.0, P1_star=20.0, P2_star=20.0).azeotrope()
    assert x_az == pytest.approx(0.5)


def test_no_azeotrope_for_ideal_solution():
    assert WilsonSolution(1.0, 1.0, 30.0, 20.0).azeotrope() is None


def test_wilson_infinite_dilution():
    m = WilsonSolution(Lambda12=0.3, Lambda21=0.8, P1_star=1.0, P2_star=1.0)
    g1, g2 = m.activity_coefficients(0.0)
    assert np.log(g1) == pytest.approx(-np.log(0.3) + 1 - 0.8)
    g1, g2 = m.activity_coefficients(1.0)
    assert np.log(g2) == pytest.approx(-np.log(0.8) + 1 - 0.3)


def test_wilson_from_energies():
    m = WilsonSolution.from_energies(V1=58.7, V2=18.07, lambda12=1000.0, lambda21=2000.0, T=350.0, P1_star=1.0, P2_star=1.0)
    RT = 8.314462618 * 350.0
    assert m.Lambda12 == pytest.approx(18.07 / 58.7 * np.exp(-1000.0 / RT))
    assert m.Lambda21 == pytest.approx(58.7 / 18.07 * np.exp(-2000.0 / RT))


def test_nrtl_infinite_dilution():
    m = NRTLSolution(tau12=0.5, tau21=1.2, P1_star=1.0, P2_star=1.0, alpha=0.3)
    _, g2 = m.activity_coefficients(1.0)
    assert np.log(g2) == pytest.approx(0.5 + 1.2 * np.exp(-0.3 * 1.2))


def test_uniquac_athermal_equal_size_is_ideal():
    m = UNIQUACSolution(r1=2.0, r2=2.0, q1=1.5, q2=1.5, tau12=1.0, tau21=1.0, P1_star=1.0, P2_star=1.0)
    g1, g2 = m.activity_coefficients(np.linspace(0.1, 0.9, 5))
    np.testing.assert_allclose(g1, 1.0)
    np.testing.assert_allclose(g2, 1.0)


def test_generic_excess_gibbs_matches_margules_closed_form():
    m = MargulesSolution(A=1.5, P1_star=1.0, P2_star=1.0)
    x = np.linspace(0.1, 0.9, 5)
    np.testing.assert_allclose(m.excess_gibbs(x, 300.0), 8.314462618 * 300.0 * 1.5 * x * (1 - x), rtol=1e-9)


def test_antoine_round_trip_and_clausius_clapeyron_limit():
    water = AntoineEquation(A=5.08354, B=1663.125, C=-45.622)
    assert water.pressure(373.15) == pytest.approx(1.01325, rel=1e-3)
    T = np.array([350.0, 360.0, 370.0])
    np.testing.assert_allclose(water.temperature(water.pressure(T)), T)
    # C = 0 is exactly the integrated Clausius-Clapeyron equation, with dH = B R ln10
    antoine = AntoineEquation(A=5.0, B=2000.0, C=0.0)
    cc = ClausiusClapeyron(delta_h_vap=2000.0 * 8.314462618 * np.log(10), T_ref=350.0, P_ref=float(antoine.pressure(350.0)))
    np.testing.assert_allclose(antoine.pressure(T), cc.pressure(T), rtol=1e-6)
