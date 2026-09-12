"""Tests for chemistrykit.analytical.systems.titration: redox and EDTA titration curves."""

import numpy as np
import pytest

from chemistrykit.analytical.systems.titration import EDTATitration, RedoxTitration


def test_redox_titration_equivalence_potential_is_weighted_average():
    titration = RedoxTitration(E1_standard=0.771, n1=1, E2_standard=1.72, n2=1, C_analyte=0.10, V_analyte=0.050, C_titrant=0.10)
    V_eq = titration.equivalence_volume()
    E_eq = titration.response_at(np.array([V_eq]))[0]
    assert E_eq == pytest.approx((0.771 + 1.72) / 2.0, abs=1e-6)


def test_redox_titration_equivalence_potential_general_n1_n2_weighted_average():
    titration = RedoxTitration(E1_standard=0.5, n1=1, E2_standard=1.0, n2=2, C_analyte=0.10, V_analyte=0.050, C_titrant=0.10)
    V_eq = titration.equivalence_volume()
    E_eq = titration.response_at(np.array([V_eq]))[0]
    expected = (1 * 0.5 + 2 * 1.0) / (1 + 2)
    assert E_eq == pytest.approx(expected, abs=1e-6)


def test_redox_titration_curve_is_monotonic_increasing():
    titration = RedoxTitration(E1_standard=0.771, n1=1, E2_standard=1.72, n2=1, C_analyte=0.10, V_analyte=0.050, C_titrant=0.10)
    V = np.linspace(1e-6, 0.09, 3000)
    curve = titration.curve(V)
    assert np.all(np.diff(curve.response) > 0)


def test_redox_titration_equivalence_point_detection_matches_stoichiometric():
    titration = RedoxTitration(E1_standard=0.5, n1=1, E2_standard=1.0, n2=2, C_analyte=0.10, V_analyte=0.050, C_titrant=0.10)
    V = np.linspace(1e-6, 0.09, 20000)
    V_eq_numeric = titration.find_equivalence_point(V)
    assert V_eq_numeric == pytest.approx(titration.equivalence_volume(), abs=1e-3)


def test_edta_titration_equivalence_volume_is_stoichiometric():
    titration = EDTATitration(C_metal=0.010, V_metal=0.050, K_conditional=1e10, C_edta=0.020)
    assert titration.equivalence_volume() == pytest.approx(0.025)


def test_edta_titration_pM_curve_is_monotonic_increasing():
    titration = EDTATitration(C_metal=0.010, V_metal=0.050, K_conditional=1e10, C_edta=0.010)
    V = np.linspace(1e-6, 0.09, 3000)
    curve = titration.curve(V)
    assert np.all(np.diff(curve.response) > 0)


def test_edta_titration_equivalence_point_approaches_large_K_approximation():
    K = 1e14
    titration = EDTATitration(C_metal=0.010, V_metal=0.050, K_conditional=K, C_edta=0.010)
    V_eq = titration.equivalence_volume()
    pM_exact = titration.response_at(np.array([V_eq]))[0]
    C_M_eq = titration.C_metal * titration.V_metal / (titration.V_metal + V_eq)
    pM_approx = 0.5 * np.log10(K / C_M_eq)
    assert pM_exact == pytest.approx(pM_approx, abs=1e-3)


def test_edta_titration_higher_K_gives_sharper_endpoint():
    V = np.linspace(1e-6, 0.09, 20000)
    low_K = EDTATitration(C_metal=0.010, V_metal=0.050, K_conditional=1e6, C_edta=0.010)
    high_K = EDTATitration(C_metal=0.010, V_metal=0.050, K_conditional=1e12, C_edta=0.010)
    low_curve = low_K.curve(V)
    high_curve = high_K.curve(V)
    low_jump = np.max(np.abs(np.gradient(low_curve.response, low_curve.V)))
    high_jump = np.max(np.abs(np.gradient(high_curve.response, high_curve.V)))
    assert high_jump > low_jump
