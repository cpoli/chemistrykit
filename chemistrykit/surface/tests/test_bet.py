"""Tests for chemistrykit.surface.systems.bet against closed-form/known results.

The Langmuir-limit test is the important cross-check flagged in the
module docstring: BET must reduce numerically to the Langmuir isotherm
as P0 -> infinity at fixed K = C/P0.
"""

import numpy as np
import pytest

from chemistrykit.surface.systems.bet import BETIsotherm, bet_loading, fit_bet
from chemistrykit.surface.systems.langmuir import langmuir_coverage


def test_bet_loading_is_zero_at_zero_pressure():
    assert bet_loading(Vm=2.0, C=50.0, P=0.0, P0=10.0) == pytest.approx(0.0)


def test_bet_reduces_to_langmuir_as_P0_to_infinity():
    K = 2.0
    P0 = 1.0e7
    C = K * P0
    Vm = 1.0
    P = np.linspace(0.01, 5.0, 50)
    V_bet = bet_loading(Vm, C, P, P0)
    theta_langmuir = langmuir_coverage(K, P)
    np.testing.assert_allclose(V_bet, Vm * theta_langmuir, rtol=1e-4)


def test_bet_langmuir_limit_improves_as_P0_grows():
    K = 2.0
    Vm = 1.0
    P = np.linspace(0.01, 5.0, 50)
    theta_langmuir = langmuir_coverage(K, P)
    errors = []
    for P0 in (1.0e3, 1.0e5, 1.0e7):
        C = K * P0
        V_bet = bet_loading(Vm, C, P, P0)
        errors.append(np.max(np.abs(V_bet - Vm * theta_langmuir)))
    assert errors[0] > errors[1] > errors[2]


def test_bet_isotherm_matches_function():
    iso = BETIsotherm(Vm=3.0, C=60.0, P0=15.0)
    P = np.array([1.0, 2.0, 3.0])
    np.testing.assert_allclose(iso.loading(P), bet_loading(3.0, 60.0, P, 15.0))


def test_bet_isotherm_rejects_nonpositive_parameters():
    with pytest.raises(ValueError):
        BETIsotherm(Vm=-1.0, C=1.0, P0=1.0)
    with pytest.raises(ValueError):
        BETIsotherm(Vm=1.0, C=0.0, P0=1.0)
    with pytest.raises(ValueError):
        BETIsotherm(Vm=1.0, C=1.0, P0=0.0)


def test_fit_bet_recovers_known_parameters_from_exact_data():
    Vm_true, C_true, P0 = 5.0, 80.0, 10.0
    P = np.array([0.5, 1.0, 1.5, 2.0, 2.5, 3.0])
    V = bet_loading(Vm_true, C_true, P, P0)
    fit = fit_bet(P, V, P0=P0)
    assert fit.Vm == pytest.approx(Vm_true, rel=1e-4)
    assert fit.C == pytest.approx(C_true, rel=1e-4)
    assert fit.r_squared == pytest.approx(1.0, abs=1e-9)
