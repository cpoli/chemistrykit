"""Tests for chemistrykit.surface.systems.tpd against closed-form results."""

import numpy as np
import pytest

from chemistrykit.constants import R
from chemistrykit.surface.systems.tpd import first_order_peak_temperature, redhead_desorption_energy, simulate_tpd


@pytest.mark.parametrize("Ed,nu,beta", [(80e3, 1e13, 5.0), (120e3, 1e13, 10.0), (100e3, 1e11, 2.0)])
def test_simulated_first_order_peak_matches_exact_condition(Ed, nu, beta):
    res = simulate_tpd(Ed, nu, beta, T_start=100.0, T_end=900.0, n_points=40001)
    assert res.peak_temperature == pytest.approx(first_order_peak_temperature(Ed, nu, beta), abs=0.05)


def test_peak_condition_is_satisfied():
    Ed, nu, beta = 90e3, 1e13, 3.0
    Tp = first_order_peak_temperature(Ed, nu, beta)
    assert Ed / (R * Tp**2) == pytest.approx(nu / beta * np.exp(-Ed / (R * Tp)), rel=1e-9)


def test_first_order_peak_independent_of_initial_coverage():
    a = simulate_tpd(100e3, 1e13, 10.0, theta0=1.0).peak_temperature
    b = simulate_tpd(100e3, 1e13, 10.0, theta0=0.1).peak_temperature
    assert a == pytest.approx(b, abs=0.05)


def test_second_order_peak_shifts_down_with_coverage():
    high = simulate_tpd(100e3, 1e13, 10.0, theta0=1.0, order=2).peak_temperature
    low = simulate_tpd(100e3, 1e13, 10.0, theta0=0.1, order=2).peak_temperature
    assert high < low


def test_coverage_desorbs_completely_and_integrates_to_theta0():
    res = simulate_tpd(100e3, 1e13, 10.0, theta0=0.6)
    assert res.coverage[-1] == pytest.approx(0.0, abs=1e-9)
    total = np.sum(0.5 * (res.desorption_rate[1:] + res.desorption_rate[:-1]) * np.diff(res.T))
    assert total == pytest.approx(0.6, rel=1e-4)


def test_redhead_recovers_Ed_for_nu_1e13():
    for Ed in [60e3, 100e3, 150e3]:
        Tp = first_order_peak_temperature(Ed, 1e13, 10.0)
        assert redhead_desorption_energy(Tp, 1e13, 10.0) == pytest.approx(Ed, rel=0.015)


def test_invalid_order_raises():
    with pytest.raises(ValueError):
        simulate_tpd(100e3, 1e13, 10.0, order=3)
