"""Tests for chemistrykit.electrochem.systems.battery against closed-form results."""

import numpy as np
import pytest

from chemistrykit.electrochem.systems.battery import ConstantCurrentBattery, effective_capacity, peukert_discharge_time


def test_ideal_peukert_discharge_time_inversely_proportional_to_current():
    t1 = peukert_discharge_time(capacity_peukert=10.0, current=1.0, k=1.0)
    t2 = peukert_discharge_time(capacity_peukert=10.0, current=2.0, k=1.0)
    assert t1 == pytest.approx(2.0 * t2)


def test_ideal_effective_capacity_independent_of_current():
    """At k=1, effective delivered capacity should be exactly the Peukert
    capacity constant, regardless of discharge current."""
    for current in (0.5, 1.0, 5.0, 20.0):
        C_eff = effective_capacity(capacity_peukert=10.0, current=current, k=1.0)
        assert C_eff == pytest.approx(10.0)


def test_effective_capacity_decreases_with_rate_when_k_greater_than_one():
    currents = np.array([0.5, 1.0, 2.0, 5.0, 10.0])
    C_eff = effective_capacity(capacity_peukert=10.0, current=currents, k=1.2)
    assert np.all(np.diff(C_eff) < 0.0)


def test_effective_capacity_reduces_to_peukert_time_times_current():
    C_p, I, k = 20.0, 3.0, 1.15
    t = peukert_discharge_time(C_p, I, k)
    C_eff = effective_capacity(C_p, I, k)
    assert C_eff == pytest.approx(I * t)


def test_constant_current_battery_discharge_time_matches_peukert_formula():
    battery = ConstantCurrentBattery(capacity_peukert=5.0, current=2.5, v_nominal=3.7, k=1.1)
    assert battery.discharge_time() == pytest.approx(peukert_discharge_time(5.0, 2.5, 1.1))


def test_constant_current_battery_voltage_sags_with_internal_resistance():
    battery = ConstantCurrentBattery(capacity_peukert=10.0, current=1.0, v_nominal=3.7, internal_resistance=0.1)
    v = battery.terminal_voltage(0.5)
    assert v == pytest.approx(3.7 - 1.0 * 0.1)


def test_constant_current_battery_voltage_drops_to_zero_after_cutoff():
    battery = ConstantCurrentBattery(capacity_peukert=2.0, current=1.0, v_nominal=3.7)
    assert battery.terminal_voltage(1.0) > 0.0
    assert battery.terminal_voltage(3.0) == pytest.approx(0.0)


def test_constant_current_battery_state_of_charge_bounds():
    battery = ConstantCurrentBattery(capacity_peukert=4.0, current=1.0, v_nominal=3.7)
    t = np.linspace(0.0, 8.0, 50)
    soc = battery.state_of_charge(t)
    assert np.all(soc >= 0.0) and np.all(soc <= 1.0)
    assert soc[0] == pytest.approx(1.0)
    assert soc[-1] == pytest.approx(0.0)


def test_discharge_curve_result_shapes():
    battery = ConstantCurrentBattery(capacity_peukert=4.0, current=1.0, v_nominal=3.7)
    t = np.linspace(0.0, 5.0, 20)
    result = battery.discharge_curve(t)
    assert result.t.shape == (20,)
    assert result.voltage.shape == (20,)
    assert result.state_of_charge.shape == (20,)
