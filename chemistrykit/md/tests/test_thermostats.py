"""Tests for chemistrykit.md.systems.thermostats."""

import pytest

from chemistrykit.md.systems.lj_fluid import LJFluid
from chemistrykit.md.systems.thermostats import NoseHooverThermostat, VelocityRescalingThermostat


def test_velocity_rescaling_thermostat_hits_target_temperature_immediately():
    fluid = LJFluid.from_lattice(n_per_side=3, cutoff=1.5, density=0.6, temperature=2.0, rng=0)
    fluid.velocities *= 3.0  # perturb well away from the target
    thermostat = VelocityRescalingThermostat(target_temperature=2.0)
    thermostat.apply(fluid, dt=0.001)
    assert fluid.temperature() == pytest.approx(2.0)


def test_velocity_rescaling_thermostat_respects_interval():
    fluid = LJFluid.from_lattice(n_per_side=3, cutoff=1.5, density=0.6, temperature=2.0, rng=1)
    fluid.velocities *= 3.0
    T_perturbed = fluid.temperature()
    thermostat = VelocityRescalingThermostat(target_temperature=2.0, interval=3)
    thermostat.apply(fluid, dt=0.001)
    thermostat.apply(fluid, dt=0.001)
    assert fluid.temperature() == pytest.approx(T_perturbed)  # no rescale yet
    thermostat.apply(fluid, dt=0.001)
    assert fluid.temperature() == pytest.approx(2.0)


def test_velocity_rescaling_thermostat_run_maintains_temperature():
    fluid = LJFluid.from_lattice(n_per_side=3, cutoff=1.5, density=0.6, temperature=1.5, rng=2)
    thermostat = VelocityRescalingThermostat(target_temperature=1.5, interval=5)
    result = fluid.run(dt=0.001, n_steps=200, thermostat=thermostat, sample_every=20)
    assert result.temperature[-1] == pytest.approx(1.5, abs=0.3)


def test_velocity_rescaling_thermostat_rejects_nonpositive_temperature():
    with pytest.raises(ValueError):
        VelocityRescalingThermostat(target_temperature=0.0)


def test_nose_hoover_thermostat_drives_temperature_toward_target():
    fluid = LJFluid.from_lattice(n_per_side=3, cutoff=1.5, density=0.6, temperature=2.0, rng=3)
    fluid.velocities *= 2.0
    initial_deviation = abs(fluid.temperature() - 2.0)
    thermostat = NoseHooverThermostat(target_temperature=2.0, Q=5.0)
    result = fluid.run(dt=0.001, n_steps=500, thermostat=thermostat, sample_every=50)
    final_deviation = abs(result.temperature[-1] - 2.0)
    assert final_deviation < initial_deviation


def test_nose_hoover_thermostat_rejects_nonpositive_parameters():
    with pytest.raises(ValueError):
        NoseHooverThermostat(target_temperature=0.0, Q=1.0)
    with pytest.raises(ValueError):
        NoseHooverThermostat(target_temperature=1.0, Q=0.0)
