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


def test_berendsen_relaxes_temperature_geometrically_without_dynamics():
    from chemistrykit.md.systems.thermostats import BerendsenThermostat

    fluid = LJFluid.from_lattice(n_per_side=3, cutoff=1.5, density=0.6, temperature=2.0, rng=3)
    thermostat = BerendsenThermostat(target_temperature=1.0, tau=0.1)
    for _ in range(10):
        thermostat.apply(fluid, dt=0.01)
    # T_n = T0 + (T_i - T0) (1 - dt/tau)^n
    assert fluid.temperature() == pytest.approx(1.0 + 1.0 * 0.9**10, rel=1e-12)


def test_berendsen_with_tau_equal_dt_is_full_rescaling():
    from chemistrykit.md.systems.thermostats import BerendsenThermostat

    fluid = LJFluid.from_lattice(n_per_side=3, cutoff=1.5, density=0.6, temperature=2.5, rng=4)
    BerendsenThermostat(target_temperature=1.2, tau=0.005).apply(fluid, dt=0.005)
    assert fluid.temperature() == pytest.approx(1.2)


def test_stochastic_velocity_rescaling_samples_canonical_kinetic_energy():
    import numpy as np

    from chemistrykit.md.systems.thermostats import StochasticVelocityRescalingThermostat

    fluid = LJFluid.from_lattice(n_per_side=3, cutoff=1.5, density=0.6, temperature=2.0, rng=5)
    T0 = 1.3
    thermostat = StochasticVelocityRescalingThermostat(target_temperature=T0, tau=0.05, rng=0)
    samples = []
    for step in range(40000):
        thermostat.apply(fluid, dt=0.01)
        if step >= 1000:
            samples.append(fluid.kinetic_energy())
    samples = np.array(samples)
    nf = fluid.degrees_of_freedom()
    k_bar = 0.5 * nf * T0
    assert samples.mean() == pytest.approx(k_bar, rel=0.02)
    assert samples.var() == pytest.approx(2.0 * k_bar**2 / nf, rel=0.1)


@pytest.mark.parametrize("cls_name", ["BerendsenThermostat", "StochasticVelocityRescalingThermostat"])
def test_new_thermostats_reject_nonpositive_tau(cls_name):
    from chemistrykit.md.systems import thermostats

    with pytest.raises(ValueError):
        getattr(thermostats, cls_name)(target_temperature=1.0, tau=0.0)
