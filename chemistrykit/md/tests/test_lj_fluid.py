"""Tests for chemistrykit.md.systems.lj_fluid against closed-form/known results."""

import numpy as np
import pytest

from chemistrykit.md.systems.lj_fluid import LennardJones, LJFluid


def test_lj_minimum_is_at_two_to_one_sixth_sigma():
    """The textbook LJ result: minimum at r=2^(1/6)*sigma, U=-epsilon, F=0 there."""
    lj = LennardJones(epsilon=2.0, sigma=1.5)
    r_min = 2.0 ** (1.0 / 6.0) * 1.5
    assert lj.r_min == pytest.approx(r_min)
    assert float(lj.energy(r_min)) == pytest.approx(-2.0)
    assert float(lj.force_scalar(r_min)) == pytest.approx(0.0, abs=1e-10)


def test_lj_zero_at_sigma():
    lj = LennardJones(epsilon=1.0, sigma=1.0)
    assert float(lj.energy(1.0)) == pytest.approx(0.0, abs=1e-12)


def test_lj_repulsive_inside_minimum_attractive_outside():
    lj = LennardJones()
    assert float(lj.force_scalar(0.9 * lj.r_min)) > 0.0  # repulsive
    assert float(lj.force_scalar(2.0 * lj.r_min)) < 0.0  # attractive


def test_lj_rejects_nonpositive_parameters():
    with pytest.raises(ValueError):
        LennardJones(epsilon=0.0, sigma=1.0)
    with pytest.raises(ValueError):
        LennardJones(epsilon=1.0, sigma=-1.0)


def test_lj_fluid_conserves_energy_in_nve():
    fluid = LJFluid.from_lattice(n_per_side=4, density=0.5, temperature=1.0, rng=1)
    result = fluid.run(dt=0.001, n_steps=500, sample_every=50)
    total = result.total_energy
    drift = np.max(np.abs(total - total[0])) / abs(total[0])
    assert drift < 0.02


def test_lj_fluid_from_lattice_matches_target_temperature():
    fluid = LJFluid.from_lattice(n_per_side=3, cutoff=1.5, density=0.7, temperature=1.5, rng=2)
    assert fluid.temperature() == pytest.approx(1.5)


def test_lj_fluid_zero_net_momentum_conserved_during_run():
    fluid = LJFluid.from_lattice(n_per_side=3, cutoff=1.5, density=0.6, temperature=1.0, rng=3)
    total_momentum0 = np.sum(fluid.masses[:, None] * fluid.velocities, axis=0)
    fluid.run(dt=0.001, n_steps=200, sample_every=50)
    total_momentum1 = np.sum(fluid.masses[:, None] * fluid.velocities, axis=0)
    np.testing.assert_allclose(total_momentum0, 0.0, atol=1e-10)
    np.testing.assert_allclose(total_momentum1, 0.0, atol=1e-8)


def test_lj_fluid_rejects_box_smaller_than_twice_cutoff():
    positions = np.zeros((2, 3))
    velocities = np.zeros((2, 3))
    with pytest.raises(ValueError):
        LJFluid(positions, velocities, box_length=1.0, cutoff=2.5)


def test_lj_fluid_pressure_approaches_ideal_gas_at_low_density():
    rng = np.random.default_rng(4)
    n, box_length, T = 30, 50.0, 1.0
    positions = rng.uniform(0.0, box_length, size=(n, 3))
    velocities = rng.normal(0.0, np.sqrt(T), size=(n, 3))
    velocities -= velocities.mean(axis=0)
    fluid = LJFluid(positions, velocities, box_length=box_length, cutoff=2.5)
    rho = n / box_length**3
    ideal_pressure = rho * T
    assert fluid.pressure() == pytest.approx(ideal_pressure, rel=0.1, abs=1e-6)


def test_radial_distribution_function_approaches_one_for_uniform_random_gas():
    """A dilute, non-interacting (uniformly random) configuration has g(r) ~ 1 away from small r."""
    rng = np.random.default_rng(5)
    n, box_length = 400, 20.0
    g_accum = None
    n_samples = 8
    for _ in range(n_samples):
        positions = rng.uniform(0.0, box_length, size=(n, 3))
        velocities = np.zeros((n, 3))
        fluid = LJFluid(positions, velocities, box_length=box_length, epsilon=1e-12, cutoff=2.5)
        r, g = fluid.radial_distribution_function(r_max=8.0, n_bins=40)
        g_accum = g.copy() if g_accum is None else g_accum + g
    g_mean = g_accum / n_samples
    mid_range = (r > 2.0) & (r < 6.0)
    assert np.mean(g_mean[mid_range]) == pytest.approx(1.0, abs=0.15)


def test_radial_distribution_function_shape():
    fluid = LJFluid.from_lattice(n_per_side=3, cutoff=1.5, density=0.8, temperature=1.0, rng=6)
    r, g = fluid.radial_distribution_function(n_bins=50)
    assert r.shape == (50,)
    assert g.shape == (50,)
    assert np.all(g >= 0.0)


def test_lj_fluid_step_and_forces_and_potential_agree_on_energy_conservation_short_run():
    fluid = LJFluid.from_lattice(n_per_side=3, cutoff=1.5, density=0.7, temperature=1.0, rng=7)
    e0 = fluid.kinetic_energy() + fluid.potential_energy()
    for _ in range(20):
        fluid.step(dt=0.0005)
    e1 = fluid.kinetic_energy() + fluid.potential_energy()
    assert e1 == pytest.approx(e0, rel=1e-2)
