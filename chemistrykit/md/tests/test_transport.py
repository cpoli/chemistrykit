"""Tests for chemistrykit.md.systems.transport."""

import numpy as np
import pytest

from chemistrykit.md.systems.transport import (
    einstein_diffusion_coefficient,
    green_kubo_diffusion_coefficient,
    mean_squared_displacement,
    unwrap_trajectory,
    velocity_autocorrelation,
)


def test_unwrap_recovers_free_flight_through_periodic_box():
    t = np.arange(50) * 0.4
    true = np.stack([1.0 + 0.7 * t, 2.0 - 0.3 * t, 0.5 * t], axis=-1)[:, None, :]
    wrapped = np.mod(true, 3.0)
    np.testing.assert_allclose(unwrap_trajectory(wrapped, 3.0), true, atol=1e-12)


def test_msd_of_ballistic_motion_is_v_squared_t_squared():
    rng = np.random.default_rng(0)
    v = rng.normal(size=(20, 3))
    t = np.arange(30.0)
    positions = t[:, None, None] * v[None, :, :]
    msd = mean_squared_displacement(positions, max_lag=10)
    np.testing.assert_allclose(msd, np.mean(np.sum(v * v, axis=1)) * np.arange(11) ** 2)


def test_einstein_recovers_random_walk_diffusion_coefficient():
    rng = np.random.default_rng(1)
    D, dt, n_frames = 0.3, 0.1, 400
    steps = rng.normal(0.0, np.sqrt(2.0 * D * dt), size=(n_frames - 1, 500, 3))
    positions = np.concatenate([np.zeros((1, 500, 3)), np.cumsum(steps, axis=0)])
    msd = mean_squared_displacement(positions, max_lag=100)
    lags = np.arange(101) * dt
    assert einstein_diffusion_coefficient(lags, msd) == pytest.approx(D, rel=0.05)


def test_green_kubo_matches_ornstein_uhlenbeck_diffusion():
    # OU velocities: C(t) = (kT/m) e^{-gamma t} per component, D = kT/(m gamma).
    rng = np.random.default_rng(2)
    kT_over_m, gamma, dt, n_frames, n_particles = 1.0, 2.0, 0.02, 20000, 200
    a = np.exp(-gamma * dt)
    v = np.empty((n_frames, n_particles, 3))
    v[0] = rng.normal(0.0, np.sqrt(kT_over_m), size=(n_particles, 3))
    noise = rng.normal(0.0, np.sqrt(kT_over_m * (1 - a * a)), size=(n_frames - 1, n_particles, 3))
    for k in range(1, n_frames):
        v[k] = a * v[k - 1] + noise[k - 1]
    vacf = velocity_autocorrelation(v, max_lag=200)
    lags = np.arange(201) * dt
    np.testing.assert_allclose(vacf[:50], 3 * kT_over_m * np.exp(-gamma * lags[:50]), atol=0.03)
    assert green_kubo_diffusion_coefficient(lags, vacf) == pytest.approx(kT_over_m / gamma, rel=0.05)


def test_vacf_normalization_and_bad_lag():
    v = np.random.default_rng(3).normal(size=(10, 4, 3))
    assert velocity_autocorrelation(v, normalize=True)[0] == pytest.approx(1.0)
    with pytest.raises(ValueError):
        velocity_autocorrelation(v, max_lag=10)
