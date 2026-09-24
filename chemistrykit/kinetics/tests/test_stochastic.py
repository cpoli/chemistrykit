"""Tests for chemistrykit.kinetics.systems.stochastic (Gillespie SSA)."""

import numpy as np
import pytest

from chemistrykit.kinetics.systems.stochastic import gillespie_ssa


def test_decay_ensemble_mean_matches_exponential():
    """Pure death A -> 0: E[n(t)] = n0 exp(-c t) exactly."""
    rng = np.random.default_rng(42)
    n0, c, t = 200, 0.5, 2.0
    finals = [gillespie_ssa([[-1]], [c], [[1]], [n0], t_max=5.0, seed=rng).sample([t])[0, 0] for _ in range(300)]
    # binomial std of a single run is ~6.9, so the mean of 300 runs has std ~0.4
    assert np.mean(finals) == pytest.approx(n0 * np.exp(-c * t), abs=2.0)


def test_immigration_death_stationary_distribution_is_poisson():
    """0 -> A (k), A -> 0 (gamma): stationary mean = variance = k/gamma."""
    k, gamma = 20.0, 1.0
    traj = gillespie_ssa([[1, -1]], [k, gamma], [[0, 1]], [0], t_max=2000.0, seed=7)
    samples = traj.sample(np.linspace(50.0, 2000.0, 4000))[:, 0]
    assert samples.mean() == pytest.approx(k / gamma, rel=0.05)
    assert samples.var() == pytest.approx(k / gamma, rel=0.15)


def test_dimerization_propensity_uses_pair_count():
    """2A -> B with c: first waiting time has rate c*n(n-1)/2; with n=2 the
    only event converts both A into one B and stops."""
    traj = gillespie_ssa([[-2], [1]], [1.0], [[2], [0]], [2, 0], t_max=1e6, species=("A", "B"), seed=0)
    assert traj.count("A").tolist() == [2, 0]
    assert traj.count("B").tolist() == [0, 1]


def test_counts_stay_nonnegative_and_seed_is_reproducible():
    args = ([[-1, 0], [1, -1]], [1.0, 2.0], [[1, 0], [0, 1]], [100, 0], 10.0)
    a = gillespie_ssa(*args, seed=11)
    b = gillespie_ssa(*args, seed=11)
    np.testing.assert_array_equal(a.counts, b.counts)
    assert a.counts.min() >= 0
