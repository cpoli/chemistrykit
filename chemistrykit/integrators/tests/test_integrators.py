"""Tests for chemistrykit.integrators itself.

Ported from physicskit.integrators.tests.test_integrators (same numerics,
same calling convention -- see chemistrykit.integrators.fixed_step for why
the port keeps rk4/leapfrog/yoshida4/dopri5's implementations unchanged and
only rewords the docstrings' physics-specific framing).
"""

import numpy as np
from numba import njit

from chemistrykit.integrators import (
    dopri5_integrate,
    leapfrog_step,
    rk4_integrate,
    velocity_verlet_step,
)


@njit
def _harmonic_force(pos, t, params):
    omega = params[0]
    return -(omega**2) * pos


@njit
def _harmonic_rhs(state, t, params):
    omega = params[0]
    out = np.empty(2)
    out[0] = state[1]
    out[1] = -(omega**2) * state[0]
    return out


def _harmonic_energy(pos, vel, omega):
    return 0.5 * vel**2 + 0.5 * omega**2 * pos**2


def test_velocity_verlet_is_an_alias_for_leapfrog():
    """The two names refer to the same kick-drift-kick scheme -- 'leapfrog'
    in orbital-mechanics literature, 'velocity Verlet' in molecular
    dynamics -- so they must be the identical function object, not merely
    numerically equivalent copies."""
    assert velocity_verlet_step is leapfrog_step


def test_dopri5_matches_analytic_harmonic_oscillator():
    omega = 2.0
    params = np.array([omega])
    state0 = np.array([1.0, 0.0])
    t_end = 10.0

    ts, ys = dopri5_integrate(_harmonic_rhs, state0, 0.0, t_end, 0.01, params, rtol=1e-10, atol=1e-12)

    assert ts[0] == 0.0
    assert ts[-1] == t_end
    expected = np.cos(omega * ts)
    np.testing.assert_allclose(ys[:, 0], expected, atol=1e-6)


def test_dopri5_energy_drift_shrinks_with_tolerance():
    omega = 2.0
    params = np.array([omega])
    state0 = np.array([1.0, 0.0])
    t_end = 20.0

    drifts = []
    for rtol in (1e-4, 1e-9):
        _, ys = dopri5_integrate(_harmonic_rhs, state0, 0.0, t_end, 0.01, params, rtol=rtol, atol=rtol * 1e-3)
        energy = _harmonic_energy(ys[:, 0], ys[:, 1], omega)
        drifts.append(np.max(np.abs(energy - energy[0]) / energy[0]))

    assert drifts[1] < drifts[0]


def test_dopri5_takes_fewer_steps_for_looser_tolerance():
    omega = 2.0
    params = np.array([omega])
    state0 = np.array([1.0, 0.0])
    t_end = 20.0

    ts_tight, _ = dopri5_integrate(_harmonic_rhs, state0, 0.0, t_end, 0.01, params, rtol=1e-10, atol=1e-12)
    ts_loose, _ = dopri5_integrate(_harmonic_rhs, state0, 0.0, t_end, 0.01, params, rtol=1e-3, atol=1e-6)

    assert len(ts_loose) < len(ts_tight)


def test_dopri5_agrees_with_fixed_step_rk4_at_tight_tolerance():
    omega = 2.0
    params = np.array([omega])
    state0 = np.array([1.0, 0.0])
    t_end = 5.0

    ts, ys = dopri5_integrate(_harmonic_rhs, state0, 0.0, t_end, 0.01, params, rtol=1e-10, atol=1e-12)
    _, ys_rk4 = rk4_integrate(_harmonic_rhs, state0, 0.0, 1e-4, 50000, params)

    np.testing.assert_allclose(ys[-1], ys_rk4[-1], atol=1e-6)


def test_rk4_matches_analytic_first_order_decay():
    """A first-order-kinetics sanity check ahead of chemistrykit.kinetics:
    dC/dt = -k*C has the exact solution C(t) = C0 * exp(-k*t); RK4 with a
    small fixed step should reproduce it to a tight tolerance."""

    @njit
    def _decay_rhs(state, t, params):
        k = params[0]
        return -k * state

    k = 0.7
    c0 = np.array([2.0])
    ts, cs = rk4_integrate(_decay_rhs, c0, 0.0, 1e-3, 5000, np.array([k]))
    expected = 2.0 * np.exp(-k * ts)
    np.testing.assert_allclose(cs[:, 0], expected, atol=1e-6)
