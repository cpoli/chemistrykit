"""Numba-accelerated fixed-step time integrators, shared across subpackages.

Right-hand-side functions passed to these integrators must be module-level
``@njit`` functions with the signature ``rhs(state, t, params) -> ndarray``
(for :func:`rk4_step` / :func:`rk4_integrate`) or ``force(pos, t, params) ->
ndarray`` (for :func:`leapfrog_step` / :func:`leapfrog_integrate` /
:func:`yoshida4_step` / :func:`yoshida4_integrate`). Keeping ``params`` as an
explicit ``float64`` array (rather than closing over Python scalars) lets
Numba compile these integrators once and reuse them, as first-class
functions, across many differently-parameterized systems -- e.g. sweeping a
grid of rate constants without recompiling per value.

:func:`rk4_integrate` is the workhorse for :mod:`chemistrykit.kinetics`'s
reaction-network ODEs ``dC/dt = f(C, t)``, which have no notion of a
"force" and are not, in general, separable or symplectic -- the two
symplectic integrators below (:func:`leapfrog_step`, 2nd order, and
:func:`yoshida4_step`, 4th order) exist for :mod:`chemistrykit.md`'s
molecular-dynamics equations of motion ``pos'' = force(pos, t)`` instead,
where the force may not depend on velocity.
"""

from __future__ import annotations

from typing import Callable

import numpy as np
from numba import njit
from numpy.typing import NDArray

__all__ = [
    "RHSFunc",
    "rk4_step",
    "rk4_integrate",
    "leapfrog_step",
    "leapfrog_integrate",
    "velocity_verlet_step",
    "velocity_verlet_integrate",
    "yoshida4_step",
    "yoshida4_integrate",
]

#: Signature shared by every right-hand-side / force function accepted by
#: the integrators below: ``f(state, t, params) -> derivative``.
RHSFunc = Callable[[NDArray[np.float64], float, NDArray[np.float64]], NDArray[np.float64]]


@njit(cache=True)
def rk4_step(
    rhs: RHSFunc,
    state: NDArray[np.float64],
    t: float,
    dt: float,
    params: NDArray[np.float64],
) -> NDArray[np.float64]:
    """Single classical 4th-order Runge-Kutta step.

    Parameters
    ----------
    rhs : callable
        Numba-jitted right-hand-side function ``rhs(state, t, params) ->
        ndarray``.
    state : ndarray of float, shape (dim,)
        Current state vector.
    t : float
        Current time.
    dt : float
        Step size.
    params : ndarray of float
        Parameter vector passed through to `rhs`.

    Returns
    -------
    ndarray of float, shape (dim,)
        The state advanced by one step of size `dt`.
    """
    k1 = rhs(state, t, params)
    k2 = rhs(state + 0.5 * dt * k1, t + 0.5 * dt, params)
    k3 = rhs(state + 0.5 * dt * k2, t + 0.5 * dt, params)
    k4 = rhs(state + dt * k3, t + dt, params)
    return state + (dt / 6.0) * (k1 + 2.0 * k2 + 2.0 * k3 + k4)


@njit(cache=True)
def rk4_integrate(
    rhs: RHSFunc,
    state0: NDArray[np.float64],
    t0: float,
    dt: float,
    n_steps: int,
    params: NDArray[np.float64],
) -> tuple[NDArray[np.float64], NDArray[np.float64]]:
    """Integrate ``n_steps`` of RK4 starting from ``state0``.

    Parameters
    ----------
    rhs : callable
        Numba-jitted right-hand-side function ``rhs(state, t, params) ->
        ndarray``.
    state0 : ndarray of float, shape (dim,)
        Initial state vector.
    t0 : float
        Initial time.
    dt : float
        Step size.
    n_steps : int
        Number of integration steps.
    params : ndarray of float
        Parameter vector passed through to `rhs`.

    Returns
    -------
    times : ndarray of float, shape (n_steps + 1,)
        Time at each step, starting at `t0`.
    states : ndarray of float, shape (n_steps + 1, dim)
        State at each step, starting at `state0`.
    """
    dim = state0.shape[0]
    states = np.empty((n_steps + 1, dim))
    times = np.empty(n_steps + 1)
    states[0] = state0
    times[0] = t0
    state = state0.copy()
    t = t0
    for i in range(n_steps):
        state = rk4_step(rhs, state, t, dt, params)
        t = t0 + (i + 1) * dt
        states[i + 1] = state
        times[i + 1] = t
    return times, states


@njit(cache=True)
def leapfrog_step(
    force: RHSFunc,
    pos: NDArray[np.float64],
    vel: NDArray[np.float64],
    t: float,
    dt: float,
    params: NDArray[np.float64],
) -> tuple[NDArray[np.float64], NDArray[np.float64]]:
    """Single velocity-Verlet (symplectic) step for ``pos'' = force(pos, t)``.

    Parameters
    ----------
    force : callable
        Numba-jitted acceleration function ``force(pos, t, params) ->
        ndarray``.
    pos : ndarray of float, shape (dim,)
        Current position vector.
    vel : ndarray of float, shape (dim,)
        Current velocity vector.
    t : float
        Current time.
    dt : float
        Step size.
    params : ndarray of float
        Parameter vector passed through to `force`.

    Returns
    -------
    pos_new : ndarray of float, shape (dim,)
        Position advanced by one step of size `dt`.
    vel_new : ndarray of float, shape (dim,)
        Velocity advanced by one step of size `dt`.
    """
    a0 = force(pos, t, params)
    pos_new = pos + vel * dt + 0.5 * a0 * dt * dt
    a1 = force(pos_new, t + dt, params)
    vel_new = vel + 0.5 * (a0 + a1) * dt
    return pos_new, vel_new


@njit(cache=True)
def leapfrog_integrate(
    force: RHSFunc,
    pos0: NDArray[np.float64],
    vel0: NDArray[np.float64],
    t0: float,
    dt: float,
    n_steps: int,
    params: NDArray[np.float64],
) -> tuple[NDArray[np.float64], NDArray[np.float64], NDArray[np.float64]]:
    """Integrate ``n_steps`` of the symplectic leapfrog scheme.

    Parameters
    ----------
    force : callable
        Numba-jitted acceleration function ``force(pos, t, params) ->
        ndarray``.
    pos0 : ndarray of float, shape (dim,)
        Initial position vector.
    vel0 : ndarray of float, shape (dim,)
        Initial velocity vector.
    t0 : float
        Initial time.
    dt : float
        Step size.
    n_steps : int
        Number of integration steps.
    params : ndarray of float
        Parameter vector passed through to `force`.

    Returns
    -------
    times : ndarray of float, shape (n_steps + 1,)
        Time at each step, starting at `t0`.
    positions : ndarray of float, shape (n_steps + 1, dim)
        Position at each step, starting at `pos0`.
    velocities : ndarray of float, shape (n_steps + 1, dim)
        Velocity at each step, starting at `vel0`.
    """
    dim = pos0.shape[0]
    positions = np.empty((n_steps + 1, dim))
    velocities = np.empty((n_steps + 1, dim))
    times = np.empty(n_steps + 1)
    positions[0] = pos0
    velocities[0] = vel0
    times[0] = t0
    pos = pos0.copy()
    vel = vel0.copy()
    t = t0
    for i in range(n_steps):
        pos, vel = leapfrog_step(force, pos, vel, t, dt, params)
        t = t0 + (i + 1) * dt
        positions[i + 1] = pos
        velocities[i + 1] = vel
        times[i + 1] = t
    return times, positions, velocities


#: ``leapfrog_step``/``leapfrog_integrate`` under the name used by
#: molecular-dynamics literature, where this same kick-drift-kick scheme is
#: usually called "velocity Verlet".
velocity_verlet_step = leapfrog_step
velocity_verlet_integrate = leapfrog_integrate


#: Yoshida (1990) 4th-order composition coefficients: the 4th-order step is
#: three 2nd-order (leapfrog) sub-steps of sizes ``w1*dt``, ``w0*dt``,
#: ``w1*dt`` (which sum to ``dt``, since ``2*w1 + w0 == 1``).
_YOSHIDA_CBRT2 = 2.0 ** (1.0 / 3.0)
_YOSHIDA_W0 = -_YOSHIDA_CBRT2 / (2.0 - _YOSHIDA_CBRT2)
_YOSHIDA_W1 = 1.0 / (2.0 - _YOSHIDA_CBRT2)


@njit(cache=True)
def yoshida4_step(
    force: RHSFunc,
    pos: NDArray[np.float64],
    vel: NDArray[np.float64],
    t: float,
    dt: float,
    params: NDArray[np.float64],
) -> tuple[NDArray[np.float64], NDArray[np.float64]]:
    """Single 4th-order symplectic (Yoshida) step for ``pos'' = force(pos, t)``.

    Composes three :func:`leapfrog_step` sub-steps with Yoshida's (1990)
    coefficients to raise the (still symplectic) accuracy from 2nd to 4th
    order, at roughly 3x the cost per step of plain leapfrog. Because it
    remains symplectic, it -- like leapfrog -- conserves energy far better
    than RK4 over long integrations, now with much smaller local truncation
    error too.

    Parameters
    ----------
    force : callable
        Numba-jitted acceleration function ``force(pos, t, params) ->
        ndarray``.
    pos : ndarray of float, shape (dim,)
        Current position vector.
    vel : ndarray of float, shape (dim,)
        Current velocity vector.
    t : float
        Current time.
    dt : float
        Step size.
    params : ndarray of float
        Parameter vector passed through to `force`.

    Returns
    -------
    pos_new : ndarray of float, shape (dim,)
        Position advanced by one step of size `dt`.
    vel_new : ndarray of float, shape (dim,)
        Velocity advanced by one step of size `dt`.
    """
    h1 = _YOSHIDA_W1 * dt
    h0 = _YOSHIDA_W0 * dt
    pos, vel = leapfrog_step(force, pos, vel, t, h1, params)
    t = t + h1
    pos, vel = leapfrog_step(force, pos, vel, t, h0, params)
    t = t + h0
    pos, vel = leapfrog_step(force, pos, vel, t, h1, params)
    return pos, vel


@njit(cache=True)
def yoshida4_integrate(
    force: RHSFunc,
    pos0: NDArray[np.float64],
    vel0: NDArray[np.float64],
    t0: float,
    dt: float,
    n_steps: int,
    params: NDArray[np.float64],
) -> tuple[NDArray[np.float64], NDArray[np.float64], NDArray[np.float64]]:
    """Integrate ``n_steps`` of the 4th-order symplectic Yoshida scheme.

    Parameters
    ----------
    force : callable
        Numba-jitted acceleration function ``force(pos, t, params) ->
        ndarray``.
    pos0 : ndarray of float, shape (dim,)
        Initial position vector.
    vel0 : ndarray of float, shape (dim,)
        Initial velocity vector.
    t0 : float
        Initial time.
    dt : float
        Step size.
    n_steps : int
        Number of integration steps.
    params : ndarray of float
        Parameter vector passed through to `force`.

    Returns
    -------
    times : ndarray of float, shape (n_steps + 1,)
        Time at each step, starting at `t0`.
    positions : ndarray of float, shape (n_steps + 1, dim)
        Position at each step, starting at `pos0`.
    velocities : ndarray of float, shape (n_steps + 1, dim)
        Velocity at each step, starting at `vel0`.
    """
    dim = pos0.shape[0]
    positions = np.empty((n_steps + 1, dim))
    velocities = np.empty((n_steps + 1, dim))
    times = np.empty(n_steps + 1)
    positions[0] = pos0
    velocities[0] = vel0
    times[0] = t0
    pos = pos0.copy()
    vel = vel0.copy()
    t = t0
    for i in range(n_steps):
        pos, vel = yoshida4_step(force, pos, vel, t, dt, params)
        t = t0 + (i + 1) * dt
        positions[i + 1] = pos
        velocities[i + 1] = vel
        times[i + 1] = t
    return times, positions, velocities
