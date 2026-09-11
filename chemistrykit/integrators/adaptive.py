"""Adaptive-step-size integration: embedded Dormand-Prince RK5(4).

This is deliberately offered only for the general first-order path
(``rhs(state, t, params) -> ndarray``), not for the symplectic
:func:`~chemistrykit.integrators.fixed_step.leapfrog_step` /
:func:`~chemistrykit.integrators.fixed_step.yoshida4_step` family. Those
integrators conserve the Hamiltonian over long integrations *because* they
take a fixed step; naively varying ``dt`` from an embedded error estimate
breaks that guarantee; a genuinely adaptive symplectic scheme needs a
reversible step-size controller (see Hairer, Lubich & Wanner, *Geometric
Numerical Integration*), which is a different algorithm, not a drop-in
extension of this one. Use :func:`dopri5_integrate` for reaction-network
kinetics with widely separated rate constants (a common source of numerical
stiffness for chain and steady-state-approximation mechanisms) or any other
non-conservative/driven system where accuracy control matters more than
exact invariant conservation; use the fixed-step symplectic integrators for
long-time Hamiltonian/molecular-dynamics evolution.
"""

from __future__ import annotations

import numpy as np
from numba import njit
from numpy.typing import NDArray

from chemistrykit.integrators.fixed_step import RHSFunc

__all__ = ["dopri5_step", "dopri5_integrate"]

# Dormand & Prince (1980) 5(4) coefficients (the classic "ode45" tableau).
_C2, _C3, _C4, _C5 = 1.0 / 5.0, 3.0 / 10.0, 4.0 / 5.0, 8.0 / 9.0

_A21 = 1.0 / 5.0
_A31, _A32 = 3.0 / 40.0, 9.0 / 40.0
_A41, _A42, _A43 = 44.0 / 45.0, -56.0 / 15.0, 32.0 / 9.0
_A51, _A52, _A53, _A54 = 19372.0 / 6561.0, -25360.0 / 2187.0, 64448.0 / 6561.0, -212.0 / 729.0
_A61, _A62, _A63, _A64, _A65 = 9017.0 / 3168.0, -355.0 / 33.0, 46732.0 / 5247.0, 49.0 / 176.0, -5103.0 / 18656.0

# 5th-order solution weights (= 7th-stage row, i.e. this method is FSAL).
_B1, _B3, _B4, _B5, _B6 = 35.0 / 384.0, 500.0 / 1113.0, 125.0 / 192.0, -2187.0 / 6784.0, 11.0 / 84.0
# Embedded 4th-order weights, for the error estimate y5 - y4.
_B1S, _B3S, _B4S, _B5S, _B6S, _B7S = (
    5179.0 / 57600.0,
    7571.0 / 16695.0,
    393.0 / 640.0,
    -92097.0 / 339200.0,
    187.0 / 2100.0,
    1.0 / 40.0,
)
_E1 = _B1 - _B1S
_E3 = _B3 - _B3S
_E4 = _B4 - _B4S
_E5 = _B5 - _B5S
_E6 = _B6 - _B6S
_E7 = -_B7S


@njit(cache=True)
def dopri5_step(
    rhs: RHSFunc,
    state: NDArray[np.float64],
    t: float,
    dt: float,
    params: NDArray[np.float64],
) -> tuple[NDArray[np.float64], NDArray[np.float64]]:
    """Single embedded Dormand-Prince RK5(4) step, with an error estimate.

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
        Step size to attempt.
    params : ndarray of float
        Parameter vector passed through to `rhs`.

    Returns
    -------
    state_new : ndarray of float, shape (dim,)
        The 5th-order state estimate at ``t + dt``.
    error : ndarray of float, shape (dim,)
        Difference between the 5th- and embedded 4th-order estimates, for
        step-size control (see :func:`dopri5_integrate`).
    """
    k1 = rhs(state, t, params)
    k2 = rhs(state + dt * (_A21 * k1), t + _C2 * dt, params)
    k3 = rhs(state + dt * (_A31 * k1 + _A32 * k2), t + _C3 * dt, params)
    k4 = rhs(state + dt * (_A41 * k1 + _A42 * k2 + _A43 * k3), t + _C4 * dt, params)
    k5 = rhs(state + dt * (_A51 * k1 + _A52 * k2 + _A53 * k3 + _A54 * k4), t + _C5 * dt, params)
    k6 = rhs(state + dt * (_A61 * k1 + _A62 * k2 + _A63 * k3 + _A64 * k4 + _A65 * k5), t + dt, params)
    state_new = state + dt * (_B1 * k1 + _B3 * k3 + _B4 * k4 + _B5 * k5 + _B6 * k6)
    k7 = rhs(state_new, t + dt, params)  # FSAL: also the first stage of the next step
    error = dt * (_E1 * k1 + _E3 * k3 + _E4 * k4 + _E5 * k5 + _E6 * k6 + _E7 * k7)
    return state_new, error


@njit(cache=True)
def dopri5_integrate(
    rhs: RHSFunc,
    state0: NDArray[np.float64],
    t0: float,
    t_end: float,
    dt0: float,
    params: NDArray[np.float64],
    rtol: float = 1e-6,
    atol: float = 1e-9,
    dt_min: float = 1e-12,
    dt_max: float = 1e6,
    safety: float = 0.9,
    max_steps: int = 100_000,
) -> tuple[NDArray[np.float64], NDArray[np.float64]]:
    """Integrate forward from ``t0`` to ``t_end`` with adaptive step-size control.

    Steps are accepted or rejected from the embedded RK5(4) error estimate
    against an ``atol + rtol * |state|`` tolerance (the same convention as
    :func:`scipy.integrate.solve_ivp`), with the step size adjusted after
    every attempt. Only forward integration (``t_end > t0``, ``dt0 > 0``) is
    supported.

    Parameters
    ----------
    rhs : callable
        Numba-jitted right-hand-side function ``rhs(state, t, params) ->
        ndarray``.
    state0 : ndarray of float, shape (dim,)
        Initial state vector.
    t0, t_end : float
        Integration interval, with ``t_end > t0``.
    dt0 : float
        Initial step size to attempt.
    params : ndarray of float
        Parameter vector passed through to `rhs`.
    rtol, atol : float
        Relative and absolute error tolerances.
    dt_min, dt_max : float
        Step-size bounds; a step is accepted once ``dt`` shrinks to
        `dt_min` regardless of its error estimate, to guarantee progress.
    safety : float, default=0.9
        Safety factor applied to the step-size update.
    max_steps : int
        Upper bound on the number of attempted steps (accepted or
        rejected), to guarantee termination.

    Returns
    -------
    times : ndarray of float, shape (n_accepted + 1,)
        Times of the accepted steps, starting at `t0` and ending at
        `t_end`.
    states : ndarray of float, shape (n_accepted + 1, dim)
        State at each accepted step.
    """
    dim = state0.shape[0]
    cap = 1024
    times = np.empty(cap)
    states = np.empty((cap, dim))
    times[0] = t0
    states[0] = state0

    t = t0
    state = state0.copy()
    dt = dt0
    count = 0
    for _ in range(max_steps):
        if t >= t_end:
            break
        if t + dt > t_end:
            dt = t_end - t

        state_new, error = dopri5_step(rhs, state, t, dt, params)
        scale = atol + rtol * np.maximum(np.abs(state), np.abs(state_new))
        err_norm = np.sqrt(np.mean((error / scale) ** 2))

        if err_norm <= 1.0 or dt <= dt_min:
            t = t + dt
            state = state_new
            count += 1
            if count >= cap:
                new_cap = cap * 2
                new_times = np.empty(new_cap)
                new_states = np.empty((new_cap, dim))
                new_times[:cap] = times
                new_states[:cap] = states
                times = new_times
                states = new_states
                cap = new_cap
            times[count] = t
            states[count] = state

        if err_norm == 0.0:
            factor = 5.0
        else:
            factor = safety * err_norm ** (-0.2)
            factor = min(max(factor, 0.2), 5.0)
        dt = min(max(dt * factor, dt_min), dt_max)

    return times[: count + 1], states[: count + 1]
