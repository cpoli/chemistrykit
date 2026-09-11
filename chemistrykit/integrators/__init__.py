"""Numba-accelerated ODE integrators shared across chemistrykit subpackages.

Right-hand-side / force callbacks use the ``f(state_or_pos, t, params) ->
ndarray`` convention throughout, so a single compiled callback can be reused
across systems with different parameter values without recompiling (see
:mod:`chemistrykit.integrators.fixed_step`).

* :func:`rk4_integrate` -- classical 4th-order Runge-Kutta (not
  norm-preserving; the standard workhorse for reaction-rate ODEs and
  general first-order kinetics).
* :func:`leapfrog_integrate` (alias :func:`velocity_verlet_integrate`) --
  2nd-order symplectic Stormer-Verlet, for separable systems
  ``pos'' = force(pos, t)`` (used by :mod:`chemistrykit.md` for molecular
  dynamics).
* :func:`yoshida4_integrate` -- 4th-order symplectic, built from three
  leapfrog sub-steps.
* :func:`dopri5_integrate` -- adaptive-step-size embedded Dormand-Prince
  RK5(4), for stiff-ish or accuracy-sensitive reaction networks (e.g. ones
  with widely separated rate constants, as in a steady-state-approximation
  chain).

This module is intentionally domain-agnostic: it was ported unchanged from
physicskit's shared ``integrators`` package (same numerics, same calling
convention), since integrating a system of ODEs forward in time is the same
numerical problem whether the state vector holds phase-space coordinates or
chemical concentrations. :mod:`chemistrykit.kinetics.core.integrators`
re-exports these functions rather than reimplementing them.
"""

from chemistrykit.integrators.adaptive import dopri5_integrate, dopri5_step
from chemistrykit.integrators.fixed_step import (
    RHSFunc,
    leapfrog_integrate,
    leapfrog_step,
    rk4_integrate,
    rk4_step,
    velocity_verlet_integrate,
    velocity_verlet_step,
    yoshida4_integrate,
    yoshida4_step,
)

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
    "dopri5_step",
    "dopri5_integrate",
]
