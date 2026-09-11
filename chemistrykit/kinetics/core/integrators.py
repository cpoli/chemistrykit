"""Numba-accelerated numerical integrators for kinetics.

Re-exports the shared, params-array-convention integrators from
:mod:`chemistrykit.integrators`, which is where the implementations live.
See that module for the ``rhs(state, t, params) -> ndarray`` calling
convention required by every function here -- kept explicit rather than
closing over Python scalars so Numba compiles these integrators once and
reuses them, as first-class functions, across every reaction network
defined in :mod:`chemistrykit.kinetics.systems`.

:mod:`chemistrykit.kinetics.core.base_system`'s :class:`ReactionNetwork`
ABC is the intended entry point for most users; this module exists for
code that wants the low-level integrator functions directly (e.g. to
sweep many rate-constant values without going through a class instance
per value).
"""

from __future__ import annotations

from chemistrykit.integrators import RHSFunc, dopri5_integrate, dopri5_step, rk4_integrate, rk4_step

__all__ = [
    "RHSFunc",
    "rk4_step",
    "rk4_integrate",
    "dopri5_step",
    "dopri5_integrate",
]
