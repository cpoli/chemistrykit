"""Abstract base classes for kinetics models, and the result container.

Reaction-network kinetics (mass-action rate laws integrated over time) is
a first-order ODE system ``dC/dt = f(C, t)`` -- the same calling
convention as :mod:`chemistrykit.integrators`, and with no notion of a
conserved "energy" or "momentum" the way there is for
``chemistrykit.classical``'s (as-yet-unbuilt) Hamiltonian systems. The
abstract interface below therefore mirrors physicskit's
``physicskit.chaos.core.base_system.DynamicalSystem`` (a bare vector field
driven by the shared RK4/adaptive integrators) rather than
``physicskit.classical.core.base_system``'s Hamiltonian/Lagrangian
hierarchy, which doesn't apply here.

Two ABCs cover the two shapes of model in this subpackage:

* :class:`RateLaw` -- a single elementary reaction with a closed-form
  integrated concentration profile (:mod:`chemistrykit.kinetics.systems.rate_laws`),
  where no numerical integration is needed at all.
* :class:`ReactionNetwork` -- a system of coupled species whose
  concentrations are only available by integrating the ODEs numerically
  (:mod:`chemistrykit.kinetics.systems.networks`,
  :mod:`chemistrykit.kinetics.systems.oscillators`, and the progress-curve
  model in :mod:`chemistrykit.kinetics.systems.enzyme`).

Notes
-----
Concrete :class:`ReactionNetwork` subclasses build a standalone
module-level ``@njit`` right-hand-side function (conventionally stored as
``self._rhs_njit``), following the first-class-function pattern documented
in ``physicskit.classical.core.base_system``: Numba's nopython mode can
only call a genuine njit dispatcher from inside another njit function, not
a bound Python method.
"""

from __future__ import annotations

from abc import ABC, abstractmethod
from collections.abc import Sequence
from dataclasses import dataclass, field

import numpy as np

from chemistrykit.integrators import dopri5_integrate, rk4_integrate

__all__ = ["KineticsResult", "RateLaw", "ReactionNetwork"]


@dataclass
class KineticsResult:
    """Container for the output of a :meth:`ReactionNetwork.integrate` call.

    Mirrors physicskit's ``SimulationResult`` (a stable, dataclass return
    type consumed by visualizers), specialized to a named-species
    concentration trajectory instead of a phase-space one.
    """

    t: np.ndarray
    """ndarray, shape (n_steps + 1,): Time samples."""

    y: np.ndarray
    """ndarray, shape (n_steps + 1, n_species): Concentration of every
    species (columns, in the order given by ``species``) at each time
    sample (rows)."""

    species: Sequence[str] = field(default_factory=tuple)
    """tuple of str: Species names, in column order matching ``y``."""

    method: str = ""
    """str: Name of the integrator used (e.g. ``"rk4"``)."""

    extra: dict = field(default_factory=dict)
    """dict: Free-form slot for any additional diagnostics a system
    chooses to attach."""

    def concentration(self, name: str) -> np.ndarray:
        """Return the concentration trajectory of a single named species.

        Parameters
        ----------
        name : str
            A species name from ``self.species``.

        Returns
        -------
        ndarray, shape (n_steps + 1,)

        Examples
        --------
        >>> import numpy as np
        >>> result = KineticsResult(t=np.array([0.0, 1.0]), y=np.array([[1.0, 0.0], [0.5, 0.5]]), species=("A", "B"))
        >>> result.concentration("B")
        array([0. , 0.5])
        """
        idx = list(self.species).index(name)
        return self.y[:, idx]


class RateLaw(ABC):
    """Base class for a single elementary reaction's integrated rate law.

    Subclasses implement the closed-form solution of ``dC/dt = -rate(C)``
    (or, for zero order, ``+rate`` with a sign convention appropriate to
    that order) directly -- there is nothing to numerically integrate for
    these textbook cases, unlike :class:`ReactionNetwork`.
    """

    @abstractmethod
    def concentration(self, t):
        """Return the reactant concentration at time(s) `t`.

        Parameters
        ----------
        t : float or array-like of float

        Returns
        -------
        float or ndarray
        """

    @abstractmethod
    def rate(self, t=None):
        """Return the instantaneous reaction rate at time(s) `t`.

        Parameters
        ----------
        t : float or array-like of float, optional
            Defaults to 0 (the initial rate) where the rate law does not
            need `t` explicitly.

        Returns
        -------
        float or ndarray
        """

    @abstractmethod
    def half_life(self):
        """Return the time for the concentration to fall to half its initial value.

        Returns
        -------
        float
        """


class ReactionNetwork(ABC):
    """Common base for a system of species evolving as ``dC/dt = f(C, t)``.

    Concrete subclasses must set, in ``__init__``:

    - ``self.species`` : tuple of str -- ordered species names.
    - ``self._rhs_njit`` : an ``@njit`` dispatcher with signature
      ``(state, t, params) -> dstate`` (the
      :data:`chemistrykit.integrators.RHSFunc` convention).
    - ``self.params`` : ndarray -- the parameter vector read by
      ``self._rhs_njit`` (or ``np.empty(0)`` if every rate constant is
      instead baked into the njit closure itself, as
      :mod:`chemistrykit.kinetics.systems.networks`'s mass-action engine
      does with its stoichiometry/rate-constant arrays).

    and call ``super().__init__(state0)``.
    """

    species: Sequence[str] = ()
    _rhs_njit = None
    params: np.ndarray = np.empty(0)

    def __init__(self, state0):
        self.state = np.asarray(state0, dtype=np.float64)
        self.t = 0.0

    @abstractmethod
    def rhs(self, state: np.ndarray, t: float) -> np.ndarray:
        """Evaluate ``dC/dt`` at ``(state, t)``.

        A thin, plain-Python wrapper around ``self._rhs_njit`` (for
        interactive use / plotting outside a numba context).

        Parameters
        ----------
        state : ndarray
            Current concentration vector.
        t : float
            Current time.

        Returns
        -------
        ndarray
        """

    def reset(self, state0=None, t0: float = 0.0):
        """Reset the system's state and clock.

        Parameters
        ----------
        state0 : array-like, optional
            New state; if omitted, the current state is kept.
        t0 : float, default 0.0
            New time.

        Returns
        -------
        ndarray
            The (possibly updated) current state.
        """
        if state0 is not None:
            self.state = np.asarray(state0, dtype=np.float64)
        self.t = t0
        return self.state

    def integrate(self, t_span, dt=None, method: str = "rk4", **kwargs) -> KineticsResult:
        """Integrate the reaction network forward in time.

        Parameters
        ----------
        t_span : tuple of float
            ``(t0, t1)``, start and end time.
        dt : float, optional
            Fixed step size, required for ``method="rk4"``. For
            ``method="dopri5"`` it is instead used as the *initial* step
            size attempt (defaulting to ``(t1 - t0) / 1000``), since the
            adaptive integrator adjusts it automatically.
        method : {"rk4", "dopri5"}
            Integrator to use. ``"dopri5"`` (adaptive Dormand-Prince) is
            recommended for stiff-ish networks with widely separated rate
            constants (e.g. a steady-state-approximation chain).
        **kwargs
            Forwarded to :func:`chemistrykit.integrators.dopri5_integrate`
            (e.g. ``rtol``, ``atol``) when ``method="dopri5"``.

        Returns
        -------
        KineticsResult
        """
        t0, t1 = t_span
        if method == "rk4":
            if dt is None:
                raise ValueError("dt is required for method='rk4'")
            n_steps = int(round((t1 - t0) / dt))
            ts, ys = rk4_integrate(self._rhs_njit, self.state, t0, dt, n_steps, self.params)
        elif method == "dopri5":
            dt0 = dt if dt is not None else (t1 - t0) / 1000.0
            ts, ys = dopri5_integrate(self._rhs_njit, self.state, t0, t1, dt0, self.params, **kwargs)
        else:
            raise ValueError(f"Unknown method '{method}' for ReactionNetwork")
        self.state = ys[-1]
        self.t = ts[-1]
        return KineticsResult(t=ts, y=ys, species=tuple(self.species), method=method)
