r"""Abstract base classes for molecular-dynamics models, and the result container.

Every model in :mod:`chemistrykit.md` is a system of point particles
evolving under Newton's second law, ``m_i d^2r_i/dt^2 = F_i(\{r\})``, which
-- like :mod:`chemistrykit.kinetics`'s reaction networks -- has exactly one
natural shape (a state that is integrated forward in time), so a single
ABC, :class:`MolecularDynamicsSystem`, covers every concrete system in
:mod:`chemistrykit.md.systems` (:class:`~chemistrykit.md.systems.lj_fluid.LJFluid`,
:class:`~chemistrykit.md.systems.pair_potentials.DiatomicOscillator`,
:class:`~chemistrykit.md.systems.pair_potentials.HarmonicMolecule`) --
mirroring :class:`chemistrykit.kinetics.core.base_system.ReactionNetwork`.
Two smaller ABCs, :class:`PairPotential` and :class:`AnglePotential`,
capture the two shapes of *interaction* a force field is built from (a
function of one interatomic distance, or of one bond angle).

Concrete :class:`MolecularDynamicsSystem` subclasses build a standalone
module-level ``@njit`` acceleration function (conventionally stored as
``self._accel_njit``), following the same factory-closure pattern
documented in :mod:`chemistrykit.kinetics.core.base_system` and used by
physicskit's N-body/classical systems: Numba's nopython mode can only call
a genuine njit dispatcher from inside another njit function (here, the
shared :func:`chemistrykit.integrators.velocity_verlet_step`), not a bound
Python method.
"""

from __future__ import annotations

from abc import ABC, abstractmethod
from dataclasses import dataclass, field
from typing import Optional

import numpy as np

from chemistrykit.integrators import velocity_verlet_step

__all__ = ["MDResult", "PairPotential", "AnglePotential", "MolecularDynamicsSystem"]


@dataclass
class MDResult:
    """Container for the output of a :meth:`MolecularDynamicsSystem.run` call.

    Mirrors :class:`chemistrykit.kinetics.core.base_system.KineticsResult`
    (a stable, dataclass return type consumed by visualizers), specialized
    to a particle trajectory instead of a concentration one.
    """

    t: np.ndarray
    """ndarray, shape (n_frames,): Time at each sampled frame."""

    positions: np.ndarray
    """ndarray, shape (n_frames, n_particles, n_dim): Particle positions at
    each sampled frame."""

    velocities: np.ndarray
    """ndarray, shape (n_frames, n_particles, n_dim): Particle velocities at
    each sampled frame."""

    kinetic_energy: np.ndarray
    """ndarray, shape (n_frames,): Total kinetic energy at each frame."""

    potential_energy: np.ndarray
    """ndarray, shape (n_frames,): Total potential energy at each frame."""

    temperature: np.ndarray
    """ndarray, shape (n_frames,): Instantaneous kinetic temperature (see
    :meth:`MolecularDynamicsSystem.temperature`) at each frame."""

    box_length: Optional[float] = None
    """float or None: Cubic periodic box side length, or ``None`` for a
    non-periodic system."""

    extra: dict = field(default_factory=dict)
    """dict: Free-form slot for any additional diagnostics a system
    chooses to attach (e.g. pressure)."""

    @property
    def total_energy(self) -> np.ndarray:
        """ndarray, shape (n_frames,): ``kinetic_energy + potential_energy``.

        For a system with no thermostat attached this should be constant
        to within integration error -- the standard NVE
        energy-conservation check for a symplectic integrator like
        velocity-Verlet.
        """
        return self.kinetic_energy + self.potential_energy

    def speeds(self) -> np.ndarray:
        """Return per-particle speeds at every sampled frame.

        Returns
        -------
        ndarray, shape (n_frames, n_particles)

        Examples
        --------
        >>> import numpy as np
        >>> result = MDResult(
        ...     t=np.array([0.0]),
        ...     positions=np.zeros((1, 2, 3)),
        ...     velocities=np.array([[[3.0, 4.0, 0.0], [0.0, 0.0, 0.0]]]),
        ...     kinetic_energy=np.array([0.0]),
        ...     potential_energy=np.array([0.0]),
        ...     temperature=np.array([0.0]),
        ... )
        >>> result.speeds()
        array([[5., 0.]])
        """
        return np.linalg.norm(self.velocities, axis=-1)


class PairPotential(ABC):
    """Common interface for a radially symmetric two-body potential U(r).

    Concrete subclasses (:class:`~chemistrykit.md.systems.lj_fluid.LennardJones`,
    :class:`~chemistrykit.md.systems.pair_potentials.Morse`,
    :class:`~chemistrykit.md.systems.pair_potentials.Buckingham`,
    :class:`~chemistrykit.md.systems.pair_potentials.HarmonicBond`)
    implement :meth:`energy` and :meth:`force_scalar` so they can be
    swapped in and compared directly, mirroring the role
    :class:`chemistrykit.thermo.core.base_system.EquationOfState` plays
    for equations of state.
    """

    @abstractmethod
    def energy(self, r):
        """Return the potential energy U(r) at separation(s) `r`.

        Parameters
        ----------
        r : float or array-like of float

        Returns
        -------
        float or ndarray
        """

    @abstractmethod
    def force_scalar(self, r):
        r"""Return the radial force magnitude :math:`f(r) = -dU/dr`.

        By convention, positive is repulsive (pushes the two particles
        apart along their separation vector :math:`\vec r_i - \vec r_j`)
        and negative is attractive -- so the force on particle `i` due to
        particle `j` is ``f(r) * (r_i - r_j) / r``.

        Parameters
        ----------
        r : float or array-like of float

        Returns
        -------
        float or ndarray
        """


class AnglePotential(ABC):
    r"""Common interface for a three-body bond-angle potential :math:`U(\theta)`.

    Used for angle-bending terms in :class:`~chemistrykit.md.systems.pair_potentials.HarmonicMolecule`,
    where :math:`\theta` is the angle at a vertex atom `j` between the two
    bonds `i-j` and `k-j`.
    """

    @abstractmethod
    def energy(self, theta):
        r"""Return :math:`U(\theta)` for angle(s) `theta`, in radians."""

    @abstractmethod
    def torque(self, theta):
        r"""Return :math:`-dU/d\theta` for angle(s) `theta`, in radians."""


class MolecularDynamicsSystem(ABC):
    """Common engine for a system of particles integrated by velocity-Verlet.

    Concrete subclasses must set, in ``__init__`` (before calling
    ``super().__init__(...)``):

    - ``self._accel_njit`` : an ``@njit`` dispatcher with signature
      ``(positions, t, params) -> acceleration`` (the
      :data:`chemistrykit.integrators.RHSFunc` convention, specialized so
      the "state" is a position array and the return value is
      *acceleration*, i.e. force already divided by mass -- exactly the
      convention :func:`chemistrykit.integrators.velocity_verlet_step`
      expects).
    - ``self.params`` : ndarray -- the parameter vector read by
      ``self._accel_njit`` (``np.empty(0)`` if every numeric parameter is
      instead baked into the njit closure itself).
    - ``self.box_length`` : float or ``None`` -- cubic periodic box side
      length, for systems that use periodic boundary conditions
      (``None`` for a finite, non-periodic cluster).

    and call ``super().__init__(positions, velocities, masses)``.

    Parameters
    ----------
    positions, velocities : array-like, shape (n_particles, n_dim)
    masses : array-like, shape (n_particles,)
    """

    box_length: Optional[float] = None
    _accel_njit = None
    params: np.ndarray = np.empty(0)

    def __init__(self, positions, velocities, masses):
        self.positions = np.array(positions, dtype=np.float64, copy=True)
        self.velocities = np.array(velocities, dtype=np.float64, copy=True)
        self.masses = np.asarray(masses, dtype=np.float64)
        self.t = 0.0

    def degrees_of_freedom(self) -> int:
        """Number of kinetic degrees of freedom, assuming zero net momentum.

        ``n_dim * n_particles - n_dim``: total Cartesian degrees of
        freedom minus the `n_dim` removed by fixing the center-of-mass
        velocity -- standard practice whenever the total momentum is
        conserved/zeroed at initialization (Allen & Tildesley, *Computer
        Simulation of Liquids*, 2nd ed., Ch. 2.6). Subclasses with
        different constraints (e.g. a fixed number of internal
        vibrational modes) should override this.

        Returns
        -------
        int
        """
        n, d = self.positions.shape
        return d * n - d

    def kinetic_energy(self) -> float:
        r"""Total kinetic energy, :math:`K = \sum_i \frac{1}{2} m_i v_i^2`.

        Returns
        -------
        float
        """
        return float(0.5 * np.sum(self.masses[:, None] * self.velocities**2))

    @abstractmethod
    def forces_and_potential(self, positions):
        """Return the force on every particle and the total potential energy.

        Parameters
        ----------
        positions : ndarray, shape (n_particles, n_dim)

        Returns
        -------
        forces : ndarray, shape (n_particles, n_dim)
        potential_energy : float
        """

    def potential_energy(self) -> float:
        """Total potential energy at the system's current configuration.

        Returns
        -------
        float
        """
        _, potential = self.forces_and_potential(self.positions)
        return float(potential)

    def temperature(self, k_b: float = 1.0) -> float:
        r"""Instantaneous kinetic temperature via the equipartition theorem.

        .. math::

            T = \frac{2K}{k_B \times \text{dof}}

        (Allen & Tildesley, *Computer Simulation of Liquids*, 2nd ed.,
        eq. 2.60.) In the reduced units used by
        :class:`~chemistrykit.md.systems.lj_fluid.LJFluid`, ``k_b=1.0``
        (the default) is the correct choice.

        Parameters
        ----------
        k_b : float, default 1.0
            Boltzmann constant in the system's unit system.

        Returns
        -------
        float
        """
        dof = self.degrees_of_freedom()
        return 2.0 * self.kinetic_energy() / (dof * k_b)

    def step(self, dt: float) -> None:
        """Advance the stored state by one velocity-Verlet step, in place.

        Delegates the actual kick-drift-kick arithmetic to the shared
        :func:`chemistrykit.integrators.velocity_verlet_step` -- this
        method exists only to supply that function with this system's
        ``_accel_njit``/``params``, exactly as
        :meth:`chemistrykit.kinetics.core.base_system.ReactionNetwork.integrate`
        wraps :func:`chemistrykit.integrators.rk4_integrate`.

        Parameters
        ----------
        dt : float
            Time step.
        """
        self.positions, self.velocities = velocity_verlet_step(self._accel_njit, self.positions, self.velocities, self.t, dt, self.params)
        self.t += dt

    def run(self, dt: float, n_steps: int, thermostat=None, sample_every: int = 1) -> MDResult:
        """Integrate forward `n_steps` velocity-Verlet steps, sampling as it goes.

        Parameters
        ----------
        dt : float
            Time step.
        n_steps : int
            Number of integration steps (must be a multiple of `sample_every`).
        thermostat : object, optional
            An object with an ``apply(system, dt)`` method, called once
            after every step (e.g.
            :class:`~chemistrykit.md.systems.thermostats.VelocityRescalingThermostat`,
            :class:`~chemistrykit.md.systems.thermostats.NoseHooverThermostat`).
            Omit for a microcanonical (NVE) run.
        sample_every : int, default 1
            Record a frame every this many steps (plus the initial state).

        Returns
        -------
        MDResult
        """
        if n_steps % sample_every != 0:
            raise ValueError("n_steps must be a multiple of sample_every")
        n_frames = n_steps // sample_every + 1
        n, d = self.positions.shape
        times = np.empty(n_frames)
        positions = np.empty((n_frames, n, d))
        velocities = np.empty((n_frames, n, d))
        kinetic = np.empty(n_frames)
        potential = np.empty(n_frames)
        temperature = np.empty(n_frames)

        def _record(idx):
            times[idx] = self.t
            positions[idx] = self.positions
            velocities[idx] = self.velocities
            kinetic[idx] = self.kinetic_energy()
            potential[idx] = self.potential_energy()
            temperature[idx] = self.temperature()

        _record(0)
        frame = 1
        for step_idx in range(1, n_steps + 1):
            self.step(dt)
            if thermostat is not None:
                thermostat.apply(self, dt)
            if step_idx % sample_every == 0:
                _record(frame)
                frame += 1
        return MDResult(
            t=times,
            positions=positions,
            velocities=velocities,
            kinetic_energy=kinetic,
            potential_energy=potential,
            temperature=temperature,
            box_length=self.box_length,
        )
