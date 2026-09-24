r"""Thermostats for canonical-ensemble (NVT) molecular dynamics.

A :class:`~chemistrykit.md.core.base_system.MolecularDynamicsSystem`
integrated by plain velocity-Verlet alone conserves total energy (an NVE,
microcanonical run). A thermostat modifies the velocities after each step
to instead sample (approximately, for the simplified schemes below) a
canonical (NVT, fixed-temperature) ensemble. Both classes here implement
``apply(system, dt)``, called once by
:meth:`chemistrykit.md.core.base_system.MolecularDynamicsSystem.run` after
every velocity-Verlet step.
"""

from __future__ import annotations

import numpy as np

__all__ = [
    "VelocityRescalingThermostat",
    "NoseHooverThermostat",
    "BerendsenThermostat",
    "StochasticVelocityRescalingThermostat",
]


class VelocityRescalingThermostat:
    r"""Simple (deterministic) velocity-rescaling thermostat.

    Rescales every particle's velocity by a common factor
    :math:`\lambda=\sqrt{T_{\text{target}}/T_{\text{current}}}` every
    `interval` MD steps, so the instantaneous kinetic temperature matches
    `target_temperature` exactly at each rescaling (Allen & Tildesley,
    *Computer Simulation of Liquids*, 2nd ed., Ch. 2.3; Frenkel & Smit,
    *Understanding Molecular Simulation*, 2nd ed., Ch. 6.1.1).

    This is the crude, deterministic member of the "velocity rescaling"
    thermostat family: it does *not* itself sample the canonical (NVT)
    ensemble correctly (the kinetic-energy fluctuations of a true
    canonical ensemble are suppressed to zero at every rescaling step),
    unlike the stochastic Bussi-Donadio-Parrinello (2007) variant, or (to
    leading order, with a finite relaxation time) the Berendsen (1984)
    weak-coupling thermostat. It is nonetheless the standard first
    introduction to temperature control in MD and is exactly what is
    meant by "velocity rescaling" here.

    Parameters
    ----------
    target_temperature : float
        Must be positive.
    interval : int, default 1
        Rescale every this many calls to :meth:`apply` (i.e. MD steps).
    k_b : float, default 1.0
        Boltzmann constant in the system's unit system (``1.0`` for the
        reduced units used by :class:`~chemistrykit.md.systems.lj_fluid.LJFluid`).

    Examples
    --------
    >>> import numpy as np
    >>> from chemistrykit.md.systems.lj_fluid import LJFluid
    >>> fluid = LJFluid.from_lattice(n_per_side=3, cutoff=1.5, density=0.6, temperature=2.0, rng=0)
    >>> fluid.velocities *= 2.0  # perturb away from the target temperature
    >>> thermostat = VelocityRescalingThermostat(target_temperature=2.0)
    >>> result = fluid.run(dt=0.001, n_steps=50, thermostat=thermostat, sample_every=10)
    >>> bool(abs(result.temperature[-1] - 2.0) < 0.5)
    True
    """

    def __init__(self, target_temperature: float, interval: int = 1, k_b: float = 1.0):
        if target_temperature <= 0:
            raise ValueError("target_temperature must be positive")
        self.target_temperature = float(target_temperature)
        self.interval = int(interval)
        self.k_b = float(k_b)
        self._calls = 0

    def apply(self, system, dt: float) -> None:
        """Rescale ``system.velocities`` in place if a rescaling is due.

        Parameters
        ----------
        system : chemistrykit.md.core.base_system.MolecularDynamicsSystem
        dt : float
            Unused (present for interface symmetry with
            :class:`NoseHooverThermostat`).
        """
        self._calls += 1
        if self._calls % self.interval != 0:
            return
        current_temperature = system.temperature(self.k_b)
        if current_temperature > 0:
            system.velocities *= np.sqrt(self.target_temperature / current_temperature)


class NoseHooverThermostat:
    r"""The Nose-Hoover extended-system thermostat.

    Couples the physical system to a fictitious "friction" degree of
    freedom :math:`\xi` whose equation of motion drives the
    time-averaged kinetic temperature to `target_temperature` while
    (in the exact, continuous-time extended Lagrangian) still generating
    genuine canonical-ensemble fluctuations -- S. Nose, *J. Chem. Phys.*
    81, 511 (1984); W. G. Hoover, *Phys. Rev. A* 31, 1695 (1985).

    .. math::

        \dot{\vec v}_i = \frac{\vec F_i}{m_i} - \xi \vec v_i, \qquad
        Q\dot\xi = \left(\sum_i m_i v_i^2\right) - \text{dof}\cdot k_BT

    This implementation integrates the friction/velocity-scaling part
    with the explicit, first-order update

    .. math::

        \xi \mathrel{+}= \frac{2K - \text{dof}\cdot k_BT}{Q}\,dt, \qquad
        \vec v_i \mathrel{*}= e^{-\xi\,dt}

    applied once after each ordinary velocity-Verlet step -- an
    approximation to the reference operator-split ("Trotter
    factorization") integrator of Martyna, Tuckerman, Tobias & Klein,
    *Mol. Phys.* 87, 1117 (1996) (see also Frenkel & Smit, *Understanding
    Molecular Simulation*, 2nd ed., Ch. 6.1.2, Algorithm 15), which
    instead splits the thermostat update into two half-steps
    symmetrically wrapped around the position/velocity update for
    improved (time-reversible) long-time behavior. The simplified,
    single-step coupling used here still relaxes the time-averaged
    kinetic temperature to `target_temperature`, but without that
    scheme's more favorable long-time energy-conservation properties --
    flagged explicitly here as the approximation it is.

    Parameters
    ----------
    target_temperature : float
        Must be positive.
    Q : float
        Thermostat "mass" (coupling strength); a larger `Q` couples more
        weakly, giving a slower thermostat response. Must be positive.
    dof : int, optional
        Degrees of freedom used for the target kinetic energy
        :math:`\text{dof}\cdot k_BT/2`; defaults to the system's own
        :meth:`~chemistrykit.md.core.base_system.MolecularDynamicsSystem.degrees_of_freedom`
        at the first call to :meth:`apply`.
    k_b : float, default 1.0

    Examples
    --------
    >>> import numpy as np
    >>> from chemistrykit.md.systems.lj_fluid import LJFluid
    >>> fluid = LJFluid.from_lattice(n_per_side=3, cutoff=1.5, density=0.6, temperature=2.0, rng=0)
    >>> fluid.velocities *= 2.0  # perturb away from the target temperature
    >>> thermostat = NoseHooverThermostat(target_temperature=2.0, Q=10.0)
    >>> result = fluid.run(dt=0.001, n_steps=200, thermostat=thermostat, sample_every=20)
    >>> bool(abs(result.temperature[-1] - 2.0) < abs(result.temperature[0] - 2.0))
    True
    """

    def __init__(self, target_temperature: float, Q: float, dof: int = None, k_b: float = 1.0):
        if target_temperature <= 0 or Q <= 0:
            raise ValueError("target_temperature and Q must be positive")
        self.target_temperature = float(target_temperature)
        self.Q = float(Q)
        self.dof = dof
        self.k_b = float(k_b)
        self.xi = 0.0

    def apply(self, system, dt: float) -> None:
        """Advance the friction coefficient and rescale ``system.velocities`` in place.

        Parameters
        ----------
        system : chemistrykit.md.core.base_system.MolecularDynamicsSystem
        dt : float
        """
        if self.dof is None:
            self.dof = system.degrees_of_freedom()
        kinetic_energy = system.kinetic_energy()
        g = (2.0 * kinetic_energy - self.dof * self.k_b * self.target_temperature) / self.Q
        self.xi += g * dt
        system.velocities *= np.exp(-self.xi * dt)


class BerendsenThermostat:
    r"""The Berendsen weak-coupling thermostat.

    Rescales every velocity by a common factor after each step so that the
    instantaneous kinetic temperature relaxes exponentially toward
    `target_temperature` with time constant :math:`\tau` -- H. J. C.
    Berendsen, J. P. M. Postma, W. F. van Gunsteren, A. DiNola & J. R. Haak,
    *J. Chem. Phys.* 81, 3684 (1984), eq. 11:

    .. math::

        \lambda^2 = 1 + \frac{\Delta t}{\tau}\left(\frac{T_0}{T} - 1\right),
        \qquad \frac{dT}{dt} = \frac{T_0 - T}{\tau}.

    For :math:`\tau = \Delta t` this reduces to the immediate, complete
    rescaling of :class:`VelocityRescalingThermostat`. Like that scheme it
    suppresses canonical kinetic-energy fluctuations, so it does not sample
    the NVT ensemble exactly (see :class:`StochasticVelocityRescalingThermostat`).

    Parameters
    ----------
    target_temperature : float
        :math:`T_0`; must be positive.
    tau : float
        Coupling time constant; must be positive.
    k_b : float, default 1.0

    Examples
    --------
    With no forces acting, each call moves the temperature a fraction
    :math:`\Delta t/\tau` of the way to the target:

    >>> from chemistrykit.md.systems.lj_fluid import LJFluid
    >>> fluid = LJFluid.from_lattice(n_per_side=3, cutoff=1.5, density=0.6, temperature=2.0, rng=0)
    >>> thermostat = BerendsenThermostat(target_temperature=1.0, tau=0.1)
    >>> thermostat.apply(fluid, dt=0.01)
    >>> round(fluid.temperature(), 10)  # 2.0 + 0.1 * (1.0 - 2.0)
    1.9
    """

    def __init__(self, target_temperature: float, tau: float, k_b: float = 1.0):
        if target_temperature <= 0 or tau <= 0:
            raise ValueError("target_temperature and tau must be positive")
        self.target_temperature = float(target_temperature)
        self.tau = float(tau)
        self.k_b = float(k_b)

    def apply(self, system, dt: float) -> None:
        """Rescale ``system.velocities`` in place by the Berendsen factor.

        Parameters
        ----------
        system : chemistrykit.md.core.base_system.MolecularDynamicsSystem
        dt : float
        """
        current_temperature = system.temperature(self.k_b)
        if current_temperature <= 0:
            return
        lam2 = 1.0 + dt / self.tau * (self.target_temperature / current_temperature - 1.0)
        system.velocities *= np.sqrt(max(lam2, 0.0))


class StochasticVelocityRescalingThermostat:
    r"""The Bussi-Donadio-Parrinello stochastic velocity-rescaling thermostat.

    Berendsen's exponential relaxation of the kinetic energy :math:`K`
    plus a correctly scaled stochastic term, chosen so that the kinetic
    energy samples its exact canonical distribution -- G. Bussi, D.
    Donadio & M. Parrinello, *J. Chem. Phys.* 126, 014101 (2007), eq. 7
    and Appendix (eq. A7):

    .. math::

        \alpha^2 = c + (1-c)\frac{\bar K}{N_f K}\sum_{i=1}^{N_f} R_i^2
        + 2 R_1 \sqrt{c(1-c)\frac{\bar K}{N_f K}},
        \qquad c = e^{-\Delta t/\tau},

    where :math:`\bar K = N_f k_B T_0/2`, the :math:`R_i` are independent
    standard normal numbers (the sum over :math:`i \ge 2` is drawn as a
    single :math:`\chi^2_{N_f-1}` variate), and every velocity is scaled
    by :math:`\alpha`. The stationary distribution of :math:`K` is the
    canonical gamma distribution with mean :math:`\bar K` and variance
    :math:`2\bar K^2/N_f`.

    Parameters
    ----------
    target_temperature : float
        :math:`T_0`; must be positive.
    tau : float
        Relaxation time; must be positive.
    rng : int, numpy.random.Generator, or None
    dof : int, optional
        :math:`N_f`; defaults to the system's
        :meth:`~chemistrykit.md.core.base_system.MolecularDynamicsSystem.degrees_of_freedom`.
    k_b : float, default 1.0

    Examples
    --------
    >>> from chemistrykit.md.systems.lj_fluid import LJFluid
    >>> fluid = LJFluid.from_lattice(n_per_side=3, cutoff=1.5, density=0.6, temperature=3.0, rng=0)
    >>> thermostat = StochasticVelocityRescalingThermostat(target_temperature=1.0, tau=0.01, rng=0)
    >>> for _ in range(200):
    ...     thermostat.apply(fluid, dt=0.01)
    >>> bool(0.5 < fluid.temperature() < 1.5)
    True
    """

    def __init__(self, target_temperature: float, tau: float, rng=None, dof: int = None, k_b: float = 1.0):
        if target_temperature <= 0 or tau <= 0:
            raise ValueError("target_temperature and tau must be positive")
        self.target_temperature = float(target_temperature)
        self.tau = float(tau)
        self.rng = np.random.default_rng(rng)
        self.dof = dof
        self.k_b = float(k_b)

    def apply(self, system, dt: float) -> None:
        """Rescale ``system.velocities`` in place by a stochastic factor.

        Parameters
        ----------
        system : chemistrykit.md.core.base_system.MolecularDynamicsSystem
        dt : float
        """
        if self.dof is None:
            self.dof = system.degrees_of_freedom()
        nf = self.dof
        kinetic_energy = system.kinetic_energy()
        if kinetic_energy <= 0:
            return
        k_target = 0.5 * nf * self.k_b * self.target_temperature
        c = np.exp(-dt / self.tau)
        r1 = self.rng.standard_normal()
        sum_r2 = r1 * r1 + (self.rng.chisquare(nf - 1) if nf > 1 else 0.0)
        ratio = k_target / (nf * kinetic_energy)
        alpha2 = c + (1.0 - c) * ratio * sum_r2 + 2.0 * r1 * np.sqrt(c * (1.0 - c) * ratio)
        sign = 1.0 if r1 + np.sqrt(c / ((1.0 - c) * ratio)) >= 0 else -1.0
        system.velocities *= sign * np.sqrt(max(alpha2, 0.0))
