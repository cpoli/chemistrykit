r"""The Lennard-Jones fluid in reduced units.

The workhorse test system of molecular simulation: identical particles
interacting through the 12-6 Lennard-Jones potential
(:class:`LennardJones`), integrated in a cubic periodic box
(:class:`LJFluid`). Following universal convention (Allen & Tildesley,
*Computer Simulation of Liquids*, 2nd ed., Ch. 1.4 and Appendix B; Frenkel
& Smit, *Understanding Molecular Simulation*, 2nd ed., Ch. 3), everything
here is expressed in *reduced units*: lengths in units of :math:`\sigma`,
energies in units of :math:`\epsilon`, masses in units of the particle
mass `m`, and :math:`k_B=1` -- so temperature, pressure, and time are all
dimensionless combinations of :math:`\epsilon,\sigma,m` too (e.g. time in
units of :math:`\sigma\sqrt{m/\epsilon}`). A real substance's properties
are recovered by multiplying back through its own :math:`\epsilon,\sigma,m`.
"""

from __future__ import annotations

from typing import Optional

import numpy as np
from numba import njit

from chemistrykit.md.core.base_system import MolecularDynamicsSystem, PairPotential
from chemistrykit.md.utils.neighbor_list import VerletNeighborList, build_neighbor_list
from chemistrykit.md.utils.pbc import minimum_image_displacement, wrap_positions

__all__ = ["LennardJones", "LJFluid"]


class LennardJones(PairPotential):
    r"""The 12-6 Lennard-Jones pair potential.

    .. math::

        U(r) = 4\epsilon\left[\left(\frac{\sigma}{r}\right)^{12} - \left(\frac{\sigma}{r}\right)^{6}\right]

    The :math:`r^{-12}` term models short-range Pauli repulsion (steeply,
    for computational convenience -- it has no deeper physical
    justification), and the :math:`r^{-6}` term is the genuine
    London-dispersion attraction (J. E. Lennard-Jones, *Proc. R. Soc. A*
    106, 463 (1924); Allen & Tildesley, *Computer Simulation of Liquids*,
    2nd ed., Ch. 1.4).

    Parameters
    ----------
    epsilon : float, default 1.0
        Well depth.
    sigma : float, default 1.0
        Finite distance at which ``U(sigma) == 0``.

    Examples
    --------
    The potential has its minimum exactly at :math:`r=2^{1/6}\sigma`
    (:attr:`r_min`), where :math:`U=-\epsilon` and the force vanishes --
    the mechanical-equilibrium separation of an isolated pair:

    >>> lj = LennardJones(epsilon=1.0, sigma=1.0)
    >>> round(float(lj.energy(lj.r_min)), 10)
    -1.0
    >>> round(float(lj.force_scalar(lj.r_min)), 10)
    -0.0

    At ``r=sigma`` the potential is exactly zero by construction:

    >>> round(float(lj.energy(1.0)), 10)
    0.0
    """

    def __init__(self, epsilon: float = 1.0, sigma: float = 1.0):
        if epsilon <= 0 or sigma <= 0:
            raise ValueError("epsilon and sigma must be positive")
        self.epsilon = float(epsilon)
        self.sigma = float(sigma)

    def energy(self, r):
        r = np.asarray(r, dtype=np.float64)
        sr6 = (self.sigma / r) ** 6
        return 4.0 * self.epsilon * (sr6**2 - sr6)

    def force_scalar(self, r):
        r = np.asarray(r, dtype=np.float64)
        sr6 = (self.sigma / r) ** 6
        return 24.0 * self.epsilon / r * (2.0 * sr6**2 - sr6)

    @property
    def r_min(self) -> float:
        r"""Separation at the potential minimum, :math:`2^{1/6}\sigma`."""
        return 2.0 ** (1.0 / 6.0) * self.sigma


def _pack_neighbor_params(pairs_i: np.ndarray, pairs_j: np.ndarray) -> np.ndarray:
    """Pack a variable-length pair list into the fixed-signature ``params`` vector.

    ``params = [n_pairs, i_0, i_1, ..., i_{n-1}, j_0, j_1, ..., j_{n-1}]``,
    all stored as ``float64`` (indices are small enough to round-trip
    exactly). Keeping the pair *count* variable while the njit
    acceleration function's signature stays fixed is what lets
    :meth:`LJFluid.step` rebuild the neighbor list every few steps
    without forcing Numba to recompile a new dispatcher each time (see
    :func:`_make_lj_accel_njit`).

    Parameters
    ----------
    pairs_i, pairs_j : ndarray of int

    Returns
    -------
    ndarray of float64, shape (2 * n_pairs + 1,)
    """
    n_pairs = pairs_i.shape[0]
    head = np.array([float(n_pairs)])
    return np.concatenate([head, pairs_i.astype(np.float64), pairs_j.astype(np.float64)])


def _make_lj_accel_njit(epsilon: float, sigma: float, mass: float, box_length: float, cutoff: float):
    """Build a standalone ``@njit`` LJ acceleration function over a neighbor-pair list.

    Factory that closes over the fluid's *fixed* scalar parameters
    (``epsilon``, ``sigma``, ``mass``, ``box_length``, ``cutoff``) and
    returns a genuine ``@njit`` dispatcher with the
    :data:`chemistrykit.integrators.RHSFunc` signature
    ``accel(positions, t, params) -> acceleration``, usable as
    :attr:`~chemistrykit.md.core.base_system.MolecularDynamicsSystem._accel_njit`
    -- the same factory-closure pattern
    :mod:`chemistrykit.kinetics.systems.networks` uses for its mass-action
    right-hand side. The *variable-length* neighbor-pair list is instead
    threaded through the ``params`` argument at call time (see
    :func:`_pack_neighbor_params`), so rebuilding the neighbor list never
    requires recompiling this closure.

    Parameters
    ----------
    epsilon, sigma : float
        Lennard-Jones parameters.
    mass : float
        Per-particle mass (uniform -- a monatomic fluid of identical
        particles).
    box_length : float
        Cubic periodic box side length.
    cutoff : float
        Interaction cutoff.

    Returns
    -------
    callable
        ``@njit`` dispatcher ``accel(positions, t, params) -> acceleration``.
    """
    sigma2 = sigma * sigma
    cutoff2 = cutoff * cutoff

    @njit(cache=True)
    def accel(positions, t, params):
        n_pairs = int(params[0])
        n, d = positions.shape
        forces = np.zeros((n, d))
        diff = np.empty(d)
        for k in range(n_pairs):
            i = int(params[1 + k])
            j = int(params[1 + n_pairs + k])
            r2 = 0.0
            for m in range(d):
                dm = positions[i, m] - positions[j, m]
                if box_length > 0.0:
                    dm -= box_length * round(dm / box_length)
                diff[m] = dm
                r2 += dm * dm
            if r2 < cutoff2:
                sr2 = sigma2 / r2
                sr6 = sr2 * sr2 * sr2
                sr12 = sr6 * sr6
                f_over_r = 24.0 * epsilon * (2.0 * sr12 - sr6) / r2
                for m in range(d):
                    fm = f_over_r * diff[m]
                    forces[i, m] += fm
                    forces[j, m] -= fm
        return forces / mass

    return accel


class LJFluid(MolecularDynamicsSystem):
    r"""A monatomic Lennard-Jones fluid in a cubic periodic box, in reduced units.

    Positions evolve under :func:`chemistrykit.integrators.velocity_verlet_step`
    driven by a ``@njit`` pairwise-force loop (:func:`_make_lj_accel_njit`)
    that reads its interacting pairs from a periodically-rebuilt
    :class:`~chemistrykit.md.utils.neighbor_list.VerletNeighborList`, under
    the minimum-image convention
    (:func:`chemistrykit.md.utils.pbc.minimum_image_displacement`). See the
    module docstring for the reduced-units convention.

    Parameters
    ----------
    positions : array-like, shape (N, 3)
    velocities : array-like, shape (N, 3)
    box_length : float
        Cubic periodic box side length; must be at least ``2 * cutoff``
        for the minimum-image convention to be unambiguous.
    epsilon, sigma, mass : float, default 1.0
        Lennard-Jones/particle parameters (all particles identical).
    cutoff : float, default 2.5
        Interaction cutoff, in units of `sigma` (2.5 is the conventional
        standard truncation for a liquid-density LJ fluid; Allen &
        Tildesley, *Computer Simulation of Liquids*, 2nd ed., Ch. 1.4).
    skin : float, default 0.3
        Neighbor-list skin distance, see :class:`~chemistrykit.md.utils.neighbor_list.VerletNeighborList`.
    rebuild_every : int, default 20
        Rebuild the neighbor list after this many :meth:`step` calls.
    """

    def __init__(
        self,
        positions,
        velocities,
        box_length: float,
        epsilon: float = 1.0,
        sigma: float = 1.0,
        mass: float = 1.0,
        cutoff: float = 2.5,
        skin: float = 0.3,
        rebuild_every: int = 20,
    ):
        if box_length <= 2.0 * cutoff:
            raise ValueError("box_length must be at least twice the cutoff for the minimum-image convention to be valid")
        self.epsilon = float(epsilon)
        self.sigma = float(sigma)
        self.cutoff = float(cutoff)
        self.box_length = float(box_length)
        n = np.asarray(positions).shape[0]
        masses = np.full(n, float(mass))
        self._neighbor_list = VerletNeighborList(cutoff=self.cutoff, skin=skin, rebuild_every=rebuild_every)
        self._accel_njit = _make_lj_accel_njit(self.epsilon, self.sigma, float(mass), self.box_length, self.cutoff)
        # Truncated-and-shifted potential (Allen & Tildesley Ch. 2.4): subtracting
        # U(cutoff) makes the *reported* energy continuous at the cutoff, without
        # changing the force (a constant shift has zero derivative). The force
        # itself still has a small discontinuity at r=cutoff (this is the
        # "shifted potential", not the (stiffer) "shifted force" scheme), which is
        # the standard, small source of energy drift in a plain truncated LJ MD run.
        self._u_shift = 4.0 * self.epsilon * ((self.sigma / self.cutoff) ** 12 - (self.sigma / self.cutoff) ** 6)
        super().__init__(positions, velocities, masses)
        self.positions = wrap_positions(self.positions, self.box_length)
        pairs_i, pairs_j = self._neighbor_list.pairs(self.positions, self.box_length)
        self.params = _pack_neighbor_params(pairs_i, pairs_j)

    def step(self, dt: float) -> None:
        pairs_i, pairs_j = self._neighbor_list.pairs(self.positions, self.box_length)
        self.params = _pack_neighbor_params(pairs_i, pairs_j)
        super().step(dt)
        self.positions = wrap_positions(self.positions, self.box_length)

    def forces_and_potential(self, positions):
        """Exact (full, non-skin-list) forces, energy, and (cached) virial.

        Unlike :meth:`step`'s force evaluation (which reuses a
        skin-enlarged neighbor list for performance), this always rebuilds
        an exact-cutoff pair list from the given `positions`, so energy,
        pressure, and g(r) queries are never affected by neighbor-list
        staleness.

        Parameters
        ----------
        positions : ndarray, shape (N, 3)

        Returns
        -------
        forces : ndarray, shape (N, 3)
        potential_energy : float
        """
        positions = np.asarray(positions, dtype=np.float64)
        forces = np.zeros_like(positions)
        pairs_i, pairs_j = build_neighbor_list(positions, self.box_length, self.cutoff)
        if pairs_i.size == 0:
            self._last_virial = 0.0
            return forces, 0.0
        diff = minimum_image_displacement(positions[pairs_i], positions[pairs_j], self.box_length)
        r2 = np.sum(diff * diff, axis=1)
        sr2 = (self.sigma**2) / r2
        sr6 = sr2**3
        sr12 = sr6**2
        energies = 4.0 * self.epsilon * (sr12 - sr6) - self._u_shift
        f_over_r = 24.0 * self.epsilon * (2.0 * sr12 - sr6) / r2
        contributions = f_over_r[:, None] * diff
        np.add.at(forces, pairs_i, contributions)
        np.add.at(forces, pairs_j, -contributions)
        self._last_virial = float(np.sum(f_over_r * r2))
        return forces, float(np.sum(energies))

    def virial(self) -> float:
        r"""The virial sum :math:`W=\sum_{i<j}\vec r_{ij}\cdot\vec F_{ij}` at the current configuration.

        Returns
        -------
        float
        """
        self.forces_and_potential(self.positions)
        return self._last_virial

    def pressure(self, k_b: float = 1.0) -> float:
        r"""Instantaneous pressure via the virial theorem.

        .. math::

            P = \frac{N k_B T}{V} + \frac{W}{dV}

        (Allen & Tildesley, *Computer Simulation of Liquids*, 2nd ed., eq.
        2.60, generalized from :math:`d=3` to general dimension.)

        Parameters
        ----------
        k_b : float, default 1.0

        Returns
        -------
        float

        Examples
        --------
        A very dilute LJ gas (large box, few particles) has pressure
        close to the ideal-gas value :math:`P=\rho k_B T`, since the
        virial correction from interactions vanishes as the density does:

        >>> import numpy as np
        >>> rng = np.random.default_rng(0)
        >>> n, L = 20, 40.0
        >>> positions = rng.uniform(0.0, L, size=(n, 3))
        >>> velocities = rng.normal(0.0, 1.0, size=(n, 3))
        >>> velocities -= velocities.mean(axis=0)
        >>> fluid = LJFluid(positions, velocities, box_length=L, cutoff=2.5)
        >>> fluid.virial()  # no pairs within cutoff at this density
        0.0
        >>> rho = n / L**3
        >>> ideal_pressure = rho * fluid.temperature()
        >>> bool(abs(fluid.pressure() - ideal_pressure) < 1e-12)
        True
        """
        n, d = self.positions.shape
        volume = self.box_length**d
        w = self.virial()
        return (n * k_b * self.temperature(k_b) + w / d) / volume

    def radial_distribution_function(self, r_max: Optional[float] = None, n_bins: int = 100, positions=None):
        r"""Compute the radial distribution function g(r) for a configuration.

        .. math::

            g(r) = \frac{2\,\langle n_{\text{pairs}}(r, r+dr)\rangle}{N\,\rho\,4\pi r^2 dr}

        the (3D) local density of particles at distance `r` from a
        reference particle, relative to the bulk density :math:`\rho`
        (Allen & Tildesley, *Computer Simulation of Liquids*, 2nd ed., Ch.
        2.6; Frenkel & Smit, *Understanding Molecular Simulation*, 2nd
        ed., Ch. 3.6). :math:`g(r) \to 1` for an ideal (non-interacting,
        spatially uniform) gas, and develops structure (a first peak near
        the LJ minimum, decaying oscillations beyond it) in a liquid or
        solid.

        Parameters
        ----------
        r_max : float, optional
            Maximum separation to histogram; defaults to half the box
            length (beyond which the minimum-image convention is
            ambiguous).
        n_bins : int, default 100
        positions : ndarray, shape (N, 3), optional
            Configuration to analyze; defaults to the system's current
            ``self.positions``.

        Returns
        -------
        r : ndarray, shape (n_bins,)
            Bin-center separations.
        g : ndarray, shape (n_bins,)
        """
        positions = self.positions if positions is None else np.asarray(positions, dtype=np.float64)
        n = positions.shape[0]
        box_length = self.box_length
        r_max = 0.5 * box_length if r_max is None else r_max
        pairs_i, pairs_j = build_neighbor_list(positions, box_length, r_max)
        counts = np.zeros(n_bins)
        edges = np.linspace(0.0, r_max, n_bins + 1)
        if pairs_i.size > 0:
            diff = minimum_image_displacement(positions[pairs_i], positions[pairs_j], box_length)
            r = np.sqrt(np.sum(diff * diff, axis=1))
            counts, edges = np.histogram(r, bins=n_bins, range=(0.0, r_max))
        r_mid = 0.5 * (edges[:-1] + edges[1:])
        dr = edges[1] - edges[0]
        volume = box_length**3
        rho = n / volume
        shell_volume = 4.0 * np.pi * r_mid**2 * dr
        n_ideal = rho * shell_volume
        with np.errstate(divide="ignore", invalid="ignore"):
            g = np.where(n_ideal > 0, (2.0 * counts) / (n * n_ideal), 0.0)
        return r_mid, g

    @classmethod
    def from_lattice(
        cls,
        n_per_side: int,
        density: float,
        temperature: float,
        epsilon: float = 1.0,
        sigma: float = 1.0,
        mass: float = 1.0,
        cutoff: float = 2.5,
        rng=None,
        **kwargs,
    ) -> LJFluid:
        r"""Build an :class:`LJFluid` on a simple-cubic lattice with Maxwell-Boltzmann velocities.

        A standard way to initialize an MD run: a regular lattice avoids
        any risk of initial particle overlap (which would otherwise blow
        up the steeply repulsive :math:`r^{-12}` term), and velocity
        components drawn independently from a Gaussian with variance
        :math:`k_BT/m` reproduce the Maxwell-Boltzmann distribution
        exactly (see
        :class:`chemistrykit.statmech.MaxwellBoltzmannSpeedDistribution`);
        the total momentum is then zeroed and the velocities rescaled so
        the realized instantaneous temperature matches `temperature`
        exactly, rather than only on average.

        Parameters
        ----------
        n_per_side : int
            Particles per lattice edge; total particle count is
            ``n_per_side**3``.
        density : float
            Number density :math:`\rho=N/V`, in units of :math:`\sigma^{-3}`.
        temperature : float
            Target temperature, in units of :math:`\epsilon/k_B`.
        epsilon, sigma, mass : float, default 1.0
        cutoff : float, default 2.5
        rng : int, numpy.random.Generator, or None
            Seed or generator for the initial velocities.
        **kwargs
            Forwarded to the :class:`LJFluid` constructor (e.g. ``skin``,
            ``rebuild_every``).

        Returns
        -------
        LJFluid

        Examples
        --------
        >>> fluid = LJFluid.from_lattice(n_per_side=4, density=0.5, temperature=1.0, rng=0)
        >>> fluid.positions.shape
        (64, 3)
        >>> round(fluid.temperature(), 6)
        1.0
        """
        rng = np.random.default_rng(rng)
        n = n_per_side**3
        volume = n / density
        box_length = volume ** (1.0 / 3.0)
        spacing = box_length / n_per_side
        idx = np.arange(n_per_side)
        xx, yy, zz = np.meshgrid(idx, idx, idx, indexing="ij")
        positions = np.stack([xx.ravel(), yy.ravel(), zz.ravel()], axis=1).astype(np.float64) * spacing
        v_scale = np.sqrt(temperature / mass)
        velocities = rng.normal(0.0, v_scale, size=(n, 3))
        velocities -= velocities.mean(axis=0)
        system = cls(positions, velocities, box_length, epsilon=epsilon, sigma=sigma, mass=mass, cutoff=cutoff, **kwargs)
        current_temperature = system.temperature()
        if current_temperature > 0:
            system.velocities *= np.sqrt(temperature / current_temperature)
        return system
