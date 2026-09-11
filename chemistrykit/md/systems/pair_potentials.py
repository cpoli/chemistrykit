r"""Additional pairwise/angle potentials, and small bonded-force-field systems.

:class:`Morse` and :class:`Buckingham` are alternative *non-bonded*
radial potentials (like :class:`~chemistrykit.md.systems.lj_fluid.LennardJones`,
but with different repulsive/attractive functional forms);
:class:`HarmonicBond` and :class:`HarmonicAngle` are the *bonded* terms
used to hold a molecule's internal geometry together in a classical force
field. :class:`DiatomicOscillator` and :class:`HarmonicMolecule` are
small, non-periodic :class:`~chemistrykit.md.core.base_system.MolecularDynamicsSystem`
demonstrations built from these terms.
"""

from __future__ import annotations

from collections.abc import Sequence

import numpy as np
from numba import njit

from chemistrykit.md.core.base_system import AnglePotential, MolecularDynamicsSystem, PairPotential

__all__ = [
    "Morse",
    "Buckingham",
    "HarmonicBond",
    "HarmonicAngle",
    "DiatomicOscillator",
    "HarmonicMolecule",
]


class Morse(PairPotential):
    r"""The Morse interatomic potential.

    .. math::

        U(r) = D_e\left[1 - e^{-a(r-r_e)}\right]^2 - D_e

    A three-parameter anharmonic bond potential: dissociation energy
    :math:`D_e`, equilibrium separation :math:`r_e`, and a width
    parameter `a` controlling the curvature at the minimum (related to
    the harmonic force constant there by :math:`k = 2D_ea^2`, see
    :attr:`force_constant`) -- P. M. Morse, *Phys. Rev.* 34, 57 (1929);
    Atkins & de Paula, *Physical Chemistry*, 11th ed., Ch. 12 (used there
    for vibrational anharmonicity; see also
    :class:`~chemistrykit.md.systems.pair_potentials.DiatomicOscillator`,
    whose classical vibration this potential drives).

    Parameters
    ----------
    De : float
        Dissociation energy (well depth), must be positive.
    a : float
        Width parameter (inverse length), must be positive.
    re : float
        Equilibrium separation.

    Examples
    --------
    Exactly reproduces the well depth and zero force at the minimum:

    >>> morse = Morse(De=4.0, a=1.5, re=1.0)
    >>> round(float(morse.energy(morse.re)), 10)
    -4.0
    >>> round(float(morse.force_scalar(morse.re)), 10)
    0.0

    And the dissociation limit :math:`U(r\to\infty) \to 0`:

    >>> round(float(morse.energy(20.0)), 6)
    -0.0
    """

    def __init__(self, De: float, a: float, re: float):
        if De <= 0 or a <= 0:
            raise ValueError("De and a must be positive")
        self.De = float(De)
        self.a = float(a)
        self.re = float(re)

    def energy(self, r):
        r = np.asarray(r, dtype=np.float64)
        x = np.exp(-self.a * (r - self.re))
        return self.De * (1.0 - x) ** 2 - self.De

    def force_scalar(self, r):
        r = np.asarray(r, dtype=np.float64)
        x = np.exp(-self.a * (r - self.re))
        return 2.0 * self.a * self.De * x * (x - 1.0)

    @property
    def force_constant(self) -> float:
        r"""Harmonic force constant at the minimum, :math:`k=2D_ea^2`."""
        return 2.0 * self.De * self.a**2

    @classmethod
    def from_force_constant(cls, De: float, k: float, re: float) -> Morse:
        r"""Build a :class:`Morse` potential from :math:`D_e`, a harmonic force constant `k`, and `re`.

        Inverts :math:`k=2D_ea^2` for :math:`a=\sqrt{k/(2D_e)}`, so the
        curvature at the minimum matches a given harmonic bond exactly
        while still capturing anharmonicity/dissociation away from it.

        Parameters
        ----------
        De : float
        k : float
        re : float

        Returns
        -------
        Morse

        Examples
        --------
        >>> morse = Morse.from_force_constant(De=4.0, k=10.0, re=1.0)
        >>> round(morse.force_constant, 6)
        10.0
        """
        a = np.sqrt(k / (2.0 * De))
        return cls(De=De, a=a, re=re)


class Buckingham(PairPotential):
    r"""The Buckingham (exp-6) pair potential.

    .. math::

        U(r) = A e^{-Br} - \frac{C}{r^6}

    An alternative to Lennard-Jones for the short-range repulsion,
    replacing the :math:`r^{-12}` term with an exponential (closer to the
    true quantum-mechanical electron-overlap repulsion) at the cost of an
    unphysical turnover to :math:`U\to-\infty` as :math:`r\to0` that pure
    LJ does not have -- R. A. Buckingham, *Proc. R. Soc. A* 168, 264
    (1938); used extensively for rare-gas potentials and ionic-crystal
    lattice energies (see the future ``chemistrykit.crystal``).

    Parameters
    ----------
    A : float
        Repulsive prefactor (energy units), must be positive.
    B : float
        Repulsive exponential decay rate (inverse length), must be positive.
    C : float
        Dispersion (attractive :math:`r^{-6}`) coefficient (energy *
        length\ :sup:`6`), must be non-negative.

    Examples
    --------
    At large separation the potential is dominated by the (attractive,
    slower-decaying) dispersion term, so both the energy and the force
    are negative:

    >>> pot = Buckingham(A=1.0e4, B=3.0, C=1.0)
    >>> bool(pot.energy(10.0) < 0.0)
    True
    >>> bool(pot.force_scalar(10.0) < 0.0)
    True
    """

    def __init__(self, A: float, B: float, C: float):
        if A <= 0 or B <= 0:
            raise ValueError("A and B must be positive")
        if C < 0:
            raise ValueError("C must be non-negative")
        self.A = float(A)
        self.B = float(B)
        self.C = float(C)

    def energy(self, r):
        r = np.asarray(r, dtype=np.float64)
        return self.A * np.exp(-self.B * r) - self.C / r**6

    def force_scalar(self, r):
        r = np.asarray(r, dtype=np.float64)
        return self.A * self.B * np.exp(-self.B * r) - 6.0 * self.C / r**7


class HarmonicBond(PairPotential):
    r"""A harmonic (Hookean) bond-stretch potential.

    .. math::

        U(r) = \frac{1}{2}k(r-r_0)^2

    The simplest bonded term in a classical force field -- the small
    -oscillation limit of any real bond potential, including
    :class:`Morse` (Atkins & de Paula, *Physical Chemistry*, 11th ed.,
    Ch. 8, the classical harmonic oscillator).

    Parameters
    ----------
    k : float
        Force constant, must be positive.
    r0 : float
        Equilibrium bond length.

    Examples
    --------
    >>> bond = HarmonicBond(k=50.0, r0=1.0)
    >>> float(bond.energy(1.0))
    0.0
    >>> round(float(bond.force_scalar(1.1)), 10)
    -5.0
    """

    def __init__(self, k: float, r0: float):
        if k <= 0:
            raise ValueError("k must be positive")
        self.k = float(k)
        self.r0 = float(r0)

    def energy(self, r):
        r = np.asarray(r, dtype=np.float64)
        return 0.5 * self.k * (r - self.r0) ** 2

    def force_scalar(self, r):
        r = np.asarray(r, dtype=np.float64)
        return -self.k * (r - self.r0)


class HarmonicAngle(AnglePotential):
    r"""A harmonic bond-angle-bending potential.

    .. math::

        U(\theta) = \frac{1}{2}k_\theta(\theta-\theta_0)^2

    Parameters
    ----------
    k_theta : float
        Force constant (energy per radian squared), must be positive.
    theta0 : float
        Equilibrium angle, in radians.

    Examples
    --------
    >>> angle = HarmonicAngle(k_theta=100.0, theta0=1.9106)  # ~109.47 deg, tetrahedral
    >>> float(angle.energy(angle.theta0))
    0.0
    >>> float(angle.torque(angle.theta0))
    -0.0
    """

    def __init__(self, k_theta: float, theta0: float):
        if k_theta <= 0:
            raise ValueError("k_theta must be positive")
        self.k_theta = float(k_theta)
        self.theta0 = float(theta0)

    def energy(self, theta):
        theta = np.asarray(theta, dtype=np.float64)
        return 0.5 * self.k_theta * (theta - self.theta0) ** 2

    def torque(self, theta):
        theta = np.asarray(theta, dtype=np.float64)
        return -self.k_theta * (theta - self.theta0)


def _make_bond_accel_njit(kind_code: int, p0: float, p1: float, p2: float, m1: float, m2: float):
    """Build a standalone ``@njit`` acceleration function for a 2-body bond potential.

    Factory closing over the fixed bond parameters and a ``kind_code``
    (0 = :class:`HarmonicBond`, 1 = :class:`Morse`) selecting which
    closed-form radial force to evaluate, following the same
    factory-closure pattern as
    :func:`chemistrykit.md.systems.lj_fluid._make_lj_accel_njit`.

    Parameters
    ----------
    kind_code : int
    p0, p1, p2 : float
        Bond potential parameters: ``(k, r0, unused)`` for harmonic,
        ``(De, a, re)`` for Morse.
    m1, m2 : float
        Masses of the two atoms.

    Returns
    -------
    callable
    """

    @njit(cache=True)
    def accel(positions, t, params):
        dx = positions[0, 0] - positions[1, 0]
        dy = positions[0, 1] - positions[1, 1]
        dz = positions[0, 2] - positions[1, 2]
        r = np.sqrt(dx * dx + dy * dy + dz * dz)
        if kind_code == 0:
            f = -p0 * (r - p1)
        else:
            x = np.exp(-p1 * (r - p2))
            f = 2.0 * p1 * p0 * x * (x - 1.0)
        fx, fy, fz = f * dx / r, f * dy / r, f * dz / r
        acc = np.empty((2, 3))
        acc[0, 0] = fx / m1
        acc[0, 1] = fy / m1
        acc[0, 2] = fz / m1
        acc[1, 0] = -fx / m2
        acc[1, 1] = -fy / m2
        acc[1, 2] = -fz / m2
        return acc

    return accel


class DiatomicOscillator(MolecularDynamicsSystem):
    r"""Two atoms bound by a radial potential: classical bond-vibration MD.

    A textbook two-body problem: with no external forces, the
    center-of-mass moves uniformly and the relative coordinate oscillates
    as if a single particle of the reduced mass :math:`\mu=m_1m_2/(m_1+m_2)`
    were moving in the bond potential itself. For a
    :class:`~chemistrykit.md.systems.pair_potentials.HarmonicBond`, this
    gives the exact classical vibration period
    :math:`T=2\pi\sqrt{\mu/k}` (see :meth:`harmonic_period`), which the
    simulated trajectory is checked against in the test suite.

    Parameters
    ----------
    potential : HarmonicBond or Morse
        The bond potential.
    m1, m2 : float
        Atomic masses.
    r0 : float
        Initial bond length.
    v_rel0 : float, default 0.0
        Initial relative (radial) velocity along `axis`.
    axis : array-like, shape (3,), default (1, 0, 0)
        Bond direction.
    """

    def __init__(self, potential, m1: float, m2: float, r0: float, v_rel0: float = 0.0, axis=(1.0, 0.0, 0.0)):
        axis = np.asarray(axis, dtype=np.float64)
        axis = axis / np.linalg.norm(axis)
        total_mass = m1 + m2
        r1 = (m2 / total_mass) * r0 * axis
        r2 = -(m1 / total_mass) * r0 * axis
        v1 = (m2 / total_mass) * v_rel0 * axis
        v2 = -(m1 / total_mass) * v_rel0 * axis
        positions = np.array([r1, r2])
        velocities = np.array([v1, v2])
        if isinstance(potential, HarmonicBond):
            self._accel_njit = _make_bond_accel_njit(0, potential.k, potential.r0, 0.0, float(m1), float(m2))
        elif isinstance(potential, Morse):
            self._accel_njit = _make_bond_accel_njit(1, potential.De, potential.a, potential.re, float(m1), float(m2))
        else:
            raise TypeError("DiatomicOscillator supports HarmonicBond or Morse potentials")
        self.potential = potential
        self.params = np.empty(0)
        super().__init__(positions, velocities, [m1, m2])

    def forces_and_potential(self, positions):
        positions = np.asarray(positions, dtype=np.float64)
        diff = positions[0] - positions[1]
        r = float(np.linalg.norm(diff))
        f = float(self.potential.force_scalar(r))
        force_on_0 = f * diff / r
        forces = np.array([force_on_0, -force_on_0])
        return forces, float(self.potential.energy(r))

    def degrees_of_freedom(self) -> int:
        """A single relative-coordinate (bond-stretch) vibrational degree of freedom."""
        return 1

    def bond_length(self) -> float:
        """Current separation between the two atoms.

        Returns
        -------
        float
        """
        return float(np.linalg.norm(self.positions[0] - self.positions[1]))

    @staticmethod
    def harmonic_period(k: float, m1: float, m2: float) -> float:
        r"""Classical vibration period for a harmonic bond, :math:`T=2\pi\sqrt{\mu/k}`.

        Parameters
        ----------
        k : float
            Harmonic force constant.
        m1, m2 : float
            Atomic masses.

        Returns
        -------
        float

        Examples
        --------
        >>> round(float(DiatomicOscillator.harmonic_period(k=1.0, m1=1.0, m2=1.0)), 6)
        4.442883
        """
        mu = m1 * m2 / (m1 + m2)
        return 2.0 * np.pi * np.sqrt(mu / k)


def _make_harmonic_molecule_accel_njit(bonds: np.ndarray, angles: np.ndarray, masses: np.ndarray):
    """Build a standalone ``@njit`` acceleration function for harmonic bonds + angles.

    Factory closing over fixed structural arrays (bond and angle
    topology/parameters, and per-atom masses), following the same
    pattern as :func:`chemistrykit.kinetics.systems.networks._make_mass_action_rhs`.

    Parameters
    ----------
    bonds : ndarray, shape (n_bonds, 4)
        Rows ``(i, j, k, r0)``.
    angles : ndarray, shape (n_angles, 5)
        Rows ``(i, j, k, k_theta, theta0)``, angle at vertex `j`.
    masses : ndarray, shape (n_atoms,)

    Returns
    -------
    callable
    """
    bonds = np.ascontiguousarray(bonds, dtype=np.float64)
    angles = np.ascontiguousarray(angles, dtype=np.float64)
    masses = np.ascontiguousarray(masses, dtype=np.float64)
    n_bonds = bonds.shape[0]
    n_angles = angles.shape[0]

    @njit(cache=True)
    def accel(positions, t, params):
        n = positions.shape[0]
        forces = np.zeros((n, 3))
        for b in range(n_bonds):
            i = int(bonds[b, 0])
            j = int(bonds[b, 1])
            k = bonds[b, 2]
            r0 = bonds[b, 3]
            dx = positions[i, 0] - positions[j, 0]
            dy = positions[i, 1] - positions[j, 1]
            dz = positions[i, 2] - positions[j, 2]
            r = np.sqrt(dx * dx + dy * dy + dz * dz)
            f = -k * (r - r0)
            forces[i, 0] += f * dx / r
            forces[i, 1] += f * dy / r
            forces[i, 2] += f * dz / r
            forces[j, 0] -= f * dx / r
            forces[j, 1] -= f * dy / r
            forces[j, 2] -= f * dz / r
        for a in range(n_angles):
            i = int(angles[a, 0])
            j = int(angles[a, 1])
            kk = int(angles[a, 2])
            k_theta = angles[a, 3]
            theta0 = angles[a, 4]
            v1x = positions[i, 0] - positions[j, 0]
            v1y = positions[i, 1] - positions[j, 1]
            v1z = positions[i, 2] - positions[j, 2]
            v2x = positions[kk, 0] - positions[j, 0]
            v2y = positions[kk, 1] - positions[j, 1]
            v2z = positions[kk, 2] - positions[j, 2]
            r1 = np.sqrt(v1x * v1x + v1y * v1y + v1z * v1z)
            r2 = np.sqrt(v2x * v2x + v2y * v2y + v2z * v2z)
            cos_theta = (v1x * v2x + v1y * v2y + v1z * v2z) / (r1 * r2)
            cos_theta = min(1.0, max(-1.0, cos_theta))
            theta = np.arccos(cos_theta)
            sin_theta = np.sqrt(max(1.0 - cos_theta * cos_theta, 1.0e-12))
            coef = k_theta * (theta - theta0) / sin_theta
            fix = coef * (v2x / (r1 * r2) - cos_theta * v1x / (r1 * r1))
            fiy = coef * (v2y / (r1 * r2) - cos_theta * v1y / (r1 * r1))
            fiz = coef * (v2z / (r1 * r2) - cos_theta * v1z / (r1 * r1))
            fkx = coef * (v1x / (r1 * r2) - cos_theta * v2x / (r2 * r2))
            fky = coef * (v1y / (r1 * r2) - cos_theta * v2y / (r2 * r2))
            fkz = coef * (v1z / (r1 * r2) - cos_theta * v2z / (r2 * r2))
            forces[i, 0] += fix
            forces[i, 1] += fiy
            forces[i, 2] += fiz
            forces[kk, 0] += fkx
            forces[kk, 1] += fky
            forces[kk, 2] += fkz
            forces[j, 0] -= fix + fkx
            forces[j, 1] -= fiy + fky
            forces[j, 2] -= fiz + fkz
        for i in range(n):
            forces[i, 0] /= masses[i]
            forces[i, 1] /= masses[i]
            forces[i, 2] /= masses[i]
        return forces

    return accel


class HarmonicMolecule(MolecularDynamicsSystem):
    r"""A small, non-periodic cluster of atoms held together by harmonic bonds and angles.

    A minimal molecular "force field": mirrors
    :class:`chemistrykit.kinetics.systems.networks.StoichiometricNetwork`'s
    structural-array-and-factory-closure pattern, but for bonded-force
    evaluation instead of mass-action reaction rates. The angle-bending
    force is derived from the chain rule
    :math:`\vec F_i=-dU/d\theta\cdot d\theta/d\vec r_i` applied to
    :math:`\theta=\arccos(\hat v_1\cdot\hat v_2)`, a standard force-field
    result (e.g. Leach, *Molecular Modelling: Principles and
    Applications*, 2nd ed., Ch. 4).

    Parameters
    ----------
    positions, velocities : array-like, shape (n_atoms, 3)
    masses : array-like, shape (n_atoms,)
    bonds : sequence of (i, j, k, r0), optional
        Harmonic bonds between atom indices `i`, `j` with force constant
        `k` and equilibrium length `r0`.
    angles : sequence of (i, j, k, k_theta, theta0), optional
        Harmonic angle terms at vertex atom `j` (between bonds `i-j` and
        `k-j`) with force constant `k_theta` and equilibrium angle
        `theta0` (radians).
    """

    def __init__(self, positions, velocities, masses, bonds: Sequence = (), angles: Sequence = ()):
        masses = np.asarray(masses, dtype=np.float64)
        self.bonds = np.asarray(bonds, dtype=np.float64).reshape(-1, 4) if len(bonds) else np.zeros((0, 4))
        self.angles = np.asarray(angles, dtype=np.float64).reshape(-1, 5) if len(angles) else np.zeros((0, 5))
        self._accel_njit = _make_harmonic_molecule_accel_njit(self.bonds, self.angles, masses)
        self.params = np.empty(0)
        super().__init__(positions, velocities, masses)

    def forces_and_potential(self, positions):
        positions = np.asarray(positions, dtype=np.float64)
        forces = np.zeros_like(positions)
        potential = 0.0
        for i, j, k, r0 in self.bonds:
            i, j = int(i), int(j)
            diff = positions[i] - positions[j]
            r = float(np.linalg.norm(diff))
            f = -k * (r - r0)
            forces[i] += f * diff / r
            forces[j] -= f * diff / r
            potential += 0.5 * k * (r - r0) ** 2
        for i, j, k, k_theta, theta0 in self.angles:
            i, j, k = int(i), int(j), int(k)
            v1 = positions[i] - positions[j]
            v2 = positions[k] - positions[j]
            r1 = float(np.linalg.norm(v1))
            r2 = float(np.linalg.norm(v2))
            cos_theta = np.clip(np.dot(v1, v2) / (r1 * r2), -1.0, 1.0)
            theta = float(np.arccos(cos_theta))
            sin_theta = max(np.sqrt(1.0 - cos_theta**2), 1e-12)
            coef = k_theta * (theta - theta0) / sin_theta
            f_i = coef * (v2 / (r1 * r2) - cos_theta * v1 / r1**2)
            f_k = coef * (v1 / (r1 * r2) - cos_theta * v2 / r2**2)
            forces[i] += f_i
            forces[k] += f_k
            forces[j] -= f_i + f_k
            potential += 0.5 * k_theta * (theta - theta0) ** 2
        return forces, potential
