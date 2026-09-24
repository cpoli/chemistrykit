r"""Holonomic bond-length constraints: the SHAKE algorithm.

Fast bond-stretch vibrations limit an MD time step far more than the slow
motions of interest. Ryckaert, Ciccotti & Berendsen replaced such bonds by
rigid constraints :math:`|\vec r_i - \vec r_j| = d_{ij}`, enforced after
each unconstrained Verlet position update by iteratively adding
displacements along the *old* bond vectors (the constraint-force
directions) until every constraint holds to a tolerance -- J.-P. Ryckaert,
G. Ciccotti & H. J. C. Berendsen, *J. Comput. Phys.* 23, 327 (1977). For
each violated constraint in turn,

.. math::

    g = \frac{d_{ij}^2 - |\vec s_{ij}|^2}{2\,(\vec s_{ij}\cdot\vec r^{\,\text{old}}_{ij})\,(1/m_i + 1/m_j)},
    \qquad \vec s_i \mathrel{+}= \frac{g}{m_i}\vec r^{\,\text{old}}_{ij},
    \quad \vec s_j \mathrel{-}= \frac{g}{m_j}\vec r^{\,\text{old}}_{ij},

where :math:`\vec s` are the updated positions. :class:`ShakeMolecule` uses
this inside velocity-Verlet, followed by the velocity stage of Andersen's
RATTLE (*J. Comput. Phys.* 52, 24 (1983)) that removes velocity components
along each constrained bond.
"""

from __future__ import annotations

from collections.abc import Sequence

import numpy as np

from chemistrykit.md.systems.pair_potentials import HarmonicMolecule

__all__ = ["shake", "ShakeMolecule"]


def shake(positions_new, positions_old, constraints, masses, tol: float = 1e-10, max_iter: int = 1000):
    r"""Correct unconstrained positions so every bond-length constraint holds (SHAKE).

    Parameters
    ----------
    positions_new : array-like, shape (n_atoms, n_dim)
        Positions after an unconstrained update.
    positions_old : array-like, shape (n_atoms, n_dim)
        Positions at the previous step (which satisfy the constraints);
        their bond vectors set the correction directions.
    constraints : sequence of (i, j, d)
        Atom pairs and their fixed separations.
    masses : array-like, shape (n_atoms,)
    tol : float, default 1e-10
        Relative tolerance on :math:`|s^2 - d^2|/d^2`.
    max_iter : int, default 1000

    Returns
    -------
    positions : ndarray, shape (n_atoms, n_dim)
        Corrected positions.
    n_iter : int
        Number of sweeps over the constraints needed to converge.

    Examples
    --------
    Two equal masses pulled apart to separation 1.2 are pulled back
    symmetrically to the constrained length 1.0:

    >>> import numpy as np
    >>> old = np.array([[0.0, 0.0, 0.0], [1.0, 0.0, 0.0]])
    >>> new = np.array([[-0.1, 0.0, 0.0], [1.1, 0.0, 0.0]])
    >>> fixed, _ = shake(new, old, [(0, 1, 1.0)], masses=[1.0, 1.0])
    >>> np.round(fixed[:, 0], 10)
    array([0., 1.])
    """
    s = np.array(positions_new, dtype=np.float64, copy=True)
    old = np.asarray(positions_old, dtype=np.float64)
    inv_m = 1.0 / np.asarray(masses, dtype=np.float64)
    for n_iter in range(1, max_iter + 1):
        converged = True
        for i, j, d in constraints:
            i, j = int(i), int(j)
            s_ij = s[i] - s[j]
            diff = d * d - float(s_ij @ s_ij)
            if abs(diff) > tol * d * d:
                converged = False
                r_ij = old[i] - old[j]
                g = diff / (2.0 * float(s_ij @ r_ij) * (inv_m[i] + inv_m[j]))
                s[i] += g * inv_m[i] * r_ij
                s[j] -= g * inv_m[j] * r_ij
        if converged:
            return s, n_iter
    raise RuntimeError("SHAKE did not converge")


class ShakeMolecule(HarmonicMolecule):
    r"""A :class:`~chemistrykit.md.systems.pair_potentials.HarmonicMolecule` with rigid (SHAKE) bonds.

    Each velocity-Verlet step makes an unconstrained half-kick and drift,
    corrects the positions with :func:`shake`, recovers the half-step
    velocities from the constrained displacement, completes the second
    half-kick, and finally removes each constrained bond's relative
    velocity component (RATTLE's velocity stage).

    Parameters
    ----------
    positions, velocities : array-like, shape (n_atoms, 3)
    masses : array-like, shape (n_atoms,)
    constraints : sequence of (i, j, d)
        Rigid bonds; `positions` should already satisfy them.
    bonds, angles : sequence, optional
        Flexible harmonic terms, as for
        :class:`~chemistrykit.md.systems.pair_potentials.HarmonicMolecule`.
    tol : float, default 1e-10

    Examples
    --------
    A rigid rotating diatomic keeps its bond length exactly:

    >>> import numpy as np
    >>> pos = np.array([[0.0, 0.0, 0.0], [1.0, 0.0, 0.0]])
    >>> vel = np.array([[0.0, -1.0, 0.0], [0.0, 1.0, 0.0]])
    >>> rotor = ShakeMolecule(pos, vel, [1.0, 1.0], constraints=[(0, 1, 1.0)])
    >>> for _ in range(500):
    ...     rotor.step(0.01)
    >>> round(float(np.linalg.norm(rotor.positions[0] - rotor.positions[1])), 8)
    1.0
    """

    def __init__(self, positions, velocities, masses, constraints: Sequence, bonds: Sequence = (), angles: Sequence = (), tol: float = 1e-10):
        super().__init__(positions, velocities, masses, bonds=bonds, angles=angles)
        self.constraints = [(int(i), int(j), float(d)) for i, j, d in constraints]
        self.tol = float(tol)
        self._remove_bond_velocities()

    def degrees_of_freedom(self) -> int:
        """Cartesian degrees of freedom minus center-of-mass motion and one per constraint.

        Returns
        -------
        int
        """
        return super().degrees_of_freedom() - len(self.constraints)

    def bond_lengths(self) -> np.ndarray:
        """Current length of every constrained bond.

        Returns
        -------
        ndarray, shape (n_constraints,)
        """
        p = self.positions
        return np.array([np.linalg.norm(p[i] - p[j]) for i, j, _ in self.constraints])

    def _remove_bond_velocities(self) -> None:
        inv_m = 1.0 / self.masses
        v, p = self.velocities, self.positions
        for _ in range(1000):
            converged = True
            for i, j, d in self.constraints:
                r_ij = p[i] - p[j]
                rv = float(r_ij @ (v[i] - v[j]))
                if abs(rv) > self.tol * d * d:
                    converged = False
                    k = rv / (float(r_ij @ r_ij) * (inv_m[i] + inv_m[j]))
                    v[i] -= k * inv_m[i] * r_ij
                    v[j] += k * inv_m[j] * r_ij
            if converged:
                return
        raise RuntimeError("RATTLE velocity stage did not converge")

    def step(self, dt: float) -> None:
        """Advance one constrained velocity-Verlet (SHAKE/RATTLE) step, in place.

        Parameters
        ----------
        dt : float
        """
        accel = self._accel_njit(self.positions, self.t, self.params)
        v_half = self.velocities + 0.5 * dt * accel
        unconstrained = self.positions + dt * v_half
        new_positions, _ = shake(unconstrained, self.positions, self.constraints, self.masses, tol=self.tol)
        v_half = (new_positions - self.positions) / dt
        self.positions = new_positions
        self.t += dt
        accel = self._accel_njit(self.positions, self.t, self.params)
        self.velocities = v_half + 0.5 * dt * accel
        self._remove_bond_velocities()
