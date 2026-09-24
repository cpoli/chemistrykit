r"""Molecular dipole moments from point charges or bond dipoles.

Peter Debye showed in 1912 that many molecules carry a permanent electric
dipole moment, and that measuring it (from the temperature dependence of
a gas's dielectric constant) tests a proposed molecular geometry (P.
Debye, *Phys. Z.* 13, 97 (1912); *Polar Molecules* (Chemical Catalog Co.,
New York, 1929)). For a set of point charges :math:`q_i` at positions
:math:`\mathbf r_i`,

.. math::

    \boldsymbol\mu = \sum_i q_i\,\mathbf r_i

which, for a neutral molecule, does not depend on the choice of origin.
Equivalently, a molecule's dipole is the vector sum of its bond dipoles,
so a symmetric arrangement (linear CO2, trigonal BF3, tetrahedral CH4)
cancels exactly while a bent or pyramidal one (H2O, NH3) does not. The
unit named after Debye, 1 D = :math:`10^{-21}/c` C m, is about 0.2082
e Å.
"""

from __future__ import annotations

import numpy as np

from chemistrykit.constants import ELEMENTARY_CHARGE, C

__all__ = ["DEBYE", "E_ANGSTROM_IN_DEBYE", "dipole_moment", "bond_dipole_sum"]

#: One debye, in coulomb metres (:math:`10^{-21}/c`).
DEBYE = 1e-21 / C

#: One elementary charge times one angstrom, expressed in debye (about 4.803).
E_ANGSTROM_IN_DEBYE = ELEMENTARY_CHARGE * 1e-10 / DEBYE


def dipole_moment(coordinates, charges) -> np.ndarray:
    r"""Dipole-moment vector of a set of point charges, in debye.

    .. math::

        \boldsymbol\mu = \sum_i q_i\,\mathbf r_i

    Parameters
    ----------
    coordinates : array-like, shape (n, 3)
        Positions in angstrom (e.g. :attr:`Molecule.coordinates
        <chemistrykit.structure.core.base_system.Molecule.coordinates>`).
    charges : array-like, shape (n,)
        Partial charges in units of the elementary charge. If they do not
        sum to zero the result depends on the origin.

    Returns
    -------
    ndarray, shape (3,)
        Dipole vector in debye (pointing from negative to positive charge,
        the physics convention).

    Examples
    --------
    Charges of +e and -e one angstrom apart make a dipole of 4.803 D:

    >>> mu = dipole_moment([[0, 0, 0], [0, 0, 1.0]], [-1.0, 1.0])
    >>> round(float(mu[2]), 3)
    4.803

    Linear CO2 with partial charges has no net dipole:

    >>> mu = dipole_moment([[0, 0, -1.16], [0, 0, 0], [0, 0, 1.16]], [-0.4, 0.8, -0.4])
    >>> float(np.linalg.norm(mu))
    0.0
    """
    x = np.asarray(coordinates, dtype=np.float64)
    q = np.asarray(charges, dtype=np.float64)
    if x.ndim != 2 or x.shape[1] != 3 or q.shape != (x.shape[0],):
        raise ValueError("coordinates must have shape (n, 3) and charges shape (n,)")
    return (q[:, None] * x).sum(axis=0) * E_ANGSTROM_IN_DEBYE


def bond_dipole_sum(bond_vectors, bond_moments) -> np.ndarray:
    r"""Molecular dipole as the vector sum of bond dipoles.

    .. math::

        \boldsymbol\mu = \sum_b \mu_b\,\hat{\mathbf u}_b

    Parameters
    ----------
    bond_vectors : array-like, shape (n_bonds, 3)
        Direction of each bond dipole (any length; normalized here).
    bond_moments : array-like, shape (n_bonds,)
        Magnitude of each bond dipole, in debye.

    Returns
    -------
    ndarray, shape (3,)
        Molecular dipole vector, in debye.

    Examples
    --------
    Two O-H bond dipoles of 1.51 D at water's 104.5 degree angle give
    :math:`2\mu_{OH}\cos(\theta/2) \approx 1.85` D, the measured value:

    >>> half = np.radians(104.5 / 2)
    >>> u = [[np.sin(half), 0, np.cos(half)], [-np.sin(half), 0, np.cos(half)]]
    >>> round(float(np.linalg.norm(bond_dipole_sum(u, [1.51, 1.51]))), 2)
    1.85
    """
    u = np.asarray(bond_vectors, dtype=np.float64)
    m = np.asarray(bond_moments, dtype=np.float64)
    if u.ndim != 2 or u.shape[1] != 3 or m.shape != (u.shape[0],):
        raise ValueError("bond_vectors must have shape (n, 3) and bond_moments shape (n,)")
    norms = np.linalg.norm(u, axis=1)
    if np.any(norms < 1e-14):
        raise ValueError("bond vectors must be nonzero")
    return ((m / norms)[:, None] * u).sum(axis=0)
