r"""The shared :class:`Molecule` container, and small geometry helpers.

Unlike :mod:`chemistrykit.quantum` (where every ``systems/`` module
implements one of exactly two shared ABCs, :class:`~chemistrykit.quantum.core.base_system.QuantumSystem`
or :class:`~chemistrykit.quantum.core.base_system.VariationalSolver`) or
:mod:`chemistrykit.solutions` (:class:`~chemistrykit.solutions.core.base_system.WeakElectrolyte`
and :class:`~chemistrykit.solutions.core.base_system.Titration`),
:mod:`chemistrykit.structure`'s four ``systems/`` modules -- VSEPR geometry
prediction, point-group determination, bond order, and formal
charge/oxidation state -- do not share a common polymorphic interface: a
VSEPR prediction *builds* a molecular geometry, a point-group
determination *classifies* one, and bond order/formal-charge are plain
numeric functions over a Lewis structure. Following
:mod:`chemistrykit.solutions.systems.solubility` and
:mod:`chemistrykit.solutions.systems.activity`'s precedent for
non-polymorphic model families, each is implemented directly as
classes/functions in its own ``systems/`` module rather than forced under
an artificial shared ABC.

What genuinely is shared -- the way :class:`~chemistrykit.solutions.core.base_system.TitrationResult`
is shared machinery for every :class:`~chemistrykit.solutions.core.base_system.Titration`
subclass -- is the lightweight :class:`Molecule` container defined here:
a plain numpy-backed atoms + 3D-coordinates + bond-list representation
(no external file-format parsing, per the spec), used as the common
currency every ``systems/`` module in this subpackage builds, consumes, or
both.
"""

from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np

from chemistrykit.periodic_table import get_element

__all__ = ["Molecule"]


@dataclass
class Molecule:
    """A minimal molecular-geometry container: element symbols, 3D coordinates, and a bond list.

    No file-format parsing (PDB/XYZ/SMILES/...) is implemented or
    depended on anywhere in chemistrykit, per the package's "no
    cheminformatics dependencies" scope decision -- a :class:`Molecule` is
    built directly from plain Python/numpy data, e.g. by
    :func:`chemistrykit.structure.systems.vsepr.build_vsepr_molecule` or by
    hand for a worked example.

    Parameters
    ----------
    symbols : sequence of str
        Element symbol of each atom, length `n_atoms`.
    coordinates : array-like, shape (n_atoms, 3)
        Cartesian coordinates, in angstrom (Å) by convention throughout
        this subpackage (chosen because bond lengths and van der Waals
        radii are most naturally tabulated in Å; nothing here depends on
        that unit beyond consistent internal use).
    bonds : sequence of tuple(int, int), optional
        0-indexed atom-index pairs describing the bond connectivity
        (an adjacency/bond list, not a bond-order-aware structure --
        see :class:`chemistrykit.structure.systems.lewis.LewisStructure`
        for a representation that also carries bond orders and lone
        pairs).

    Examples
    --------
    A bent water molecule (experimental geometry: r(O-H) = 0.958 Å,
    angle(H-O-H) = 104.5 degrees):

    >>> import numpy as np
    >>> angle = np.radians(104.5)
    >>> r = 0.958
    >>> water = Molecule(
    ...     symbols=["O", "H", "H"],
    ...     coordinates=[
    ...         [0.0, 0.0, 0.0],
    ...         [r * np.sin(angle / 2), 0.0, r * np.cos(angle / 2)],
    ...         [-r * np.sin(angle / 2), 0.0, r * np.cos(angle / 2)],
    ...     ],
    ...     bonds=[(0, 1), (0, 2)],
    ... )
    >>> round(water.bond_length(0, 1), 3)
    0.958
    >>> round(water.bond_angle(1, 0, 2), 1)
    104.5
    """

    symbols: list
    coordinates: np.ndarray
    bonds: list = field(default_factory=list)

    def __post_init__(self):
        self.symbols = list(self.symbols)
        self.coordinates = np.asarray(self.coordinates, dtype=np.float64)
        if self.coordinates.shape != (len(self.symbols), 3):
            raise ValueError(f"coordinates must have shape (n_atoms, 3) = ({len(self.symbols)}, 3), got {self.coordinates.shape}")
        self.bonds = [(int(i), int(j)) for i, j in self.bonds]
        for i, j in self.bonds:
            if not (0 <= i < self.n_atoms and 0 <= j < self.n_atoms):
                raise ValueError(f"bond ({i}, {j}) references an atom index out of range for {self.n_atoms} atoms")

    @property
    def n_atoms(self) -> int:
        """int: Number of atoms."""
        return len(self.symbols)

    def bond_length(self, i: int, j: int) -> float:
        """Distance between atoms `i` and `j`, in the same units as `coordinates`.

        Parameters
        ----------
        i, j : int
            0-indexed atom indices (need not be an entry in `bonds`).

        Returns
        -------
        float
        """
        return float(np.linalg.norm(self.coordinates[i] - self.coordinates[j]))

    def bond_angle(self, i: int, j: int, k: int) -> float:
        r"""Angle :math:`\angle ijk` at vertex atom `j`, in degrees.

        .. math::

            \theta = \arccos\left(\frac{\vec{u}_{ji}\cdot\vec{u}_{jk}}{|\vec{u}_{ji}||\vec{u}_{jk}|}\right)

        Parameters
        ----------
        i, j, k : int
            0-indexed atom indices; `j` is the vertex (central atom).

        Returns
        -------
        float
            Angle, in degrees, in :math:`[0, 180]`.
        """
        v1 = self.coordinates[i] - self.coordinates[j]
        v2 = self.coordinates[k] - self.coordinates[j]
        cos_theta = np.dot(v1, v2) / (np.linalg.norm(v1) * np.linalg.norm(v2))
        cos_theta = float(np.clip(cos_theta, -1.0, 1.0))
        return float(np.degrees(np.arccos(cos_theta)))

    def center_of_mass(self) -> np.ndarray:
        """Mass-weighted centroid of the atoms, using standard atomic weights.

        Returns
        -------
        ndarray, shape (3,)
        """
        masses = np.array([get_element(s).atomic_mass for s in self.symbols])
        return np.average(self.coordinates, axis=0, weights=masses)

    def centroid(self) -> np.ndarray:
        """Unweighted (geometric) centroid of the atoms.

        Returns
        -------
        ndarray, shape (3,)
        """
        return self.coordinates.mean(axis=0)

    def centered_coordinates(self, weighted: bool = False) -> np.ndarray:
        """Coordinates translated so the centroid (or center of mass) sits at the origin.

        Parameters
        ----------
        weighted : bool, default False
            If True, center on :meth:`center_of_mass`; otherwise on the
            unweighted :meth:`centroid`. For a molecule with genuine point-
            group symmetry, every symmetry element passes through both
            (equivalent atoms share a mass, so the two centers coincide;
            see :mod:`chemistrykit.structure.systems.point_group`), so
            either choice is valid for symmetry-element detection.

        Returns
        -------
        ndarray, shape (n_atoms, 3)
        """
        center = self.center_of_mass() if weighted else self.centroid()
        return self.coordinates - center

    def adjacency(self) -> dict:
        """Build a neighbor list from `bonds`.

        Returns
        -------
        dict of int to list of int
            Maps each atom index with at least one bond to the sorted
            list of its bonded-neighbor indices.
        """
        neighbors: dict = {i: [] for i in range(self.n_atoms)}
        for i, j in self.bonds:
            neighbors[i].append(j)
            neighbors[j].append(i)
        return {i: sorted(js) for i, js in neighbors.items()}

    def is_linear(self, tol: float = 1e-4) -> bool:
        """Whether all atoms lie on a single straight line.

        Parameters
        ----------
        tol : float, default 1e-4
            Tolerance on the normalized cross product used to detect
            collinearity.

        Returns
        -------
        bool

        Examples
        --------
        >>> co2 = Molecule(symbols=["O", "C", "O"], coordinates=[[0, 0, -1.16], [0, 0, 0], [0, 0, 1.16]])
        >>> co2.is_linear()
        True
        """
        if self.n_atoms < 3:
            return True
        centered = self.centered_coordinates()
        reference = None
        for v in centered:
            if np.linalg.norm(v) > 1e-8:
                reference = v / np.linalg.norm(v)
                break
        if reference is None:
            return True
        for v in centered:
            norm = np.linalg.norm(v)
            if norm < 1e-8:
                continue
            u = v / norm
            if np.linalg.norm(np.cross(u, reference)) > tol:
                return False
        return True


def unit_vector(v: np.ndarray) -> np.ndarray:
    """Normalize a vector to unit length.

    Parameters
    ----------
    v : array-like of float

    Returns
    -------
    ndarray
    """
    v = np.asarray(v, dtype=np.float64)
    norm = np.linalg.norm(v)
    if norm < 1e-14:
        raise ValueError("cannot normalize a (near-)zero vector")
    return v / norm


def angle_between(v1: np.ndarray, v2: np.ndarray) -> float:
    """Angle between two vectors, in degrees.

    Parameters
    ----------
    v1, v2 : array-like of float

    Returns
    -------
    float

    Examples
    --------
    >>> round(angle_between([1, 0, 0], [0, 1, 0]), 6)
    90.0
    """
    cos_theta = float(np.dot(v1, v2) / (np.linalg.norm(v1) * np.linalg.norm(v2)))
    cos_theta = float(np.clip(cos_theta, -1.0, 1.0))
    return float(np.degrees(np.arccos(cos_theta)))


__all__.extend(["unit_vector", "angle_between"])
