r"""VSEPR (valence-shell electron-pair repulsion) geometry prediction.

See Gillespie & Nyholm, *Q. Rev. Chem. Soc.* 11, 339 (1957) (the original
statement of the model) and R. J. Gillespie, *J. Chem. Educ.* 47, 18
(1970), or Atkins & de Paula, *Physical Chemistry*, 11th ed., Ch. 9.1(b),
for the qualitative theory throughout.

VSEPR treats a central atom's bonding and lone electron pairs as point
charges that arrange themselves on a sphere to minimize mutual repulsion.
For `N` electron domains (the *steric number*, :math:`N` = number of
sigma bonds + number of lone pairs -- multiple bonds count as one domain,
same as a single bond, since VSEPR only cares about the number of
electron-dense regions, not the electron count within one), the
lowest-repulsion arrangement is one of five idealized polyhedra (`N`
= 2 through 6): linear, trigonal planar, tetrahedral, trigonal
bipyramidal, and octahedral. This module generates those polyhedra's
vertex positions as genuine 3D coordinates (not a shape-name lookup
table), then -- when some domains are lone pairs rather than bonds --
selects *which* vertices the lone pairs occupy using the standard
lone-pair-placement rule (lone pairs preferentially occupy the least
crowded positions), reproducing the textbook real-molecule shapes (e.g.
water's bent AX2E2 geometry, xenon tetrafluoride's square-planar AX4E2
geometry) from that placement.

**Approximation flagged explicitly**: the bond angles produced here are
the *idealized*, unperturbed polyhedron angles (e.g. exactly 109.47
degrees between any two bonding directions of a tetrahedral AX4E0
molecule). Real lone-pair-containing molecules deviate from these ideal
angles because a lone pair, being more diffuse and closer to the central
atom, repels other domains more strongly than a bonding pair does (the
lone-pair > bonding-pair repulsion-strength rule) -- e.g. real ammonia's
H-N-H angle is compressed to about 106.7 degrees from the idealized
109.47 degrees, and real water's H-O-H angle to about 104.5 degrees. This
module reproduces the correct idealized *shape* (which vertices carry
bonds vs. lone pairs) but not that further angle-compression refinement,
which VSEPR only predicts qualitatively (Atkins & de Paula, *Physical
Chemistry*, 11th ed., Ch. 9.1(b)).
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from chemistrykit.structure.core.base_system import Molecule

__all__ = [
    "AXE_SHAPE_NAMES",
    "IDEAL_BOND_ANGLES",
    "domain_positions",
    "VSEPRGeometry",
    "build_vsepr_molecule",
]

#: For each steric number (electron-domain count), the characteristic
#: idealized bond angle(s) between adjacent bonding directions, in
#: degrees, when no lone pairs are present (Gillespie & Nyholm, *Q. Rev.
#: Chem. Soc.* 11, 339 (1957)). Steric number 5 and 6 have more than one
#: distinct angle (axial-equatorial vs. equatorial-equatorial for the
#: trigonal bipyramid; adjacent vs. opposite for the octahedron); given as
#: a tuple, smallest first.
IDEAL_BOND_ANGLES: dict = {
    2: (180.0,),
    3: (120.0,),
    4: (109.4712206,),
    5: (90.0, 120.0, 180.0),
    6: (90.0, 180.0),
}

#: AXE-method shape name for each (steric_number, lone_pairs) combination
#: (A = central atom, X = bonded ligand, E = lone pair; e.g. AX2E2 for
#: water). This is reference nomenclature data, the same spirit as
#: :mod:`chemistrykit.periodic_table`'s element table -- the actual
#: *geometry* for each entry is still generated as real 3D coordinates by
#: :func:`domain_positions` / :func:`build_vsepr_molecule`, not looked up
#: here.
AXE_SHAPE_NAMES: dict = {
    (2, 0): "linear",
    (3, 0): "trigonal planar",
    (3, 1): "bent",
    (4, 0): "tetrahedral",
    (4, 1): "trigonal pyramidal",
    (4, 2): "bent",
    (5, 0): "trigonal bipyramidal",
    (5, 1): "seesaw",
    (5, 2): "T-shaped",
    (5, 3): "linear",
    (6, 0): "octahedral",
    (6, 1): "square pyramidal",
    (6, 2): "square planar",
}


def domain_positions(steric_number: int) -> np.ndarray:
    r"""Return unit-vector electron-domain positions for the ideal `N`-domain polyhedron.

    Real 3D coordinate generation (not a lookup table): each polyhedron's
    vertices are placed from their defining symmetry, then normalized to
    the unit sphere.

    Parameters
    ----------
    steric_number : int
        Number of electron domains (sigma bonds + lone pairs), 2 through 6.

    Returns
    -------
    ndarray, shape (steric_number, 3)
        Unit vectors from the central atom, ordered so that -- for a
        steric number with structurally inequivalent sites -- the sites
        *least* favorable to a lone pair (see :func:`build_vsepr_molecule`)
        come *last*: equatorial-before-axial for `steric_number=5`
        (equatorial positions have only two 90-degree neighbors vs. an
        axial position's three, so a lone pair placed equatorially incurs
        less strong lone-pair/bonding-pair repulsion -- Gillespie &
        Nyholm, *Q. Rev. Chem. Soc.* 11, 339 (1957)).

    Examples
    --------
    Every pair of tetrahedral vertices subtends exactly the tetrahedral
    angle, 109.47 degrees:

    >>> import numpy as np
    >>> from chemistrykit.structure.core.base_system import angle_between
    >>> verts = domain_positions(4)
    >>> angles = [angle_between(verts[i], verts[j]) for i in range(4) for j in range(i + 1, 4)]
    >>> bool(np.allclose(angles, 109.4712206))
    True
    """
    if steric_number == 2:
        return np.array([[0.0, 0.0, 1.0], [0.0, 0.0, -1.0]])
    if steric_number == 3:
        angles = np.radians([0.0, 120.0, 240.0])
        return np.column_stack([np.cos(angles), np.sin(angles), np.zeros(3)])
    if steric_number == 4:
        verts = np.array(
            [
                [1.0, 1.0, 1.0],
                [1.0, -1.0, -1.0],
                [-1.0, 1.0, -1.0],
                [-1.0, -1.0, 1.0],
            ]
        )
        return verts / np.linalg.norm(verts[0])
    if steric_number == 5:
        eq_angles = np.radians([0.0, 120.0, 240.0])
        equatorial = np.column_stack([np.cos(eq_angles), np.sin(eq_angles), np.zeros(3)])
        axial = np.array([[0.0, 0.0, 1.0], [0.0, 0.0, -1.0]])
        return np.vstack([equatorial, axial])
    if steric_number == 6:
        return np.array(
            [
                [1.0, 0.0, 0.0],
                [-1.0, 0.0, 0.0],
                [0.0, 1.0, 0.0],
                [0.0, -1.0, 0.0],
                [0.0, 0.0, 1.0],
                [0.0, 0.0, -1.0],
            ]
        )
    raise ValueError("steric_number must be 2, 3, 4, 5, or 6")


def _lone_pair_indices(steric_number: int, lone_pairs: int) -> list:
    """Indices (into :func:`domain_positions`'s output) that lone pairs occupy.

    For `steric_number` 5, lone pairs fill the equatorial positions
    first (:func:`domain_positions` orders equatorial before axial), the
    least sterically crowded sites (Gillespie & Nyholm, *Q. Rev. Chem.
    Soc.* 11, 339 (1957)) -- reproducing e.g. seesaw SF4 (1 lone pair,
    equatorial), T-shaped ClF3 (2 lone pairs, both equatorial), and
    linear XeF2 (3 lone pairs, all three equatorial positions).

    For `steric_number` 6, all six octahedral positions are equivalent,
    so the first lone pair may take any of them; a second lone pair then
    goes *trans* (opposite) to the first rather than adjacent to it,
    minimizing lone-pair/lone-pair repulsion -- reproducing square-planar
    XeF4's real geometry (2 lone pairs, trans to each other).
    """
    if steric_number == 6:
        # Position 5 first, then its trans partner (position 4), then the
        # remaining two trans pairs -- order is arbitrary beyond ensuring
        # the second lone pair lands trans to the first.
        order = [5, 4, 0, 1, 2, 3]
        return order[:lone_pairs]
    return list(range(0, lone_pairs))


@dataclass
class VSEPRGeometry:
    """The predicted shape for a given steric number and lone-pair count.

    Parameters
    ----------
    steric_number : int
        Number of electron domains (sigma bonds + lone pairs).
    lone_pairs : int
        Number of lone pairs on the central atom, ``<= steric_number``.

    Examples
    --------
    >>> VSEPRGeometry(steric_number=4, lone_pairs=2).shape_name
    'bent'
    >>> VSEPRGeometry(steric_number=4, lone_pairs=0).shape_name
    'tetrahedral'
    """

    steric_number: int
    lone_pairs: int

    def __post_init__(self):
        if not (2 <= self.steric_number <= 6):
            raise ValueError("steric_number must be 2, 3, 4, 5, or 6")
        if not (0 <= self.lone_pairs <= self.steric_number):
            raise ValueError("lone_pairs must be between 0 and steric_number")

    @property
    def n_bonding_domains(self) -> int:
        """int: Number of sigma-bond (ligand) positions, ``steric_number - lone_pairs``."""
        return self.steric_number - self.lone_pairs

    @property
    def shape_name(self) -> str:
        """str: The AXE-method molecular-shape name (Gillespie's nomenclature)."""
        return AXE_SHAPE_NAMES[(self.steric_number, self.lone_pairs)]

    @property
    def electron_domain_positions(self) -> np.ndarray:
        """ndarray, shape (steric_number, 3): All domain unit-vector directions, bonds and lone pairs alike."""
        return domain_positions(self.steric_number)

    @property
    def bonding_positions(self) -> np.ndarray:
        """ndarray, shape (n_bonding_domains, 3): Unit-vector directions of just the bonded ligands.

        Lone pairs preferentially occupy the least sterically crowded
        domain positions (see :func:`domain_positions`), so these are the
        *remaining* positions after that assignment -- e.g. for
        `steric_number=5`, `lone_pairs=1` (seesaw), one equatorial
        position is given up to the lone pair and the four bonding
        positions are the other two equatorial plus both axial sites.
        """
        all_positions = self.electron_domain_positions
        lp_indices = set(_lone_pair_indices(self.steric_number, self.lone_pairs))
        return np.array([p for i, p in enumerate(all_positions) if i not in lp_indices])


def build_vsepr_molecule(steric_number: int, lone_pairs: int, bond_length: float = 1.0, central_symbol: str = "C", ligand_symbol: str = "H") -> Molecule:
    """Build a :class:`~chemistrykit.structure.core.base_system.Molecule` from a VSEPR prediction.

    Parameters
    ----------
    steric_number : int
        Number of electron domains (sigma bonds + lone pairs).
    lone_pairs : int
        Number of lone pairs on the central atom.
    bond_length : float, default 1.0
        Central-atom-to-ligand distance, in angstrom.
    central_symbol : str, default "C"
        Element symbol placed at the central atom.
    ligand_symbol : str, default "H"
        Element symbol placed at every ligand position.

    Returns
    -------
    Molecule
        Central atom at the origin plus `n_bonding_domains` identical
        ligand atoms; lone pairs are not represented as atoms (they carry
        no nuclear position), only via their effect on which directions
        the ligands occupy.

    Examples
    --------
    Methane (AX4E0): a perfect tetrahedron, every H-C-H angle exactly
    109.47 degrees:

    >>> methane = build_vsepr_molecule(steric_number=4, lone_pairs=0, bond_length=1.09, central_symbol="C", ligand_symbol="H")
    >>> round(methane.bond_angle(1, 0, 2), 4)
    109.4712

    Water (AX2E2): the idealized VSEPR angle is the parent tetrahedron's
    109.47 degrees (see the module docstring for why this is an
    idealization -- the real H-O-H angle is compressed to 104.5 degrees
    by extra lone-pair repulsion):

    >>> water = build_vsepr_molecule(steric_number=4, lone_pairs=2, bond_length=0.96, central_symbol="O", ligand_symbol="H")
    >>> round(water.bond_angle(1, 0, 2), 4)
    109.4712
    """
    geometry = VSEPRGeometry(steric_number=steric_number, lone_pairs=lone_pairs)
    ligand_positions = geometry.bonding_positions * bond_length
    symbols = [central_symbol] + [ligand_symbol] * geometry.n_bonding_domains
    coordinates = np.vstack([np.zeros((1, 3)), ligand_positions])
    bonds = [(0, i + 1) for i in range(geometry.n_bonding_domains)]
    return Molecule(symbols=symbols, coordinates=coordinates, bonds=bonds)
