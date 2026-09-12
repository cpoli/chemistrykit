"""chemistrykit.structure: molecular structure and bonding.

A lightweight :class:`~chemistrykit.structure.core.base_system.Molecule`
container (atoms + 3D coordinates + bond list, no external file-format
parsing); VSEPR geometry prediction from steric number, generating real
3D coordinates for the five idealized electron-domain polyhedra; point-
group determination from 3D coordinates (genuine symmetry-element
detection via distance/angle checks, not a formula lookup) plus
character tables for the common point groups; bond order from the
Pauling bond-length correlation and from Huckel-theory MO coefficients
(the Coulson bond order); and formal-charge/oxidation-state assignment
from a Lewis structure, using electronegativities from
:mod:`chemistrykit.periodic_table`.
"""

__version__ = "0.1.0"

from chemistrykit.structure.core.base_system import Molecule, angle_between, unit_vector
from chemistrykit.structure.systems.bonding import bond_length_from_order, bond_order_from_length, coulson_pi_bond_order
from chemistrykit.structure.systems.lewis import LewisStructure
from chemistrykit.structure.systems.point_group import (
    CHARACTER_TABLES,
    PointGroupCharacterTable,
    PointGroupResult,
    determine_point_group,
    get_character_table,
)
from chemistrykit.structure.systems.vsepr import AXE_SHAPE_NAMES, IDEAL_BOND_ANGLES, VSEPRGeometry, build_vsepr_molecule, domain_positions

__all__ = [
    "__version__",
    "Molecule",
    "unit_vector",
    "angle_between",
    "IDEAL_BOND_ANGLES",
    "AXE_SHAPE_NAMES",
    "domain_positions",
    "VSEPRGeometry",
    "build_vsepr_molecule",
    "PointGroupResult",
    "determine_point_group",
    "PointGroupCharacterTable",
    "CHARACTER_TABLES",
    "get_character_table",
    "bond_order_from_length",
    "bond_length_from_order",
    "coulson_pi_bond_order",
    "LewisStructure",
]
