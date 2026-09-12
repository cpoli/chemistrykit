"""chemistrykit.crystal: crystallography and solid-state chemistry.

The 7 crystal systems and the general unit-cell-volume formula; hard-
sphere packing (packing fraction, coordination number, atoms per cell)
for the simple cubic, body-centered cubic, face-centered cubic, and
(ideal) hexagonal close-packed lattices; ionic-crystal lattice energy via
the Born-Lande and Kapustinskii equations, and the Madelung constant of
the NaCl structure from a genuinely converging (Evjen-method) lattice
summation; Bragg's law and powder-XRD peak positions with structure
factors (including systematic absences) for cubic lattices; and
Schottky/Frenkel point-defect equilibrium.
"""

__version__ = "0.1.0"

from chemistrykit.crystal.core.base_system import LatticeEnergyModel, LatticePacking
from chemistrykit.crystal.systems.crystal_systems import classify_crystal_system, unit_cell_volume
from chemistrykit.crystal.systems.defects import frenkel_defect_concentration, schottky_defect_concentration
from chemistrykit.crystal.systems.lattice_energy import BornLandeLatticeEnergy, KapustinskiiLatticeEnergy
from chemistrykit.crystal.systems.madelung import MADELUNG_CONSTANT_NACL_LITERATURE, madelung_constant_nacl
from chemistrykit.crystal.systems.packing import (
    BodyCenteredCubicPacking,
    FaceCenteredCubicPacking,
    HexagonalClosePacking,
    SimpleCubicPacking,
)
from chemistrykit.crystal.systems.xrd import XRDPeak, bragg_angle, d_spacing_cubic, powder_xrd_peaks, structure_factor

__all__ = [
    "__version__",
    "LatticePacking",
    "LatticeEnergyModel",
    "classify_crystal_system",
    "unit_cell_volume",
    "SimpleCubicPacking",
    "BodyCenteredCubicPacking",
    "FaceCenteredCubicPacking",
    "HexagonalClosePacking",
    "BornLandeLatticeEnergy",
    "KapustinskiiLatticeEnergy",
    "MADELUNG_CONSTANT_NACL_LITERATURE",
    "madelung_constant_nacl",
    "bragg_angle",
    "d_spacing_cubic",
    "structure_factor",
    "XRDPeak",
    "powder_xrd_peaks",
    "schottky_defect_concentration",
    "frenkel_defect_concentration",
]
