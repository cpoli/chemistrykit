r"""
Musher's hypervalent molecules: main-group centres beyond the octet
===================================================================

James Musher (1969) grouped molecules such as :math:`\mathrm{PCl_5}`,
:math:`\mathrm{SF_6}`, :math:`\mathrm{ClF_3}` and the xenon fluorides into
one class: *hypervalent* molecules, whose central main-group atom appears
to hold more than the eight electrons of Lewis's octet. Counting the
Lewis electrons around the central atom (two per bond plus two per lone
pair) gives 10 or 12 for all of them.

Whatever the bonding explanation (d-orbital participation, or the
three-centre four-electron bond favoured today), VSEPR domain counting
predicts their shapes. This example builds six hypervalent molecules
with
:func:`~chemistrykit.structure.systems.vsepr.build_vsepr_molecule`,
counts the electrons around each central atom, and confirms each shape's
symmetry with
:func:`~chemistrykit.structure.systems.point_group.determine_point_group`.
"""

# %%
import matplotlib.pyplot as plt

from chemistrykit.structure.systems.point_group import determine_point_group
from chemistrykit.structure.systems.vsepr import VSEPRGeometry, build_vsepr_molecule
from chemistrykit.structure.visualizers.structure_plots import plot_molecule_3d

# name: (steric number, lone pairs, central atom, ligand, bond length in A, expected point group)
hypervalent = {
    "PCl5": (5, 0, "P", "Cl", 2.07, "D3h"),
    "SF4": (5, 1, "S", "F", 1.60, "C2v"),
    "ClF3": (5, 2, "Cl", "F", 1.65, "C2v"),
    "XeF2": (5, 3, "Xe", "F", 1.98, "D_inf_h"),
    "SF6": (6, 0, "S", "F", 1.56, "Oh"),
    "XeF4": (6, 2, "Xe", "F", 1.95, "D4h"),
}

molecules = {}
for name, (sn, lp, central, ligand, length, expected) in hypervalent.items():
    molecule = build_vsepr_molecule(sn, lp, bond_length=length, central_symbol=central, ligand_symbol=ligand)
    molecules[name] = molecule
    n_bonds = sn - lp
    electrons = 2 * n_bonds + 2 * lp  # single bonds: one shared pair each
    group = determine_point_group(molecule).group_name
    shape = VSEPRGeometry(sn, lp).shape_name
    print(f"{name:5s}: {electrons:2d} electrons around {central:2s} (octet exceeded by {electrons - 8}), AX{n_bonds}E{lp} {shape}, point group {group}")
    assert electrons > 8
    assert group == expected, (name, group)

# %%
# The six shapes, all derived from the five- and six-domain polyhedra:

fig = plt.figure(figsize=(12, 8))
for i, (name, molecule) in enumerate(molecules.items(), start=1):
    ax = fig.add_subplot(2, 3, i, projection="3d")
    plot_molecule_3d(molecule, ax=ax)
    ax.set_title(f"{name} ({2 * (hypervalent[name][0])} e around {hypervalent[name][2]})", fontsize=10)
fig.suptitle("Hypervalent molecules: more than an octet, shapes still from VSEPR")
fig.tight_layout()
plt.show()
