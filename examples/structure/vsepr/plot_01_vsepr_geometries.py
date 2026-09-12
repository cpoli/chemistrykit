r"""
VSEPR geometry prediction: from steric number to real 3D coordinates
======================================================================

Builds four molecules spanning the five idealized electron-domain
polyhedra -- methane (AX4E0, tetrahedral), water (AX2E2, bent from a
tetrahedral parent), sulfur tetrafluoride (AX4E1, seesaw from a trigonal
bipyramidal parent), and xenon tetrafluoride (AX4E2, square planar from
an octahedral parent) -- via
:func:`~chemistrykit.structure.systems.vsepr.build_vsepr_molecule`, and
verifies the exact idealized bond angles each shape predicts.
"""

# %%
import matplotlib.pyplot as plt

from chemistrykit.structure.systems.vsepr import VSEPRGeometry, build_vsepr_molecule
from chemistrykit.structure.visualizers.structure_plots import plot_molecule_3d

shapes = {
    "methane (AX4E0)": (4, 0, "C", "H", 1.09),
    "water (AX2E2)": (4, 2, "O", "H", 0.96),
    "sulfur tetrafluoride (AX4E1)": (5, 1, "S", "F", 1.60),
    "xenon tetrafluoride (AX4E2)": (6, 2, "Xe", "F", 1.95),
}
molecules = {}
for name, (steric_number, lone_pairs, central, ligand, length) in shapes.items():
    molecules[name] = build_vsepr_molecule(steric_number, lone_pairs, bond_length=length, central_symbol=central, ligand_symbol=ligand)
    print(f"{name}: {VSEPRGeometry(steric_number, lone_pairs).shape_name}")

# %%
# Each shape's characteristic bond angle is the *exact* idealized
# polyhedron angle -- 109.47 degrees for every tetrahedral-parent shape,
# 90/180 degrees for the octahedral-parent square-planar XeF4:

methane = molecules["methane (AX4E0)"]
water = molecules["water (AX2E2)"]
xef4 = molecules["xenon tetrafluoride (AX4E2)"]

print(f"CH4 H-C-H angle: {methane.bond_angle(1, 0, 2):.4f} degrees")
print(f"H2O H-O-H angle: {water.bond_angle(1, 0, 2):.4f} degrees")
print(f"XeF4 F-Xe-F (trans) angle: {xef4.bond_angle(1, 0, 2):.4f} degrees")
print(f"XeF4 F-Xe-F (cis) angle: {xef4.bond_angle(1, 0, 3):.4f} degrees")

# %%
# Plot all four geometries side by side:

fig = plt.figure(figsize=(11, 9))
for i, (name, molecule) in enumerate(molecules.items(), start=1):
    ax = fig.add_subplot(2, 2, i, projection="3d")
    plot_molecule_3d(molecule, ax=ax)
    ax.set_title(name, fontsize=10)
fig.tight_layout()
plt.show()
