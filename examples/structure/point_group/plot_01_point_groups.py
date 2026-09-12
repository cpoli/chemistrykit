r"""
Point-group determination for water, ammonia, methane, and carbon dioxide
==========================================================================

Builds each molecule's actual 3D geometry, then runs
:func:`~chemistrykit.structure.systems.point_group.determine_point_group`
-- which finds rotation axes, mirror planes, and an inversion center by
direct geometric testing, not by recognizing the molecular formula -- to
recover the textbook point-group assignments C2v, C3v, Td, and D_inf_h.
"""

# %%
import matplotlib.pyplot as plt
import numpy as np

from chemistrykit.structure.core.base_system import Molecule
from chemistrykit.structure.systems.point_group import determine_point_group, get_character_table
from chemistrykit.structure.visualizers.structure_plots import plot_molecule_3d

# Water: experimental geometry, r(O-H) = 0.958 A, angle(H-O-H) = 104.5 deg.
angle = np.radians(104.5)
r = 0.958
water = Molecule(
    symbols=["O", "H", "H"],
    coordinates=[
        [0.0, 0.0, 0.0],
        [r * np.sin(angle / 2), 0.0, r * np.cos(angle / 2)],
        [-r * np.sin(angle / 2), 0.0, r * np.cos(angle / 2)],
    ],
    bonds=[(0, 1), (0, 2)],
)

# Methane: perfect tetrahedron, r(C-H) = 1.09 A.
verts = np.array([[1.0, 1.0, 1.0], [1.0, -1.0, -1.0], [-1.0, 1.0, -1.0], [-1.0, -1.0, 1.0]])
verts = verts / np.linalg.norm(verts[0]) * 1.09
methane = Molecule(symbols=["C", "H", "H", "H", "H"], coordinates=np.vstack([[0.0, 0.0, 0.0], verts]))

# Carbon dioxide: linear, r(C=O) = 1.16 A.
co2 = Molecule(symbols=["O", "C", "O"], coordinates=[[0.0, 0.0, -1.16], [0.0, 0.0, 0.0], [0.0, 0.0, 1.16]], bonds=[(0, 1), (1, 2)])

molecules = {"water": water, "methane": methane, "carbon dioxide": co2}

for name, molecule in molecules.items():
    result = determine_point_group(molecule)
    print(
        f"{name}: {result.group_name}  "
        f"(principal axis order={result.principal_axis_order}, "
        f"sigma_v={result.has_sigma_v}, sigma_h={result.has_sigma_h}, "
        f"i={result.has_inversion_center})"
    )

# %%
# The C2v character table -- 4 one-dimensional irreducible
# representations, since C2v is abelian:

c2v = get_character_table("C2v")
print(f"\nC2v operations: {c2v.operations}")
for irrep in c2v.irreps:
    print(f"  {irrep}: {c2v.characters[c2v.irreps.index(irrep)]}")

# %%
# Plot each molecule:

fig = plt.figure(figsize=(12, 4))
for i, (name, molecule) in enumerate(molecules.items(), start=1):
    ax = fig.add_subplot(1, 3, i, projection="3d")
    plot_molecule_3d(molecule, ax=ax)
    ax.set_title(f"{name}: {determine_point_group(molecule).group_name}", fontsize=10)
fig.tight_layout()
plt.show()
