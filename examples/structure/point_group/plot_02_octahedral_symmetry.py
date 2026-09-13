r"""
Octahedral symmetry: SF6, Werner's six-coordinate complexes, and the Oh character table
=============================================================================================

Alfred Werner's 1893 coordination theory proposed that a six-coordinate
metal complex's ligands sit at the vertices of an octahedron around the
central atom -- the structural model he supported entirely through
painstaking isomer-counting arguments, decades before any direct
structural technique existed to see it. Sulfur hexafluoride is a
main-group molecule built on exactly that same six-coordinate octahedral
skeleton, and its symmetry -- four independent three-fold (C3) body-
diagonal axes on top of the three obvious four-fold (C4) axes along the
S-F bonds -- makes it the standard worked example for the cubic Oh point
group. This also exercises the specific candidate-axis machinery
(:func:`~chemistrykit.structure.utils.symmetry_ops.candidate_axes`'s
three-way position-vector sums) needed to find those body-diagonal C3
axes at all: with the six ligands placed exactly on the cardinal axes, no
*pairwise* combination of two ligand vectors can ever point along a body
diagonal such as (1, 1, 1).
"""

# %%
import matplotlib.pyplot as plt

from chemistrykit.structure.core.base_system import Molecule
from chemistrykit.structure.systems.point_group import determine_point_group, get_character_table
from chemistrykit.structure.visualizers.structure_plots import plot_molecule_3d

bond_length = 1.56  # angstrom, typical S-F distance
sf6 = Molecule(
    symbols=["S", "F", "F", "F", "F", "F", "F"],
    coordinates=[
        [0.0, 0.0, 0.0],
        [bond_length, 0.0, 0.0],
        [-bond_length, 0.0, 0.0],
        [0.0, bond_length, 0.0],
        [0.0, -bond_length, 0.0],
        [0.0, 0.0, bond_length],
        [0.0, 0.0, -bond_length],
    ],
    bonds=[(0, i) for i in range(1, 7)],
)

result = determine_point_group(sf6)
print(f"SF6 point group: {result.group_name}")
print(f"Distinct C3 axes found (body diagonals): {result.n_c3_axes}  (Werner's octahedron needs exactly 4)")
print(f"Inversion center present: {result.has_inversion_center}  (every Oh complex is centrosymmetric)")
assert result.group_name == "Oh"
assert result.n_c3_axes >= 4

# %%
# The Oh character table -- looked up by Robert Mulliken's now-universal
# Mulliken-symbol notation (A/B/E/T + g/u + numeric subscripts), the same
# labeling convention every irrep in this module uses. Its Eg and T2g
# rows are, not coincidentally, the exact labels Hans Bethe's 1929
# crystal-field theory attaches to the two ways a transition-metal
# d-orbital set splits in an octahedral ligand field:

oh = get_character_table("Oh")
print(f"\nOh irreducible representations: {oh.irreps}")
print(f"Dimension of each irrep (character under E): {[oh.character(irrep, 'E') for irrep in oh.irreps]}")
print(f"T2g and Eg (the crystal-field d-orbital splitting labels) have dimensions {oh.character('T2g', 'E'):.0f} and {oh.character('Eg', 'E'):.0f}")

# %%
ax = plot_molecule_3d(sf6)
ax.set_title(f"SF6: {result.group_name} (Werner's octahedral six-coordination)")
plt.show()
