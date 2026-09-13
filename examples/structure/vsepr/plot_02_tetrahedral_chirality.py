r"""
Van't Hoff and Le Bel's tetrahedral carbon: chirality from real 3D coordinates
==================================================================================

Jacobus van't Hoff and Joseph Le Bel independently proposed in 1874 that
a carbon atom's four bonds point toward the vertices of a regular
tetrahedron -- and that a carbon bonded to four *different* substituents
therefore has no internal mirror symmetry: its mirror image cannot be
superimposed on itself by any rotation, giving two distinct, non-
interconvertible molecules (enantiomers). This builds a genuinely
tetrahedral AX4 center -- bromochlorofluoromethane, CHFClBr, the textbook
minimal chiral carbon -- from
:func:`~chemistrykit.structure.systems.vsepr.domain_positions`'s real 3D
vertex coordinates, confirms the exact 109.47-degree bond angles van't
Hoff and Le Bel's model predicts, and then demonstrates the chirality
itself numerically: a scalar-triple-product handedness invariant flips
sign under reflection but is exactly preserved under any proper rotation
(built with :func:`~chemistrykit.structure.utils.symmetry_ops.rotation_matrix`),
so no rotation can ever turn the molecule into its own mirror image.
"""

# %%
import matplotlib.pyplot as plt
import numpy as np

from chemistrykit.structure.core.base_system import Molecule
from chemistrykit.structure.systems.vsepr import domain_positions
from chemistrykit.structure.utils.symmetry_ops import rotation_matrix
from chemistrykit.structure.visualizers.structure_plots import plot_molecule_3d

# Four different substituents (F, Cl, Br, H) placed at the four vertices
# of a regular tetrahedron -- exactly van't Hoff and Le Bel's proposed
# arrangement of a carbon atom's four bonds.
bond_length = 1.5  # angstrom, illustrative
vertices = domain_positions(4) * bond_length
symbols = ["C", "F", "Cl", "Br", "H"]
coordinates = np.vstack([np.zeros((1, 3)), vertices])
molecule = Molecule(symbols=symbols, coordinates=coordinates, bonds=[(0, 1), (0, 2), (0, 3), (0, 4)])

print("Bond angles at the central carbon (van't Hoff/Le Bel predict exactly 109.4712 degrees for every pair):")
for i in range(1, 5):
    for j in range(i + 1, 5):
        print(f"  {symbols[i]}-C-{symbols[j]}: {molecule.bond_angle(i, 0, j):.4f} degrees")

# %%
# A handedness invariant: the signed volume (scalar triple product) of
# the vectors from the central carbon to F, Cl, and Br, in that fixed
# order. This changes sign under a reflection (which inverts handedness)
# but -- being built entirely from dot and cross products of vectors
# attached to the rigid molecule -- is exactly unchanged by any *proper*
# rotation of the whole structure, since rotations preserve both lengths
# and orientation:


def handedness(coords: np.ndarray) -> float:
    """Signed volume of the (F, Cl, Br) vectors from the central atom (index 0)."""
    v_f, v_cl, v_br = coords[1] - coords[0], coords[2] - coords[0], coords[3] - coords[0]
    return float(np.dot(v_f, np.cross(v_cl, v_br)))


original_handedness = handedness(molecule.coordinates)
print(f"\nHandedness of the original molecule: {original_handedness:+.4f}")

# Mirror image: reflect every coordinate through the yz-plane (x -> -x).
mirror_coordinates = molecule.coordinates.copy()
mirror_coordinates[:, 0] *= -1.0
mirror_molecule = Molecule(symbols=symbols, coordinates=mirror_coordinates, bonds=molecule.bonds)
mirror_handedness = handedness(mirror_molecule.coordinates)
print(f"Handedness of the mirror image:      {mirror_handedness:+.4f}  (sign flipped by the reflection)")

# A proper rotation -- no matter which axis or angle -- leaves the
# handedness invariant unchanged, confirming it genuinely measures
# handedness rather than some artifact of the coordinate choice:
R = rotation_matrix(axis=[0.3, 0.7, 1.0], angle=1.9)
rotated_coordinates = molecule.coordinates @ R.T
rotated_handedness = handedness(rotated_coordinates)
print(f"Handedness after an arbitrary proper rotation: {rotated_handedness:+.4f}  (unchanged)")

assert np.sign(original_handedness) == np.sign(rotated_handedness)
assert np.sign(original_handedness) != np.sign(mirror_handedness)
print("\n=> no rotation can superimpose the molecule on its mirror image: CHFClBr is chiral.")

# %%
# Plot both enantiomers side by side:

fig = plt.figure(figsize=(10, 5))
ax1 = fig.add_subplot(1, 2, 1, projection="3d")
plot_molecule_3d(molecule, ax=ax1)
ax1.set_title("CHFClBr (R or S)")
ax2 = fig.add_subplot(1, 2, 2, projection="3d")
plot_molecule_3d(mirror_molecule, ax=ax2)
ax2.set_title("Mirror image (the other enantiomer)")
fig.tight_layout()
plt.show()
