r"""
Werner's coordination theory: counting isomers to prove the octahedron
=======================================================================

In 1893 Alfred Werner proposed that a six-coordinate metal complex places
its ligands at the corners of an octahedron. No experiment could see the
shape directly, so he argued from isomer counts. A complex
:math:`\mathrm{MA_4B_2}` such as :math:`[\mathrm{Co(NH_3)_4Cl_2}]^+` has only
two forms (green *trans* "praseo" and violet *cis* "violeo"). An octahedron
gives exactly two ways to place the two B ligands; a flat hexagon or a
trigonal prism gives three. For :math:`\mathrm{MA_3B_3}` the octahedron
again gives two (*fac* and *mer*).

This example counts those isomers by brute force. For each candidate shape
it finds every permutation of the six vertices that a rotation or
reflection of the shape can produce, then counts the distinct ways to
place the B ligands up to those symmetry operations. Mirror images are
counted as one here, because Werner compared these geometric isomers.
Finally it builds the octahedral complex and checks that it has
:math:`O_h` symmetry.
"""

# %%
import itertools

import matplotlib.pyplot as plt
import numpy as np

from chemistrykit.structure.core.base_system import Molecule
from chemistrykit.structure.systems.point_group import determine_point_group
from chemistrykit.structure.systems.vsepr import domain_positions
from chemistrykit.structure.visualizers.structure_plots import plot_molecule_3d

hexagon = np.array([[np.cos(a), np.sin(a), 0.0] for a in np.radians(np.arange(0, 360, 60))])
prism = np.array([[np.cos(a), np.sin(a), z] for z in (0.8, -0.8) for a in np.radians([90, 210, 330])])
shapes = {
    "octahedron (Werner)": domain_positions(6),
    "hexagonal plane": hexagon,
    "trigonal prism": prism,
}


def symmetry_permutations(vertices, tol=1e-6):
    """Vertex permutations that an orthogonal transform (rotation or reflection) of the shape achieves."""
    x = vertices - vertices.mean(axis=0)
    perms = []
    for perm in itertools.permutations(range(len(x))):
        y = x[list(perm)]
        # Best orthogonal map x -> y (orthogonal Procrustes); exact symmetries map with zero residual.
        u, _, vt = np.linalg.svd(y.T @ x)
        q = u @ vt
        if np.allclose(x @ q.T, y, atol=tol):
            perms.append(perm)
    return perms


def count_isomers(vertices, n_b):
    """Number of distinct placements of n_b B ligands on the vertices, up to symmetry."""
    perms = symmetry_permutations(vertices)
    seen, n_classes = set(), 0
    for combo in itertools.combinations(range(len(vertices)), n_b):
        key = frozenset(combo)
        if key in seen:
            continue
        n_classes += 1
        for perm in perms:
            seen.add(frozenset(perm[i] for i in combo))
    return n_classes


counts = {name: (count_isomers(v, 2), count_isomers(v, 3)) for name, v in shapes.items()}
for name, (n_ma4b2, n_ma3b3) in counts.items():
    print(f"{name:20s}: MA4B2 -> {n_ma4b2} isomers, MA3B3 -> {n_ma3b3} isomers")
print("observed for cobalt(III) ammines: 2 and 2 -> only the octahedron fits")
assert counts["octahedron (Werner)"] == (2, 2)

# %%
# The octahedral complex itself: an :math:`\mathrm{MA_6}` ion such as
# :math:`[\mathrm{Co(NH_3)_6}]^{3+}` (ligands drawn as single atoms). Its
# symmetry elements add up to :math:`O_h`, including the four
# three-fold axes through opposite faces.

bond_length = 1.97  # angstrom, typical Co-N distance
ligands = domain_positions(6) * bond_length
complex_ion = Molecule(symbols=["Co"] + ["N"] * 6, coordinates=np.vstack([[0.0, 0.0, 0.0], ligands]), bonds=[(0, i) for i in range(1, 7)])
result = determine_point_group(complex_ion)
print(f"\nMA6 point group: {result.group_name}, C3 axes found: {result.n_c3_axes}, inversion centre: {result.has_inversion_center}")
assert result.group_name == "Oh"

# %%
# Draw the two MA4B2 isomers Werner separated, cis and trans:

first = int(np.argmax(ligands[:, 0]))  # ligand on +x
partners = {"cis-[MA4B2]": int(np.argmax(ligands[:, 1])), "trans-[MA4B2]": int(np.argmin(ligands[:, 0]))}

fig = plt.figure(figsize=(10, 5))
for k, (label, partner) in enumerate(partners.items(), start=1):
    symbols = ["Co"] + ["Cl" if i in (first, partner) else "N" for i in range(6)]
    mol = Molecule(symbols=symbols, coordinates=complex_ion.coordinates, bonds=complex_ion.bonds)
    ax = fig.add_subplot(1, 2, k, projection="3d")
    plot_molecule_3d(mol, ax=ax)
    ax.set_title(label)
fig.tight_layout()
plt.show()
