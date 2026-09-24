r"""
Schoenflies point-group symbols for water, ammonia, methane, and carbon dioxide
================================================================================

Arthur Schoenflies named every point group after the symmetry elements that
generate it: a principal axis :math:`C_n`, extra :math:`C_2` axes
(:math:`D_n`), mirror planes :math:`\sigma_v` or :math:`\sigma_h`, an
inversion centre :math:`i`, and the special cubic groups :math:`T_d` and
:math:`O_h`. This example builds four molecules from their real 3D
geometry, lets
:func:`~chemistrykit.structure.systems.point_group.determine_point_group`
find the symmetry elements by testing the coordinates directly, and shows
how the elements it finds spell out each Schoenflies symbol: :math:`C_{2v}`,
:math:`C_{3v}`, :math:`T_d` and :math:`D_{\infty h}`.
"""

# %%
import matplotlib.pyplot as plt
import numpy as np

from chemistrykit.structure.core.base_system import Molecule
from chemistrykit.structure.systems.point_group import determine_point_group
from chemistrykit.structure.visualizers.structure_plots import plot_molecule_3d

# Water: r(O-H) = 0.958 A, angle(H-O-H) = 104.5 deg.
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

# Ammonia: trigonal pyramid, r(N-H) = 1.012 A, the three H atoms 0.38 A
# below the nitrogen.
r_nh, drop = 1.012, 0.38
rho = np.sqrt(r_nh**2 - drop**2)
phis = np.radians([90.0, 210.0, 330.0])
ammonia = Molecule(
    symbols=["N", "H", "H", "H"],
    coordinates=np.vstack([[0.0, 0.0, 0.0], np.column_stack([rho * np.cos(phis), rho * np.sin(phis), -drop * np.ones(3)])]),
    bonds=[(0, 1), (0, 2), (0, 3)],
)

# Methane: regular tetrahedron, r(C-H) = 1.09 A.
verts = np.array([[1.0, 1.0, 1.0], [1.0, -1.0, -1.0], [-1.0, 1.0, -1.0], [-1.0, -1.0, 1.0]])
verts = verts / np.linalg.norm(verts[0]) * 1.09
methane = Molecule(symbols=["C", "H", "H", "H", "H"], coordinates=np.vstack([[0.0, 0.0, 0.0], verts]), bonds=[(0, k) for k in range(1, 5)])

# Carbon dioxide: linear, r(C=O) = 1.16 A.
co2 = Molecule(symbols=["O", "C", "O"], coordinates=[[0.0, 0.0, -1.16], [0.0, 0.0, 0.0], [0.0, 0.0, 1.16]], bonds=[(0, 1), (1, 2)])

molecules = {"water": water, "ammonia": ammonia, "methane": methane, "carbon dioxide": co2}

# %%
# Read the Schoenflies symbol off the elements that were found. The
# subscript number is the order of the principal axis, ``v`` means mirror
# planes that contain that axis, and ``d`` / ``h`` mark the diagonal and
# horizontal planes of the higher groups:

results = {}
for name, molecule in molecules.items():
    res = determine_point_group(molecule)
    results[name] = res
    print(
        f"{name:15s} -> {res.group_name:8s} principal C_n: n={res.principal_axis_order}, "
        f"mirror planes: {res.n_mirror_planes}, sigma_v: {res.has_sigma_v}, "
        f"sigma_h: {res.has_sigma_h}, inversion centre: {res.has_inversion_center}, C3 axes: {res.n_c3_axes}"
    )

assert [results[k].group_name for k in molecules] == ["C2v", "C3v", "Td", "D_inf_h"]

# %%
# Plot each molecule with the symbol :func:`determine_point_group` assigned:

fig = plt.figure(figsize=(14, 4))
for i, (name, molecule) in enumerate(molecules.items(), start=1):
    ax = fig.add_subplot(1, 4, i, projection="3d")
    plot_molecule_3d(molecule, ax=ax)
    ax.set_title(f"{name}: {results[name].group_name}", fontsize=10)
fig.suptitle("Schoenflies symbols found from the geometry")
fig.tight_layout()
plt.show()
