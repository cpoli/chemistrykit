r"""
Steno's law: constancy of interfacial angles
=============================================

Nicolas Steno (1669) found that quartz crystals, however unevenly grown,
always meet at the same angles between corresponding faces. The angle
between two faces depends only on their Miller indices -- not on how big
the faces are -- which
:func:`~chemistrykit.crystal.systems.crystal_systems.interplanar_angle_cubic`
computes for a cubic crystal from the face normals.
"""

# %%
import matplotlib.pyplot as plt
import numpy as np

from chemistrykit.crystal.systems.crystal_systems import interplanar_angle_cubic

pairs = [
    ((1, 0, 0), (0, 1, 0), "cube / cube"),
    ((1, 0, 0), (1, 1, 1), "cube / octahedron"),
    ((1, 1, 1), (1, 1, -1), "octahedron / octahedron"),
    ((1, 0, 0), (1, 1, 0), "cube / dodecahedron"),
    ((1, 1, 0), (1, 1, 1), "dodecahedron / octahedron"),
]
for hkl_1, hkl_2, label in pairs:
    print(f"{str(hkl_1):12s} {str(hkl_2):12s} {label:28s} {interplanar_angle_cubic(hkl_1, hkl_2):7.2f} deg")

# %%
# Steno's point: growing one face bigger than another (a distorted crystal)
# changes the face *sizes*, never the angles. Scaling the face normals --
# the same faces pushed outward by different amounts -- leaves every
# interfacial angle unchanged:
rng = np.random.default_rng(0)
for _ in range(3):
    s1, s2 = rng.uniform(0.5, 3.0, size=2)
    angle = np.degrees(np.arccos(np.dot(s1 * np.array([1, 0, 0]), s2 * np.array([1, 1, 1])) / (s1 * s2 * np.sqrt(3.0))))
    print(f"face distances {s1:.2f}, {s2:.2f}: cube/octahedron angle = {angle:.4f} deg")
    assert np.isclose(angle, interplanar_angle_cubic((1, 0, 0), (1, 1, 1)))

# %%
# A regular and a distorted cuboctahedron-like cross section (the {100}
# and {110} faces cut in the (001) plane): same angles, different shapes.
fig, axes = plt.subplots(1, 2, figsize=(9, 4.5))
normals = [np.array([np.cos(t), np.sin(t)]) for t in np.radians(np.arange(0, 360, 45))]
for ax, distances, title in (
    (axes[0], np.ones(8), "regular crystal"),
    (axes[1], np.array([1.0, 0.8, 1.6, 1.2, 1.3, 0.9, 1.1, 1.4]), "distorted crystal"),
):
    corners = []
    for i in range(8):
        n1, n2 = normals[i], normals[(i + 1) % 8]
        corners.append(np.linalg.solve(np.array([n1, n2]), np.array([distances[i], distances[(i + 1) % 8]])))
    corners = np.array(corners + [corners[0]])
    ax.plot(corners[:, 0], corners[:, 1], "k-")
    ax.fill(corners[:, 0], corners[:, 1], alpha=0.2)
    ax.set_aspect("equal")
    ax.set_title(f"{title}: every interior angle = 135 deg")
    ax.axis("off")
plt.tight_layout()
plt.show()
