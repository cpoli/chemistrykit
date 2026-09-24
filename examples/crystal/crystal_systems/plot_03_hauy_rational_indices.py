r"""
Hauy's law of rational indices: Miller indices from face intercepts
====================================================================

Hauy (1784-1801) argued that every natural crystal face cuts the three
axes at intercepts in small whole-number ratios.
:func:`~chemistrykit.crystal.systems.crystal_systems.miller_indices_from_intercepts`
turns those rational intercepts into the integer Miller indices
:math:`(hkl)` by taking reciprocals and clearing fractions; the
:math:`(hkl)` triple then fixes each face's interplanar spacing.
"""

# %%
import math

import matplotlib.pyplot as plt
import numpy as np

from chemistrykit.crystal.systems.crystal_systems import miller_indices_from_intercepts
from chemistrykit.crystal.systems.xrd import d_spacing_cubic

faces = {
    "cube face": (1, math.inf, math.inf),
    "dodecahedron face": (1, 1, math.inf),
    "octahedron face": (1, 1, 1),
    "intercepts 1, 2, 3": (1, 2, 3),
    "intercepts 1/2, 1, inf": (0.5, 1, math.inf),
    "intercepts 2, 3, 6": (2, 3, 6),
}
a = 5.64  # NaCl, Angstrom
for name, intercepts in faces.items():
    hkl = miller_indices_from_intercepts(*intercepts)
    print(f"{name:24s} intercepts={str(intercepts):22s} -> (hkl)={hkl}   d={d_spacing_cubic(a, *hkl):.3f} A")

# %%
# The intercept ratios are rational, so the reciprocals clear to small
# integers: intercepts (2, 3, 6) give reciprocals (1/2, 1/3, 1/6) = (3, 2, 1)/6.
assert miller_indices_from_intercepts(2, 3, 6) == (3, 2, 1)

# %%
# Draw a few faces as planes cutting the unit cube at their intercepts.
fig = plt.figure(figsize=(9, 3.5))
for i, (hkl, color) in enumerate((((1, 0, 0), "C0"), ((1, 1, 0), "C1"), ((1, 1, 1), "C2"))):
    ax = fig.add_subplot(1, 3, i + 1, projection="3d")
    for s, e in [((0, 0, 0), (1, 0, 0)), ((0, 0, 0), (0, 1, 0)), ((0, 0, 0), (0, 0, 1))]:
        ax.plot(*zip(s, e, strict=True), "k-", lw=0.8)
    h, k, l = hkl
    g = np.linspace(0, 1, 2)
    if l != 0:
        # the (111) face is the triangle through the three unit intercepts
        ax.plot_trisurf(np.array([1, 0, 0]), np.array([0, 1, 0]), np.array([0, 0, 1]), color=color, alpha=0.5)
    elif k != 0:
        X, Z = np.meshgrid(g, g)
        ax.plot_surface(X, (1 - h * X) / k, Z, color=color, alpha=0.5)
    else:
        Y, Z = np.meshgrid(g, g)
        ax.plot_surface(np.ones_like(Y) / h, Y, Z, color=color, alpha=0.5)
    ax.set_title(f"({h}{k}{l})")
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.set_zlim(0, 1)
    ax.set_xlabel("a")
    ax.set_ylabel("b")
    ax.set_zlabel("c")
plt.tight_layout()
plt.show()
