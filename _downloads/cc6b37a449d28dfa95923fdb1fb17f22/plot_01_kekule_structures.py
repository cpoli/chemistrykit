r"""
Kekulé's benzene: enumerating the Kekulé structures of aromatic hydrocarbons
============================================================================

August Kekulé (1865) proposed that benzene's six carbons form a ring with
alternating single and double bonds, and soon after that the two
alternating arrangements interconvert, which explained why benzene has
only one ortho-disubstituted isomer. A Kekulé structure is a choice of
double bonds that gives every carbon exactly one; in graph language, a
perfect matching of the carbon skeleton.

This example uses
:func:`~chemistrykit.structure.systems.kekule.kekule_structures` to list
benzene's two Kekulé structures and naphthalene's three, and counts them
for larger benzenoid hydrocarbons: a linear acene with :math:`r` rings has
:math:`r+1`, while angular phenanthrene has five, one reason it is more
stable than its linear isomer anthracene (four).
"""

# %%
import matplotlib.pyplot as plt
import numpy as np

from chemistrykit.structure.systems.kekule import count_kekule_structures, kekule_structures


def hexagon(cx, cy, start=90.0):
    """Six vertex positions of a unit hexagon centred at (cx, cy)."""
    a = np.radians(start + 60.0 * np.arange(6))
    return np.column_stack([cx + np.cos(a), cy + np.sin(a)])


# Benzene: atoms 0-5 around the ring.
benzene_xy = hexagon(0.0, 0.0)
benzene_bonds = [(i, (i + 1) % 6) for i in range(6)]

# Naphthalene: two hexagons sharing the bond between atoms 4 and 9.
naphthalene_bonds = [(0, 1), (1, 2), (2, 3), (3, 4), (4, 9), (9, 0), (4, 5), (5, 6), (6, 7), (7, 8), (8, 9)]
left = hexagon(-np.sqrt(3) / 2, 0.0, start=30.0)
right = hexagon(np.sqrt(3) / 2, 0.0, start=30.0)
naphthalene_xy = np.array([left[1], left[2], left[3], left[4], left[5], right[4], right[5], right[0], right[1], right[2]])

benzene_k = kekule_structures(6, benzene_bonds)
naphthalene_k = kekule_structures(10, naphthalene_bonds)
print(f"benzene: {len(benzene_k)} Kekule structures -> {benzene_k}")
print(f"naphthalene: {len(naphthalene_k)} Kekule structures")
assert len(benzene_k) == 2 and len(naphthalene_k) == 3

# %%
# Each benzene bond is double in exactly one of the two structures, so
# averaging them gives every bond the same order, 1.5. That averaging is
# how Kekulé's oscillation hypothesis makes the six bonds equivalent.

double_fraction = [sum((min(b), max(b)) in s for s in benzene_k) / len(benzene_k) for b in benzene_bonds]
print(f"fraction of Kekule structures in which each benzene bond is double: {double_fraction}")


# %%
# Kekulé structure counts of larger benzenoid hydrocarbons:


def acene(n_rings):
    """Linear acene skeleton with n_rings fused rings (4n+2 carbons)."""
    top = list(range(2 * n_rings + 1))
    bottom = list(range(2 * n_rings + 1, 4 * n_rings + 2))
    bonds = [(a, a + 1) for a in top[:-1]] + [(b, b + 1) for b in bottom[:-1]]
    bonds += [(top[2 * k], bottom[2 * k]) for k in range(n_rings + 1)]
    return 4 * n_rings + 2, bonds


phenanthrene = [(0, 1), (1, 2), (2, 3), (3, 4), (4, 5), (5, 0), (4, 6), (6, 7), (7, 8), (8, 9), (9, 3), (8, 10), (10, 11), (11, 12), (12, 13), (13, 9)]
counts = {name: count_kekule_structures(*acene(r)) for r, name in enumerate(["benzene", "naphthalene", "anthracene", "tetracene", "pentacene"], start=1)}
counts["phenanthrene"] = count_kekule_structures(14, phenanthrene)
for name, k in counts.items():
    print(f"{name:12s}: K = {k}")
assert counts["anthracene"] == 4 and counts["phenanthrene"] == 5

# %%
# Draw benzene's two and naphthalene's three Kekulé structures:


def draw(ax, xy, bonds, doubles, title):
    double_set = {tuple(sorted(b)) for b in doubles}
    for i, j in bonds:
        is_double = tuple(sorted((i, j))) in double_set
        ax.plot(*zip(xy[i], xy[j], strict=True), color="crimson" if is_double else "black", linewidth=4.0 if is_double else 1.5)
    ax.set_aspect("equal")
    ax.axis("off")
    ax.set_title(title, fontsize=10)


fig, axes = plt.subplots(1, 5, figsize=(15, 3.2))
for k, s in enumerate(benzene_k):
    draw(axes[k], benzene_xy, benzene_bonds, s, f"benzene, structure {k + 1}")
for k, s in enumerate(naphthalene_k):
    draw(axes[2 + k], naphthalene_xy, naphthalene_bonds, s, f"naphthalene, structure {k + 1}")
fig.suptitle("Kekulé structures (double bonds drawn thick red)")
fig.tight_layout()
plt.show()
