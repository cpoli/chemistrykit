r"""
Periodic boundary conditions and the minimum-image convention
================================================================

Metropolis and coworkers (1953) avoided the walls of a finite box by
surrounding it with periodic copies of itself: a particle leaving one
face re-enters through the opposite one
(:func:`~chemistrykit.md.utils.pbc.wrap_positions`), and each pair
interacts through whichever periodic image of the partner is nearest
(:func:`~chemistrykit.md.utils.pbc.minimum_image_displacement`):

.. math::

   \Delta\vec r = \vec r_i - \vec r_j - L\,\mathrm{round}\!\left(\frac{\vec r_i - \vec r_j}{L}\right).

The left panel shows a two-dimensional box with its eight neighboring
images and, for one particle, the line to the nearest image of each
partner; the right panel follows a particle drifting across the box,
wrapped back in each time it crosses an edge.
"""

# %%
import matplotlib.pyplot as plt
import numpy as np

from chemistrykit.md.utils.pbc import minimum_image_displacement, wrap_positions

L = 10.0
rng = np.random.default_rng(3)
positions = rng.uniform(0.0, L, size=(6, 2))
positions[0] = [1.0, 8.8]  # the reference particle, near a corner

fig, axes = plt.subplots(1, 2, figsize=(11, 5.2))
ax = axes[0]
for sx in (-1, 0, 1):
    for sy in (-1, 0, 1):
        shift = np.array([sx, sy]) * L
        primary = sx == 0 and sy == 0
        ax.add_patch(plt.Rectangle(shift, L, L, fill=False, edgecolor="black" if primary else "lightgray", linewidth=1.5 if primary else 0.8))
        ax.scatter(*(positions + shift).T, s=25, color="steelblue" if primary else "lightsteelblue", zorder=3)
ax.scatter(*positions[0], s=60, color="firebrick", zorder=4, label="reference particle")
for j in range(1, len(positions)):
    d = minimum_image_displacement(positions[j], positions[0], L)
    ax.plot(*np.stack([positions[0], positions[0] + d]).T, color="firebrick", linewidth=1.0)
ax.set_xlim(-L, 2 * L)
ax.set_ylim(-L, 2 * L)
ax.set_aspect("equal")
ax.set_title("Nearest images of each partner")
ax.legend(loc="lower left")

# %%
# A particle moving at constant velocity through the periodic box. The
# wrapped coordinate jumps back each time it crosses a face, while the
# motion itself is unbroken:

t = np.linspace(0.0, 30.0, 600)
unwrapped = np.array([2.0, 3.0]) + np.outer(t, [0.9, 0.4])
wrapped = wrap_positions(unwrapped, L)
axes[1].plot(t, unwrapped[:, 0], color="gray", linestyle=":", label="x, unwrapped")
axes[1].plot(t, wrapped[:, 0], color="steelblue", label="x, wrapped into [0, L)")
axes[1].set_xlabel("t")
axes[1].set_ylabel("x")
axes[1].set_title("A particle crossing the periodic boundary")
axes[1].legend()
fig.tight_layout()

# %%
# Two particles on opposite edges of the box are actually neighbors:

r_i, r_j = np.array([0.5, 5.0]), np.array([9.5, 5.0])
print("plain displacement:        ", r_i - r_j)
print("minimum-image displacement:", minimum_image_displacement(r_i, r_j, L))

plt.show()
