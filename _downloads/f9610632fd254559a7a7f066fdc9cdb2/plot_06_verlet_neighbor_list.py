r"""
Verlet's neighbor list: a skin that lets one pair list serve many steps
=========================================================================

Checking every pair of :math:`N` particles costs :math:`N(N-1)/2`
distance evaluations per step. Verlet's 1967 trick: list the pairs within
:math:`r_c + r_\text{skin}` and reuse that list for several steps; a pair
outside the list can only come inside the true cutoff :math:`r_c` after
the two particles have closed the skin distance between them.
:class:`~chemistrykit.md.utils.neighbor_list.VerletNeighborList`
implements this. Here a Lennard-Jones fluid is run with a list that is
never rebuilt, and at every step the listed pairs are checked against the
exact set of pairs inside the cutoff
(:func:`~chemistrykit.md.utils.neighbor_list.build_neighbor_list`), for
lists built with no skin and with a skin of :math:`0.3\sigma`.
"""

# %%
import matplotlib.pyplot as plt
import numpy as np

from chemistrykit.md.systems.lj_fluid import LJFluid
from chemistrykit.md.utils.neighbor_list import build_neighbor_list

cutoff, dt, n_steps = 2.5, 0.005, 60
fluid = LJFluid.from_lattice(n_per_side=7, density=0.8, temperature=1.0, cutoff=cutoff, rng=0)
fluid.run(dt=dt, n_steps=400, sample_every=400)  # melt the lattice first
n = fluid.positions.shape[0]
L = fluid.box_length

lists = {skin: set(zip(*[a.tolist() for a in build_neighbor_list(fluid.positions, L, cutoff + skin)], strict=True)) for skin in (0.0, 0.3)}
missed = {skin: [] for skin in lists}
for _ in range(n_steps):
    fluid.step(dt)
    exact = set(zip(*[a.tolist() for a in build_neighbor_list(fluid.positions, L, cutoff)], strict=True))
    for skin, pairs in lists.items():
        missed[skin].append(len(exact - pairs))

fig, ax = plt.subplots(figsize=(7, 4.5))
steps = np.arange(1, n_steps + 1)
ax.plot(steps, missed[0.0], color="firebrick", label="list built with no skin")
ax.plot(steps, missed[0.3], color="steelblue", label=r"list built with skin $0.3\sigma$")
ax.set_xlabel("steps since the list was built")
ax.set_ylabel("interacting pairs missing from the list")
ax.set_title("A skinned Verlet list stays complete for many steps")
ax.legend()
fig.tight_layout()

# %%
# Without a skin the list goes stale after a single step. With a
# :math:`0.3\sigma` skin it still holds every interacting pair many steps
# later, while holding only a fraction of all :math:`N(N-1)/2` pairs --
# the saving Verlet's list delivers at every step it is reused:

first_miss = next((k + 1 for k, m in enumerate(missed[0.3]) if m > 0), None)
print(f"all pairs: {n * (n - 1) // 2}, skinned list: {len(lists[0.3])}, exact list at build: {len(lists[0.0])}")
print(f"no-skin list first misses a pair after {next(k + 1 for k, m in enumerate(missed[0.0]) if m > 0)} step(s)")
print(f"skinned list first misses a pair after: {first_miss if first_miss else f'more than {n_steps}'} steps")

plt.show()
