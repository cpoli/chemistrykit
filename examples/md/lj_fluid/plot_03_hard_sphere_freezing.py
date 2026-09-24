r"""
Alder and Wainwright's hard spheres: packing-driven order with no attraction
===============================================================================

Alder and Wainwright's 1957 computer experiments showed that particles
with *no attraction at all* -- hard spheres -- order into a crystal-like
arrangement once they are packed densely enough: freezing driven by
excluded volume (entropy), not by energy. This package integrates smooth
forces rather than hard-sphere collision events, so the hard spheres are
approximated by steeply repulsive particles: the Lennard-Jones potential
cut off at its own minimum :math:`r=2^{1/6}\sigma`
(:attr:`~chemistrykit.md.systems.lj_fluid.LennardJones.r_min`) and shifted,
which leaves only a repulsive wall. Comparing the radial distribution
function g(r) of this purely repulsive
:class:`~chemistrykit.md.systems.lj_fluid.LJFluid` at low density and at
a density near hard-sphere freezing shows sharp positional order
appearing from packing alone.
"""

# %%
import matplotlib.pyplot as plt

from chemistrykit.md.systems.lj_fluid import LennardJones, LJFluid
from chemistrykit.md.systems.thermostats import VelocityRescalingThermostat
from chemistrykit.md.visualizers.md_plots import plot_radial_distribution_function

r_wca_cutoff = LennardJones().r_min  # 2**(1/6) sigma: purely repulsive, hard-sphere-like
T_target = 1.0
dt, n_steps = 0.003, 4000

fig, ax = plt.subplots(figsize=(7, 5))
peak_heights = {}
for density, color in [(0.30, "steelblue"), (1.05, "darkorange")]:
    fluid = LJFluid.from_lattice(
        n_per_side=6,
        density=density,
        temperature=T_target,
        cutoff=r_wca_cutoff,
        rng=0,
        skin=0.3,
        rebuild_every=10,
    )
    thermostat = VelocityRescalingThermostat(target_temperature=T_target, interval=5)
    fluid.run(dt=dt, n_steps=n_steps, thermostat=thermostat, sample_every=200)
    r, g = fluid.radial_distribution_function(n_bins=120)
    plot_radial_distribution_function(r, g, ax=ax, color=color, label=f"density = {density}")
    peak_heights[density] = float(g.max())

ax.set_title("Hard-sphere-like particles: dilute gas vs. densely packed")
ax.legend()
fig.tight_layout()

# %%
# The dilute configuration's g(r) is essentially featureless beyond the
# repulsive core -- close to the flat g(r)=1 of an ideal gas, since
# particles rarely encounter one another. At more than three times the
# density, with exactly the same purely repulsive potential (no
# attractive term added), the first-neighbor peak becomes far taller and
# sharper: particles are packed closely enough that excluded volume alone
# forces pronounced local order, exactly the entropic (packing-driven)
# ordering mechanism Alder and Wainwright found in hard-sphere computer
# experiments -- reproducing it here only qualitatively, since a
# convincing first-order fluid-solid transition needs a much larger
# system and a much longer run than this short pedagogical example:

for density, height in peak_heights.items():
    print(f"density = {density}: g(r) first-peak height = {height:.2f}")

plt.show()
