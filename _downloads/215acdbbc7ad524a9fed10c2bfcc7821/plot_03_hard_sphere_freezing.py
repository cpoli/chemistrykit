r"""
Purely repulsive (WCA) particles: packing order without attraction
=====================================================================

Truncating the Lennard-Jones potential at its own minimum,
:math:`r=2^{1/6}\sigma` (:attr:`~chemistrykit.md.systems.lj_fluid.LennardJones.r_min`),
and shifting it to be continuous there removes the attractive tail
entirely and leaves a purely repulsive potential -- the
Weeks-Chandler-Andersen (WCA) reference system, and (in reduced units) a
smoothed stand-in for a fluid of hard spheres. Comparing the radial
distribution function g(r) of this purely repulsive
:class:`~chemistrykit.md.systems.lj_fluid.LJFluid` at low and high density
shows that pronounced short-range positional order can appear from
*packing alone*, with no attractive force anywhere in the potential --
the essential physics behind Alder and Wainwright's 1957 discovery that
hard spheres alone can freeze.
"""

# %%
import matplotlib.pyplot as plt

from chemistrykit.md.systems.lj_fluid import LennardJones, LJFluid
from chemistrykit.md.systems.thermostats import VelocityRescalingThermostat
from chemistrykit.md.visualizers.md_plots import plot_radial_distribution_function

r_wca_cutoff = LennardJones().r_min  # 2**(1/6) sigma: purely repulsive LJ (WCA)
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

ax.set_title("g(r) of a purely repulsive (WCA) fluid: dilute vs. dense")
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
