r"""
An NVE Lennard-Jones fluid: energy conservation and g(r)
==========================================================

A monatomic Lennard-Jones fluid (:class:`~chemistrykit.md.systems.lj_fluid.LJFluid`)
started on a simple-cubic lattice with Maxwell-Boltzmann velocities
(:meth:`~chemistrykit.md.systems.lj_fluid.LJFluid.from_lattice`), integrated
with no thermostat (a microcanonical, NVE run) via the shared
velocity-Verlet integrator. Two standard checks: total energy should be
constant to within integration error, and the radial distribution
function g(r) should show liquid-like structure (a first peak near the
Lennard-Jones minimum) rather than the flat g(r)=1 of an ideal gas.
"""

# %%
import matplotlib.pyplot as plt

from chemistrykit.md.systems.lj_fluid import LJFluid
from chemistrykit.md.visualizers.md_plots import plot_energy_conservation, plot_radial_distribution_function

fluid = LJFluid.from_lattice(n_per_side=6, density=0.7, temperature=1.0, cutoff=2.5, rng=0)
result = fluid.run(dt=0.002, n_steps=2000, sample_every=20)

fig, axes = plt.subplots(1, 2, figsize=(11, 4.5))
plot_energy_conservation(result, ax=axes[0])
axes[0].set_title("Energy conservation (NVE)")

r, g = fluid.radial_distribution_function(n_bins=100)
plot_radial_distribution_function(r, g, ax=axes[1], color="steelblue")
axes[1].set_title("Radial distribution function g(r)")
fig.tight_layout()

# %%
# The relative energy drift over the whole run is small -- the hallmark
# of a symplectic integrator (velocity-Verlet) applied to a Hamiltonian
# system with no thermostat:

total = result.total_energy
drift = abs(total.max() - total.min()) / abs(total[0])
print(f"Relative energy drift over the run: {drift:.2e}")

# %%
# g(r) peaks close to the Lennard-Jones equilibrium separation
# :math:`r=2^{1/6}\sigma\approx1.122\sigma`, the signature of short-range
# liquid-like order (nearest-neighbor coordination shells) that a
# non-interacting ideal gas (g(r)=1 everywhere) would not show:

r_peak = r[g.argmax()]
print(f"g(r) peaks at r = {r_peak:.3f} sigma (LJ minimum at 2^(1/6) = {2.0 ** (1.0 / 6.0):.3f} sigma)")

plt.show()
