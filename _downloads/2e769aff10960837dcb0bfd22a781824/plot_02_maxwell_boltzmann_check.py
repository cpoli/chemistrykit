r"""
Cross-checking an MD trajectory against the Maxwell-Boltzmann distribution
=============================================================================

A direct link between :mod:`chemistrykit.md` and :mod:`chemistrykit.statmech`:
molecular dynamics produces an actual ensemble of particle velocities,
and statistical mechanics predicts what their distribution *should* look
like at equilibrium. After a short equilibration run, an
:class:`~chemistrykit.md.systems.lj_fluid.LJFluid`'s particle speeds are
histogrammed and compared against
:class:`chemistrykit.statmech.MaxwellBoltzmannSpeedDistribution` evaluated
at the simulation's own (instantaneous, measured) temperature -- no
fitting involved.
"""

# %%
import matplotlib.pyplot as plt
import numpy as np

from chemistrykit.md.systems.lj_fluid import LJFluid
from chemistrykit.md.visualizers.md_plots import plot_speed_distribution
from chemistrykit.statmech import MaxwellBoltzmannSpeedDistribution

fluid = LJFluid.from_lattice(n_per_side=7, density=0.5, temperature=1.5, cutoff=2.5, rng=1)
fluid.run(dt=0.002, n_steps=1000, sample_every=1000)  # brief equilibration, discard the trajectory

speeds = np.linalg.norm(fluid.velocities, axis=1)
T_measured = fluid.temperature()

fig, ax = plt.subplots(figsize=(7, 5))
plot_speed_distribution(speeds, mass=1.0, temperature=T_measured, ax=ax, n_bins=30)
ax.set_title(f"MD speed histogram vs. Maxwell-Boltzmann (T = {T_measured:.3f})")
fig.tight_layout()

# %%
# The simulated mean speed agrees with the analytic Maxwell-Boltzmann
# prediction at the same (measured) temperature to within a few percent
# -- the residual gap is ordinary finite-N sampling noise, not a
# systematic discrepancy:

distribution = MaxwellBoltzmannSpeedDistribution(mass=1.0, temperature=T_measured, k_b=1.0)
print(f"simulated mean speed:  {speeds.mean():.4f}")
print(f"predicted mean speed:  {distribution.mean_speed():.4f}")
print(f"relative difference:   {abs(speeds.mean() - distribution.mean_speed()) / distribution.mean_speed():.2%}")

plt.show()
