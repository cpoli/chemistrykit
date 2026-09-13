r"""
Verifying the speed distribution: Stern's molecular-beam experiment
=======================================================================

Otto Stern's 1920 molecular-beam apparatus gave the first direct
experimental measurement of a gas's actual molecular speeds, rather than
only their statistical consequences (pressure, diffusion, viscosity) -- a
rotating-drum time-of-flight technique later refined by Zartman and Ko
(1930-1934) into a quantitative confirmation of the full Maxwell-Boltzmann
functional form, not just its mean.
:meth:`~chemistrykit.statmech.MaxwellBoltzmannSpeedDistribution.sample`
draws exactly such a synthetic "beam" of molecular speeds; binning it into
a histogram and comparing it to the exact analytic pdf reproduces the kind
of check Stern, and later Zartman and Ko, made directly against a real gas.
"""

# %%
import matplotlib.pyplot as plt
import numpy as np

from chemistrykit.statmech import MaxwellBoltzmannSpeedDistribution

mass = 2.18e-25  # potassium-like, kg -- close to the alkali-metal beams Stern's and Zartman/Ko's apparatus actually used
temperature = 470.0  # K, a plausible oven temperature for an alkali-metal beam source

dist = MaxwellBoltzmannSpeedDistribution(mass=mass, temperature=temperature)
sampled_speeds = dist.sample(200_000, rng=0)

fig, ax = plt.subplots(figsize=(7, 5))
ax.hist(sampled_speeds, bins=80, density=True, alpha=0.5, color="steelblue", label="synthetic 'beam' (200,000 draws)")
v = np.linspace(0.0, float(sampled_speeds.max()), 400)
ax.plot(v, dist.pdf(v), color="crimson", linewidth=2, label="exact Maxwell-Boltzmann pdf")
ax.set_xlabel("speed (m/s)")
ax.set_ylabel("probability density")
ax.set_title("Sampling the Maxwell-Boltzmann distribution (a Stern-type beam measurement)")
ax.legend()
fig.tight_layout()

# %%
# The sampled mean speed converges to the exact analytic mean speed as
# the number of "molecules" in the synthetic beam grows -- exactly the
# statistical convergence a real, finite molecular-beam measurement is
# limited by:

for n in (100, 10_000, 1_000_000):
    speeds = dist.sample(n, rng=0)
    print(f"n = {n:>9,}: sampled mean speed = {speeds.mean():.2f} m/s  (exact: {dist.mean_speed():.2f} m/s)")

plt.show()
