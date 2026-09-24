r"""
Smoluchowski's diffusion-limited reaction rate
=================================================

Smoluchowski (1917) asked how fast two species meet if they react on
first contact, so that diffusion alone limits the rate. Solving the
diffusion equation around one reactant gives

.. math::

    k_D = 4\pi D R^* N_A, \qquad
    k(t) = k_D\left(1 + \frac{R^*}{\sqrt{\pi D t}}\right)

(:func:`~chemistrykit.kinetics.systems.rate_theory.smoluchowski_rate_constant`,
:func:`~chemistrykit.kinetics.systems.rate_theory.smoluchowski_transient_rate_constant`).
With the Stokes-Einstein diffusion coefficient the molecular size cancels
and :math:`k_D = 8RT/3\eta`
(:func:`~chemistrykit.kinetics.systems.rate_theory.diffusion_limited_rate_constant`):
the ceiling on any bimolecular rate constant in solution depends only on
the solvent's viscosity.
"""

# %%
import matplotlib.pyplot as plt
import numpy as np

from chemistrykit.kinetics.systems.rate_theory import (
    diffusion_limited_rate_constant,
    smoluchowski_rate_constant,
    smoluchowski_transient_rate_constant,
)

D = 4e-9  # m^2/s, two small solutes in water
R_contact = 5e-10  # m
kD = smoluchowski_rate_constant(D, R_contact)
print(f"steady-state k_D = {kD * 1000:.3e} L/(mol s)")

fig, axes = plt.subplots(1, 2, figsize=(12, 4.5))

# %%
# The transient: at first the reactants' neighbourhoods have not yet been
# depleted, so encounters are more frequent than at steady state.

t = np.logspace(-13, -7, 200)
axes[0].semilogx(t, smoluchowski_transient_rate_constant(D, R_contact, t) * 1000, color="steelblue", label="k(t)")
axes[0].axhline(kD * 1000, color="gray", linestyle="--", label="steady state $k_D = 4\\pi D R^* N_A$")
axes[0].axvline(R_contact**2 / (np.pi * D), color="crimson", linestyle=":", label=r"$t = R^{*2}/\pi D$ (k = 2 k_D)")
axes[0].set_xlabel("t after mixing (s)")
axes[0].set_ylabel("k (L mol$^{-1}$ s$^{-1}$)")
axes[0].set_title("Smoluchowski transient decays to the steady rate")
axes[0].legend(fontsize=8)

# %%
# Viscosity sets the ceiling: 8RT/(3 eta) for several solvents at 298 K
# (approximate literature viscosities).

solvents = {"diethyl ether": 2.2e-4, "acetone": 3.1e-4, "water": 8.9e-4, "ethanol": 1.07e-3, "ethylene glycol": 1.6e-2, "glycerol": 0.95}
names = list(solvents)
k_solvent = np.array([diffusion_limited_rate_constant(298.15, eta) * 1000 for eta in solvents.values()])
axes[1].barh(names, k_solvent, color="steelblue")
axes[1].set_xscale("log")
axes[1].set_xlabel("diffusion-limited k (L mol$^{-1}$ s$^{-1}$)")
axes[1].set_title("Diffusion limit 8RT/(3 eta) at 298 K")
for name, kv in zip(names, k_solvent):
    print(f"{name:16s}: {kv:.2e} L/(mol s)")

fig.tight_layout()
plt.show()
