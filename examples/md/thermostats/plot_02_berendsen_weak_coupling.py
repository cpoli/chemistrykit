r"""
Berendsen's weak coupling: exponential relaxation to the bath temperature
============================================================================

Berendsen and coworkers (1984) couple the system to a heat bath at
:math:`T_0` by rescaling all velocities each step by

.. math::

   \lambda = \sqrt{1 + \frac{\Delta t}{\tau}\left(\frac{T_0}{T} - 1\right)},

so that the temperature obeys :math:`dT/dt = (T_0 - T)/\tau` and relaxes
exponentially with time constant :math:`\tau`
(:class:`~chemistrykit.md.systems.thermostats.BerendsenThermostat`). The
left panel cools a dilute, nearly ideal Lennard-Jones gas from
:math:`T=3` toward :math:`T_0=1` for three coupling constants and
overlays :math:`T_0 + (T(0)-T_0)e^{-t/\tau}`. The right panel shows the
method's known flaw in a dense liquid: once equilibrated, weak coupling
suppresses the kinetic-energy fluctuations a canonical ensemble should
have.
"""

# %%
import matplotlib.pyplot as plt
import numpy as np

from chemistrykit.md.systems.lj_fluid import LJFluid
from chemistrykit.md.systems.thermostats import BerendsenThermostat

T0, T_start, dt = 1.0, 3.0, 0.005

fig, axes = plt.subplots(1, 2, figsize=(11, 4.5))
for tau, color in [(0.1, "firebrick"), (0.5, "darkorange"), (2.0, "steelblue")]:
    gas = LJFluid.from_lattice(n_per_side=5, density=0.01, temperature=T_start, cutoff=2.5, rng=0)
    result = gas.run(dt=dt, n_steps=1200, thermostat=BerendsenThermostat(T0, tau=tau), sample_every=10)
    axes[0].plot(result.t, result.temperature, color=color, label=rf"$\tau$ = {tau}")
    axes[0].plot(result.t, T0 + (T_start - T0) * np.exp(-result.t / tau), color="black", linestyle=":", linewidth=0.9)
axes[0].axhline(T0, color="gray", linewidth=0.8)
axes[0].set_xlabel("t")
axes[0].set_ylabel("temperature")
axes[0].set_title("Dilute gas: exponential relaxation (dotted: theory)")
axes[0].legend()

# %%
# In a dense liquid, compare the spread of the kinetic energy under
# Berendsen coupling with the canonical value
# :math:`\sigma_K = \bar K\sqrt{2/N_f}`:

liquid = LJFluid.from_lattice(n_per_side=5, density=0.8, temperature=T0, cutoff=2.4, rng=1)
thermostat = BerendsenThermostat(T0, tau=0.1)
liquid.run(dt=dt, n_steps=1000, thermostat=thermostat, sample_every=1000)
result = liquid.run(dt=dt, n_steps=3000, thermostat=thermostat, sample_every=5)
nf = liquid.degrees_of_freedom()
k_bar = 0.5 * nf * T0
axes[1].plot(result.t, result.kinetic_energy, color="firebrick", label=r"Berendsen, $\tau$ = 0.1")
axes[1].axhspan(k_bar * (1 - np.sqrt(2 / nf)), k_bar * (1 + np.sqrt(2 / nf)), color="gray", alpha=0.25, label=r"canonical $\bar K \pm \sigma_K$")
axes[1].axhline(k_bar, color="black", linestyle=":", linewidth=0.8)
axes[1].set_xlabel("t")
axes[1].set_ylabel("kinetic energy")
axes[1].set_title("Dense liquid: suppressed fluctuations")
axes[1].legend()
fig.tight_layout()

print(f"kinetic-energy standard deviation: Berendsen {result.kinetic_energy.std():.2f}, canonical {k_bar * np.sqrt(2 / nf):.2f}")

plt.show()
