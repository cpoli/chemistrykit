r"""
Green-Kubo: the diffusion coefficient as the integral of the velocity autocorrelation
========================================================================================

Green (1954) and Kubo (1957) showed that transport coefficients are time
integrals of equilibrium correlation functions. For self-diffusion the
relevant function is the velocity autocorrelation :math:`C(t)=\langle\vec
v(0)\cdot\vec v(t)\rangle`, and

.. math::

   D = \frac{1}{3}\int_0^\infty \langle \vec v(0)\cdot\vec v(t)\rangle\, dt .

This example computes :math:`C(t)` for a dense Lennard-Jones liquid
(:func:`~chemistrykit.md.systems.transport.velocity_autocorrelation`)
and the running Green-Kubo integral
(:func:`~chemistrykit.md.systems.transport.green_kubo_diffusion_coefficient`),
whose plateau is :math:`D`, and checks it against the Einstein
mean-squared-displacement estimate from the same trajectory.
"""

# %%
import matplotlib.pyplot as plt
import numpy as np

from chemistrykit.md.systems.lj_fluid import LJFluid
from chemistrykit.md.systems.thermostats import VelocityRescalingThermostat
from chemistrykit.md.systems.transport import (
    einstein_diffusion_coefficient,
    green_kubo_diffusion_coefficient,
    mean_squared_displacement,
    velocity_autocorrelation,
)

T, dt, sample_every = 1.0, 0.004, 2
fluid = LJFluid.from_lattice(n_per_side=6, density=0.8, temperature=T, cutoff=2.5, rng=0)
fluid.run(dt=dt, n_steps=2500, thermostat=VelocityRescalingThermostat(T, interval=10), sample_every=2500)
result = fluid.run(dt=dt, n_steps=6000, sample_every=sample_every)

max_lag = 400
lag_t = np.arange(max_lag + 1) * dt * sample_every
vacf = velocity_autocorrelation(result.velocities, max_lag=max_lag)
running_D = np.array([green_kubo_diffusion_coefficient(lag_t[: k + 1], vacf[: k + 1]) for k in range(1, max_lag + 1)])
D_gk = float(np.mean(running_D[-100:]))

msd = mean_squared_displacement(result.positions, box_length=result.box_length, max_lag=max_lag)
D_einstein = einstein_diffusion_coefficient(lag_t, msd, fit_from=1.0)

fig, axes = plt.subplots(1, 2, figsize=(11, 4.5))
axes[0].plot(lag_t, vacf / vacf[0], color="steelblue")
axes[0].axhline(0.0, color="gray", linestyle=":", linewidth=0.8)
axes[0].set_xlabel("t")
axes[0].set_ylabel(r"$C(t)/C(0)$")
axes[0].set_title("Velocity autocorrelation function")

axes[1].plot(lag_t[1:], running_D, color="steelblue", label="running Green-Kubo integral")
axes[1].axhline(D_einstein, color="black", linestyle="--", label=f"Einstein MSD estimate ({D_einstein:.3f})")
axes[1].set_xlabel("upper integration limit t")
axes[1].set_ylabel("D(t)")
axes[1].set_title("Green-Kubo integral converges to D")
axes[1].legend()
fig.tight_layout()

# %%
# In the dense liquid the VACF dips below zero: an atom rebounds off
# the cage of neighbors surrounding it. The negative lobe lowers the
# integral, and the plateau agrees with the independent Einstein
# estimate to within the statistical noise of a short run:

print(f"Green-Kubo D = {D_gk:.4f}, Einstein D = {D_einstein:.4f}")
print(f"most negative C(t)/C(0) = {np.min(vacf / vacf[0]):.3f}")

plt.show()
