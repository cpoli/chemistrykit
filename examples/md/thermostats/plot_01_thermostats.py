r"""
Thermostats: NVE vs. velocity-rescaling vs. Nose-Hoover
==========================================================

Starting the same Lennard-Jones fluid configuration well away from its
target temperature, three runs compare how the instantaneous kinetic
temperature evolves: an unthermostatted (NVE, microcanonical) run, where
it simply stays wherever it started; the crude, deterministic
:class:`~chemistrykit.md.systems.thermostats.VelocityRescalingThermostat`,
which snaps it back to the target at every rescaling; and the
:class:`~chemistrykit.md.systems.thermostats.NoseHooverThermostat`, which
relaxes it smoothly toward the target via an extended-system friction
variable.
"""

# %%
import matplotlib.pyplot as plt

from chemistrykit.md.systems.lj_fluid import LJFluid
from chemistrykit.md.systems.thermostats import NoseHooverThermostat, VelocityRescalingThermostat

T_target = 1.0
dt, n_steps, sample_every = 0.002, 1000, 10


def make_perturbed_fluid():
    fluid = LJFluid.from_lattice(n_per_side=6, density=0.6, temperature=T_target, cutoff=2.5, rng=0)
    fluid.velocities *= 1.8  # perturb well away from the target temperature
    return fluid


runs = {
    "NVE (no thermostat)": (make_perturbed_fluid(), None),
    "velocity rescaling": (make_perturbed_fluid(), VelocityRescalingThermostat(target_temperature=T_target, interval=5)),
    "Nose-Hoover": (make_perturbed_fluid(), NoseHooverThermostat(target_temperature=T_target, Q=5.0)),
}

fig, ax = plt.subplots(figsize=(7, 5))
for label, (fluid, thermostat) in runs.items():
    result = fluid.run(dt=dt, n_steps=n_steps, thermostat=thermostat, sample_every=sample_every)
    ax.plot(result.t, result.temperature, label=label)

ax.axhline(T_target, color="gray", linestyle=":", linewidth=0.8, label="target T")
ax.set_xlabel("t")
ax.set_ylabel("instantaneous temperature")
ax.set_title("Thermostat comparison, starting well above the target temperature")
ax.legend()
fig.tight_layout()

# %%
# The NVE run's temperature stays near its (perturbed) starting value
# throughout -- with no thermostat, there is nothing to drive it back to
# `T_target`; both thermostats instead pull it there, the
# velocity-rescaling one abruptly and the Nose-Hoover one smoothly.

plt.show()
