r"""
Verlet's integrator: bounded energy error instead of drift
============================================================

Verlet's 1967 algorithm advances positions from the two previous
positions and the current force,

.. math::

   \vec r(t+\Delta t) = 2\vec r(t) - \vec r(t-\Delta t) + \frac{\vec F(t)}{m}\Delta t^2,

which is algebraically the same trajectory as the kick-drift-kick
"velocity Verlet" form used by
:meth:`~chemistrykit.md.core.base_system.MolecularDynamicsSystem.step`.
Because the scheme is time-reversible and symplectic, its energy error
stays bounded however long the run. This example integrates the same
Lennard-Jones fluid (:class:`~chemistrykit.md.systems.lj_fluid.LJFluid`)
with the same time step using velocity Verlet and using the naive
forward-Euler update :math:`\vec r \mathrel{+}= \vec v\,\Delta t`,
:math:`\vec v \mathrel{+}= (\vec F/m)\,\Delta t`, and follows the
total energy.
"""

# %%
import matplotlib.pyplot as plt
import numpy as np

from chemistrykit.md.systems.lj_fluid import LJFluid
from chemistrykit.md.utils.pbc import wrap_positions

dt, n_steps = 0.002, 1500


def make_fluid():
    return LJFluid.from_lattice(n_per_side=5, density=0.6, temperature=1.0, cutoff=2.4, rng=0)


verlet = make_fluid()
verlet_result = verlet.run(dt=dt, n_steps=n_steps, sample_every=1)
e_verlet = verlet_result.total_energy / verlet.positions.shape[0]

euler = make_fluid()
e_euler = []
for _ in range(n_steps + 1):
    forces, potential = euler.forces_and_potential(euler.positions)
    e_euler.append((euler.kinetic_energy() + potential) / euler.positions.shape[0])
    euler.positions = wrap_positions(euler.positions + dt * euler.velocities, euler.box_length)
    euler.velocities = euler.velocities + dt * forces / euler.masses[:, None]
e_euler = np.array(e_euler)

t = np.arange(n_steps + 1) * dt
fig, ax = plt.subplots(figsize=(7, 4.5))
ax.plot(t, e_euler - e_euler[0], color="firebrick", label="forward Euler")
ax.plot(t, e_verlet - e_verlet[0], color="steelblue", label="Verlet (velocity form)")
ax.set_xlabel("t")
ax.set_ylabel("change in total energy per particle")
ax.set_title("Same fluid, same time step: Verlet vs. forward Euler")
ax.legend()
fig.tight_layout()

# %%
# Forward Euler pumps energy into the fluid at a steady rate -- it heats
# up for no physical reason -- while Verlet's energy only jitters about
# its starting value:

print(f"energy change per particle after t = {t[-1]:.1f}:")
print(f"  forward Euler: {e_euler[-1] - e_euler[0]:+.4f}")
print(f"  Verlet:        {e_verlet[-1] - e_verlet[0]:+.4f} (max excursion {np.max(np.abs(e_verlet - e_verlet[0])):.4f})")

plt.show()
