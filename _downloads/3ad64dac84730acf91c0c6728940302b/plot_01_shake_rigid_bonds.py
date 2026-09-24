r"""
SHAKE: rigid bonds let a bending molecule take larger time steps
===================================================================

Ryckaert, Ciccotti, and Berendsen's SHAKE algorithm (1977) replaces stiff
bond springs by exact length constraints, restored after every Verlet
position update by iterative corrections along the old bond vectors
(:func:`~chemistrykit.md.systems.constraints.shake`). A water-like
triatomic is simulated twice: with stiff harmonic O-H springs
(:class:`~chemistrykit.md.systems.pair_potentials.HarmonicMolecule`),
and with the two bonds held rigid by SHAKE
(:class:`~chemistrykit.md.systems.constraints.ShakeMolecule`), keeping
the same flexible angle-bending term. Both runs use the same time step, chosen comfortably for the slow
bending motion but too long for the fast bond stretch.
"""

# %%
import matplotlib.pyplot as plt
import numpy as np

from chemistrykit.md.systems.constraints import ShakeMolecule
from chemistrykit.md.systems.pair_potentials import HarmonicMolecule

theta0, r0 = np.deg2rad(104.5), 1.0
masses = [1.0, 16.0, 1.0]
k_bond, k_angle = 3000.0, 50.0
positions = np.array([[np.sin(theta0 / 2), np.cos(theta0 / 2), 0.0], [0.0, 0.0, 0.0], [-np.sin(theta0 / 2), np.cos(theta0 / 2), 0.0]])
velocities = np.zeros_like(positions)
velocities[0] = [np.cos(theta0 / 2), -np.sin(theta0 / 2), 0.0]  # tangential kicks that open the angle
velocities[2] = [-np.cos(theta0 / 2), -np.sin(theta0 / 2), 0.0]
velocities -= np.average(velocities, axis=0, weights=masses)

angles = [(0, 1, 2, k_angle, theta0)]
flexible = HarmonicMolecule(positions, velocities, masses, bonds=[(0, 1, k_bond, r0), (1, 2, k_bond, r0)], angles=angles)
rigid = ShakeMolecule(positions, velocities, masses, constraints=[(0, 1, r0), (1, 2, r0)], angles=angles)

dt, n_steps = 0.04, 250


def angle_deg(p):
    v1, v2 = p[0] - p[1], p[2] - p[1]
    return np.degrees(np.arccos(v1 @ v2 / np.linalg.norm(v1) / np.linalg.norm(v2)))


records = {"flexible (stiff springs)": flexible, "rigid bonds (SHAKE)": rigid}
traces = {name: {"bond": [], "angle": [], "energy": []} for name in records}
alive = {name: True for name in records}
for _ in range(n_steps + 1):
    for name, mol in records.items():
        p = mol.positions
        bond = np.linalg.norm(p[0] - p[1])
        alive[name] = alive[name] and abs(bond - r0) < 0.5  # stop following a run once it has blown up
        tr = traces[name]
        tr["bond"].append(bond if alive[name] else np.nan)
        tr["angle"].append(angle_deg(p) if alive[name] else np.nan)
        tr["energy"].append(mol.kinetic_energy() + mol.potential_energy() if alive[name] else np.nan)
        if alive[name]:
            mol.step(dt)

t = np.arange(n_steps + 1) * dt
fig, axes = plt.subplots(1, 3, figsize=(14, 4.2))
for name, color in zip(records, ["firebrick", "steelblue"]):
    tr = traces[name]
    axes[0].plot(t, tr["bond"], color=color, label=name)
    axes[1].plot(t, tr["angle"], color=color, label=name)
    axes[2].plot(t, np.array(tr["energy"]) - tr["energy"][0], color=color, label=name)
axes[0].set_ylabel("O-H bond length")
axes[0].set_ylim(0.4, 1.6)
axes[1].set_ylabel("H-O-H angle (degrees)")
axes[2].set_ylabel("change in total energy")
axes[2].set_yscale("symlog", linthresh=1e-3)
axes[2].set_ylim(-10, 10)
for ax, title in zip(axes, ["Bond length", "Bond angle", "Energy conservation"]):
    ax.set_xlabel("t")
    ax.set_title(title)
axes[0].legend()
fig.tight_layout()

# %%
# The stiff O-H springs vibrate far faster than the bend: Verlet is only
# stable for :math:`\omega\Delta t < 2`, and at this time step the
# flexible bonds exceed that limit, so their energy runs away within a
# few steps. With SHAKE the bond lengths stay fixed to the solver
# tolerance, only the slow bend is left to integrate, and the energy only
# fluctuates within a few percent, with no drift:

rigid_bonds = np.array(traces["rigid bonds (SHAKE)"]["bond"])
bond_period = 2 * np.pi / np.sqrt(k_bond * (1 / 1.0 + 1 / 16.0))
print(f"stiff-bond vibration period {bond_period:.3f}: omega * dt = {2 * np.pi / bond_period * dt:.2f} (> 2)")
print(f"flexible run blew up after {np.argmax(np.isnan(traces['flexible (stiff springs)']['bond']))} steps")
print(f"SHAKE bond-length deviation: {np.max(np.abs(rigid_bonds - r0)):.1e}")
energy = np.array(traces["rigid bonds (SHAKE)"]["energy"])
print(f"SHAKE run: largest energy change {np.max(np.abs(energy - energy[0])):.1e} (total energy {energy[0]:.3f})")

plt.show()
