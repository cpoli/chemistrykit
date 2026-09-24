r"""
Pauling's resonance energy as Huckel delocalization energy
============================================================

Linus Pauling's *The Nature of the Chemical Bond* made "resonance energy"
the chemist's measure of the extra stability a molecule gains because its
true structure is a superposition of several Lewis structures (benzene's
two Kekule structures foremost). In Huckel theory that stabilization is
the delocalization energy
(:meth:`~chemistrykit.quantum.systems.huckel.HuckelSystem.delocalization_energy`):
the computed pi energy minus that of the same electrons held in isolated,
ethylene-like double bonds, :math:`2(\alpha+\beta)` per bond.
"""

# %%
import matplotlib.pyplot as plt
import numpy as np

from chemistrykit.quantum.systems.huckel import HuckelSystem

molecules = {
    "ethylene": (HuckelSystem.linear_polyene(2), 2),
    "butadiene": (HuckelSystem.linear_polyene(4), 4),
    "hexatriene": (HuckelSystem.linear_polyene(6), 6),
    "cyclobutadiene": (HuckelSystem.cyclic_polyene(4), 4),
    "benzene": (HuckelSystem.cyclic_polyene(6), 6),
}
resonance = {}
for name, (system, n_pi) in molecules.items():
    resonance[name] = system.delocalization_energy(n_pi)
    print(f"{name:15s} pi energy = {system.pi_electron_energy(n_pi):7.4f} |beta|, resonance energy = {resonance[name]:7.4f}")

# %%
# With beta = -1, benzene's resonance energy is exactly 2|beta|, far more
# than open-chain hexatriene with the same six pi electrons, and
# cyclobutadiene gets none at all. The thermochemical benzene
# benzene resonance energy of about 36 kcal/mol (from heats of
# hydrogenation) corresponds to |beta| near 18 kcal/mol on this scale.

fig, ax = plt.subplots(figsize=(7, 4))
names = list(resonance)
values = -np.array([resonance[n] for n in names])
ax.bar(names, values, color=["gray", "gray", "gray", "lightcoral", "seagreen"])
ax.set_ylabel(r"resonance stabilization $(-E_{deloc})/|\beta|$")
ax.set_title("Resonance (delocalization) energy: benzene stands out")
fig.tight_layout()

# %%
# Resonance energy per pi electron in linear polyenes creeps up with chain
# length but stays well below benzene's 1/3 |beta| per electron.

n_values = np.arange(2, 21, 2)
per_electron = [-HuckelSystem.linear_polyene(n).delocalization_energy(n) / n for n in n_values]
fig2, ax2 = plt.subplots(figsize=(6, 4))
ax2.plot(n_values, per_electron, "o-", label="linear polyenes")
ax2.axhline(-resonance["benzene"] / 6, color="seagreen", linestyle="--", label="benzene")
ax2.set_xlabel("number of pi electrons")
ax2.set_ylabel(r"resonance energy per electron $/|\beta|$")
ax2.set_title("Resonance energy per pi electron")
ax2.legend()
fig2.tight_layout()

plt.show()
