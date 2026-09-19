r"""
Particle-in-a-box wavefunctions, cubic-box degeneracy, and dye color
======================================================================

Plots the first few 1D particle-in-a-box wavefunctions
(:class:`~chemistrykit.quantum.systems.particle_in_box.ParticleInBox1D`),
each offset by its own energy level, then checks the textbook 3D
cubic-box degeneracy pattern
(:class:`~chemistrykit.quantum.systems.particle_in_box.ParticleInBox3D`),
and finally applies Kuhn's free-electron model
(:func:`~chemistrykit.quantum.systems.particle_in_box.conjugated_dye_absorption_wavelength`)
to predict how a cyanine dye's absorption color shifts with conjugation
length.
"""

# %%
import matplotlib.pyplot as plt
import numpy as np

from chemistrykit.quantum.systems.particle_in_box import ParticleInBox1D, ParticleInBox3D, conjugated_dye_absorption_wavelength
from chemistrykit.quantum.visualizers.quantum_plots import plot_particle_in_box_wavefunctions

box = ParticleInBox1D(length=1.0e-9)
for n in (1, 2, 3, 4):
    print(f"E_{n} = {box.energy(n) / 1.602176634e-19:.4f} eV")

fig, ax = plt.subplots(figsize=(7, 5))
plot_particle_in_box_wavefunctions(box, [1, 2, 3, 4], ax=ax, scale=2.0e-20)
ax.set_title("1D particle-in-a-box wavefunctions (1 nm box)")
fig.tight_layout()

# %%
# A cubic 3D box shows genuine degeneracy from symmetry alone: the three
# permutations of (2, 1, 1) share exactly the same energy.

cubic_box = ParticleInBox3D(Lx=1.0e-9, Ly=1.0e-9, Lz=1.0e-9)
e211 = cubic_box.energy(2, 1, 1)
e121 = cubic_box.energy(1, 2, 1)
e112 = cubic_box.energy(1, 1, 2)
print(f"E(2,1,1) = E(1,2,1) = E(1,1,2)? {np.allclose([e211, e121, e112], e211)}")

degeneracies = cubic_box.degeneracy(n_max=3)
fig2, ax2 = plt.subplots(figsize=(6, 4))
energies_ev = np.array(sorted(degeneracies.keys())) / 1.602176634e-19
counts = [degeneracies[e] for e in sorted(degeneracies.keys())]
ax2.bar(range(len(energies_ev)), counts, tick_label=[f"{e:.2f}" for e in energies_ev])
ax2.set_xlabel("energy (eV)")
ax2.set_ylabel("degeneracy")
ax2.set_title("Cubic-box level degeneracies")
fig2.tight_layout()

# %%
# Kuhn's free-electron model: treating a linear conjugated dye's pi
# electrons as particles in a 1D box predicts that longer conjugation
# (more pi electrons, longer effective box) red-shifts the absorption
# -- the qualitative basis of cyanine-dye color tuning.

chain_lengths_nm = np.linspace(0.6, 2.0, 8)
electron_counts = [6, 8, 10, 12, 14, 16, 18, 20]
wavelengths_nm = [
    conjugated_dye_absorption_wavelength(box_length=L * 1.0e-9, n_pi_electrons=n) * 1.0e9 for L, n in zip(chain_lengths_nm, electron_counts, strict=True)
]
for L, n, wl in zip(chain_lengths_nm, electron_counts, wavelengths_nm, strict=True):
    print(f"L = {L:.2f} nm, {n} pi electrons -> lambda = {wl:.0f} nm")

fig3, ax3 = plt.subplots(figsize=(6, 4))
ax3.plot(electron_counts, wavelengths_nm, "o-", color="darkorange")
ax3.set_xlabel("number of pi electrons")
ax3.set_ylabel("predicted absorption wavelength (nm)")
ax3.set_title("Kuhn free-electron model: longer conjugation red-shifts absorption")
fig3.tight_layout()

plt.show()
