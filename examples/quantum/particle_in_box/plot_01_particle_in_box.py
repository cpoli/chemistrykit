r"""
Solving Schrodinger's equation for a particle in a box
========================================================

Plots the first few 1D particle-in-a-box wavefunctions
(:class:`~chemistrykit.quantum.systems.particle_in_box.ParticleInBox1D`),
each offset by its own energy level, then checks the textbook 3D
cubic-box degeneracy pattern
(:class:`~chemistrykit.quantum.systems.particle_in_box.ParticleInBox3D`).
These are the simplest exact solutions of Schrodinger's 1926
time-independent equation :math:`-\frac{\hbar^2}{2m}\nabla^2\psi=E\psi`
inside hard walls: quantized energies :math:`E_n=n^2h^2/8mL^2` appear
from the boundary conditions alone, with no quantization postulate.
"""

# %%
import matplotlib.pyplot as plt
import numpy as np

from chemistrykit.quantum.systems.particle_in_box import ParticleInBox1D, ParticleInBox3D
from chemistrykit.quantum.visualizers.quantum_plots import plot_particle_in_box_wavefunctions

box = ParticleInBox1D(length=1.0e-9)
for n in (1, 2, 3, 4):
    print(f"E_{n} = {box.energy(n) / 1.602176634e-19:.4f} eV, E_{n}/E_1 = {box.energy(n) / box.energy(1):.4f} (n^2 = {n**2})")

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

plt.show()
