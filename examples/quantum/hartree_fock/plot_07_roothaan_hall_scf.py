r"""
Roothaan-Hall self-consistent field: RHF for HeH+ and H2
========================================================

With two or more electrons the Fock matrix depends on the orbitals being
solved for, so the Roothaan-Hall equations :math:`FC=SC\varepsilon` must
be iterated. :class:`~chemistrykit.quantum.RestrictedHartreeFock` starts
from the core Hamiltonian (no electron repulsion), builds the density
matrix from the occupied orbital, rebuilds :math:`F`, and repeats until
the density stops changing. Left: the SCF energy of HeH+ (STO-3G, 1.4632
bohr) converging to Szabo & Ostlund's -2.860662 hartree. Right: the RHF
potential curve of H2. It has a good minimum near 1.35 bohr but rises
far above two separate hydrogen atoms at large :math:`R`. That is the
well-known failure of a single doubly occupied orbital to describe bond
breaking, since it keeps an ionic H+ H- component at every distance.
"""

# %%
import matplotlib.pyplot as plt
import numpy as np

from chemistrykit.quantum import RestrictedHartreeFock
from chemistrykit.quantum.systems.helium import HARTREE_ENERGY
from chemistrykit.quantum.utils.basis_sets import BOHR_RADIUS, contracted_one_electron, kinetic_integral, nuclear_attraction_integral, sto3g_1s

heh = RestrictedHartreeFock.heh_plus().scf()
print(f"HeH+ converged in {heh.n_iterations} iterations: E = {heh.total_energy / HARTREE_ENERGY:.6f} hartree")
print("orbital energies (hartree):", np.round(heh.orbital_energies / HARTREE_ENERGY, 6))

R = np.linspace(0.7, 6.0, 40)
E_h2 = np.array([RestrictedHartreeFock.h2(bond_length=r * BOHR_RADIUS).scf().total_energy for r in R]) / HARTREE_ENERGY

# One STO-3G hydrogen atom: a single electron, so its energy is just <T + V>.
phi = sto3g_1s(1.24, [0.0, 0.0, 0.0])
T_H = contracted_one_electron(kinetic_integral, phi, phi)
V_H = contracted_one_electron(nuclear_attraction_integral, phi, phi, 1.0, [0.0, 0.0, 0.0])
E_H = (T_H + V_H) / HARTREE_ENERGY

fig, axes = plt.subplots(1, 2, figsize=(12, 4.8))
axes[0].plot(np.arange(1, heh.n_iterations + 1), heh.energy_history / HARTREE_ENERGY, "o-", color="darkorange")
axes[0].axhline(-2.860662, color="gray", linestyle="--", label="Szabo & Ostlund")
axes[0].set_xlabel("SCF iteration")
axes[0].set_ylabel("total energy (hartree)")
axes[0].set_title("HeH+ SCF convergence (STO-3G)")
axes[0].legend()

axes[1].plot(R, E_h2, color="steelblue", label="RHF / STO-3G")
axes[1].axhline(2 * E_H, color="gray", linestyle="--", label="two H atoms, same basis")
axes[1].set_xlabel("H-H distance (bohr)")
axes[1].set_ylabel("total energy (hartree)")
axes[1].set_title("RHF H2: good near equilibrium, wrong at dissociation")
axes[1].legend()
fig.tight_layout()

# %%
i = int(np.argmin(E_h2))
print(f"H2 minimum near R = {R[i]:.2f} bohr, E = {E_h2[i]:.4f} hartree")
print(f"at R = {R[-1]:.1f} bohr RHF lies {E_h2[-1] - 2 * E_H:.3f} hartree above two H atoms")

plt.show()
