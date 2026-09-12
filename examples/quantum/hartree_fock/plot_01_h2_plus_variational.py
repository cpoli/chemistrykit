r"""
H2+ in a minimal 2-Gaussian basis: the variational principle at work
========================================================================

:class:`~chemistrykit.quantum.systems.hartree_fock.H2PlusVariational`
builds the 2x2 Hamiltonian and overlap matrices for a minimal one-Gaussian-
per-proton basis, solves the secular equation ``HC=SCE``, and then
variationally optimizes the shared Gaussian exponent -- a genuine
application of the Rayleigh-Ritz variational theorem: minimizing over a
real parameter can only lower (never raise) the computed energy toward
the true ground state.
"""

# %%
import matplotlib.pyplot as plt
import numpy as np

from chemistrykit.constants import ELECTRONVOLT
from chemistrykit.quantum.systems.hartree_fock import H2PlusVariational

BOHR_RADIUS = 5.29177e-11
BOND_LENGTH = 106.0e-12  # H2+'s experimental equilibrium bond length

h2plus = H2PlusVariational(bond_length=BOND_LENGTH)
alpha_naive = 1.0 / BOHR_RADIUS**2  # a hydrogen-atom-1s-sized guess

result = h2plus.solve(alpha_naive)
print(f"Bonding MO energy:     {result.energies[0] / ELECTRONVOLT:.3f} eV")
print(f"Antibonding MO energy: {result.energies[1] / ELECTRONVOLT:.3f} eV")
print(f"Overlap S_AB at the naive exponent: {result.overlap[0, 1]:.4f}")

# %%
# Scanning the exponent shows the total energy has a genuine minimum --
# exactly what the Rayleigh-Ritz variational theorem promises for a
# well-posed trial wavefunction:

alphas = np.linspace(0.3, 3.0, 60) / BOHR_RADIUS**2
total_energies_ev = np.array([h2plus.total_energy(a) for a in alphas]) / ELECTRONVOLT

optimization = h2plus.optimize_exponent(alpha_naive)
print(f"\nNaive-guess total energy:     {optimization.naive_energy / ELECTRONVOLT:.3f} eV")
print(f"Optimized total energy:       {optimization.optimized_energy / ELECTRONVOLT:.3f} eV")
print(f"Improvement from optimizing:  {optimization.improvement / ELECTRONVOLT:.3f} eV")
print(f"Optimized exponent / naive exponent: {optimization.optimized_alpha / alpha_naive:.3f}")

fig, ax = plt.subplots(figsize=(7, 5))
ax.plot(alphas * BOHR_RADIUS**2, total_energies_ev, color="steelblue")
ax.axvline(alpha_naive * BOHR_RADIUS**2, color="gray", linestyle=":", label="naive guess (H-atom-like)")
ax.axvline(optimization.optimized_alpha * BOHR_RADIUS**2, color="crimson", linestyle="--", label="variationally optimized")
ax.set_xlabel(r"exponent $\alpha$ (units of $1/a_0^2$)")
ax.set_ylabel("total energy (eV)")
ax.set_title("H2+ total energy vs. Gaussian orbital exponent")
ax.legend()
fig.tight_layout()

# %%
# Repeating the exponent optimization at each bond length traces out a
# genuine molecular potential-energy curve, with a real minimum -- a
# bound H2+ molecular ion, not a hardcoded number:

bond_lengths_pm = np.linspace(60.0, 300.0, 25)
energies_ev = []
for R_pm in bond_lengths_pm:
    system = H2PlusVariational(bond_length=R_pm * 1.0e-12)
    opt = system.optimize_exponent(alpha_naive)
    energies_ev.append(opt.optimized_energy / ELECTRONVOLT)
energies_ev = np.array(energies_ev)

equilibrium_idx = np.argmin(energies_ev)
print(f"\nPredicted equilibrium bond length: {bond_lengths_pm[equilibrium_idx]:.0f} pm (experimental: 106 pm)")
print(f"Predicted total energy at equilibrium: {energies_ev[equilibrium_idx]:.2f} eV")

fig2, ax2 = plt.subplots(figsize=(6, 5))
ax2.plot(bond_lengths_pm, energies_ev, "o-", color="darkgreen")
ax2.axvline(bond_lengths_pm[equilibrium_idx], color="crimson", linestyle="--", label="predicted minimum")
ax2.set_xlabel("bond length (pm)")
ax2.set_ylabel("optimized total energy (eV)")
ax2.set_title("H2+ potential-energy curve (minimal single-Gaussian basis)")
ax2.legend()
fig2.tight_layout()

plt.show()
