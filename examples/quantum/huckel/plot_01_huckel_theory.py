r"""
Huckel theory: butadiene, benzene, and the 4n+2 aromaticity rule
====================================================================

Builds the Huckel secular matrix for butadiene and benzene
(:class:`~chemistrykit.quantum.systems.huckel.HuckelSystem`), diagonalizes
it via :func:`numpy.linalg.eigh` (a genuine eigenvalue problem, not a
lookup table), cross-checks the result against the closed-form
Coulson/Frost eigenvalue formulas, and checks Huckel's 4n+2 rule
(:func:`~chemistrykit.quantum.systems.huckel.is_aromatic_by_huckel_rule`)
against the *computed* spectrum's degeneracies.
"""

# %%
import matplotlib.pyplot as plt
import numpy as np

from chemistrykit.quantum.systems.huckel import (
    HuckelSystem,
    cyclic_polyene_eigenvalues,
    is_aromatic_by_huckel_rule,
    linear_polyene_eigenvalues,
)
from chemistrykit.quantum.visualizers.quantum_plots import plot_huckel_levels

butadiene = HuckelSystem.linear_polyene(4)
result_butadiene = butadiene.solve()
print("Butadiene Huckel eigenvalues (alpha=0, beta=-1):")
print(f"  numerical:  {np.sort(result_butadiene.energies)}")
print(f"  closed form: {linear_polyene_eigenvalues(4)}")
print(f"  pi-electron energy (4 pi electrons): {butadiene.pi_electron_energy(4):.4f}")

# %%
benzene = HuckelSystem.cyclic_polyene(6)
result_benzene = benzene.solve()
print("\nBenzene Huckel eigenvalues (alpha=0, beta=-1):")
print(f"  numerical:   {np.sort(result_benzene.energies)}")
print(f"  closed form: {cyclic_polyene_eigenvalues(6)}")
print("  (the textbook pattern: alpha+2beta, alpha+beta x2, alpha-beta x2, alpha-2beta)")

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(9, 5))
plot_huckel_levels(result_butadiene, n_pi_electrons=4, ax=ax1)
ax1.set_title("Butadiene")
plot_huckel_levels(result_benzene, n_pi_electrons=6, ax=ax2)
ax2.set_title("Benzene")
fig.tight_layout()

# %%
# Huckel's 4n+2 rule, checked against the *computed* spectrum -- not just
# electron counting: benzene (6 pi electrons) is aromatic, cyclobutadiene
# (4 pi electrons) is not, and cyclooctatetraene (8 pi electrons) fails
# the electron-count test outright.

for name, n_atoms, n_pi in [("benzene", 6, 6), ("cyclobutadiene", 4, 4), ("cyclooctatetraene", 8, 8)]:
    energies = HuckelSystem.cyclic_polyene(n_atoms).solve().energies
    aromatic = is_aromatic_by_huckel_rule(n_pi, energies)
    print(f"{name}: {n_pi} pi electrons, aromatic by Huckel's rule? {aromatic}")

# %%
# The delocalization (resonance) energy quantifies benzene's extra
# stability relative to three isolated (localized) double bonds:

delocalization = benzene.delocalization_energy(6)
print(f"\nBenzene delocalization energy: {delocalization:.4f} beta (negative = stabilized)")

plt.show()
