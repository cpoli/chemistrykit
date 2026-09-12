r"""
Bond order: the Pauling length correlation, and the Huckel Coulson bond order
================================================================================

Two independent routes to a (generally non-integer) bond order: Pauling's
empirical bond-length correlation applied to benzene's intermediate C-C
bond length, and the Coulson bond order computed directly from Huckel
molecular-orbital coefficients for the same molecule -- the two
independent methods agree that benzene's C-C bond order sits between 1
and 2.
"""

# %%
import matplotlib.pyplot as plt
import numpy as np

from chemistrykit.quantum.systems.huckel import HuckelSystem
from chemistrykit.structure.systems.bonding import bond_order_from_length, coulson_pi_bond_order
from chemistrykit.structure.visualizers.structure_plots import plot_bond_order_correlation

# Route 1: Pauling's empirical correlation from the observed bond length.
single_bond_length = 1.54  # C-C single bond (ethane), angstrom
benzene_length = 1.397  # angstrom
pauling_order = bond_order_from_length(single_bond_length, benzene_length)
print(f"Pauling bond order from benzene's C-C length: {pauling_order:.3f}")

# %%
# Route 2: the Coulson pi bond order from Huckel theory (a totally
# independent, from-scratch molecular-orbital calculation) -- the
# textbook 2/3 pi bond order for every C-C bond in benzene, on top of the
# sigma-bond framework Huckel theory doesn't model at all:

benzene_huckel = HuckelSystem.cyclic_polyene(n_atoms=6)
result = benzene_huckel.solve()
order = np.argsort(result.energies)
occupations = np.zeros(6)
occupations[order[:3]] = 2.0  # 6 pi electrons fill the 3 bonding MOs

pi_bond_order = coulson_pi_bond_order(result.coefficients, occupations, 0, 1)
total_bond_order = 1.0 + pi_bond_order  # 1 sigma bond + the Huckel pi contribution
print(f"Huckel pi bond order (C1-C2): {pi_bond_order:.4f}")
print(f"Total bond order (sigma + pi): {total_bond_order:.4f}")

# %%
# Plot the Pauling correlation curve, marking benzene's position on it:

ax = plot_bond_order_correlation(single_bond_length=single_bond_length, label="Pauling correlation")
ax.scatter([pauling_order], [benzene_length], color="crimson", zorder=3, label="benzene C-C")
ax.legend()

plt.show()
