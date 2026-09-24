r"""
Staudinger's macromolecules: solution viscosity grows with chain length
==========================================================================

Staudinger's case that rubber and cellulose are genuine covalent
macromolecules rested heavily on dilute-solution viscosity: he found the
specific viscosity per unit concentration rising in proportion to the
chain length, :math:`\eta_\text{sp}/c=K_mM`
(:func:`~chemistrykit.polymer.systems.solution_viscosity.staudinger_specific_viscosity`).
If the chains were colloidal aggregates held together by weak
association, the "particle size" would fall apart on dilution or change
of solvent; covalent chains keep their length, so a polymer-analogous
reaction (e.g. hydrogenating rubber) leaves the degree of polymerization
-- and hence :math:`\eta_\text{sp}/c` -- unchanged. This example contrasts
the two pictures.
"""

# %%
import matplotlib.pyplot as plt
import numpy as np

from chemistrykit.polymer.systems.solution_viscosity import staudinger_specific_viscosity

Km = 1.0e-4  # Staudinger constant (per molar-mass unit, per concentration unit)
M = np.array([1e4, 3e4, 1e5, 3e5])
c = np.linspace(0.001, 0.01, 10)

for Mi in M:
    ratio = staudinger_specific_viscosity(c, Mi, Km) / c
    print(f"M = {Mi:8.0f}:  eta_sp/c = {ratio[0]:.2f} at every concentration (constant: {np.allclose(ratio, ratio[0])})")

# %%
# Colloidal-aggregate picture: particles of 1e5 units dissociate on
# dilution (a simple association equilibrium, apparent size ~ sqrt(c)),
# so eta_sp/c would *fall* as the solution is diluted, unlike a real
# macromolecule.
M_apparent = 1e5 * np.sqrt(c / c[-1])
eta_colloid = staudinger_specific_viscosity(c, M_apparent, Km) / c
eta_macro = staudinger_specific_viscosity(c, 1e5, Km) / c

# %%
fig, axes = plt.subplots(1, 2, figsize=(11, 4.2))
M_scan = np.logspace(3.5, 6, 50)
axes[0].loglog(M_scan, staudinger_specific_viscosity(1.0, M_scan, Km), "k-")
axes[0].loglog(M, staudinger_specific_viscosity(1.0, M, Km), "o")
axes[0].set_xlabel("molar mass M")
axes[0].set_ylabel(r"$\eta_{sp}/c$")
axes[0].set_title(r"Staudinger's rule: $\eta_{sp}/c = K_m M$")

axes[1].plot(c, eta_macro, label="covalent macromolecule (Staudinger)")
axes[1].plot(c, eta_colloid, "--", label="dissociating colloidal aggregate")
axes[1].set_xlabel("concentration c")
axes[1].set_ylabel(r"$\eta_{sp}/c$")
axes[1].set_title("Dilution test: real molecules keep their size")
axes[1].legend()
plt.tight_layout()
plt.show()
