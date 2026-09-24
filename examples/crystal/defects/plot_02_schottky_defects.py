r"""
Schottky defects: paired cation and anion vacancies in rock salt
==================================================================

Wagner and Schottky (1930) identified the second point-defect mechanism
of ionic crystals: a cation and an anion leave the bulk together,
keeping the crystal neutral and leaving a vacancy pair behind. The
equilibrium number is :math:`n_S=N\exp(-\Delta H_S/2k_BT)`
(:func:`~chemistrykit.crystal.systems.defects.schottky_defect_concentration`).
In NaCl, the classic Schottky-dominated solid, :math:`\Delta H_S\approx2.3` eV.
"""

# %%
import matplotlib.pyplot as plt
import numpy as np

from chemistrykit.constants import ELEMENTARY_CHARGE
from chemistrykit.crystal.systems.defects import schottky_defect_concentration
from chemistrykit.crystal.visualizers.crystal_plots import plot_defect_concentration_vs_temperature

N = 2.24e22  # NaCl formula units per cm^3
T = np.linspace(300.0, 1070.0, 200)  # NaCl melts at 1074 K

for dh_eV in (2.0, 2.3, 2.6):
    n_melt = schottky_defect_concentration(N, dh_eV * ELEMENTARY_CHARGE, 1070.0)
    print(f"Delta H_S = {dh_eV:.1f} eV: vacancy-pair fraction just below melting = {n_melt / N:.2e}")

# %%
# Each Schottky defect is one Na+ vacancy plus one Cl- vacancy, so the
# crystal stays stoichiometric -- the ratio of the two vacancy counts is 1
# at every temperature -- while the total grows by many orders of
# magnitude from room temperature to the melting point:
dh = 2.3 * ELEMENTARY_CHARGE
n_300 = schottky_defect_concentration(N, dh, 300.0)
n_1070 = schottky_defect_concentration(N, dh, 1070.0)
print(f"\nn_S(300 K) = {n_300:.2e} per cm^3,  n_S(1070 K) = {n_1070:.2e} per cm^3,  ratio {n_1070 / n_300:.1e}")

# %%
fig, ax = plt.subplots()
for dh_eV in (2.0, 2.3, 2.6):
    n = schottky_defect_concentration(N, dh_eV * ELEMENTARY_CHARGE, T)
    plot_defect_concentration_vs_temperature(T, n, ax=ax, label=rf"$\Delta H_S$ = {dh_eV:.1f} eV")
ax.set_title("Schottky vacancy pairs in NaCl")
plt.tight_layout()
plt.show()
