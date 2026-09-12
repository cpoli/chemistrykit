r"""
Schottky and Frenkel defect equilibrium
==========================================

:func:`~chemistrykit.crystal.systems.defects.schottky_defect_concentration`
and
:func:`~chemistrykit.crystal.systems.defects.frenkel_defect_concentration`
both give a Boltzmann-factor defect population, :math:`n\propto\exp(-\Delta
H/2k_BT)` -- more defects at higher temperature, and fewer for a larger
formation enthalpy.
"""

# %%
import matplotlib.pyplot as plt
import numpy as np

from chemistrykit.crystal.systems.defects import frenkel_defect_concentration, schottky_defect_concentration
from chemistrykit.crystal.visualizers.crystal_plots import plot_defect_concentration_vs_temperature

N = 1.0e22  # cation sites (illustrative)
N_interstitial = 1.0e21  # available interstitial sites
delta_h_schottky = 2.0e-19  # J, ~1.25 eV, typical alkali-halide scale
delta_h_frenkel = 3.2e-19  # J, ~2 eV, typically larger than Schottky in ionic solids

T = np.linspace(300.0, 1200.0, 200)
n_schottky = schottky_defect_concentration(N, delta_h_schottky, T)
n_frenkel = frenkel_defect_concentration(N, N_interstitial, delta_h_frenkel, T)

for T_check in (300.0, 600.0, 900.0, 1200.0):
    ns = schottky_defect_concentration(N, delta_h_schottky, T_check)
    nf = frenkel_defect_concentration(N, N_interstitial, delta_h_frenkel, T_check)
    print(f"T={T_check:6.0f} K   Schottky n/N={ns / N:.3e}   Frenkel n/N={nf / N:.3e}")

# %%
# The defect population grows by many orders of magnitude between room
# temperature and typical solid-state processing temperatures:
schottky_ratio = schottky_defect_concentration(N, delta_h_schottky, 1200.0) / schottky_defect_concentration(N, delta_h_schottky, 300.0)
print(f"\nSchottky ratio (1200 K / 300 K): {schottky_ratio:.3e}")

# %%
fig, ax = plt.subplots()
plot_defect_concentration_vs_temperature(T, n_schottky, ax=ax, label="Schottky")
plot_defect_concentration_vs_temperature(T, n_frenkel, ax=ax, label="Frenkel")
plt.tight_layout()
plt.show()
