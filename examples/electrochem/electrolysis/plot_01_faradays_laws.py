r"""
Faraday's laws of electrolysis, and galvanic vs. electrolytic cells
=======================================================================

:func:`~chemistrykit.electrochem.systems.electrolysis.faradays_law_mass`
gives the mass of a species deposited or consumed by a constant current
over time. This example verifies the linear proportionality to both
current and time (Faraday's first law), computes the minimum applied
voltage needed to force a non-spontaneous (electrolytic) reaction to
run, and contrasts a spontaneous galvanic reaction against one that
needs external electrolysis.
"""

# %%
import matplotlib.pyplot as plt
import numpy as np

from chemistrykit.electrochem.systems.electrolysis import faradays_law_mass, minimum_applied_voltage_electrolytic
from chemistrykit.electrochem.systems.standard_potentials import (
    STANDARD_REDUCTION_POTENTIALS,
    cell_potential,
    is_spontaneous,
)

# %%
# Copper electrorefining: Cu2+ + 2e- -> Cu (M = 63.546 g/mol).
molar_mass_Cu, n_Cu = 63.546, 2
currents = np.linspace(1.0, 20.0, 8)
time_s = 3600.0
mass_deposited = faradays_law_mass(currents, time_s, molar_mass_Cu, n_Cu)

fig, ax = plt.subplots()
ax.plot(currents, mass_deposited, "o-")
ax.set_xlabel("Current (A)")
ax.set_ylabel("Mass Cu deposited in 1 h (g)")
ax.set_title("Faraday's first law: mass proportional to current")
fig.tight_layout()

# %%
# A galvanic (spontaneous) reaction: the Daniell cell.
cu = STANDARD_REDUCTION_POTENTIALS["Cu2+/Cu"]
zn = STANDARD_REDUCTION_POTENTIALS["Zn2+/Zn"]
E_daniell = cell_potential(cu, zn)
print(f"Daniell cell (Cu cathode / Zn anode): E = {E_daniell:.2f} V, spontaneous = {is_spontaneous(E_daniell)}")
print("-> runs galvanically, delivering electrical work.")

# %%
# The *reverse* pairing (Zn cathode / Cu anode) is non-spontaneous and
# needs to be driven electrolytically.
E_reverse = cell_potential(zn, cu)
V_min = minimum_applied_voltage_electrolytic(E_reverse)
print(f"\nReversed pairing (Zn cathode / Cu anode): E = {E_reverse:.2f} V, spontaneous = {is_spontaneous(E_reverse)}")
print(f"-> requires at least {V_min:.2f} V applied externally (electrolysis) to run.")

plt.show()
