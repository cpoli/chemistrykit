r"""
The Nernst equation: standard cells, concentration cells, and activity corrections
=====================================================================================

:func:`~chemistrykit.electrochem.systems.nernst.nernst_potential` relates
a cell's potential under non-standard conditions to its standard
potential and reaction quotient. A concentration cell (built from
identical half-cells differing only in concentration) makes the entire
driving force come from this activity term. Finally,
:func:`~chemistrykit.electrochem.systems.nernst.nernst_potential_with_activity`
shows how Debye-Huckel activity-coefficient corrections (reused from
:mod:`chemistrykit.solutions.systems.activity`) pull the cell potential
away from the ideal, unit-activity-coefficient prediction as ionic
strength grows.
"""

# %%
import matplotlib.pyplot as plt
import numpy as np

from chemistrykit.electrochem.systems.nernst import (
    activity_corrected_reaction_quotient,
    concentration_cell_potential,
    nernst_potential,
)
from chemistrykit.electrochem.visualizers.electrochem_plots import plot_nernst_concentration_dependence

E_standard, n = 0.34, 2  # Cu2+/Cu half-reaction vs. SHE
Q_values = np.logspace(-3, 3, 7)
for Q in Q_values:
    print(f"Q={Q:9.3g}  E = {nernst_potential(E_standard, n, Q):+.4f} V")

# %%
# At Q=1 (unit activity everywhere), E reduces to E_standard exactly.
print(f"\nE_standard = {E_standard} V; Nernst at Q=1: {nernst_potential(E_standard, n, Q=1.0)} V")

ax = plot_nernst_concentration_dependence(E_standard=E_standard, n=n)
plt.tight_layout()

# %%
# A concentration cell: same electrode material on both sides, so
# E_standard = 0 identically -- the entire driving force is the Nernst
# activity term.
n_conc = 1
ratios = np.array([2.0, 5.0, 10.0, 50.0, 100.0])
E_conc = concentration_cell_potential(n=n_conc, C_cathode=ratios * 0.01, C_anode=0.01)
fig, ax2 = plt.subplots()
ax2.plot(ratios, E_conc, "o-")
ax2.set_xscale("log")
ax2.set_xlabel(r"$C_{cathode}/C_{anode}$")
ax2.set_ylabel("E (V)")
ax2.set_title("Concentration cell potential vs. concentration ratio")
fig.tight_layout()

# %%
# Activity-coefficient corrections: at higher ionic strength, a reaction
# quotient built from Debye-Huckel-corrected activities deviates from the
# raw-concentration (ideal) quotient -- but only once the two species
# carry *different* charge magnitudes (equal charge magnitude gives both
# species the same activity coefficient, which cancels exactly in the
# ratio).
concentrations = [0.1, 0.1]
charges = [1, -2]
stoich = [-1.0, 1.0]
Q_ideal = concentrations[1] / concentrations[0]
Q_corrected = activity_corrected_reaction_quotient(concentrations, charges, stoich)
print(f"\nRaw-concentration Q:        {Q_ideal:.6f}")
print(f"Activity-corrected Q:       {Q_corrected:.6f}")

plt.show()
